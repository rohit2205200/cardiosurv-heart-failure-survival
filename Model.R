
library(survival)
library(survminer)
library(tidyverse)
library(gtsummary)
library(readr)
raw_data<- read_csv("D:/project folder/heart_failure_clinical_records_dataset.csv")
# Verify dimensions (should be 299 rows and 13 columns)
dim(raw_data)
head(raw_data)
# Quick sanity check
glimpse(raw_data)
df <- raw_data %>%
  mutate(
    # Ensure binary indicators are factor variables for proper Table 1 rendering
    sex = factor(sex, levels = c(0, 1), labels = c("Female", "Male")),
    high_blood_pressure = factor(high_blood_pressure, levels = c(0, 1), labels = c("No", "Yes")),
    diabetes = factor(diabetes, levels = c(0, 1), labels = c("No", "Yes")),
    smoking = factor(smoking, levels = c(0, 1), labels = c("No", "Yes")),
    event = DEATH_EVENT
  )

# Verify Event count (Rule of thumb check: 96 events allows 6 to 9 variables comfortably)
table(df$event)

table1 <- df %>%
  select(age, sex, ejection_fraction, serum_creatinine, serum_sodium, high_blood_pressure, event) %>%
  tbl_summary(
    by = event,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      age ~ "Age (Years)",
      sex ~ "Biological Sex",
      ejection_fraction ~ "Ejection Fraction (%)",
      serum_creatinine ~ "Serum Creatinine (mg/dL)",
      serum_sodium ~ "Serum Sodium (mEq/L)",
      high_blood_pressure ~ "History of Hypertension"
    )
  ) %>%
  add_p() %>% # Automatically runs t-tests / Wilcoxon / Chi-Square
  bold_labels()

table1 <- df %>%
  select(age, sex, ejection_fraction, serum_creatinine, serum_sodium, 
         high_blood_pressure, diabetes, smoking, event) %>%
  tbl_summary(
    by = event,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      age ~ "Age (Years)",
      sex ~ "Biological Sex",
      ejection_fraction ~ "Ejection Fraction (%)",
      serum_creatinine ~ "Serum Creatinine (mg/dL)",
      serum_sodium ~ "Serum Sodium (mEq/L)",
      high_blood_pressure ~ "Hypertension",
      diabetes ~ "Diabetes",
      smoking ~ "Smoking Status"
    )
  ) %>%
  add_p() %>%
  bold_labels() %>%
  bold_p(t = 0.05)

# Display Table 1 in Viewer
table1
# 1. Create a binary clinical threshold
df <- df %>%
  mutate(ef_group = ifelse(ejection_fraction < 30, "EF < 30% (High Risk)", "EF >= 30% (Normal/Mild)"))

# 2. Fit Kaplan-Meier Model
km_fit <- survfit(Surv(time, event) ~ ef_group, data = df)

# 3. Run the Log-Rank Hypothesis Test
log_rank_test <- survdiff(Surv(time, event) ~ ef_group, data = df)
print(log_rank_test)

# 4. Generate the KM Plot with Risk Table
km_plot <- ggsurvplot(
  km_fit,
  data = df,
  pval = TRUE,
  pval.coord = c(20, 0.2),
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  linetype = "strata",
  surv.median.line = "hv",
  palette = c("#E7B800", "#2E9FDF"),
  title = "Kaplan-Meier Survival Curves Stratified by Ejection Fraction",
  xlab = "Time (Days of Follow-up)",
  ylab = "Cumulative Survival Probability",
  legend.title = "Risk Stratum",
  ggtheme = theme_minimal(base_size = 12)
)

print(km_plot)

# ==============================================================================
# PHASE 3: MULTIVARIABLE COX PROPORTIONAL HAZARDS MODEL & DIAGNOSTICS
# ==============================================================================

# 1. Create a clean modeling subset to eliminate any name-matching or factor bugs
cox_data <- df %>%
  select(
    time, 
    event, 
    age, 
    sex, 
    ejection_fraction, 
    serum_creatinine, 
    serum_sodium, 
    high_blood_pressure
  ) %>%
  drop_na()

# 2. Fit the Multivariable Cox Proportional Hazards Model
cox_fit <- coxph(
  Surv(time, event) ~ age + sex + ejection_fraction + serum_creatinine + serum_sodium + high_blood_pressure,
  data = cox_data,
  x = TRUE,
  y = TRUE
)

# 3. Print Statistical Summary (Betas, HR, SE, and C-Index)
cat("=== COX MODEL SUMMARY ===\n")
summary(cox_fit)

# 4. Generate Publication Regression Table (using gtsummary)
cox_table <- tbl_regression(cox_fit, exponentiate = TRUE) %>%
  bold_p(t = 0.05) %>%
  bold_labels()

# Display table in RStudio Viewer
cox_table

# 5. Build Bulletproof Native Forest Plot (No ggforest dependency)
hr_data <- data.frame(
  Variable = c(
    "Age (Years)", 
    "Sex: Male vs Female", 
    "Ejection Fraction (%)", 
    "Serum Creatinine (mg/dL)", 
    "Serum Sodium (mEq/L)", 
    "Hypertension: Yes vs No"
  ),
  HR = exp(coef(cox_fit)),
  Lower_CI = exp(confint(cox_fit))[, 1],
  Upper_CI = exp(confint(cox_fit))[, 2]
)

# Reorder so variables display cleanly top-to-bottom
hr_data$Variable <- factor(hr_data$Variable, levels = rev(hr_data$Variable))

forest_plot <- ggplot(hr_data, aes(x = HR, y = Variable)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "firebrick", size = 0.8) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.25, color = "#2b7489", size = 1) +
  geom_point(size = 3.5, color = "#2b7489") +
  scale_x_log10(breaks = c(0.5, 0.8, 1.0, 1.2, 1.5, 2.0)) +
  labs(
    title = "Multivariable Cox Proportional Hazards Model",
    subtitle = "Adjusted Hazard Ratios (95% Confidence Intervals)",
    x = "Hazard Ratio (log scale, Null Effect = 1.0)",
    y = ""
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold", size = 11)
  )

# Render Forest Plot in Plots Pane
print(forest_plot)

# 6. Test the Proportional Hazards Assumption (Schoenfeld Residuals)
cat("\n=== TESTING PROPORTIONAL HAZARDS ASSUMPTION (cox.zph) ===\n")
ph_test <- cox.zph(cox_fit)
print(ph_test)

# 7. Plot Schoenfeld Residuals Diagnostics
ggcoxzph(ph_test)

if(!require(shiny)) install.packages("shiny")

library(shiny)
library(survival)
library(ggplot2)

# Ensure the baseline model is stored cleanly
cox_model_app <- cox_fit

# 1. User Interface
ui <- fluidPage(
  titlePanel(
    div(
      h2("Clinical Survival Trajectory & Prognostic Risk Estimator", style = "color: #1a365d; font-weight: bold;"),
      h5("Interactive Biostatistical Demonstration | Heart Failure Inpatient Cohort", style = "color: #4a5568;")
    )
  ),
  hr(),
  
  sidebarLayout(
    sidebarPanel(
      h4("Patient Baseline Covariates", style = "color: #2b6cb0; font-weight: bold;"),
      
      sliderInput("app_age", "Age (Years):", 
                  min = 40, max = 95, value = 60, step = 1),
      
      selectInput("app_sex", "Biological Sex:", 
                  choices = c("Male", "Female"), selected = "Male"),
      
      sliderInput("app_ef", "Left Ventricular Ejection Fraction (%):", 
                  min = 14, max = 80, value = 38, step = 1),
      
      sliderInput("app_creatinine", "Serum Creatinine (mg/dL):", 
                  min = 0.5, max = 9.0, value = 1.1, step = 0.1),
      
      sliderInput("app_sodium", "Serum Sodium (mEq/L):", 
                  min = 115, max = 150, value = 137, step = 1),
      
      selectInput("app_hbp", "History of Hypertension:", 
                  choices = c("No", "Yes"), selected = "No"),
      
      hr(),
      helpText("Adjust parameters to observe real-time recalculation of the patient's individual survival curve S(t | X).")
    ),
    
    mainPanel(
      h4("Predicted Individualized Survival Trajectory", style = "color: #2d3748; font-weight: bold;"),
      plotOutput("survPlot", height = "380px"),
      hr(),
      fluidRow(
        column(4, 
               wellPanel(
                 style = "background-color: #ebf8ff; border-color: #bee3f8; text-align: center;",
                 h5("90-Day Survival", style = "color: #2b6cb0; margin: 0;"),
                 h2(textOutput("prob90"), style = "color: #2b6cb0; font-weight: bold; margin: 5px 0;")
               )
        ),
        column(4, 
               wellPanel(
                 style = "background-color: #ebf8ff; border-color: #bee3f8; text-align: center;",
                 h5("180-Day Survival", style = "color: #2b6cb0; margin: 0;"),
                 h2(textOutput("prob180"), style = "color: #2b6cb0; font-weight: bold; margin: 5px 0;")
               )
        ),
        column(4, 
               wellPanel(
                 style = "background-color: #fffaf0; border-color: #feebc8; text-align: center;",
                 h5("Prognostic Risk Score", style = "color: #c05621; margin: 0;"),
                 h2(textOutput("riskTier"), style = "font-weight: bold; margin: 5px 0;")
               )
        )
      )
    )
  )
)

# 2. Server Logic (Fixed Plot Rendering)
server <- function(input, output) {
  
  # Reactive dataframe holding user inputs
  patient_profile <- reactive({
    data.frame(
      age = input$app_age,
      sex = factor(input$app_sex, levels = c("Female", "Male")),
      ejection_fraction = input$app_ef,
      serum_creatinine = input$app_creatinine,
      serum_sodium = input$app_sodium,
      high_blood_pressure = factor(input$app_hbp, levels = c("No", "Yes"))
    )
  })
  
  # Compute individualized survival fit
  patient_surv <- reactive({
    survfit(cox_model_app, newdata = patient_profile())
  })
  
  # Fixed Plot: uses standard geom_ribbon and geom_step
  output$survPlot <- renderPlot({
    sf <- patient_surv()
    plot_df <- data.frame(
      time = sf$time,
      surv = sf$surv,
      lower = sf$lower,
      upper = sf$upper
    )
    
    ggplot(plot_df, aes(x = time, y = surv)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#2b6cb0") +
      geom_step(color = "#2b6cb0", linewidth = 1.2) +
      geom_vline(xintercept = c(90, 180), linetype = "dashed", color = "gray50") +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
      scale_x_continuous(breaks = seq(0, 280, by = 30)) +
      labs(
        x = "Time Following Admission (Days)",
        y = "Cumulative Survival Probability S(t)",
        caption = "Dashed markers at 90 and 180 days. Shaded band indicates 95% Confidence Interval."
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold")
      )
  })
  
  # Calculate 90-day survival probability
  output$prob90 <- renderText({
    sf <- patient_surv()
    idx <- which.min(abs(sf$time - 90))
    paste0(round(sf$surv[idx] * 100, 1), "%")
  })
  
  # Calculate 180-day survival probability
  output$prob180 <- renderText({
    sf <- patient_surv()
    idx <- which.min(abs(sf$time - 180))
    paste0(round(sf$surv[idx] * 100, 1), "%")
  })
  
  # Compute linear predictor score & classify risk tier
  output$riskTier <- renderText({
    lp <- predict(cox_model_app, newdata = patient_profile(), type = "lp")
    if (lp < -0.3) {
      paste0(round(lp, 2), " (Low)")
    } else if (lp <= 0.5) {
      paste0(round(lp, 2), " (Moderate)")
    } else {
      paste0(round(lp, 2), " (High)")
    }
  })
}

# 3. Launch Application
shinyApp(ui = ui, server = server)

```

Run this block now. The window will open and display the step curve with its shaded confidence band with no external errors.
