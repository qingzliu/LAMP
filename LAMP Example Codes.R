#### Notes: This script demonstrates the application of the Linear Analysis Mixed Model Phenotyper (LAMP) method 
####        to study longitudinal trajectories of mood scores in patients, using PHQ scores as an example.  
####        The example data used here is synthetic and exhibits various longitudinal patterns for different patients. 
####        In Step 4, the logistic regression analysis is illustrative and does not yield significant results due to the synthetic nature of the data. 


# Load necessary libraries
library(lme4)
library(dplyr)
library(ggplot2)
library(tidyr)

#### Generate Example Dataset
# Simulate a dataset with PHQ scores, gender, antidepressant frequency, and years from baseline.
set.seed(123)
example_data <- data.frame(
  record_id = rep(1:100, each = 10),
  years_from_baseline = rep(0:9, times = 100),
  gender = rep(sample(c("Male", "Female"), 100, replace = TRUE), each = 10),
  Antidepressants_Freq = rep(runif(100, min = 0, max = 1), each = 10),  # Constant for each ID
  PHQ = NA
)

# Assign PHQ trends with noise based on predefined patterns (Static, Worsening, Recovering).
example_data <- example_data %>%
  group_by(record_id) %>%
  mutate(
    PHQ = case_when(
      record_id %% 4 == 1 ~ round(10 + rnorm(n(), mean = 0, sd = 2)),  # Static
      record_id %% 4 == 2 ~ round(10 + 0.5 * years_from_baseline + rnorm(n(), mean = 0, sd = 2)),  # Worsening
      record_id %% 4 == 3 ~ round(10 - 0.5 * years_from_baseline + rnorm(n(), mean = 0, sd = 2)),  # Recovering
      record_id %% 4 == 0 ~ round(10 + rnorm(n(), mean = 0, sd = 2))   # Static
    )
  ) %>%
  ungroup()

#### LAMP Step 1: Fit Linear Mixed Model and Smooth Noisy Longitudinal PHQ Scores
# Smooth noisy longitudinal PHQ trajectories into linear trends using a linear mixed model.
lmm_example <- lmer(PHQ ~ years_from_baseline + Antidepressants_Freq + gender +
                      (1 + years_from_baseline | record_id),
                    data = example_data, REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))

# Display model coefficients
print(summary(lmm_example)$coefficients)

#### LAMP Step 2: Categorize Trajectories
# Create groups based on baseline severity and trajectory trends (Recovering, Static, Worsening).
example_data_summary <- example_data %>%
  group_by(record_id) %>%
  summarise(
    PHQ_intercept = mean(PHQ),
    PHQ_slope = lm(PHQ ~ years_from_baseline)$coefficients[2],
    gender = first(gender),  # Ensure gender is included
    .groups = 'drop'
  )

# Determine cutoff for slope groups
find_interval <- function(vec) {
  sorted_vec <- sort(abs(vec))
  target_length <- length(vec) / 2
  idx <- which.min(abs(seq_along(vec) - target_length))
  return(sorted_vec[idx])
}

PHQ_slope_cut <- find_interval(example_data_summary$PHQ_slope)
example_data_summary <- example_data_summary %>%
  mutate(
    PHQ_slope_group = cut(PHQ_slope, breaks = c(-Inf, -PHQ_slope_cut, PHQ_slope_cut, Inf),
                          labels = c("Recovering", "Static", "Worsening")),
    PHQ_intercept_group = cut(PHQ_intercept, breaks = c(-Inf, 5, 10, 15, Inf),
                              labels = c("Low", "Mild", "Moderate", "Severe"))
  )

# Merge categorized data back into the original dataset for visualization.
example_longitudinal_data <- merge(example_data, example_data_summary[, c("record_id", "PHQ_slope_group", "PHQ_intercept_group")], by = "record_id")

# Plot longitudinal trajectories grouped by severity and trend
ggplot(data = example_longitudinal_data, aes(x = years_from_baseline, y = PHQ, group = record_id, colour = record_id)) +
  geom_line(show.legend = FALSE) +
  facet_grid(factor(str_to_title(PHQ_intercept_group), levels = c("Severe", "Moderate", "Mild", "Low"), ordered = TRUE) ~ PHQ_slope_group) +
  ylab("PHQ-9") +
  xlab("Years") +
  ggtitle("Longitudinal PHQ Trajectories for Bipolar 1 Patients")

# Visualize Trends: Plot individual and average trajectories over time categorized by severity and trends.
time_points <- 0:15
example_trend_data <- example_data_summary %>%
  mutate(time = list(time_points)) %>%
  unnest(cols = c(time)) %>%
  mutate(value = PHQ_intercept + PHQ_slope * time)

avg_data <- example_trend_data %>%
  group_by(PHQ_intercept_group, PHQ_slope_group, time) %>%
  summarise(avg_value = mean(value), .groups = 'drop')

sample_size_data <- example_trend_data %>%
  group_by(PHQ_intercept_group, PHQ_slope_group, time) %>%
  summarise(N = n(), .groups = 'drop')

avg_data <- left_join(avg_data, sample_size_data, by = c("PHQ_intercept_group", "PHQ_slope_group", "time"))

avg_intercept_slope <- example_trend_data %>%
  group_by(PHQ_intercept_group, PHQ_slope_group) %>%
  summarise(avg_intercept = mean(PHQ_intercept),
            avg_slope = mean(PHQ_slope) %>% round(2), .groups = 'drop')

avg_data <- left_join(avg_data, avg_intercept_slope, by = c("PHQ_intercept_group", "PHQ_slope_group"))

# Plot trends
ggplot(data = example_trend_data, aes(x = time, y = value, group = record_id)) +
  geom_line(aes(color = record_id), show.legend = FALSE, color = "grey") +
  geom_line(data = avg_data, aes(y = avg_value, group = interaction(PHQ_intercept_group, PHQ_slope_group)), color = "blue") +
  facet_grid(factor(PHQ_intercept_group, levels = c("Severe", "Moderate", "Mild", "Low")) ~ PHQ_slope_group) +
  ylab("PHQ-9") +
  xlab("Years") +
  coord_cartesian(ylim = c(0, 35))

#### LAMP Step 3: Define LAMP Phenotype
# Define a binary phenotype (Stable/Recovering vs Resistant/Worsening) based on trajectory groups.
example_data_summary <- example_data_summary %>%
  mutate(
    New_Group = case_when(
      PHQ_slope_group == "Recovering" | (PHQ_slope_group == "Static" & PHQ_intercept_group == "Low") ~ "Stable/Recovering",
      PHQ_slope_group == "Worsening" | (PHQ_slope_group == "Static" & PHQ_intercept_group %in% c("Mild", "Moderate", "Severe")) ~ "Resistant/Worsening",
      TRUE ~ NA_character_
    )
  )

#### LAMP Step 4: Logistic Regression
# Conduct logistic regression to examine the associations between covariates and binary phenotype.
example_logit_data <- example_data_summary %>%
  filter(!is.na(New_Group)) %>%
  mutate(New_Group = factor(New_Group, levels = c("Stable/Recovering", "Resistant/Worsening")))

logit_model <- glm(New_Group ~ PHQ_intercept + PHQ_slope + gender, family = binomial(link = "logit"), data = example_logit_data)
summary(logit_model)


