library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(yardstick)

# load cleaned dataset
load("data/biomarker-clean.RData")

#---------------------------------------------
# 1. Partition into training and test sets first
#---------------------------------------------
set.seed(101422)
biomarker_split <- initial_split(biomarker_clean, prop = 0.8)
train <- training(biomarker_split)
test  <- testing(biomarker_split)

#---------------------------------------------
# 2. FEATURE SELECTION — training data only
#---------------------------------------------

## ---- MULTIPLE TESTING ----
test_fn <- function(.df) {
  t_test(
    .df,
    formula = level ~ group,
    order = c("ASD", "TD"),
    alternative = "two-sided",
    var.equal = FALSE
  )
}

ttests_out <- train %>%
  select(-ados) %>%
  pivot_longer(-group, names_to = "protein", values_to = "level") %>%
  nest(data = c(level, group)) %>%
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(
    m = n(),
    hm = log(m) + 1 / (2 * m) - digamma(1),
    rank = row_number(),
    p.adj = m * hm * p_value / rank
  )

# top 10 proteins by adjusted p-value
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## ---- RANDOM FOREST ----
predictors_train <- train %>%
  select(-c(group, ados))
response_train <- factor(train$group)

set.seed(101422)
rf_out <- randomForest(
  x = predictors_train,
  y = response_train,
  ntree = 1000,
  importance = TRUE
)

proteins_s2 <- rf_out$importance %>%
  as_tibble(rownames = "protein") %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

#---------------------------------------------
# 3. HARD INTERSECTION of top features
#---------------------------------------------
proteins_sstar <- intersect(proteins_s1, proteins_s2)

#---------------------------------------------
# 4. LOGISTIC REGRESSION — training only
#---------------------------------------------
biomarker_sstar <- train %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = group == "ASD") %>%
  select(-group)

fit <- glm(class ~ ., data = biomarker_sstar, family = "binomial")

#---------------------------------------------
# 5. Evaluate on test data (never seen before)
#---------------------------------------------
test_sstar <- test %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = group == "ASD") %>%
  select(-group)

# get predictions
test_sstar <- test_sstar %>%
  mutate(pred = predict(fit, newdata = test_sstar, type = "response"))

# metrics
class_metrics <- metric_set(sensitivity, specificity, accuracy, roc_auc)

# Add prediction and classification columns
test_sstar_eval <- test_sstar %>%
  mutate(
    pred = predict(fit, newdata = test_sstar, type = "response"),
    estimate = factor(pred > 0.5, levels = c(FALSE, TRUE)),
    truth = factor(class, levels = c(FALSE, TRUE))
  )

# Define metrics
class_metrics <- metric_set(sensitivity, specificity, accuracy, roc_auc)

# Evaluate
test_sstar_eval %>%
  class_metrics(
    truth = truth,
    estimate = estimate,
    pred,
    event_level = "second"
  )


## -----------------------------------------------
## BENCHMARK: Compare your "training-only" pipeline 
##             vs. in-class (full-data) pipeline
## -----------------------------------------------

#--------------------------------------------------
# 1. Load data and your existing train/test split
#--------------------------------------------------
load("data/biomarker-clean.RData")

set.seed(101422)
biomarker_split <- initial_split(biomarker_clean, prop = 0.8)
train <- training(biomarker_split)
test  <- testing(biomarker_split)

#--------------------------------------------------
# 2. YOUR PIPELINE (training-only feature selection)
#--------------------------------------------------
# (Assumes you already ran your full code up to test_sstar_eval)
# test_sstar_eval should already include: truth, estimate, pred

metrics_yours <- test_sstar_eval %>%
  summarise(
    sensitivity = sensitivity_vec(truth, estimate, event_level = "second"),
    specificity = specificity_vec(truth, estimate, event_level = "second"),
    accuracy    = accuracy_vec(truth, estimate),
    roc_auc     = roc_auc_vec(truth, pred, event_level = "second")
  ) %>%
  mutate(method = "Training-only selection")

#--------------------------------------------------
# 3. IN-CLASS PIPELINE (feature selection on full data)
#--------------------------------------------------

## ---- Multiple testing on full data ----
test_fn <- function(.df){
  t_test(.df,
         formula = level ~ group,
         order = c("ASD", "TD"),
         alternative = "two-sided",
         var.equal = FALSE)
}

ttests_out_full <- biomarker_clean %>%
  select(-ados) %>%
  pivot_longer(-group, names_to = "protein", values_to = "level") %>%
  nest(data = c(level, group)) %>%
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(
    m = n(),
    hm = log(m) + 1/(2*m) - digamma(1),
    rank = row_number(),
    p.adj = m * hm * p_value / rank
  )

proteins_s1_full <- ttests_out_full %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## ---- Random forest on full data ----
predictors_full <- biomarker_clean %>%
  select(-c(group, ados))
response_full <- factor(biomarker_clean$group)

set.seed(101422)
rf_out_full <- randomForest(
  x = predictors_full,
  y = response_full,
  ntree = 1000,
  importance = TRUE
)

proteins_s2_full <- rf_out_full$importance %>%
  as_tibble(rownames = "protein") %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

## ---- Hard intersection ----
proteins_sstar_full <- intersect(proteins_s1_full, proteins_s2_full)

## ---- Logistic regression ----
biomarker_sstar_full <- biomarker_clean %>%
  select(group, any_of(proteins_sstar_full)) %>%
  mutate(class = (group == "ASD")) %>%
  select(-group)

fit_full <- glm(class ~ ., data = biomarker_sstar_full, family = "binomial")

#--------------------------------------------------
# 4. Evaluate in-class model on YOUR held-out test set
#--------------------------------------------------
test_full <- test %>%
  select(group, any_of(proteins_sstar_full)) %>%
  mutate(class = (group == "ASD")) %>%
  select(-group)

test_full_eval <- test_full %>%
  mutate(
    pred = predict(fit_full, newdata = test_full, type = "response"),
    estimate = factor(pred > 0.5, levels = c(FALSE, TRUE)),
    truth = factor(class, levels = c(FALSE, TRUE))
  )

metrics_full <- test_full_eval %>%
  summarise(
    sensitivity = sensitivity_vec(truth, estimate, event_level = "second"),
    specificity = specificity_vec(truth, estimate, event_level = "second"),
    accuracy    = accuracy_vec(truth, estimate),
    roc_auc     = roc_auc_vec(truth, pred, event_level = "second")
  ) %>%
  mutate(method = "In-class (full-data selection)")

#--------------------------------------------------
# 5. Combine and compare
#--------------------------------------------------
bind_rows(metrics_yours, metrics_full)


