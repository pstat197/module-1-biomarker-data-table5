#code for q1/q2
set.seed(10302025)
# Change Headers to Protein Acronyms
raw.biomarkers1 <- raw.biomarkers
colnames(raw.biomarkers1) <- raw.biomarkers1[1,]
raw.data <- raw.biomarkers1[-1,]
sample.proteins <- sample(colnames(raw.data), 5)
# RS of n=4 data
cat('Our proteins sampled are: ', sample.proteins)

# Checking distributions
long.data <- raw.data %>%
  select(all_of(sample.proteins)) %>%
  pivot_longer(cols = everything(), names_to = "Protein", values_to = "Value") %>%
  mutate(Value = as.numeric(Value))

# Histogram
long.data %>% 
  drop_na() %>% 
  ggplot(aes(x = Value)) +
  geom_histogram(color = "white", fill = "black", bins = 30) +
  facet_wrap(~ Protein, scales = "free_x") +
  labs(
    title = "Distributions of Selected Protein Biomarkers",
    x = "Protein Level",
    y = "Frequency"
  ) +
  theme_minimal()


set.seed(10302025)

long.data %>% 
  drop_na() %>% 
  ggplot(aes(x = log(Value))) +
  geom_histogram(color = "white", fill = "black", bins = 30) +
  facet_wrap(~ Protein, scales = "free_x") +
  labs(
    title = "Distributions of Selected Protein Biomarkers",
    x = "Log-Transformed Protein Level",
    y = "Frequency"
  ) +
  theme_minimal()

proteins <- setdiff(names(raw.biomarkers), 
                    c("Group", "Target Full Name"))

z.scores <- raw.biomarkers %>% 
  mutate(across(all_of(proteins), ~scale(as.numeric(.x))))


outlier_summary <- z.scores %>% 
  mutate(across(all_of(proteins), ~abs(.x) > 3)) %>% 
  mutate(n_outliers = rowSums(across(all_of(proteins)), na.rm = TRUE)) %>% 
  filter(!is.na(Group) & Group != '') %>%
  group_by(Group) %>% 
  summarise(
    mean_outliers = mean(n_outliers, na.rm = TRUE),
    median_outliers = median(n_outliers, na.rm = TRUE),
    sd_outliers = sd(n_outliers, na.rm = TRUE),
    max_outliers = max(n_outliers, na.rm = TRUE),
    .groups = "drop"
  )

outlier_summary %>% kable(caption = 'Outlier Distribution Table', digits=2)

#q3a
library(tidyverse)
library(tidymodels)
library(modelr)
library(rsample)
library(yardstick)

load("../data/biomarker-clean.RData")

biomarker <- biomarker_clean

# These are  5 proteins of interest that were used in the log regression example
s_star <- c("DERM", "RELT", "IgD", "PTN", "FSTL1")

biomarker <- biomarker %>%
  select(group, all_of(s_star)) %>%
  mutate(class = (group == "ASD")) %>%
  select(-group)

set.seed(123)

# 80% used for model fitting (training) and 20% held out for evaluation (testing).
partitions <- initial_split(biomarker, prop = 0.8)
train <- training(partitions)
test <- testing(partitions)

# Predicts the probability of ASD (class == TRUE) based on protein levels.
fit <- glm(class ~ ., data = train, family = binomial(link = "logit"))
summary(fit)

# Computing predictions on the test set
pred_df <- test %>%
  add_predictions(fit, type = "response") %>%
  mutate(pred_class = (pred > 0.5),
         group = factor(class, labels = c("TD", "ASD")),
         pred_group = factor(pred_class, labels = c("TD", "ASD")))

# Checking factor order
levels(pred_df$group)

# Defining evaluation metrics
panel_fn <- metric_set(accuracy, sensitivity, specificity, roc_auc)

results <- pred_df %>%
  panel_fn(truth = group,
           estimate = pred_group,
           pred,
           event_level = "second")

results


#3b 
top_n <- 20

dat <- biomarker_clean
dat$group <- factor(dat$group)
proteins <- setdiff(names(dat), c('group','ados'))

tt_res <- sapply(proteins, function(p){
  x <- dat[[p]]
  grp <- dat$group
  ok <- !is.na(x) & !is.na(grp)
  if(sum(ok) < 3) return(NA)
  t <- try(t.test(x[ok] ~ grp[ok]), silent=TRUE)
  if(inherits(t,'try-error')) return(NA)
  t$p.value
})
tt_df <- tibble::tibble(protein = proteins,
                        pvalue = as.numeric(tt_res)) %>%
  arrange(pvalue) %>%
  slice_head(n = top_n)

rf_dat <- dat %>% dplyr::select(dplyr::all_of(proteins))
rf_resp <- dat$group
rf_fit <- randomForest::randomForest(x = rf_dat, y = rf_resp, ntree = 1000, importance = TRUE)
rf_imp_mat <- randomForest::importance(rf_fit, type = 2)
rf_imp_val <- if(is.matrix(rf_imp_mat)) rf_imp_mat[,1] else rf_imp_mat
rf_df <- tibble::tibble(protein = names(rf_imp_val), importance = as.numeric(rf_imp_val)) %>%
  arrange(desc(importance)) %>%
  slice_head(n = top_n)

X <- as.matrix(rf_dat)
y <- as.numeric(dat$group) - 1
cv <- glmnet::cv.glmnet(X, y, family = 'binomial', alpha = 1, nfolds = 5)
coef_min <- as.matrix(coef(cv, s = 'lambda.min'))
coefs <- coef_min[-1,1]
lasso_df <- tibble::tibble(protein = proteins, coef = as.numeric(coefs)) %>%
  mutate(abscoef = abs(coef)) %>%
  arrange(desc(abscoef)) %>%
  slice_head(n = top_n)


method_comparison <- bind_rows(
  mutate(tt_df %>% mutate(rank = row_number()), method = "T-test"),
  mutate(rf_df %>% mutate(rank = row_number()), method = "Random Forest"),
  mutate(lasso_df %>% mutate(rank = row_number()), method = "LASSO")
)

ggplot(method_comparison, aes(x = rank, y = protein, color = method)) +
  geom_point(size = 3) +
  theme_minimal() +
  scale_x_continuous(breaks = 1:top_n) +
  labs(title = paste("Top", top_n, "Proteins Selected by Each Method"),
       x = "Rank within Method",
       y = "Protein",
       color = "Selection Method") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 8))


tt_top <- tt_df$protein
rf_top <- rf_df$protein
lasso_top <- lasso_df$protein

all_methods <- intersect(intersect(tt_top, rf_top), lasso_top)

two_plus_methods <- unique(c(
  intersect(tt_top, rf_top),
  intersect(tt_top, lasso_top),
  intersect(rf_top, lasso_top)
))

overlap_summary <- tibble::tribble(
  ~"Overlap Type", ~"Count", ~"Percentage",
  "Selected by all methods", length(all_methods), length(all_methods)/top_n*100,
  "Selected by >=2 methods", length(two_plus_methods), length(two_plus_methods)/top_n*100,
  "Unique proteins total", length(unique(c(tt_top, rf_top, lasso_top))), 
  length(unique(c(tt_top, rf_top, lasso_top)))/top_n*100
)

kable(overlap_summary, 
      caption = paste("Overlap Analysis of Top", top_n, "Proteins"),
      digits = 1)

#3c 
library(tidyverse)
library(randomForest)
library(yardstick)
library(modelr)
library(tidymodels)

#t-tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = FALSE)
}

ttests_out <- biomarker_clean %>%
  select(-ados) %>%
  pivot_longer(-group, names_to = 'protein', values_to = 'level') %>%
  nest(data = c(level, group)) %>%
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

cat("Top proteins from t-tests:\n")
print(proteins_s1)

#rf
predictors <- biomarker_clean %>% select(-c(group, ados))
response <- factor(biomarker_clean$group)

set.seed(101422)
rf_out <- randomForest(x = predictors, y = response, ntree = 1000, importance = TRUE)
cat("\nRandom Forest confusion matrix:\n")
print(rf_out$confusion)

proteins_s2 <- rf_out$importance %>%
  as_tibble(rownames = "protein") %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

cat("\nTop proteins from Random Forest:\n")
print(proteins_s2)

#fuzzy intersection
protein_counts <- table(c(proteins_s1, proteins_s2))
fuzzy_proteins <- names(protein_counts[protein_counts >= 2])

cat("\nProteins selected by fuzzy intersection:\n")
print(fuzzy_proteins)

#log regression on fuzzy intersection

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(fuzzy_proteins)) %>%
  mutate(class = factor(ifelse(group == "ASD", "ASD", "TD"))) %>%
  select(-group)

set.seed(101422)
biomarker_split <- initial_split(biomarker_sstar, prop = 0.8, strata = class)

fit <- glm(class ~ ., data = training(biomarker_split), family = "binomial")

class_metrics <- metric_set(sensitivity, specificity, accuracy, roc_auc)

results <- testing(biomarker_split) %>%
  add_predictions(fit, type = "response") %>%
  mutate(
    class = factor(class, levels = c("TD", "ASD")),
    pred_label = factor(ifelse(pred > 0.5, "ASD", "TD"), levels = c("TD", "ASD"))
  ) %>%
  class_metrics(truth = class, estimate = pred_label, pred, event_level = "second")

#cat("\nLogistic Regression Performance (Fuzzy Intersection):\n")
#print(results, n = Inf)
kable(results, caption = "Logistic Regression Performance (Fuzzy Intersection)")


#q4
library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(yardstick)

#partition
set.seed(101422)
biomarker_split <- initial_split(biomarker_clean, prop = 0.8)
train <- training(biomarker_split)
test  <- testing(biomarker_split)

#training data only

#multiple testing
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

# rf
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

#hard intersect on top features
proteins_sstar <- intersect(proteins_s1, proteins_s2)

#training log regression
biomarker_sstar <- train %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = group == "ASD") %>%
  select(-group)

fit <- glm(class ~ ., data = biomarker_sstar, family = "binomial")

#eval on test data
test_sstar <- test %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = group == "ASD") %>%
  select(-group)

# get predictions
test_sstar <- test_sstar %>%
  mutate(pred = predict(fit, newdata = test_sstar, type = "response"))

# metrics
class_metrics <- metric_set(sensitivity, specificity, accuracy, roc_auc)

# add prediction and classification columns
test_sstar_eval <- test_sstar %>%
  mutate(
    pred = predict(fit, newdata = test_sstar, type = "response"),
    estimate = factor(pred > 0.5, levels = c(FALSE, TRUE)),
    truth = factor(class, levels = c(FALSE, TRUE))
  )

# define metrics
class_metrics <- metric_set(sensitivity, specificity, accuracy, roc_auc)

# evaluate
test_sstar_eval %>%
  class_metrics(
    truth = truth,
    estimate = estimate,
    pred,
    event_level = "second"
  )


set.seed(101422)
biomarker_split <- initial_split(biomarker_clean, prop = 0.8)
train <- training(biomarker_split)
test  <- testing(biomarker_split)

#our approach

metrics_yours <- test_sstar_eval %>%
  summarise(
    sensitivity = sensitivity_vec(truth, estimate, event_level = "second"),
    specificity = specificity_vec(truth, estimate, event_level = "second"),
    accuracy    = accuracy_vec(truth, estimate),
    roc_auc     = roc_auc_vec(truth, pred, event_level = "second")
  ) %>%
  mutate(method = "Training-only selection")



#t testing on full data (og in-class approach)
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

# RF on full data
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

#hard intersection
proteins_sstar_full <- intersect(proteins_s1_full, proteins_s2_full)

#log regression
biomarker_sstar_full <- biomarker_clean %>%
  select(group, any_of(proteins_sstar_full)) %>%
  mutate(class = (group == "ASD")) %>%
  select(-group)

fit_full <- glm(class ~ ., data = biomarker_sstar_full, family = "binomial")

#eval model on our test set
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

#combine and compare
bind_rows(metrics_yours, metrics_full)


kable(dplyr::bind_rows(metrics_yours, metrics_full),
      caption = "Training vs Full-Data Performance",
      digits = 3)





















