library(furrr)
plan(multisession, workers=12)
# --- 1. Create the Simulation Function ---
# This function runs a single instance of our DAG simulation and returns the
# estimated log-odds ratios from TND and Cohort models.

run_one_simulation <- function(n, true_beta_exp) {
  # --- Data Generation ---
  HSB <- rnorm(n, mean = 0, sd = 1)
  logit_E <- 2 + 1.5 * HSB
  Exposure <- rbinom(n, size = 1, prob = plogis(logit_E))
  logit_OthDis <- -2.0 + 1.2 * Exposure
  Other_diseases <- rbinom(n, 1, plogis(logit_OthDis))
  logit_D <- -3.5 + true_beta_exp * Exposure + 1 * HSB 
  D <- rbinom(n, 1, plogis(logit_D))
  logit_Symptoms <- -3.5 + 2 * D + 2 * Other_diseases
  Symptoms <- rbinom(n, 1, plogis(logit_Symptoms))
  Symptoms[D == 0 & Other_diseases == 0] <- 0
  logit_Testing <- -4.0 + 1 * HSB + 4 * Symptoms
  Testing <- rbinom(n, 1, plogis(logit_Testing))
  #Testing[Symptoms == 0 ] <- 0
  D_star <- D * Testing
  D_star_symp <- D * Testing * Symptoms
  sim_data <- data.frame(Exposure, HSB, Other_diseases, D, Symptoms, Testing, D_star,D_star_symp )
  tnd_data_symp <- sim_data |> filter(Testing == 1, Symptoms==1)
  tnd_data_all <- sim_data |> filter(Testing == 1)
  # --- Model Fitting ---
  true_model_D <- glm(D ~ Exposure + HSB, data = sim_data, family = binomial)
  true_model_D_OD <- glm(D ~ Exposure + HSB +Other_diseases, data = sim_data, family = binomial)
  tnd_symp <- glm(D_star ~ Exposure, data = tnd_data_symp, family = binomial)
  tnd_all <- glm(D_star ~ Exposure, data = tnd_data_all, family = binomial)
  cohort_symp <- glm(D_star_symp ~ Exposure, data = sim_data, family = binomial)
  cohort_all <- glm(D_star ~ Exposure, data = sim_data, family = binomial)
  # --- Extract and Return Coefficients ---
  coef_true_d <- coef(true_model_D)["Exposure"]
  coef_true_d_od <- coef(true_model_D_OD)["Exposure"]
  coef_tnd_symp <- coef(tnd_symp)["Exposure"]
  coef_cohort_symp <- coef(cohort_symp)["Exposure"]
  coef_cohort_all <- coef(cohort_all)["Exposure"]
  coef_tnd_all <- coef(tnd_all)["Exposure"]

  return(c(
    true_d = coef_true_d,
    true_d_od = coef_true_d_od,
    tnd_estimate_symp = coef_tnd_symp,
    tnd_estimate_all = coef_tnd_all,
    cohort_symp = coef_cohort_symp,
    cohort_all = coef_cohort_all
  ))
}



# 2. Run the Simulation
n_replications <- 1000
n_subjects_per_sim <- 100000
true_value <- 1.5   #logOR

# This may take a minute to run.
message("Running simulations...")
results_df_dag3 <- future_map_dfr(
  .x = 1:n_replications, 
  .f = ~run_one_simulation(n = n_subjects_per_sim, true_beta_exp = true_value),
  .options = furrr_options(seed = 42), # Ensures reproducibility
  .progress=T
)
message("Done.")
# 3 Results 
summary(results_df_dag3)
results_df_dag3 |> 
  reframe(
    across(everything(),sd)
  )
results_df_dag3 |>
  reframe(
    tnd_wo_asymp =tnd_estimate_symp.Exposure - true_value,
    tnd_w_asymp = tnd_estimate_all.Exposure - true_value,
    cohort_symp_bias = cohort_symp.Exposure - true_value,
    cohort_all_bias = cohort_all.Exposure - true_value,
  ) |>
  pivot_longer(cols = everything()) |>
  ggplot(aes(name, value)) +
  scale_y_continuous(limits = c(-2, 2)) +
  geom_hline(aes(yintercept = 0), linetype=2)+
  geom_boxplot() +
  labs(
    x = "Design",
    y = "Bias"
  )+
  scale_x_discrete(labels=c("Cohort with asymptomatic",
                            "Cohort - symptomatic",
                            "TND with asymptomatic",
                            "TND - symptomatic"))+
  theme_minimal()



HSB <- rnorm(n, mean = 0, sd = 1)
logit_E <- 2 + 1.5 * HSB
Exposure <- rbinom(n, size = 1, prob = plogis(logit_E))
logit_OthDis <- -2.0 + 1.2 * Exposure
Other_diseases <- rbinom(n, 1, plogis(logit_OthDis))
logit_D <- -3.5 + true_beta_exp * Exposure + 1 * HSB 
D <- rbinom(n, 1, plogis(logit_D))
logit_Symptoms <- -3.5 + 2 * D + 2 * Other_diseases
Symptoms <- rbinom(n, 1, plogis(logit_Symptoms))
Symptoms[D == 0 & Other_diseases == 0] <- 0
logit_Testing <- -4.0 + 1 * HSB + 4 * Symptoms
Testing <- rbinom(n, 1, plogis(logit_Testing))
#Testing[Symptoms == 0 ] <- 0
D_star <- D * Testing
D_star_symp <- D * Testing * Symptoms
sim_data <- data.frame(Exposure, HSB, Other_diseases, D, Symptoms, Testing, D_star,D_star_symp )
tnd_data_symp <- sim_data |> filter(Testing == 1, Symptoms==1)
tnd_data_all <- sim_data |> filter(Testing == 1)
glm(D ~ Exposure + HSB, data = sim_data, family = binomial)
glm(D ~ Exposure + HSB +Other_diseases, data = sim_data, family = binomial)
glm(D_star ~ Exposure, data = tnd_data_symp, family = binomial)
glm(D_star ~ Exposure + Symptoms, data = tnd_data_all, family = binomial)
glm(D_star_symp ~ Exposure, data = sim_data, family = binomial)
glm(D_star ~ Exposure, data = sim_data, family = binomial)
