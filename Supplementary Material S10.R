### Paper title: A method to estimate actual infrastructure-induced mortality by integrating sampling biases

### Supplementary Material S10

### R script to simulate and analyse road survey census data for the Prior Sensitivity Analysis for vertebrate
#groups only affected by carcass location bias (Mammals G5)

## R version 4.2.2

# -----------------------------------------------------------------------------
# SIMULATION WORKFLOW DESCRIPTION
# -----------------------------------------------------------------------------
# This script executes the Prior Sensitivity Analysis workflow described in 
# Section 2.4.3 of the main text for vertebrates groups only affected by carcass
# location bias(Mammals G5). We simulate data by fixing the daily mean number of roadkills
#(lambda) and carcass observation probability per survey method (po), while 
#varying carcass location probabilities (pL).

# The design follows a full factorial approach resulting in 9 distinct scenarios.
# We perform 20 simulations per scenario, totaling 180 simulations per vertebrate group.

#    BIOLOGICAL PARAMETERS VARIATION:
#    We simulate scenarios using empirically derived values from literature, 
#    as well as scenarios where pL is shifted by +/- 10% away from their 
#    values in literature. We selected a +/- 10% deviation because larger deviations 
#    would result in parameter combinations that are rarely observed in nature.
#    - pL levels:  Literature value, Literature + 10%, Literature - 10%.

#    PRIOR SPECIFICATIONS TESTED:
#    For each simulated scenario, we test the model's performance using three 
#    prior specifications for pL:
#
#    A. Accurate priors: 
#       Priors centered on the value used in the simulations.
#
#    B. Inaccurate priors: 
#       Priors deliberately selected to contradict the values used in simulation.
#       Logic: Bias the prior in the opposite direction.
#       - If true parameter < 0.5 -> Prior mean set to 0.7.
#       - If true parameter > 0.5 -> Prior mean set to 0.3.
#
#    C. Uninformative priors: 
#       Using uniform distributions from 0 to 1 (Beta(1,1)).

#Here's an overview of the different sections of the code

# 1. Global Simulation Parameters: Sets up fixed parameters (lambda, po) and 
# loads the specific vertebrate group characteristics.

# 2. Load Monthly Trend: Loads the annual tendency of roadkill abundance.

# 3. Function to calculate beta parameters: Calculates alpha and beta from Mean and SE.

# 4.  Define Factorial Design: Constructs the grid of 9 scenarios combining
# biological variations (+/- 10%) and prior types.

# 5. Main Simulation Loop - Factorial Design: Iterates through the 27 scenarios.
#Simulates census data based on the specific pL of the current scenario,
#and runs the model applying the specific Prior logic (Accurate/Inaccurate/Uninformative)


# -----------------------------------------------------------------------------
# RENV & LIBRARIES
# This project uses the 'renv' package to manage dependencies and ensure 
# reproducibility. To replicate the exact computational environment:
#
# 1. Download and open Supplementary_Material project in RStudio (preferably via
#    the .Rproj file).
# 2. If this is your first time running this code on this machine, uncomment 
#    and run the line below to install the specific package versions required:
#
#    renv::restore()
#
# Note: You do not need to manually run 'install.packages()'.
# -----------------------------------------------------------------------------
library(jagsUI)
library(dplyr)
library(truncnorm)
# -----------------------------------------------------------------------------
# DIRECTORY & PATHS
# -----------------------------------------------------------------------------
# Set up relative paths for a project-based workflow
if(!dir.exists("prior_sensitivity_results")) dir.create("prior_sensitivity_results")
save_path <- "prior_sensitivity_results/" 

# -----------------------------------------------------------------------------
# 1. Global Simulation Parameters & Config
# -----------------------------------------------------------------------------

nroad_transects <- 10
nmonth <- 12   
nmethod <- 3   
nsimulations <- 20  # Replications per scenario

# Standard Error for priors (Only pl applies here)
se_pl <- 0.05

# Mammals G5 Configuration
select_group <- function(group_name) {
  if (group_name == "mammals_g5") {
    anual_tendency <- "Supplementary material S5/anual_tendency_mammals_g5.csv"
    pl_group <- 0.5              # Carcass location probability
    pom_group <- c(1, 0.9, 0.8)  # Carcass observation probabilities
  } else {
    stop("Invalid group. Please choose 'mammals_g5'.")
  }
  return(list(anual_tendency = anual_tendency, pl_group = pl_group, pom_group = pom_group))
}

group <- "mammals_g5"
selected_group <- select_group(group)

# Create save directory if it doesn't exist
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

# Detect cores for parallel processing 
n_cores <- 3
cat("Using", n_cores, "cores for JAGS parallel processing\n")

# -----------------------------------------------------------------------------
# 2. Load Monthly Trend
# -----------------------------------------------------------------------------
if(file.exists(selected_group$anual_tendency)){
  df <- read.csv(selected_group$anual_tendency, sep = ";")
  monthly_trend <- df$Roadkill
} else {
  stop(paste("Annual trend file not found:", selected_group$anual_tendency))
}

# -----------------------------------------------------------------------------
# 3. Function to calculate beta parameters
# -----------------------------------------------------------------------------

get_beta_params <- function(mean, se) {
  se <- max(se, 1e-5)
  var <- se^2
  # Beta distribution constraint: var < mean * (1 - mean)
  if (var >= mean * (1 - mean)) var <- mean * (1 - mean) * 0.99
  
  alpha <- ((1 - mean) / var - 1 / mean) * (mean^2)
  beta  <- alpha * (1 / mean - 1)
  
  list(alpha = max(alpha, 0.01), beta = max(beta, 0.01))
}

# -----------------------------------------------------------------------------
# 4. Define Factorial Design
# -----------------------------------------------------------------------------

# Level 1: pL variations (±10%)
pl_scenarios <- list(
  pl_literature = selected_group$pl_group,
  pl_plus10 = min(0.99, selected_group$pl_group * 1.10),
  pl_minus10 = selected_group$pl_group * 0.90
)

# Level 2: Prior types
prior_types <- c("accurate", "inaccurate", "uninformative")

# Create factorial design
factorial_design <- expand.grid(
  pl_scenario = names(pl_scenarios),
  prior_type = prior_types,
  stringsAsFactors = FALSE
)

cat("\n=== Factorial Design (Mammals G5) ===\n")
cat("Total scenarios:", nrow(factorial_design), "\n")
cat("Simulations per scenario:", nsimulations, "\n")
cat("Total simulations:", nrow(factorial_design) * nsimulations, "\n\n")

# -----------------------------------------------------------------------------
# 5. Main Simulation Loop - Factorial Design
# -----------------------------------------------------------------------------

total_sim_counter <- 0

# Loop through each scenario
for (scenario_id in 1:nrow(factorial_design)) {
  
  # Extract Params for this scenario
  pl_scenario_name <- factorial_design$pl_scenario[scenario_id]
  prior_type <- factorial_design$prior_type[scenario_id]
  
  # Get actual parameter values
  pl_true <- pl_scenarios[[pl_scenario_name]]
  
  cat("\n", rep("=", 80), "\n")
  cat("SCENARIO", scenario_id, "of", nrow(factorial_design), "\n")
  cat(" Params: pL=", round(pl_true, 4), "(", pl_scenario_name, ")\n")
  cat(" Prior:", prior_type, "\n")
  cat(rep("=", 80), "\n\n")
  
  # -------------------------------------
  # Configure Priors
  # -------------------------------------
  if (prior_type == "accurate") {
    pl_prior_params <- get_beta_params(pl_true, se_pl)
  } else if (prior_type == "inaccurate") {
    # Shift mean to opposite side of 0.5
    pl_inaccurate_mean <- ifelse(pl_true < 0.5, 0.7, 0.3)
    pl_prior_params <- get_beta_params(pl_inaccurate_mean, se_pl)
  } else {  # uninformative
    pl_prior_params <- list(alpha = 1, beta = 1)
  }
  
  # Loop through simulations
  for (sim in 1:nsimulations) {
    
    total_sim_counter <- total_sim_counter + 1
    
    # Log
    cat("\n", rep("-", 70), "\n")
    cat("  Simulation", sim, "/", nsimulations, 
        "(Total Progress:", total_sim_counter, "/", nrow(factorial_design) * nsimulations, ")\n")
    cat("  Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    
    tryCatch({
      
      # -------------------------------------
      # Simulate Data 
      # -------------------------------------
      # No daily loop, single Lambda per month
      
      # A. Lambda Simulation
      lambdasim <- monthly_trend
      
      # B. Total N_it (Poisson)
      N_it <- matrix(NA, nrow = nroad_transects, ncol = nmonth)
      for (t in 1:nmonth) {
        if (length(lambdasim) >= t) {
          N_it[, t] <- rpois(nroad_transects, lambda = lambdasim[t]) 
        }
      }
      
      # C. N2_it (Location on Road) - Binomial with pL
      N2_it <- matrix(NA, nrow = nroad_transects, ncol = nmonth)
      for (i in 1:nroad_transects) {
        for (t in 1:nmonth) {
          N2_it[i, t] <- rbinom(1, N_it[i, t], pl_true)
        }
      }
      
      # D. Census Data C (Observation) - Binomial with pO
      # Note: No N3 step, N2 goes directly to Observation
      nsurveys <- 3 
      C <- array(dim = c(nroad_transects, nmonth, nmethod, nsurveys))
      
      for (m in 1:nmethod) {
        for (k in 1:nsurveys) {
          C[,,m,k] <- matrix(rbinom(nroad_transects * nmonth, N2_it, selected_group$pom_group[m]), 
                             nrow = nroad_transects, ncol = nmonth)
        }
      }
      
      # -------------------------------------
      # Save Simulated Data
      # -------------------------------------
      scenario_suffix <- paste0("_plScen_", pl_scenario_name, "_prior_", prior_type)
      
      if (!dir.exists(save_path)) stop(paste("Save path does not exist:", save_path))
      
      saveRDS(lambdasim, file = paste0(save_path, group, scenario_suffix, "_lambdasim_sim_", sim, ".RData"))
      saveRDS(N_it, file = paste0(save_path, group, scenario_suffix, "_N_it_sim_", sim, ".RData"))
      saveRDS(N2_it, file = paste0(save_path, group, scenario_suffix, "_N2_it_sim_", sim, ".RData"))
      saveRDS(C, file = paste0(save_path, group, scenario_suffix, "_census_data_sim_", sim, ".RData"))
      
      # -------------------------------------
      # JAGS Analysis
      # -------------------------------------
      bdata <- list(
        C = C, 
        nroad_transects = nroad_transects, 
        nmonth = nmonth, 
        nmet = nmethod, 
        nsurveys = nsurveys,
        alpha_prior_pl = pl_prior_params$alpha, 
        beta_prior_pl = pl_prior_params$beta
      )
      
      # Unique model filename per Sim
      model_filename <- paste0("model_g5_scen", scenario_id, "_sim", sim, ".txt")
      
      sink(model_filename)
      cat("
      model {
        # Priors
        pl ~ dbeta(alpha_prior_pl, beta_prior_pl)  # Carcass location probability
        
        for (m in 1:nmet) {
          po[m] ~ dunif(0, 1)  # Carcass observation probability
        }
        
        for (t in 1:nmonth) {
          lambda[t] ~ dunif(0, 300)  # Expected abundance
        }
        
        # Likelihood
        for (i in 1:nroad_transects) {
          for (t in 1:nmonth) {
            # 1. Total Number
            N[i, t] ~ dpois(lambda[t])
            
            # 2. Location (N -> N2)
            N2[i,t] ~ dbin(pl, N[i,t])T(0,100000)
            
            # 3. Observation (N2 -> C) - No Persistence Step
            for (m in 1:nmet) {
              for (j in 1:nsurveys) {
                C[i, t, m, j] ~ dbin(po[m], N2[i, t])T(0, 100000)
              }
            }
          }
        }
        
        # Derived quantities
        for (t in 1:nmonth) {
          totalN[t] <- sum(N[, t]) 
          totalN2[t] <- sum(N2[,t])
        }
      }
      ", fill = TRUE)
      sink()
      
      # Initial Values
      Nst <- apply(C, c(1, 2), max, na.rm = TRUE) + 500
      Nst[!is.finite(Nst)] <- 500 
      
      inits <- function() {
        list(
          lambda = runif(nmonth, 5, 100), 
          pl = rbeta(1, pl_prior_params$alpha, pl_prior_params$beta),        
          po = runif(3),                                                              
          N = Nst                                                               
        )
      }
      
      params <- c("lambda", "pl", "po", "totalN", "totalN2")
      
      # MCMC Settings
      na <- 100000; ni <- 400000; nt <- 1000; nb <- 100000
      nc <- min(n_cores, 3)
      
      jags_start_time <- Sys.time()
      cat("  Starting JAGS at:", format(jags_start_time, "%H:%M:%S"), "\n")
      
      output <- jags(
        bdata, inits, params, model_filename,
        n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
        parallel = TRUE, n.cores = nc
      )
      
      jags_end_time <- Sys.time()
      elapsed_time <- difftime(jags_end_time, jags_start_time, units = "mins")
      
      # Save Output
      saveRDS(output, file = paste0(save_path, group, scenario_suffix, "_output_analysis_sim_", sim, ".RData"))
      
      file.remove(model_filename)
      
      cat("  ✓ Completed at:", format(jags_end_time, "%H:%M:%S"), "\n")
      cat("  ⏱ JAGS runtime:", round(elapsed_time, 2), "minutes\n")
      
    }, error = function(e) {
      cat("  ✗ FATAL ERROR at:", format(Sys.time(), "%H:%M:%S"), "\n")
      cat("  >>> MESSAGE:", conditionMessage(e), "\n")
      cat("  >>> SCENARIO:", scenario_id, "SIM:", sim, "\n")
    }) # End tryCatch
    
  } # End sim loop
} # End scenario loop

cat("\n", rep("=", 80), "\n")
cat("ALL SIMULATIONS COMPLETED!\n")
cat("Finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(rep("=", 80), "\n")