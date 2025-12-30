### Paper title: A method to estimate actual infrastructure-induced mortality by integrating sampling biases

### Supplementary Material S9

### R script to simulate and analyse road survey census data for the Prior Sensitivity Analysis for vertebrate
#groups only affected by carcass persistence bias

## R version 4.2.2

# -----------------------------------------------------------------------------
# SIMULATION WORKFLOW DESCRIPTION
# -----------------------------------------------------------------------------
# This script executes the Prior Sensitivity Analysis workflow described in 
# Section 2.4.3 of the main text for vertebrates groups only affected by carcass
# persistence bias. We simulate data by fixing the daily mean number of roadkills
#(lambda) and carcass observation probability per survey method (po), while 
#varying daily carcass persistence probabilities (pPd).

# The design follows a full factorial approach resulting in 9 distinct scenarios.
# We perform 20 simulations per scenario, totaling 180 simulations per vertebrate group.

#    BIOLOGICAL PARAMETERS VARIATION:
#    We simulate scenarios using empirically derived values from literature, 
#    as well as scenarios where pPd is shifted by +/- 10% away from their 
#    values in literature. We selected a +/- 10% deviation because larger deviations 
#    would result in parameter combinations that are rarely observed in nature.
#    - pPd levels: Literature value, Literature + 10%, Literature - 10%.

#    PRIOR SPECIFICATIONS TESTED:
#    For each simulated scenario, we test the model's performance using three 
#    prior specifications for pP (average persistence):
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

# 1. Global Simulation Parameters Setup: Sets up fixed parameters (lambda, po) and 
# loads the specific vertebrate group characteristics.

# 2. Load Monthly Trend: Loads the annual tendency of roadkill abundance.

# 3. Function to calculate beta parameters: Calculates alpha and beta from Mean and SE.

# 4.  Define Factorial Design: Constructs the grid of 9 scenarios combining
# biological variations (+/- 10%) and prior types.

# 5. Main Simulation Loop - Factorial Design: Iterates through the 27 scenarios.
#Simulates census data based on the specific pPd of the current scenario,
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
# 1. Global Simulation Parameters Setup
# -----------------------------------------------------------------------------

nroad_transects <- 10
nmonth <- 12  
nmethod <- 3  
nsimulations <- 20  # Replications per scenario

# Standard Error for priors (Only pp, as pl is not applicable)
se_pp <- 0.05

# Vertebrate group selection (Groups NOT affected by carcass location bias)
# Values taken from S3 Supplementary Material
select_group <- function(group_name) {
  if (group_name == "amphibians") {
    anual_tendency <- c("Supplementary material S5/anual_tendency_amphibians.csv")
    ppd_group <- 0.392598565
    pom_group <- c(0.5, 0.3, 0.02)
  } else if (group_name == "reptiles_g1") {
    anual_tendency <- c("Supplementary material S5/anual_tendency_reptiles_g1.csv")
    ppd_group <- 0.05655527
    pom_group <- c(0.5, 0.3, 0.02)
  } else if (group_name == "mammals_g1") {
    anual_tendency <- c("Supplementary material S5/anual_tendency_mammals_g1.csv")
    ppd_group <- 0.388601036
    pom_group <- c(0.6, 0.4, 0.05)
  } else if (group_name == "mammals_g2") {
    anual_tendency <- c("Supplementary material S5/anual_tendency_mammals_g2.csv")
    ppd_group <- 0.506361323
    pom_group <- c(0.8, 0.6, 0.2)
  } else if (group_name == "mammals_g3") {
    anual_tendency <- c("Supplementary material S5/anual_tendency_mammals_g3.csv")
    ppd_group <- 0.78250591
    pom_group <- c(0.8, 0.6, 0.2)
  } else {
    stop("Invalid group. For NO-pL models, choose: 'amphibians', 'reptiles_g1', 'mammals_g1', 'mammals_g2', or 'mammals_g3'.")
  }
  
  # Removed pl_group from return
  return(list(anual_tendency = anual_tendency, ppd_group = ppd_group, pom_group = pom_group))
}

group <- "amphibians" # Change this to your desired group
selected_group <- select_group(group)

# Create save directory if it doesn't exist
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

# Detect cores for parallel processing
n_cores <- 3
cat("Using", n_cores, "cores for parallel processing\n")

# -----------------------------------------------------------------------------
# 2. Load Monthly Trend
# -----------------------------------------------------------------------------

# Ensure the csv file exists in your WD
df <- read.csv(selected_group$anual_tendency, sep = ";")
month_order <- c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", 
                 "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")
df$Month <- factor(df$Month, levels = month_order)
monthly_trend <- df$Roadkill

# -----------------------------------------------------------------------------
# 3. Function to calculate beta parameters
# -----------------------------------------------------------------------------

get_beta_params <- function(mean, se) {
  se <- max(se, 1e-5)
  alpha <- ((1 - mean) / (se^2) - 1 / mean) * (mean^2)
  beta  <- alpha * (1 / mean - 1)
  alpha <- max(alpha, 1e-5)
  beta  <- max(beta, 1e-5)
  list(alpha = alpha, beta = beta)
}

# -----------------------------------------------------------------------------
# 4. Define Factorial Design
# -----------------------------------------------------------------------------

# Level 1: pPd variations (±10%)
ppd_scenarios <- list(
  ppd_literature = selected_group$ppd_group,
  ppd_plus10 = min(0.99, selected_group$ppd_group * 1.10), # Cap at 0.99
  ppd_minus10 = selected_group$ppd_group * 0.90
)

# Level 2: Prior types
prior_types <- c("accurate", "inaccurate", "uninformative")

# Create factorial design
factorial_design <- expand.grid(
  ppd_scenario = names(ppd_scenarios),
  prior_type = prior_types,
  stringsAsFactors = FALSE
)

# Add a dummy pL column just for naming consistency if needed, or omit it.
# Here we stick to relevant parameters.

cat("\n=== Factorial Design (No pL Model) ===\n")
cat("Total scenarios:", nrow(factorial_design), "\n")
cat("Simulations per scenario:", nsimulations, "\n")
cat("Total simulations:", nrow(factorial_design) * nsimulations, "\n\n")

print(factorial_design)

# -----------------------------------------------------------------------------
# 5. Main Simulation Loop - Factorial Design
# -----------------------------------------------------------------------------

# Counter for total simulations
total_sim_counter <- 0

# Loop through each scenario in factorial design
for (scenario_id in 1:nrow(factorial_design)) {
  
  # Extract scenario parameters
  ppd_scenario_name <- factorial_design$ppd_scenario[scenario_id]
  prior_type_raw <- factorial_design$prior_type[scenario_id]
  
  # Standardize prior type naming (map 'accurate' to 'accurate' if necessary)
  if(prior_type_raw == "accurate") prior_type <- "accurate"
  else if(prior_type_raw == "inaccurate") prior_type <- "inaccurate"
  else prior_type <- "uninformative"
  
  # Get actual parameter values for this scenario
  ppd_true <- ppd_scenarios[[ppd_scenario_name]]
  
  cat("\n", rep("=", 80), "\n")
  cat("SCENARIO", scenario_id, "of", nrow(factorial_design), "\n")
  cat("pPd:", ppd_scenario_name, "(", round(ppd_true, 4), ")\n")
  cat("Prior:", prior_type, "\n")
  cat(rep("=", 80), "\n\n")
  
  # -------------------------------------
  # Calculate D and pp_true
  # -------------------------------------
  num_days <- 30
  x_values <- 0:(num_days - 1)
  y_values <- numeric(num_days)
  y_values[1] <- 1
  
  for (day in 2:num_days) {
    y_values[day] <- y_values[day - 1] * ppd_true
  }
  
  target_probability <- 0.05
  D <- round(log(target_probability) / log(ppd_true), 0)
  
  # Safety check for D
  if(is.infinite(D) || is.na(D) || D < 1) D <- 1
  
  total_area <- D * 0.95
  # Ensure we don't index out of bounds if D is large
  len_vals <- min(length(x_values), D+1)
  
  curve_function <- approxfun(x_values[1:len_vals], y_values[1:len_vals])
  curve_integral <- integrate(curve_function, lower = 0, upper = D)$value
  pp_true <- curve_integral / total_area
  
  cat("Calculated D:", D, "days\n")
  cat("Calculated pp (true):", round(pp_true, 4), "\n\n")
  
  # -------------------------------------
  # Configure Priors
  # -------------------------------------
  if (prior_type == "accurate") {
    # accurate informative: centered on true values
    pp_prior_params <- get_beta_params(pp_true, se_pp)
    
  } else if (prior_type == "inaccurate") {
    # inaccurate informative: opposite tail logic
    pp_inaccurate_mean <- ifelse(pp_true < 0.5, 0.7, 0.3)
    pp_prior_params <- get_beta_params(pp_inaccurate_mean, se_pp)
    
  } else {  # uninformative
    pp_prior_params <- list(alpha = 1, beta = 1)
  }
  
  # -----------------------------------------------------------------------------
  # Run simulations for this scenario
  # -----------------------------------------------------------------------------
  
  for (sim in 1:nsimulations) {
    
    total_sim_counter <- total_sim_counter + 1
    
    # Display timestamp
    cat("\n", rep("-", 70), "\n")
    cat("  Simulation", sim, "/", nsimulations, 
        "(Total:", total_sim_counter, "/", nrow(factorial_design) * nsimulations, ")\n")
    cat("  Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    
    # Define SDs (fixed at 0 for this sensitivity analysis)
    sd_random_values <- 0
    sd_ppd <- 0
    
    # Initialize Storage Lists
    lambdas_storage <- list()
    N_itd_storage <- list()
    N3_itd_storage <- list()
    
    # Generate random values
    random_values <- rtruncnorm(n = D, a = 0, mean = 1, sd = sd_random_values)
    
    # Generate lambdas
    for (j in 1:D) {
      lambdasim <- monthly_trend * random_values[j]
      assign(paste0("lambdasim", j), lambdasim) # Backward compatibility
      lambdas_storage[[j]] <- lambdasim         # Store in list for saving
    }
    
    # Generate N_itd matrices (Poisson)
    for (j in 1:D) {
      lambdasim_current <- lambdas_storage[[j]]
      N_itd <- matrix(NA, nrow = nroad_transects, ncol = nmonth)
      
      for (t in 1:nmonth) {
        if (length(lambdasim_current) >= t) {
          N_itd[, t] <- rpois(nroad_transects, lambda = lambdasim_current[t])
        }
      }
      assign(paste0("N_itd", j), N_itd)
      N_itd_storage[[j]] <- N_itd # Store in list
    }
    
    # Generate ppd values
    ppd_values <- rnorm(D, mean = ppd_true, sd = sd_ppd)
    
    # Initialize N3_itD (Accumulator)
    N3_itD <- matrix(0, nrow = nroad_transects, ncol = nmonth)
    
    # Generate N3_itd matrices
    for (j in 1:D) {
      ppd_t <- if (j == 1) prod(ppd_values) else prod(ppd_values[j:length(ppd_values)])
      
      # Probability bounds check
      if(ppd_t > 1) ppd_t <- 1; if(ppd_t < 0) ppd_t <- 0
      
      N3_itd <- matrix(NA, nrow = nroad_transects, ncol = nmonth)
      
      # NOTE: Using N_itd directly (instead of N2_itd because No Location Bias model)
      N_itd_current <- N_itd_storage[[j]]
      
      for (i in 1:nroad_transects) {
        for (t in 1:nmonth) {
          # Binomial: Survival from Total N directly to Persistence N3
          N3_itd[i, t] <- rbinom(1, N_itd_current[i, t], ppd_t)
        }
      }
      assign(paste0("N3_itd", j), N3_itd)
      N3_itd_storage[[j]] <- N3_itd # Store in list
      
      N3_itD <- N3_itD + N3_itd
    }
    
    # Census data simulation
    nsurveys <- 3
    C <- array(dim = c(nroad_transects, nmonth, nmethod, nsurveys))
    
    for (m in 1:3) {
      for (k in 1:nsurveys) {
        C[,,m,k] <- matrix(rbinom(nroad_transects * nmonth, N3_itD, selected_group$pom_group[m]), 
                           nrow = nroad_transects, ncol = nmonth)
      }
    }
    
    # -------------------------------------
    # Save Simulated Data & JAGS
    # -------------------------------------
    
    # Start tryCatch block
    tryCatch({
      
      scenario_suffix <- paste0("_ppdScen_", ppd_scenario_name, "_prior_", prior_type)
      
      if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
      
      saveRDS(random_values, file = paste0(save_path, group, scenario_suffix, "_random_values_sim_", sim, ".RData"))
      saveRDS(C, file = paste0(save_path, group, scenario_suffix, "_census_data_sim_", sim, ".RData"))
      
      # Save Intermediate Matrices (lambda, N_itd, N3_itd) using the lists
      for (j in 1:D) {
        saveRDS(lambdas_storage[[j]], 
                file = paste0(save_path, group, scenario_suffix, "_lambda", j, "_sim_", sim, ".RData"))
        
        saveRDS(N_itd_storage[[j]], 
                file = paste0(save_path, group, scenario_suffix, "_N_itd", j, "_sim_", sim, ".RData"))
        
        saveRDS(N3_itd_storage[[j]], 
                file = paste0(save_path, group, scenario_suffix, "_N3_itd", j, "_sim_", sim, ".RData"))
      }
      
      # Prepare data for JAGS model
      bdata <- list(
        C = C,
        nroad_transects = nroad_transects,
        nmonth = nmonth,
        nmet = nmethod,
        nsurveys = nsurveys,
        alpha_prior_pp = pp_prior_params$alpha,
        beta_prior_pp = pp_prior_params$beta
      )
      
      # Write JAGS model
      model_filename <- paste0("model_scenario", scenario_id, "_sim", sim, ".txt")
      
      sink(model_filename)
      cat("
      model {
        # Priors
        pp ~ dbeta(alpha_prior_pp, beta_prior_pp) # Carcass persistence probability
        # Note: No pL prior here
        
        for (m in 1:nmet) {
          po[m] ~ dunif(0, 1)  # Carcass observation probability
        }
      
        for (t in 1:nmonth) {
          lambda[t] ~ dunif(0, 300) # Expected abundance
        }
      
        # Likelihood
        for (i in 1:nroad_transects) {
          for (t in 1:nmonth) {
            # 1. Total Number
            N[i, t] ~ dpois(lambda[t])
            
            # 2. Persistence (Directly from N) - NO N2 step
            N3[i, t] ~ dbin(pp, N[i, t])T(0, 100000)
            
            # 3. Observation
            for (m in 1:nmet) {
              for (j in 1:nsurveys) {
                C[i, t, m, j] ~ dbin(po[m], N3[i, t])T(0, 100000)
              }
            }
          }
        }
      
        # Derived quantities
        for (t in 1:nmonth) {
          totalN[t] <- sum(N[, t]) 
          # No totalN2
          totalN3[t] <- sum(N3[, t])    
        }
      }
      ", fill = TRUE)
      sink()
      
      # Initialize state variables
      Nst <- apply(C, c(1, 2), max, na.rm = TRUE) + 500
      Nst[!is.finite(Nst)] <- 500
      
      # Initial values (Removed pl)
      inits <- function() {
        list(
          lambda = runif(nmonth, 5, 100),
          pp = rbeta(1, pp_prior_params$alpha, pp_prior_params$beta),
          po = runif(3),
          N = Nst
        )
      }
      
      # Parameters to monitor (Removed pl and totalN2)
      params <- c("lambda", "pp", "po", "totalN", "totalN3")
      
      # MCMC settings
      na <- 100000
      ni <- 400000
      nt <- 1000
      nb <- 100000
      nc <- min(n_cores, 3)
      
      # Record start time for JAGS
      jags_start_time <- Sys.time()
      cat("  Starting JAGS at:", format(jags_start_time, "%H:%M:%S"), "\n")
      
      # Run JAGS model
      output <- jags(
        bdata, inits, params, 
        model_filename,
        n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
        parallel = TRUE,
        n.cores = n_cores
      )
      
      # Calculate elapsed time
      jags_end_time <- Sys.time()
      elapsed_time <- difftime(jags_end_time, jags_start_time, units = "mins")
      
      # Save output
      saveRDS(output, file = paste0(save_path, group, scenario_suffix,
                                    "_output_analysis_sim_", sim, ".RData"))
      
      # Clean up model file
      file.remove(model_filename)
      
      cat("  ✓ Completed at:", format(jags_end_time, "%H:%M:%S"), "\n")
      cat("  ⏱ JAGS runtime:", round(elapsed_time, 2), "minutes\n")
      
    }, error = function(e) {
      cat("  ✗ Error at:", format(Sys.time(), "%H:%M:%S"), "\n")
      warning(paste("Error in scenario", scenario_id, "sim", sim, ":", e$message))
    }) # End of tryCatch
    
  } # End of sim loop
} # End of scenario loop

cat("\n", rep("=", 80), "\n")
cat("ALL NO-pL SIMULATIONS COMPLETED!\n")
cat("Total scenarios:", nrow(factorial_design), "\n")
cat("Total simulations:", total_sim_counter, "\n")
cat("Finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(rep("=", 80), "\n")