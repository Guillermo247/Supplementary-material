### Paper title: A method to estimate actual infrastructure-induced mortality by integrating sampling biases

### Supplementary Material S8

### R script to simulate and analyse road survey census data for the Prior Sensitivity Analysis for vertebrate
#groups  affected by carcass location and persistence bias

## R version 4.2.2

# -----------------------------------------------------------------------------
# SIMULATION WORKFLOW DESCRIPTION
# -----------------------------------------------------------------------------
# This script executes the Prior Sensitivity Analysis workflow described in 
# Section 2.4.3 of the main text for vertebrates groups affected by carcass
# location and persistence bias. We simulate data by fixing the daily mean number 
# of roadkills (lambda) and carcass observation probability per survey method (po), 
# while varying carcass location (pL) and daily carcass persistence (pPd) probabilities.

# The design follows a full factorial approach resulting in 27 distinct scenarios.
# We perform 20 simulations per scenario, totaling 540 simulations per vertebrate group.

#    BIOLOGICAL PARAMETERS VARIATION:
#    We simulate scenarios using empirically derived values from literature, 
#    as well as scenarios where pL and pPd are shifted by +/- 10% away from their 
#    values in literature. We selected a +/- 10% deviation because larger deviations 
#    would result in parameter combinations that are rarely observed in nature.
#    - pL levels:  Literature value, Literature + 10%, Literature - 10%.
#    - pPd levels: Literature value, Literature + 10%, Literature - 10%.

#    PRIOR SPECIFICATIONS TESTED:
#    For each simulated scenario, we test the model's performance using three 
#    prior specifications for pL and pP (average persistence):
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

# 4.  Define Factorial Design: Constructs the grid of 27 scenarios combining
# biological variations (+/- 10%) and prior types.

# 5. Main Simulation Loop - Factorial Design: Iterates through the 27 scenarios.
#Simulates census data based on the specific pL and pPd of the current scenario,
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

# Standard Error for priors
se_pl <- 0.05
se_pp <- 0.05

# Vertebrate group selection
select_group <- function(group_name) {
  if (group_name == "reptiles_g2") {
    anual_tendency <- c("Supplementary material S5/anual_tendency_reptiles_g2.csv")
    pl_group <- 0.42857
    ppd_group <- 0.398477157
    pom_group <- c(0.7, 0.5, 0.1)
  } else if (group_name == "birds_bats_g1") {
    anual_tendency <- c("Supplementary material S5/anual_tendency_birds_bats_g1.csv")
    pl_group <- 0.506329114
    ppd_group <- 0.358632169
    pom_group <- c(0.6, 0.4, 0.05)
  } else if (group_name == "birds_g2") {
    anual_tendency <- c("Supplementary material S5/anual_tendency_birds_g2.csv")
    pl_group <- 0.692307692
    ppd_group <- 0.747433082
    pom_group <- c(0.8, 0.6, 0.2)
  } else if (group_name == "mammals_g4") {
    anual_tendency <- c("Supplementary material S5/anual_tendency_mammals_g4.csv")
    pl_group <- 0.647058824
    ppd_group <- 0.805626598
    pom_group <- c(0.9, 0.7, 0.3)
  } else {
    stop("Invalid group. Please choose between 'reptiles_g2', 'birds_bats_g1', 'birds_g2' or 'mammals_g4'.")
  }
  
  return(list(anual_tendency = anual_tendency, pl_group = pl_group, ppd_group = ppd_group, pom_group = pom_group))
}

group <- "birds_bats_g1"
selected_group <- select_group(group)

# Create save directory if it doesn't exist
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

# Setting 3 cores for 3 MCMC chain
n_cores <- 3
cat("Using", n_cores, "cores for parallel processing\n")

# -----------------------------------------------------------------------------
# 2. Load Monthly Trend
# -----------------------------------------------------------------------------

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

# Level 1: pL variations (±10%)
pl_scenarios <- list(
  pl_literature = selected_group$pl_group,
  pl_plus10 = selected_group$pl_group * 1.10,
  pl_minus10 = selected_group$pl_group * 0.90
)

# Level 2: pPd variations (±10%)
ppd_scenarios <- list(
  ppd_literature = selected_group$ppd_group,
  ppd_plus10 = selected_group$ppd_group * 1.10,
  ppd_minus10 = selected_group$ppd_group * 0.90
)

# Level 3: Prior types
prior_types <- c("accurate", "inaccurate", "uninformative")

# Create factorial design
factorial_design <- expand.grid(
  pl_scenario = names(pl_scenarios),
  ppd_scenario = names(ppd_scenarios),
  prior_type = prior_types,
  stringsAsFactors = FALSE
)

cat("\n=== Factorial Design ===\n")
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
  pl_scenario_name <- factorial_design$pl_scenario[scenario_id]
  ppd_scenario_name <- factorial_design$ppd_scenario[scenario_id]
  prior_type <- factorial_design$prior_type[scenario_id]
  
  # Get actual parameter values for this scenario
  pl_true <- pl_scenarios[[pl_scenario_name]]
  ppd_true <- ppd_scenarios[[ppd_scenario_name]]
  
  cat("\n", rep("=", 80), "\n")
  cat("SCENARIO", scenario_id, "of", nrow(factorial_design), "\n")
  cat("pL:", pl_scenario_name, "(", round(pl_true, 4), ")\n")
  cat("pPd:", ppd_scenario_name, "(", round(ppd_true, 4), ")\n")
  cat("Prior:", prior_type, "\n")
  cat(rep("=", 80), "\n\n")
  
  # Calculate D for this ppd scenario
  num_days <- 30
  x_values <- 0:(num_days - 1)
  y_values <- numeric(num_days)
  y_values[1] <- 1
  
  for (day in 2:num_days) {
    y_values[day] <- y_values[day - 1] * ppd_true
  }
  
  target_probability <- 0.05
  D <- round(log(target_probability) / log(ppd_true), 0)
  
  total_area <- D * 0.95
  curve_function <- approxfun(x_values[1:min(length(x_values), D+1)], 
                              y_values[1:min(length(y_values), D+1)])
  curve_integral <- integrate(curve_function, lower = 0, upper = D)$value
  pp_true <- curve_integral / total_area
  
  cat("Calculated D:", D, "days\n")
  cat("Calculated pp:", round(pp_true, 4), "\n\n")
  
  # Configure priors for this scenario
  if (prior_type == "accurate") {
    # accurate informative: centered on true values
    pl_prior_params <- get_beta_params(pl_true, se_pl)
    pp_prior_params <- get_beta_params(pp_true, se_pp)
    
  } else if (prior_type == "inaccurate") {
    # inaccurate informative: opposite tail
    pl_inaccurate_mean <- ifelse(pl_true < 0.5, 0.7, 0.3)
    pp_inaccurate_mean <- ifelse(pp_true < 0.5, 0.7, 0.3)
    pl_prior_params <- get_beta_params(pl_inaccurate_mean, se_pl)
    pp_prior_params <- get_beta_params(pp_inaccurate_mean, se_pp)
    
  } else {  # uninformative
    pl_prior_params <- list(alpha = 1, beta = 1)
    pp_prior_params <- list(alpha = 1, beta = 1)
  }
  
  # -----------------------------------------------------------------------------
  # Data generation and model run for this scenario
  # -----------------------------------------------------------------------------
  
  for (sim in 1:nsimulations) {
    
    total_sim_counter <- total_sim_counter + 1
    
    # Display timestamp at the start of each simulation
    cat("\n", rep("-", 70), "\n")
    cat("  Simulation", sim, "/", nsimulations, 
        "(Total:", total_sim_counter, "/", nrow(factorial_design) * nsimulations, ")\n")
    cat("  Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat(rep("-", 70), "\n")
    
    # Define SDs (fixed at 0 for this analysis)
    sd_random_values <- 0
    sd_ppd <- 0
    
    # Initialize Storage Lists
    lambdas_storage <- list()
    N_itd_storage   <- list()
    N2_itd_storage  <- list()
    N3_itd_storage  <- list()
    
    # Generate random values
    random_values <- rtruncnorm(n = D, a = 0, mean = 1, sd = sd_random_values)
    
    # Generate lambdas
    for (j in 1:D) {
      lambdasim <- monthly_trend * random_values[j]
      assign(paste0("lambdasim", j), lambdasim)
      lambdas_storage[[j]] <- lambdasim # Store in list
    }
    
    # Generate N_itd matrices
    for (j in 1:D) {
      lambdasim_current <- get(paste0("lambdasim", j))
      N_itd <- matrix(NA, nrow = nroad_transects, ncol = nmonth)
      
      for (t in 1:nmonth) {
        if (length(lambdasim_current) >= t) {
          N_itd[, t] <- rpois(nroad_transects, lambda = lambdasim_current[t])
        }
      }
      assign(paste0("N_itd", j), N_itd)
      N_itd_storage[[j]] <- N_itd # Store in list
    }
    
    # Generate N2_itd matrices (using true pL value)
    for (j in 1:D) {
      N_itd_current <- get(paste0("N_itd", j))
      N2_itd <- matrix(NA, nrow = nroad_transects, ncol = nmonth)
      
      for (i in 1:nroad_transects) {
        for (t in 1:nmonth) {
          N2_itd[i, t] <- rbinom(1, N_itd_current[i, t], pl_true)
        }
      }
      assign(paste0("N2_itd", j), N2_itd)
      N2_itd_storage[[j]] <- N2_itd # Store in list
    }
    
    # Generate ppd values (using true pPd value)
    ppd_values <- rnorm(D, mean = ppd_true, sd = sd_ppd)
    
    # Initialize N3_itD
    N3_itD <- matrix(0, nrow = nroad_transects, ncol = nmonth)
    
    # Generate N3_itd matrices
    for (j in 1:D) {
      ppd_t <- if (j == 1) prod(ppd_values) else prod(ppd_values[j:length(ppd_values)])
      N3_itd <- matrix(NA, nrow = nroad_transects, ncol = nmonth)
      N2_itd_current <- get(paste0("N2_itd", j))
      
      for (i in 1:nroad_transects) {
        for (t in 1:nmonth) {
          N3_itd[i, t] <- rbinom(1, N2_itd_current[i, t], ppd_t)
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
    
    # Create scenario suffix for file naming
    scenario_suffix <- paste0("_plScen_", pl_scenario_name, 
                              "_ppdScen_", ppd_scenario_name,
                              "_prior_", prior_type)
    
    # Save data files
    tryCatch({
      
      # Save simulated data
      saveRDS(random_values, file = paste0(save_path, group, scenario_suffix,
                                           "_random_values_sim_", sim, ".RData"))
      
      # Save Intermediate Matrices
      for (j in 1:D) {
        saveRDS(lambdas_storage[[j]], file = paste0(save_path, group, scenario_suffix,
                                                    "_lambda", j, "_sim_", sim, ".RData"))
        saveRDS(N_itd_storage[[j]], file = paste0(save_path, group, scenario_suffix,
                                                  "_N_itd", j, "_sim_", sim, ".RData"))
        saveRDS(N2_itd_storage[[j]], file = paste0(save_path, group, scenario_suffix,
                                                   "_N2_itd", j, "_sim_", sim, ".RData"))
        saveRDS(N3_itd_storage[[j]], file = paste0(save_path, group, scenario_suffix,
                                                   "_N3_itd", j, "_sim_", sim, ".RData"))
      }
      
      saveRDS(C, file = paste0(save_path, group, scenario_suffix,
                               "_census_data_sim_", sim, ".RData"))
      
      # Prepare data for JAGS model
      bdata <- list(
        C = C,
        nroad_transects = nroad_transects,
        nmonth = nmonth,
        nmet = nmethod,
        nsurveys = nsurveys,
        alpha_prior_pl = pl_prior_params$alpha,
        beta_prior_pl = pl_prior_params$beta,
        alpha_prior_pp = pp_prior_params$alpha,
        beta_prior_pp = pp_prior_params$beta
      )
      
      # Write JAGS model
      model_filename <- paste0("model_scenario", scenario_id, "_sim", sim, ".txt")
      
      sink(model_filename)
      cat("
      model {
        # Priors
        pl ~ dbeta(alpha_prior_pl, beta_prior_pl)
        pp ~ dbeta(alpha_prior_pp, beta_prior_pp)
        
        for (m in 1:nmet) {
          po[m] ~ dunif(0, 1)
        }
      
        for (t in 1:nmonth) {
          lambda[t] ~ dunif(0, 300)
        }
       
        # Likelihood
        for (i in 1:nroad_transects) {
          for (t in 1:nmonth) {
            N[i, t] ~ dpois(lambda[t])
            N2[i,t] ~ dbin(pl, N[i,t])T(0,100000)
            N3[i, t] ~ dbin(pp, N2[i, t])T(0, 100000)
            
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
          totalN2[t] <- sum(N2[,t])
          totalN3[t] <- sum(N3[, t])    
        }
      }
      ", fill = TRUE)
      sink()
      
      # Initialize state variables
      Nst <- apply(C, c(1, 2), max, na.rm = TRUE) + 500
      Nst[!is.finite(Nst)] <- 500
      
      # Initial values
      inits <- function() {
        list(
          lambda = runif(nmonth, 5, 100),
          pl = rbeta(1, pl_prior_params$alpha, pl_prior_params$beta),
          pp = rbeta(1, pp_prior_params$alpha, pp_prior_params$beta),
          po = runif(3),
          N = Nst
        )
      }
      
      # Parameters to monitor
      params <- c("lambda", "pl", "pp", "po", "totalN", "totalN2", "totalN3")
      
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
    }) # Closing tryCatch properly
  }
}

cat("\n", rep("=", 80), "\n")
cat("ALL SIMULATIONS COMPLETED!\n")
cat("Total scenarios:", nrow(factorial_design), "\n")
cat("Total simulations:", total_sim_counter, "\n")
cat("Finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(rep("=", 80), "\n")

