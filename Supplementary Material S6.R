### Paper title: A novel method to estimate actual infrastructure-induced mortality by integrating sampling biases

### Supplementary Material S6

### R script to simulate and analyse road survey census data for amphibians and reptiles g1 groups only accounting for peak abundance months.

## R version 4.2.2

#In this script, we simulate census data for amphibians and reptiles g1 groups 
#only accounting for peak abundance months under various scenarios involving the
#number of road transects, as well as the variability in both the daily number 
#of roadkills (lambda_itd) and the daily carcass persistence probability (ppd), 
#as described in Section 2.4 of the main text of the paper.
#Using the modeling framework outlined in Section 2.2, this script estimate 
#the simulated values for the total number of roadkills (referred to as N_it in 
#the simulation section or totalN in the model output), as well as the carcass 
#persistence (pp), and carcass observation per survey method (pom) probabilities.
#We use in the analysis standard errors for thea pp prior set at either 0.05 or 0.1.

#In this model, carcass location bias is not accounted. As a result, Equation 2
#(see Section 2.2 of the main text) and the associated carcass location probability
#are excluded.

#Here's an overview of the different sections of the code

#1. Global Simulation Parameters Setup: This part of the code let you customise 
#each parameter in the simulation of data census and their posterior analysis

#2. Load Monthly Trend: This part of the code load the annual tendency of 
#roadkill abundance along a year

#3. Filter Only the Annual Abundance Peak: This part of the code filter only the
#peak abundance months from the monthly trend

#4. Generate Carcass Persistence Curve for pp Prior Based on Santos, Carvalho, and Mira (2011)
#Daily pp (ppd): This part of the code calculates the average probability of 
#carcass persistence probability (pp) over the days it remains on the road without
#decomposing (D-day period).

#5. Calculate alpha and beta Parameters of the Beta Distributed Prior for pp: This
#section makes a beta distributed prior from a mean value and standar error of pp

#6. Census data simulation and Analysis: This section simulates census data 
#considering the variability in the daily number of roadkills (lambda_td) and 
#daily carcass persistence (ppd). Additionally, this section applies the model 
#to estimate the total number of roadkills and pom (po1 =walking, po2 = cycling 
#and po3 = driving), while incorporating information of pl and pp as prior distributions


#7. Execution Loop: This section iterates through the defined scenarios:
# - Census data simulation on SD lambda_td 0.5 and SD ppd 0.15
# - Census data simulation on SD lambda_td 0.5 and SD ppd 0
# - Census data simulation on SD lambda_td 1.5 and SD ppd 0.05
# - Census data simulation on SD lambda_td 1.5 and SD ppd 0.15
# - Census data simulation on SD lambda_td 1.5 and SD ppd 0
# - Census data simulation on SD lambda_td 0 and SD ppd 0.05
# - Census data simulation on SD lambda_td 0 and SD ppd 0.15
# - Census data simulation on SD lambda_td 0 and SD ppd 0

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
if(!dir.exists("simulation_results")) dir.create("simulation_results")
save_path <- "simulation_results/" 

# -----------------------------------------------------------------------------
# 1. Global Simulation Parameters Setup
# -----------------------------------------------------------------------------

nroad_transects <- 10   # Scenarios: 10 or 100
nmonth <- 12  
nmethod <- 3  
nsimulations <- 20      # Number of simulations per scenario
nsurveys <- 3

# Standard Error for Carcass Persistence Probability (pp)
se_pp <- 0.05  # Scenarios: 0.05 or 0.10

# Animal Group Select Fuction

select_group <- function(group_name) {
  if (group_name == "amphibians") {
    anual_tendency <- "Supplementary material S5/anual_tendency_amphibians.csv"
    ppd_group <- 0.392598565 
    pom_group <- c(0.5, 0.3, 0.02) 
  } else if (group_name == "reptiles_g1") {
    anual_tendency <- "Supplementary material S5/anual_tendency_reptiles_g1.csv"
    ppd_group <- 0.05655527
    pom_group <- c(0.5, 0.3, 0.02)
  } else {
    stop("Invalid group. Please choose between 'amphibians' or 'reptiles_g1'.")
  }
  
  return(list(anual_tendency= anual_tendency, ppd_group = ppd_group, pom_group = pom_group))
}

# Select vertebrate group 
group <- "amphibians"  # write the group of interest
selected_group <- select_group(group)

# Define Lambda upper limit based on the group
if (group == "amphibians") {
  lambda_upper_limit <- 600
} else {
  lambda_upper_limit <- 300 # Default value for the rest of vertebrate groups
}

# -----------------------------------------------------------------------------
# 2. Load Monthly Trend
# -----------------------------------------------------------------------------

# Load Annual Trend
# Ensure the file path exists
if(file.exists(selected_group$anual_tendency)){
  df <- read.csv(selected_group$anual_tendency, sep = ";")
  monthly_trend <- df$Roadkill
} else {
  stop(paste("Annual trend file not found:", selected_group$anual_tendency))
}

# -----------------------------------------------------------------------------
# 3. Filter Only the Annual Abundance Peak
# -----------------------------------------------------------------------------

# Extract annual abundance peak data for further analysis (values greater than 5)
monthly_trend <- df$Roadkill
monthly_trend <- monthly_trend[monthly_trend > 5]
nmonth <- length(monthly_trend)

# -----------------------------------------------------------------------------
# 4. Generate Carcass Persistence Curve for pp Prior
# -----------------------------------------------------------------------------

# Carcass Persistence (pp) Curve & D calculation
num_days <- 30  
x_values <- 0:(num_days - 1)
y_values <- numeric(num_days)
y_values[1] <- 1

for (day in 2:num_days) {
  y_values[day] <- y_values[day - 1] * selected_group$ppd_group
}

target_probability <- 0.05
D <- round(log(target_probability) / log(selected_group$ppd_group), 0)

# Calculate area under curve for pp prior
total_area <- D * 0.95
curve_function <- approxfun(x_values, y_values)
curve_integral <- integrate(curve_function, lower = 0, upper = D)$value
pp <- curve_integral / total_area

# -----------------------------------------------------------------------------
# 5. Calculate alpha and beta Parameters of the Beta Distributed Prior for pp
# -----------------------------------------------------------------------------

mean_prior_pp <- pp
se_prior_pp <- se_pp
alpha_prior_pp <- ((1 - mean_prior_pp) / (se_prior_pp^2) - 1 / mean_prior_pp) * (mean_prior_pp^2)
beta_prior_pp <- alpha_prior_pp * (1 / mean_prior_pp - 1)

cat("Prior Setup Complete. D =", D, "days. No Location Bias (pl) in this model.\n")


# -----------------------------------------------------------------------------
# 6. Census data simulation and Analysis 
# -----------------------------------------------------------------------------

run_scenario_simulation <- function(sd_lambda, sd_ppd, scenario_id) {
  
  message(paste0("\n--- Starting Scenario ", scenario_id, ": SD_Lambda=", sd_lambda, ", SD_PPD=", sd_ppd, " ---"))
  

  # --------------------------------------------------------------------------
  # SIMULATION LOOP (Generation + Analysis)
  # --------------------------------------------------------------------------
  
  for (sim in 1:nsimulations) {
    
    message(paste0("   > Running Simulation ", sim, " of ", nsimulations, "..."))
    
    # --- PART A: DATA GENERATION ---
    
    # 1. Generate random values for Lambda variability
    if (sd_lambda > 0) {
      random_values <- rtruncnorm(n = D, a = 0, mean = 1, sd = sd_lambda)
    } else {
      random_values <- rep(1, D)
    }
    
    # Save random values
    saveRDS(random_values, file = paste0(save_path, group, "_SE_pp_", se_pp, 
                                         "_nroad_transects_", nroad_transects, 
                                         "_SD_N_", sd_lambda, "_SD_ppd_", sd_ppd, 
                                         "_rnd_vals_sim_", sim, ".RData"))
    
    # 2. Initialize Accumulator
    N3_itD <- matrix(0, nrow = nroad_transects, ncol = nmonth)
    
    # 3. Generate Daily Persistence (ppd)
    if (sd_ppd > 0) {
      ppd <- rtruncnorm(n = D, a = 0, b = 1, mean = selected_group$ppd_group, sd = sd_ppd)
    } else {
      ppd <- rep(selected_group$ppd_group, D)
    }
    
  
    saveRDS(ppd, file = paste0(save_path, group, "_SE_pp_", se_pp, 
                               "_nroad_transects_", nroad_transects, 
                               "_SD_N_", sd_lambda, "_SD_ppd_", sd_ppd, 
                               "_ppd_sim_", sim, ".RData"))
    
    # 4. Loop over Days (D) to calculate survival
    for (j in 1:D) {
      
      # Calculate Lambda for this day layer
      lambdasim <- monthly_trend * random_values[j]
      
      # --- N_itd: Poisson (Total Roadkill) ---
      N_itd <- matrix(NA, nrow = nroad_transects, ncol = nmonth)
      for (t in 1:nmonth) {
        if (length(lambdasim) >= t) {
          N_itd[, t] <- rpois(nroad_transects, lambda = lambdasim[t])
        }
      }
      
      # SAVE N_itd
      saveRDS(N_itd, file = paste0(save_path, group, "_SE_pp_", se_pp, 
                                   "_nroad_transects_", nroad_transects, 
                                   "_SD_N_", sd_lambda, "_SD_ppd_", sd_ppd, 
                                   "_N_itd", j, "_sim_", sim, ".RData"))
      
      # --- N3_itd: Binomial (Carcass persistence bias) ---
      ppd_t <- if (j == 1) prod(ppd) else prod(ppd[j:length(ppd)])
      N3_itd <- matrix(NA, nrow = nroad_transects, ncol = nmonth)
      for (i in 1:nroad_transects) {
        for (t in 1:nmonth) {
          N3_itd[i, t] <- rbinom(1, N_itd[i, t], ppd_t) 
        }
      }
      
      # SAVE N3_itd
      saveRDS(N3_itd, file = paste0(save_path, group, "_SE_pp_", se_pp, 
                                    "_nroad_transects_", nroad_transects, 
                                    "_SD_N_", sd_lambda, "_SD_ppd_", sd_ppd, 
                                    "_N3_itd", j, "_sim_", sim, ".RData"))
      
      # Accumulate N3
      N3_itD <- N3_itD + N3_itd
    }
    
    # SAVE N3_itD (Accumulated)
    saveRDS(N3_itD, file = paste0(save_path, group, "_SE_pp_", se_pp, 
                                  "_nroad_transects_", nroad_transects, 
                                  "_SD_N_", sd_lambda, "_SD_ppd_", sd_ppd, 
                                  "_N3_itD_sim_", sim, ".RData"))
    
    
    # 5. Census Data (Carcass observation bias)
    C <- array(dim = c(nroad_transects, nmonth, nmethod, nsurveys))
    for (m in 1:nmethod) {
      for (k in 1:nsurveys) {
        C[,,m,k] <- matrix(rbinom(nroad_transects * nmonth, N3_itD, selected_group$pom_group[m]), 
                           nrow = nroad_transects, ncol = nmonth)
      }
    }
    
    # Save Final Census Data
    saveRDS(C, file = paste0(save_path, group, "_SE_pp_", se_pp, 
                             "_nroad_transects_", nroad_transects, 
                             "_SD_N_", sd_lambda, "_SD_ppd_", sd_ppd, 
                             "_census_data_sim_", sim, ".RData"))
    
    
    # --- PART B: JAGS ANALYSIS ---
    # Running sequentially to use the data generated in this iteration
    
    tryCatch({
      current_census_data <- C
      
      # Data List 
      bdata <- list(
        C = current_census_data,
        nroad_transects = dim(current_census_data)[1],
        nmonth = dim(current_census_data)[2],
        nmet = dim(current_census_data)[3],
        nsurveys = dim(current_census_data)[4],
        alpha_prior_pp = alpha_prior_pp,
        beta_prior_pp = beta_prior_pp,
        lambda_upper_limit = lambda_upper_limit
      )
      
      # JAGS Model 
      model_file <- paste0("model_no_persistence_bias_scenario_", scenario_id, ".txt")
      sink(model_file)
      cat("
      model {
        # Priors
        pp ~ dbeta(alpha_prior_pp, beta_prior_pp)
        
        for (m in 1:nmet) {
          po[m] ~ dunif(0, 1)
        }
       
        for (t in 1:nmonth) {
          lambda[t] ~ dunif(0, lambda_upper_limit)
        }
        
        # Likelihood
        for (i in 1:nroad_transects) {
          for (t in 1:nmonth) {
            N[i, t] ~ dpois(lambda[t])
            # N3 comes directly from N 
            N3[i, t] ~ dbin(pp, N[i, t])T(0, 100000)
            
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
          totalN3[t] <- sum(N3[, t])    
        }
      }
      ", fill = TRUE)
      sink()
      
      # Initial Values
      Nst <- apply(current_census_data, c(1, 2), max, na.rm = TRUE) + 500
      Nst[!is.finite(Nst)] <- 500
      
      inits <- function() {
        list(
          lambda = runif(nmonth, 5, 100),
          pp = rbeta(1, alpha_prior_pp, beta_prior_pp),
          po = runif(3),
          N = Nst
        )
      }
      
      params <- c("lambda", "pp", "po", "totalN", "totalN3")
      
      # Run JAGS
      output <- jags(
        bdata, inits, params, model_file,
        n.adapt = 100000, n.chains = 3, n.thin = 1000, 
        n.iter = 400000, n.burnin = 100000, parallel = TRUE
      )
      
      # Save Output (Long Filename)
      saveRDS(output, file = paste0(save_path, group, "_SE_pp_", se_pp, 
                                    "_nroad_transects_", nroad_transects, 
                                    "_SD_N_", sd_lambda, "_SD_ppd_", sd_ppd, 
                                    "_output_analysis_sim_", sim, ".RData"))
      
      unlink(model_file) # Clean up text file
      
    }, error = function(e) {
      warning(paste("Error in Scenario", scenario_id, "Simulation", sim, ":", conditionMessage(e)))
    })
    
  } # End of simulation loop
}


# -----------------------------------------------------------------------------
# 7. Execution Loop: This section iterates through the defined scenarios
# -----------------------------------------------------------------------------

# Define all simulation scenarios
scenarios <- list(
  list(id = 1, sd_lambda = 0.5, sd_ppd = 0.05),
  list(id = 2, sd_lambda = 0.5, sd_ppd = 0.15),
  list(id = 3, sd_lambda = 0.5, sd_ppd = 0),
  list(id = 4, sd_lambda = 1.5, sd_ppd = 0.05),
  list(id = 5, sd_lambda = 1.5, sd_ppd = 0.15),
  list(id = 6, sd_lambda = 1.5, sd_ppd = 0),
  list(id = 7, sd_lambda = 0,   sd_ppd = 0.05),
  list(id = 8, sd_lambda = 0,   sd_ppd = 0.15),
  list(id = 9, sd_lambda = 0,   sd_ppd = 0)
)

# Run Loop
for (scn in scenarios) {
  run_scenario_simulation(sd_lambda = scn$sd_lambda, 
                          sd_ppd = scn$sd_ppd, 
                          scenario_id = scn$id)
}

print("All simulations completed successfully.")