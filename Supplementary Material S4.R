### Paper title: A novel method to estimate actual infrastructure-induced mortality by integrating sampling biases

### Supplementary Material S4s

### R script to simulate and analyse road survey census data for vertebrate groups only affected by carcass persistence bias (mammals g5)

## R version 4.2.2

#In this script, we simulate census data for mammals g5 group.
#Using the modeling framework outlined in Section 2.2, this script estimate 
#the simulated values for the total number of roadkills, as well as the carcass 
#location (pl) and carcass observation probabilities (pom). 
#We use in the analysis standard errors for the pl prior set at either 0.05 or 0.1.

#Specifics for Mammals G5:
#This section simulates census data assuming that mammals g5 carcasses are not 
#affected by carcass persistence probability (they remain on the road all month).
#Therefore, we simulate a single lambda_td value for the entire month rather than 
#a daily loop.

#Here's an overview of the different sections of the code

#1. Global Simulation Parameters Setup: This part of the code let you customise 
#each parameter in the simulation of data census and their posterior analysis

#2. Load and Monthly Trend: This part of the code load and shows in a 
#plot the annual tendency of roadkill abundance along a year

#3. Calculate alpha and beta Parameters of the Beta Distributed Prior for pl 
#Based on Roman et al. (2024): This section makes a beta distributed prior from
#a mean value and standar error of pl

#4. Census data Simulation and Analysis: this section applies the model 
#to estimate the total number of roadkills and pom (po1 =walking, po2 = cycling 
#and po3 = driving), while incorporating information of pl as prior distributions


#5. Execution

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

# Number of survey road_transects
nroad_transects <- 10  # Scenarios: 10 or 100

# Number of months included in the simulation
nmonth <- 12  

# Number of survey methods considered in the simulation
nmethod <- 3  

# Number of simulations to run
nsimulations <- 20  

# Standard Error for Carcass Location Probability (pl) prior distribution
se_pl <- 0.05  # Scenarios: 0.05 or 0.10  

# Animal Group Select Fuction

select_group <- function(group_name) {
  if (group_name == "mammals_g5") {
    anual_tendency <- "Supplementary material S5/anual_tendency_mammals_g5.csv"
    pl_group <- 0.5              # Carcass location probability
    pom_group <- c(1, 0.9, 0.8)  # Carcass observation probabilities
  } else {
    stop("Invalid group. Please choose 'mammals_g5'.")
  }
  
  return(list(anual_tendency= anual_tendency, pl_group = pl_group, pom_group = pom_group))
}

# Select vertebrate group 
group <- "mammals_g5"  
selected_group <- select_group(group)

# -----------------------------------------------------------------------------
# 2. Load Monthly Trend
# -----------------------------------------------------------------------------

# Load Annual Trend
if(file.exists(selected_group$anual_tendency)){
  df <- read.csv(selected_group$anual_tendency, sep = ";")
  monthly_trend <- df$Roadkill
} else {
  stop(paste("Annual trend file not found:", selected_group$anual_tendency))
}

# -----------------------------------------------------------------------------
# 3. Calculate alpha and beta Parameters of the Beta Distributed Prior for pl
# -----------------------------------------------------------------------------

mean_prior_pl <- selected_group$pl_group  # Prior Distribution Mean
se_prior_pl <- se_pl  # Prior Distribution Standard error

# Calculate alpha and beta parameters of the beta distribution
alpha_prior_pl <- ((1 - mean_prior_pl) / (se_prior_pl^2) - 1 / mean_prior_pl) * (mean_prior_pl^2)
beta_prior_pl <- alpha_prior_pl * (1 / mean_prior_pl - 1)

cat("Prior Setup Complete. Alpha_pl:", alpha_prior_pl, "Beta_pl:", beta_prior_pl, "\n")


# -----------------------------------------------------------------------------
# 4. Census data Simulation and Analysis 
# -----------------------------------------------------------------------------

run_simulation_mammals_g5 <- function() {
  
  message(paste0("\n--- Starting Simulation for Mammals G5 (Single Scenario) ---"))
  
  for (sim in 1:nsimulations) {
    
    message(paste0("   > Running Simulation ", sim, " of ", nsimulations, "..."))
    
    # --- PART A: DATA GENERATION ---
    
    # 1. Lambda Simulation (Single value per month, no daily variation)
    lambdasim <- monthly_trend  
    
    # Save lambdasim
    # NOTE: No loop 'i' here because it's a single step per month
    saveRDS(lambdasim, file = paste0(save_path, group, "SE_pl_", se_pl, 
                                     "_nroad_transects_", nroad_transects, 
                                     "_lambdasim_sim_", sim, ".RData"))
    
    # 2. Total number of roadkills (N_it) based on lambdasim (Poisson)
    N_it <- matrix(NA, nrow = nroad_transects, ncol = nmonth)
    
    for (t in 1:nmonth) {
      # Use length check to stay safe
      if (length(lambdasim) >= t) {
        N_it[, t] <- rpois(nroad_transects, lambda = lambdasim[t]) 
      }
    }
    
    # Save N_it
    saveRDS(N_it, file = paste0(save_path, group, "SE_pl_", se_pl, 
                                "_nroad_transects_", nroad_transects, 
                                "_N_it_sim_", sim, ".RData"))
    
    # 3. N2_itd: Binomial (Carcass location bias) 
    N2_it <- matrix(NA, nrow = nroad_transects, ncol = nmonth)
    
    for (i in 1:nroad_transects) {
      for (t in 1:nmonth) {
        N2_it[i, t] <- rbinom(1, N_it[i, t], selected_group$pl_group)
      }
    }
    
    # Save N2_it
    saveRDS(N2_it, file = paste0(save_path, group, "SE_pl_", se_pl, 
                                 "_nroad_transects_", nroad_transects, 
                                 "_N2_it_sim_", sim, ".RData"))
    
    # 4. Census Data Simulation (Carcass observation Bias)
    nsurveys <- 3 
    C <- array(dim = c(nroad_transects, nmonth, nmethod, nsurveys))
    
    for (m in 1:nmethod) {
      for (k in 1:nsurveys) {
        C[,,m,k] <- matrix(rbinom(nroad_transects * nmonth, N2_it, selected_group$pom_group[m]), 
                           nrow = nroad_transects, ncol = nmonth)
      }
    }
    
    # Save Census Data
    saveRDS(C, file = paste0(save_path, group, "SE_pl_", se_pl, 
                             "_nroad_transects_", nroad_transects, 
                             "_census_data_sim_", sim, ".RData"))
    
    
    # --- PART B: JAGS ANALYSIS ---
    
    tryCatch({
      current_census_data <- C
      
      # Data List
      bdata <- list(
        C = current_census_data,
        nroad_transects = dim(current_census_data)[1],
        nmonth = dim(current_census_data)[2],
        nmet = dim(current_census_data)[3],
        nsurveys = dim(current_census_data)[4],
        alpha_prior_pl = alpha_prior_pl,
        beta_prior_pl = beta_prior_pl
      )
      
      # JAGS Model 
      model_file <- paste0("modelplpo(m)_census_data_sim_", sim, ".txt")
      sink(model_file)
      cat("
      model {
        # Priors
        pl ~ dbeta(alpha_prior_pl, beta_prior_pl)  # Carcass location probability
        
        for (m in 1:nmet) {
          po[m] ~ dunif(0, 1)  # Carcass observation probability
        }
       
        for (t in 1:nmonth) {  # Loop over months
          lambda[t] ~ dunif(0, 300)  # Expected abundance
        }
        
        # Likelihood
        for (i in 1:nroad_transects) {
          for (t in 1:nmonth) {
            N[i, t] ~ dpois(lambda[t])                    # Total number of roadkill
            N2[i,t] ~ dbin(pl, N[i,t])T(0,100000)         # Roadkills inside road
            
            for (m in 1:nmet) {
              for (j in 1:nsurveys) {
                C[i, t, m, j] ~ dbin(po[m], N2[i, t])T(0, 100000)  # Census data
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
      Nst <- apply(current_census_data, c(1, 2), max, na.rm = TRUE) + 500
      Nst[!is.finite(Nst)] <- 500 
      
      inits <- function() {
        list(
          lambda = runif(dim(current_census_data)[2], 5, 100), 
          pl = rbeta(1, alpha_prior_pl, beta_prior_pl),        
          po = runif(3),                                       
          N = Nst                                              
        )
      }
      
      params <- c("lambda", "pl", "po", "totalN", "totalN2")
      
      # MCMC Settings
      na <- 100000 
      ni <- 400000 
      nt <- 1000   
      nb <- 100000 
      nc <- 3      
      
      # Run JAGS
      output <- jags(bdata, inits, params, model_file,
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE
      )
      
      # Save Output
      saveRDS(output, file = paste0(save_path, group, "SE_pl_", se_pl, 
                                    "_nroad_transects_", nroad_transects, 
                                    "_output_analysis_sim_", sim, ".RData"))
      
      unlink(model_file) 
      
    }, error = function(e) {
      warning(paste("Error in Simulation", sim, ":", conditionMessage(e)))
    })
    
  } # End of simulation loop
}


# -----------------------------------------------------------------------------
# 5. EXECUTION
# -----------------------------------------------------------------------------

# Run the simulation function
run_simulation_mammals_g5()

print("All Mammals G5 simulations completed successfully.")