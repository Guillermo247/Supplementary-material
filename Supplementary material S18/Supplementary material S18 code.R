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
library(dplyr)

#Load case study dataset
df <- read.csv("Supplementary material S18/Supplementary material S18 data.csv", 
               sep = ";")

# Count how many collisions with confirmed observation (OBSERVED == 1) have a certain number of unique methods
df %>%
  filter(OBSERVED == 1) %>%  # Filter only seen collisions
  group_by(ID_ROADKILL) %>%  # Group by collision ID
  summarise(unique_methods = n_distinct(METHODOLOGY)) %>%  # Count unique methods per collision
  count(unique_methods)  # Count how many collisions have N unique methods

# Filter only the collisions with OBSERVED == 1
df_seen_1 <- df %>% 
  filter(OBSERVED == 1)

# View the first records to confirm the filtering
head(df_seen_1)

# Filter collisions that have been seen with only one method
df_single_method <- df_seen_1 %>%
  group_by(ID_ROADKILL) %>%
  filter(n_distinct(METHODOLOGY) == 1) %>%
  ungroup()

# View the first records to confirm correct filtering
head(df_single_method)

#Check the number of roadkills recorded by a single method without counting repeated entries of the same roadkill ID by Personnel (P) or Volunteers (V)
n_distinct(df_single_method$ID_ROADKILL)

# Step 1: Filter collision events that were observed either exclusively by one type
#of observer (Personnel (P) or Volunteer (V)), or by both (P_V), but retains only 
#the records reported by Personnel when both types of observers are present.
df_filtered <- df_single_method %>%
  group_by(ID_ROADKILL) %>%
  filter(n_distinct(P_V) == 1 | (n_distinct(P_V) == 2 & P_V == "P")) %>%
  ungroup()

# View the result
head(df_filtered)

# Filter collisions seen only while walking survey method (WALKING)
df_filtered_walking <- df_single_method %>%
  filter(METHODOLOGY == "WALKING")

roadkill_only_walking <- length(unique(df_filtered_walking$ID_ROADKILL))
print(roadkill_only_walking)

# Filter collisions seen only while cycling survey method (CYCLING)
df_filtered_cycling <- df_single_method %>%
  filter(METHODOLOGY == "CYCLING")

roadkill_only_cycling <- length(unique(df_filtered_cycling$ID_ROADKILL))
print(roadkill_only_cycling)

# Filter collisions seen only while driving s¡survey method (DRIVING)
df_filtered_driving <- df_single_method %>%
  filter(METHODOLOGY == "DRIVING")

roadkill_only_driving <- length(unique(df_filtered_driving$ID_ROADKILL))
print(roadkill_only_driving)
