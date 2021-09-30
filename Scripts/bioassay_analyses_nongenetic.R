###############################
# bioassay_analyses_nongentic #
###############################

# This script creates and interprets the generalized mixed effects models 
# of colonization 

# Load required packages
library(dplyr)
library(lme4)

#===========
# Load Data
#============

# Load percent colonization data
colonization <- read.csv("Data/bioassay/raw_data/pct_colonization_bioassay.csv")

# Load biomass data
biomass <- read.csv("Data/bioassay/raw_data/biomass.csv")

# Load treatment data
treatment <- read.csv("Data/bioassay/raw_data/treatment_bioassay.csv")

#============================================
# Format Data to Verify Colonization Results
#============================================

# Extract rows containing colonization results from both Lauren and Stuart in order to verify 
# colonization accuracy 
colonization_overlap <- colonization[c(1,18,44,55,56,57,64), c("num_myco_stuart", "num_nonmyco_stuart", "num_myco_lauren", "num_nonmyco_lauren")] 

# Add a column of the sum of root tips for Stuart and Lauren using rowSums
colonization_overlap["sum_stuart"] <- rowSums(colonization_overlap[,1:2]) 
colonization_overlap["sum_lauren"] <- rowSums(colonization_overlap[,3:4])

# Create a column for percent colonization of mycorrhizae with Stuart's results
colonization_overlap["pct_myco_stuart"] <- colonization_overlap[,1]/colonization_overlap[,5] *100

# Create a column for percent colonization of mycorrhizae with Lauren's results
colonization_overlap["pct_myco_lauren"] <- colonization_overlap[,3]/colonization_overlap[,6] *100

#==============
# Explore Data
#==============

# Test for significant difference and correlation between the two groups of results
t.test(x = colonization_overlap$pct_myco_stuart, y = colonization_overlap$pct_myco_lauren, paired = TRUE)
cor.test(x = colonization_overlap$pct_myco_stuart, y = colonization_overlap$pct_myco_lauren)

# Plot data
plot(x = colonization_overlap$pct_myco_stuart, y = colonization_overlap$pct_myco_lauren)
abline(a = 0, b = 1)

#=======================
# Format Treatment Data
#=======================

col_treat <- colonization %>% left_join(treatment)#, by = seedling_id)

# Add a column of the sum of root tips for Stuart 
col_treat["sum_stuart"] <- rowSums(col_treat[,3:4])

# Create a column for percent colonization of mycorrhizae
col_treat["pct_myco_staurt"] <- col_treat$num_myco_stuart / col_treat$sum_stuart *100

#==============
# Explore Data
#==============

# Run a linear regression model to predict percent colonization by treatment type
control_test <- lm(pct_myco_staurt~treatment, data = col_treat)
summary(control_test)

# Plot data
boxplot(pct_myco_staurt ~ treatment, data = col_treat, xlab = "Treatment Type", ylab = "Percent Colonization")

#==================================
# Filter by Experimental Treatment
#==================================

exp_only <- col_treat %>% 
  filter(treatment == "Experimental") 

# Run a linear regression model to predict percent colonization by soil type
exp_test <- lm(pct_myco_staurt~soil_type, data = exp_only)
summary(exp_test)

# Plot Data
boxplot(pct_myco_staurt ~ soil_type, data = exp_only, ylab ="Percent Colinzation", xlab = "Soil Type")

# Check the mean for each soil type
tapply(exp_only$pct_myco_staurt, exp_only$soil_type, FUN = mean)

#=========================================
# Explore Biomass Data Between Soil Types
#=========================================

#============
#Format Data
#============

# Join biomass data with colonization data 
col_biomass <- col_treat %>% left_join(biomass)#, by = seedling_id)

# Filter out experimental treatment data
exp_biomass_only <- col_biomass %>%
  filter(treatment == "Experimental")

#==============
# Explore Data
#==============

# Run a linear regression model to predict biomass by soil type
biomass_test <- lm(biomass_g~soil_type, data = exp_biomass_only)
summary(biomass_test)

#Plot Data
boxplot(biomass_g ~ soil_type, data = exp_biomass_only, xlab = "Soil Type", ylab = "Biomass (grams)")

# Run a linear regression model to predict biomass by treatment type
biomass_control_test <- lm(biomass_g~treatment, data = col_biomass)
summary(biomass_control_test)


fix1 <- lmer(pct_myco_staurt ~ soil_type + (1 | pair_id), REML = FALSE, data = exp_only)
summary(fix1)
fix2 <- lmer(pct_myco_staurt ~ (1 | pair_id), REML = FALSE, data = exp_only)
summary(fix2)

# Compare both models (fix1 and fix2) using likelihood ration test
anova(fix1, fix2)
anova(fix2, fix1)
