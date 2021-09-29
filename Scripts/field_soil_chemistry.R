###############################
# Field soil chemistry analysis
###############################

# Load required packages
library(ggpubr)

# Create standard error function
se <- function(x){sd(x, na.rm = T) / sqrt(length(na.omit(x)))}

# Load soil chemistry data
chem <- read.csv("../Data/field/soil_chemistry.csv")

# Calculate C:N ratio
chem$cn <- chem$carbon_pct / chem$nitrogen_pct

# Water holding capacity
tapply(chem$whc_pct, chem$location, mean)
tapply(chem$whc_pct, chem$location, se)
t.test(whc_pct ~ location, data = chem, paired = T)
ggpaired(chem, x = "location", y = "whc_pct", fill = "location",
         line.color = "gray", palette = "jco",
         xlab = "Environment", ylab = "whc_pct") +
  stat_compare_means(method = "t.test", paired = TRUE)

# Plant-available phosphorus
tapply(chem$phosphate_mgkg, chem$treatment, mean)
tapply(chem$phosphate_mgkg, chem$treatment, se)
t.test(phosphate_mgkg ~ location, data = chem, paired = T)
ggpaired(chem, x = "location", y = "phosphate_mgkg", fill = "location",
         line.color = "gray", palette = "jco",
         xlab = "Environment", ylab = "phosphate_mgkg") +
  stat_compare_means(method = "t.test", paired = TRUE)

# C:N ratio
tapply(chem$cn, chem$treatment, mean)
tapply(chem$cn, chem$treatment, se)
t.test(cn ~ location, data = chem, paired = T)
ggpaired(chem, x = "location", y = "cn", fill = "location",
         line.color = "gray", palette = "jco",
         xlab = "Environment", ylab = "cn") +
  stat_compare_means(method = "t.test", paired = TRUE)

# Plant-available nitrogen
tapply(chem$ammonium_mg_kg, chem$treatment, mean, na.rm = T)
tapply(chem$ammonium_mg_kg, chem$treatment, se)
t.test(ammonium_mg_kg ~ location, data = chem, paired = T)
tapply(chem$nitrate_mg_kg, chem$treatment, mean, na.rm = T)
tapply(chem$nitrate_mg_kg, chem$treatment, se)
t.test(nitrate_mg_kg ~ location, data = chem, paired = T)
ggpaired(chem, x = "location", y = "nitrate_mg_kg", fill = "location",
         line.color = "gray", palette = "jco",
         xlab = "Environment", ylab = "nitrate_mg_kg") +
  stat_compare_means(method = "t.test", paired = TRUE)

# pH
tapply(chem$ph, chem$treatment, mean)
tapply(chem$ph, chem$treatment, se)
t.test(ph ~ location, data = chem, paired = T)
ggpaired(chem, x = "location", y = "ph", fill = "location",
         line.color = "gray", palette = "jco",
         xlab = "Environment", ylab = "ph") +
  stat_compare_means(method = "t.test", paired = TRUE)
