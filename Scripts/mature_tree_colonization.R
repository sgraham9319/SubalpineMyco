#####################################
# Mature tree % colonization analysis
#####################################

# Load required packages
library(dplyr)
library(ggpubr)

# Create standard error function
se <- function(x){sd(x, na.rm = T) / sqrt(length(na.omit(x)))}

# Load mature root data
roots <- read.csv("Data/field/pct_colonization_mature_tree.csv")

# Calculate % colonization
roots <- roots %>%
  mutate(pct_col = 100 * (myco_tips / total_tips))

# Did forest and meadow samples differ in % colonization?
tapply(roots$pct_col, roots$location, mean)
tapply(roots$pct_col, roots$location, se)
t.test(pct_col ~ location, data = roots, paired = T)

# Format data for plotting
plot_data <- roots %>%
  rename("Tree location" = location)

# Plot data
ggpaired(plot_data, x = "Tree location", y = "pct_col", fill = "Tree location",
         line.color = "gray", palette = "jco",
         xlab = "Tree location", ylab = "Root-tips colonized (%)") +
  stat_compare_means(method = "t.test", paired = TRUE, label.y.npc = "bottom")

# Did forest and meadow samples differ in number morphotypes?
tapply(roots$num_morphotypes, roots$location, mean)
tapply(roots$num_morphotypes, roots$location, se)
t.test(num_morphotypes ~ location, data = roots, paired = T)
ggpaired(roots, x = "location", y = "num_morphotypes", fill = "location",
         line.color = "gray", palette = "jco",
         xlab = "Environment", ylab = "# morphotypes") +
  stat_compare_means(paired = TRUE, label.y.npc = "bottom")



