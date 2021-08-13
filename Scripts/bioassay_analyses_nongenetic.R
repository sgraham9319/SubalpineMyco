# bioassay_analyses_nongentic
library(dplyr)


colonization <- read.csv("Data/bioassay/raw_data/pct_colonization_bioassay.csv")

colonization_overlap <- colonization[c(1,18,44,55,56,57,64), c("num_myco_stuart", "num_nonmyco_stuart", "num_myco_lauren", "num_nonmyco_lauren")] 
#extracting rows with Lauren and Stuart data

colonization_overlap["sum_stuart"] <- rowSums(colonization_overlap[,1:2]) 
#adding a column of the sum of root tips for Stuart using rowSums
colonization_overlap["sum_lauren"] <- rowSums(colonization_overlap[,3:4])
#adding a column of the sum of roots tips for Lauren using rowSums

colonization_overlap["pct_myco_stuart"] <- colonization_overlap[,1]/colonization_overlap[,5] *100
#creating a percent colonization of myco for Stuarts data
colonization_overlap["pct_myco_lauren"] <- colonization_overlap[,3]/colonization_overlap[,6] *100

t.test(x = colonization_overlap$pct_myco_stuart, y = colonization_overlap$pct_myco_lauren, paired = TRUE)
cor.test(x = colonization_overlap$pct_myco_stuart, y = colonization_overlap$pct_myco_lauren)
plot(x = colonization_overlap$pct_myco_stuart, y = colonization_overlap$pct_myco_lauren)
?abline()
abline(a = 0, b = 1)
treatment <- read.csv("Data/bioassay/raw_data/treatment_bioassay.csv")
col_treat <- colonization %>% left_join(treatment)#, by = seedling_id)
