## This code calculates distance similarity measures and then organizes pairwise comparisons 
library(reshape2)
in_out_beta <- function(phylo, dist) {
     pull out the sample data 
  cent.subsite <- sample_data(phylo)
  
  # set up the distance matrix.  
  dist.mat <- as.matrix(phyloseq::distance(phylo, method = dist, binary = TRUE)) 
  
  cent.dist <- dist.mat
  distance.name = dist
  
  # change the row names to the sample ID
  rownames(cent.subsite) -> cent.subsite$Sample 

  ##------------
  ## Comparison 1.  
  ##------------
  #subset
  ## This gives us the PCoA distances between each point between the pairwise comparision
  cent.dist.1 <- cent.dist[rownames(cent.dist) %in% subset(cent.subsite, Date=="25_Apr")$Sample, ] 
  cent.dist.1 <- cent.dist.1[, colnames(cent.dist) %in% subset(cent.subsite, Date == '05_May')$Sample] 
  
  ## melt it
  cent.melt.1 <- melt(cent.dist.1) 
  cent.melt.1$tag.1 <- cent.subsite$Plot[match(cent.melt.1$Var1, cent.subsite$Sample)] 
  cent.melt.1$tag.2 <- cent.subsite$Plot[match(cent.melt.1$Var2, cent.subsite$Sample)] 
  
  ## type factor
  cent.melt.1$same <- rep(FALSE)
  cent.melt.1$same[which(cent.melt.1$tag.1 == cent.melt.1$tag.2)] <- rep(TRUE)
  
  ## relevent factors
  cent.melt.1$Farming <- cent.subsite$Farming[match(cent.melt.1$Var1, cent.subsite$Sample)] 
  cent.melt.1$Climate <- cent.subsite$Climate[match(cent.melt.1$Var1, cent.subsite$Sample)]
  cent.melt.1$Virus <- cent.subsite$Virus[match(cent.melt.1$Var1, cent.subsite$Sample)] 
  cent.melt.1$Treatment <- cent.subsite$Treatment[match(cent.melt.1$Var1, cent.subsite$Sample)] 
  
  
  ## treatment and Time columns
  cent.1.same <- subset(cent.melt.1, same == TRUE)
  cent.1.same$Date_comparison <- rep('Apr_to_May') 

  
  ##--------------------
  ## merge arrays if you made more than one comparison
  ##--------------------
  #cent.same <- rbind(cent.1.same, cent.2.same, cent.3.same, cent.4.same) 
  cent.same
}
  
## calculations
paired_uJ <- in_out_beta(phylo = EX_ps_clean.rar, dist = 'jaccard')

paired_BC <- in_out_beta(phylo = EX_ps_clean.rar, dist = 'bray') 



#####stats 
summary(aov(value~Farming*Virus*Climate, data=paired_uJ)) 


summary(aov(value~Date_comparison, data=paired_uJ))  


lm <- (lmer(value ~ Currently_have_pets + (1/Home_ID|Cohort) + (1|Cohort) + (1|Daily_window_use_past_month) + (1|Pick_up_month), data=paired_uJ_bac))  
anova(lm)
   


# plot it into violin plots 
theme_minimal() + 
geom_violin(trim = TRUE, aes(color= Climate, fill=Climate)) + 
geom_boxplot(width = 0.1, aes(group=Climate)) + 
ylab("Unweighted Jaccard distance") + xlab("") + 
facet_grid(Date_comparison~Farming) + 
theme(legend.position = "none")  


# plot linear correlation 
ggplot(paired_uJ_bac, aes(x=log(parks_Park1Distance), y= value)) + 
  geom_point(aes(color=Daily_window_use_past_month, shape=Indoor_living_plants), size=3) + 
  geom_smooth(alpha=1/6, method="lm") + 
  theme_minimal() + 
  ylab("Unweighted Jaccard distance between In/Out") + 
  xlab("Log Distance to Nearest Park") + ylim(0.8,1) + 
  labs(color="Daily window use\nin the past month", shape="Indoor living plants")
