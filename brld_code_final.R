#######Data setup###

####call in packages###
library(e1071)
library(tidyverse)
library(ggplot2)
library(emmeans)
library(doBy)
library(multcomp)
library(DHARMa)
library(DataExplorer)
library(glmmTMB)

options(scipen=999)


######################################### WATER CHEMISTRY CODE ########################################################################

###tests on water data
water_data = data.frame(read.csv(file="water_data_BRLD.csv", header=TRUE))
water_data


###Descriptive statistics
####User specified functions###
se =  function(x) sd(x)/(sqrt(length(x)))


##########Basic descriptive statistics################
percentile = c(0.05,0.95, 0.25, 0.75)
Descriptive_stats_water = summarise(group_by(water_data,river),
                              Mean = mean(srca), 
                              N = length(srca), 
                              SD = sd(srca), 
                              SE = se(srca), 
                              Median = median(srca), 
                              Skewness = skewness(srca, type = 2),
                              Kurtosis = kurtosis(srca, type = 2),
                              Fifth_percentile= quantile(srca, probs = percentile[1], type = 2),
                              Ninetyfifth_percentile= quantile(srca, probs = percentile[2], type = 2), 
                              Twenty_five_percentile = quantile(srca, probs = percentile[3], type = 2),
                              Seventy_five_percentile = quantile(srca, probs = percentile[4], type = 2), 
                              IQR = Seventy_five_percentile - Twenty_five_percentile, 
                              Lower_IQR_outliers = Twenty_five_percentile - 1.5*(IQR), 
                              Upper_IQR_outliers = Seventy_five_percentile + 1.5*(IQR))

write.table(Descriptive_stats_water,"descriptive_stats_water_final.csv",sep=",", row.names = T)

#########glm to compare water chemistry among rivers/reaches ######

water_comparison= glm(srca ~ river, family = Gamma(link="log"), data = water_data)
summary(water_comparison)


##### compare water chemistry#####

temp_emmeans_water = emmeans(water_comparison, "river")
contrast(temp_emmeans_water, method = "pairwise", adjust="tukey")
cld(temp_emmeans_water)
plot(temp_emmeans_water, comparisons = T)
plot(residuals(water_comparison))

res_water =  simulateResiduals(water_comparison)
plot(res_water)

###compare water chemistry in upper Des Plaines River#####
DPRdata = read.table(header=T,colClasses=c("factor","numeric"),text="
site SrCa
A   1.95
A   2.95
A   1.75
A   1.84
A   1.61
A   2.01
A   1.91
B   1.91
B   2.12
B   2
B   1.97
B   1.77
")

wilcox.test(SrCa ~ site, data=DPRdata)

#######create water box plot to visualize water chemistry###########
water_box = data.frame(read.csv(file="water_data_BRLD.csv", header=TRUE))

jpeg("water_boxplot_final.jpeg",width = 6.95, height = 4.89,units = 'in', res = 1080)
water_plot = ggplot(water_box, aes(x=river, y=srca, fill=river))+
  stat_boxplot(coef=1.5, geom='errorbar', width=0.5) + 
  geom_boxplot(width=0.5) +
  xlab("\nRiver") + 
  ylab("Water Sr:Ca (mmol/mol)") + 
  ylim(0, 3.5) +
  scale_x_discrete(labels=str_wrap(c("a_DPR_upper" = "Upper Des Plaines", "b_DPR_lower" = "Lower Des Plaines", "CAWS" = "Chicago Area Waterway System",  "ILR" = "Illinois", "KKR" = "Kankakee"), width=15)) +
  scale_fill_manual(labels = c("Upper Des Plaines", "Lower Des Plaines", "Chicago Area Waterway System", "Illinois", "Kankakee"), values = c("white", "white", "gray80", "gray60", "gray30")) 
  water_plot + theme_classic() + 
  theme(text=element_text(family="serif", size=14), axis.text=element_text(size=14, color="black"), legend.position="none", axis.title.y = element_text(vjust=2)) 
dev.off()


####################################### FIN RAY EDGE MICROCHEMISTRY CODE ###########################################################

species_mix = data.frame(read.csv(file="species_combo_final_BRLD.csv", header=TRUE)) 
species_mix

edge_values = summaryBy(srca ~ id, data = species_mix, 
          FUN = list(mean))

write.table(edge_values,"average_edge.csv",sep=",", row.names = T)


species_mix_stats = data.frame(read.csv(file="average_edge_all_values.csv", header=TRUE))

#########User specified functions########## 
se =  function(x) sd(x)/(sqrt(length(x)))
##########Basic descriptive statistics ######### 
percentile = c(0.25,0.75) # Percentiles of interest 
Descriptive_stats_fish = summarise(group_by(species_mix_stats,river, lump), 
                              Mean = mean(srca_edge), 
                              N = length(srca_edge), 
                              SD = sd(srca_edge),
                              SE = se(srca_edge),
                              Median = median(srca_edge), 
                              Skewness = skewness(srca_edge, type = 2),
                              Kurtosis = kurtosis(srca_edge, type = 2),
                              Min_value = min(srca_edge),
                              Max_value = max(srca_edge),
                              Twenty_five_percentile = quantile(srca_edge, probs = percentile[1], type = 2),
                              Seventy_five_percentile = quantile(srca_edge, probs = percentile[2], type = 2), 
                              IQR = Seventy_five_percentile - Twenty_five_percentile, 
                              Lower_IQR_outliers = Twenty_five_percentile - 1.5*(IQR),  ## use to identify upper bound for outliers
                              Upper_IQR_outliers = Seventy_five_percentile + 1.5*(IQR)) ## use to identify lower bound for outliers

write.table(Descriptive_stats_fish,"descriptive_stats_average_edge.csv",sep=",", row.names = T)

####Outlier removal###########
##outliers identified using descriptive statistics run above

##centrarchids removed (7) 
species_mix = droplevels(species_mix[-which(species_mix$id == "D1-19"),])
species_mix = droplevels(species_mix[-which(species_mix$id == "D4-37"),])
species_mix = droplevels(species_mix[-which(species_mix$id == "I5-17"),])
species_mix = droplevels(species_mix[-which(species_mix$id == "I1-64"),])
species_mix = droplevels(species_mix[-which(species_mix$id == "K8-28"),])
species_mix = droplevels(species_mix[-which(species_mix$id == "K1-43"),])
species_mix = droplevels(species_mix[-which(species_mix$id == "K5-7"),])


###catostomids removed (4)
species_mix = droplevels(species_mix[-which(species_mix$id == "I6-6"),])
species_mix = droplevels(species_mix[-which(species_mix$id == "K6-17"),])
species_mix = droplevels(species_mix[-which(species_mix$id == "K6-10"),])
species_mix = droplevels(species_mix[-which(species_mix$id == "K6-13"),])

## no ictalurids removed

###lepisosteids removed (4)
species_mix = droplevels(species_mix[-which(species_mix$id == "D11-7"),])
species_mix = droplevels(species_mix[-which(species_mix$id == "I3-14"),])
species_mix = droplevels(species_mix[-which(species_mix$id == "I3-19"),])
species_mix = droplevels(species_mix[-which(species_mix$id == "I2-12"),])

species_mix 
write.table(species_mix,"edge_values_outliers_removed.csv",sep=",", row.names = T)

##################GLMMTMB with outliers removed ########################
### this  models the mean for the final 25um of each fin ray transect (with outliers removed)
### this will be used for comparison among taxa/rivers using the modelled edge means


edge_values_outliers_removed= glmmTMB(srca ~ lump + river + lump*river + (1|id),
                     data = species_mix,
                     family = Gamma(link="log"),
                     dispformula = ~lump/river + river,
                     na.action = na.omit)
summary(edge_values_outliers_removed)

### compare fin ray edge microchemistry among taxa/rivers

temp_emmeans = emmeans(edge_values_outliers_removed, ~ lump*river)
contrast(temp_emmeans,method = "pairwise", adjust="tukey")
cld(temp_emmeans)
plot(temp_emmeans, comparisons = T)

res =  simulateResiduals(edge_values_outliers_removed)
plot(res)

##### create violin plot to visualize comparisons among taxa/rivers

species_mix_violin = data.frame(read.csv(file="average_edge_outliers_removed.csv"))

jpeg("BRLD_Fin_Ray_violin_average_edge.jpeg",width = 6.95, height = 4.89,units = 'in', res = 1080)
fin_plot = ggplot(species_mix_violin, aes(x=lump, y=srca_edge, fill=river))+
  geom_violin() + 
  xlab("\nTaxonomic Group") + 
  ylab("Fin Ray Sr:Ca (mmol/mol)") + 
  labs(fill = "River") +
  ylim(0, 0.8) +
  scale_x_discrete(labels=c("bass" = "Centrarchids",  "buff" = "Catostomids", "cat" = "Ictalurids", "gar" = "Lepisosteids")) +
  scale_fill_manual(labels = c("Des Plaines", "Illinois", "Kankakee"), values = c("white", "gray70", "gray30"))

fin_plot + theme_classic() + 
  theme(text=element_text(family="serif", size=14), axis.text=element_text(size=14, color="black"), axis.title.y = element_text(vjust=2)) 
dev.off()


#### Run descriptive statistics on mean edge values to find max/min values and overlap zones

Descriptive_stats_outliers_removed = summarise(group_by(species_mix_violin,river, lump), 
                              Mean = mean(srca_edge), 
                              N = length(srca_edge), 
                              SD = sd(srca_edge),
                              SE = se(srca_edge),
                              Median = median(srca_edge), 
                              Min_value = min(srca_edge), #### use to identify overlap ranges for reclassification as "uncertain"
                              Max_value = max(srca_edge)) #### use to identify overlap ranges for reclassification as "uncertain"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 Max_value = max(scra_edge)) #### use to identify overlap ranges for reclassification as "uncertain"
              

write.table(Descriptive_stats_outliers_removed,"descriptive_stats_outliers_removed.csv",sep=",", row.names = T)

### Use the max and min values in table "descriptive_stats_outliers_removed.csv" 
### to classify Sr:Ca data points along each transect as "Des Plaines", "Illinois", "Kankakee", or within overlap zones
### after examining all taxonomic groups for evidence of passage, no evidence of passage, or inteterminate downstream residency, 
### create a histogram showing percentages of each category

fish_hist = as.data.frame(read.table(header=T,colClasses=c("factor", "numeric", "factor"),
                                     text="
Passage Proportion Taxa
aYes 6 bass
aYes 26 buff
aYes 37 cat
aYes 23 gar
bUncertain 91 bass
bUncertain 37 buff
bUncertain 19 cat
bUncertain 50 gar
cNo 3 bass
cNo 37 buff
cNo 44 cat
cNo 27 gar
"))

fish_hist

jpeg("Passage_histogram_BRLD.jpeg",width = 6.95, height = 4.89,units = 'in', res = 1080)
fish_plot = ggplot(data=fish_hist, aes(x=Passage, y=Proportion, fill=Taxa)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_text(aes(label=Proportion), vjust=-1, color="black", position = position_dodge(0.9), size=3.5, family="serif") +
  xlab("\nPassage designation") + 
  ylab("Percent of individuals in taxonomic group") + 
  labs(fill = "Taxonomic group") +
  scale_y_continuous(expand = c(0,0),limits = c(0,100))+
  scale_x_discrete(labels=str_wrap(c("aYes" = "Evidence of passage",  "bUncertain" = "Indeterminate downstream residency", "cNo" = "No evidence of passage"), width=15)) +
  scale_fill_manual(labels = c("Centrarchids (n=114)", "Catostomids (n=27)", "Ictalurids (n=41)", "Lepisosteids (n=22)"), values = c("black","gray70", "gray30", "white"))

fish_plot + theme_classic() + 
  theme(text=element_text(family="serif", size=14), axis.text=element_text(size=14, color="black"), axis.title.y = element_text(vjust=2)) 
dev.off()
