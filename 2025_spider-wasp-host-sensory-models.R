#load required packages; if they aren't already installed you will need to do that step first
library(ape)
library(geiger)
library(nlme)
library(phytools)
library(picante)
library(caper)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(stringr)
library(multcomp)
library(table1)

#Open the data for analysing sensory organ morphology between wasps that speciliase in spiders with different types of morphology, and call this object "waspHostData"
waspHostData<-read.csv("pompilid_host_sensory_data.csv", header = TRUE)
#Read the phylogenetic tree for the speices of spider wasp in this dataset, and create an object from it called "waspHostTree"
waspHostTree<-read.tree("pompilinae_pruned.tre")
#Plot this phylogenetic tree
plot(waspHostTree, cex=0.7) # cex adjusts text size
d<-matrix(nrow=1, ncol=24)
colnames(d)<-waspHostData$species
#creates a list based on the species listed in datasheet with phenotypic traits
prunedTree<-prune.sample(phylo=waspHostTree,samp=d)
#creates tree where all branches are pruned that lack entry with host phenotypic traits in dataset 
plot(prunedTree)
#the first row of the dataset are headings for the columns
head(waspHostData)
#Check the names between the tree and dataset match
name.check(prunedTree, waspHostData)
#create an object called "behaviour.simplenames" from column 3 of the waspHostData dataset, which will indicate if the wasp species is a host generalist or specialist
behaviour.simplenames<-waspHostData[,3]
names(behaviour.simplenames)<-rownames(waspHostData)
# assign colors to whether wasp is a host generalist or specialist
cols<-setNames(palette()[1:length(unique(behaviour.simplenames))],sort(unique(behaviour.simplenames)))
# plot tree
plotTree(prunedTree, fsize=0.8,ftype="i");
# add tip labels to indicate whether a host specialist or generalist
tiplabels(pie=to.matrix(behaviour.simplenames,sort(unique(behaviour.simplenames))),piecol=cols,cex=0.4)
# add legend to the labelled tree
add.simmap.legend(colors=cols,vertical=TRUE,prompt=TRUE,fsize=0.8)

#create a generalised linear model called "modelOmmat"; this model has ommatidia diamater ("ommatidia" column in the dataset) as the dependent variable with the variables (column names) spider.specialism (generalist or specialist), hosts_webs (host spider is web building or not), hosts_activity (diurnal & nocturnal, nocturnal), host_vert_stratification(ground only, vegetation only, ground and vegetation), and wasp head width (mm) as fixed factors
modelOmmat<-gls(ommatidia ~ spider.specialism + hosts_webs + hosts_activity + host_vert_stratification + head_width, data=waspHostData, correlation=corBrownian(1, phy=prunedTree), method="ML", verbose = TRUE)
#run an ANOVA on this model and produce a summary, model coefficients, and confidence intervals
anova(modelOmmat, verbose = TRUE)
coef(modelOmmat)
summary(modelOmmat)
intervals(modelOmmat)

#define host_vert_stratification as a factor with three levels: "ground", "ground and vegetation", "vegetation"; this is the vertical stratification of the habitat of the host spiders parasitised by a wasp species
waspHostData$host_vert_stratification <- factor(waspHostData$host_vert_stratification, levels=c("Ground", "Ground and vegetation", "Vegetation"))
class(waspHostData$host_vert_stratification)

#create a generalised linear model for the compound eye surface area, using fixed factors as per line 43
modelEyeSA<-gls(eye_sa ~ spider.specialism + hosts_webs + hosts_activity + host_vert_stratification + head_width, data=waspHostData, correlation=corBrownian(1, phy=prunedTree), method="ML")
#create an oboject called anovaEyeSA that will be an ANOVA on that model, then run the ANOVA and provide the outputs
anovaEyeSA <- anova(modelEyeSA)
anova(modelEyeSA)
coef(modelEyeSA)
summary(modelEyeSA)
intervals(modelEyeSA)

#the effect of host activity time on wasp eye surface area was p < 0.05 in the model, so create a boxplot to visualise that pattern 
eyeSA_activity_plot <- ggplot(waspHostData, aes(x = hosts_activity, 
                                         y = eye_sa,
                                         fill = hosts_activity)) +
  geom_boxplot(width=0.6, position = position_dodge(0.7)) + 
  labs(y = "Eye surface area", x = "Host spider activity time") + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black", fill="black") +
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.1) +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line = element_line(size = 0.5, linetype = "solid"), 
        panel.grid.major = element_line(colour = NA), 
        axis.text = element_text(colour = "black", hjust = 0, face = "plain"), 
        axis.text.x = element_text(size = 10, colour = "black", face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10, colour = "black", face = "bold",  vjust = 1),
        panel.background = element_rect(fill = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10, face = "bold"), 
        legend.background = element_rect(fill = "white"), 
        legend.position = "none", legend.direction = "horizontal",
        plot.title = element_text(size = 10, face = "bold"))
eyeSA_activity_plot

#the effect of host vertical stratification on wasp eye surface area was p < 0.05 in the model, so create a boxplot to visualise that pattern
eyeSA_strat_plot <- ggplot(waspHostData, aes(x = host_vert_stratification, 
                                                y = eye_sa,
                                                fill = host_vert_stratification)) +
  geom_boxplot(width=0.6, position = position_dodge(0.7)) + 
  labs(y = "Eye surface area", x = "Host spider vertical stratification") + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black", fill="black") +
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.1) +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line = element_line(size = 0.5, linetype = "solid"), 
        panel.grid.major = element_line(colour = NA), 
        axis.text = element_text(colour = "black", hjust = 0, face = "plain"), 
        axis.text.x = element_text(size = 10, colour = "black", face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10, colour = "black", face = "bold",  vjust = 1),
        panel.background = element_rect(fill = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10, face = "bold"), 
        legend.background = element_rect(fill = "white"), 
        legend.position = "none", legend.direction = "horizontal",
        plot.title = element_text(size = 10, face = "bold"))
eyeSA_strat_plot

#create a generalised linear model for the density of antennal sensilla trichodea, using fixed factors as per line 43, then run an ANOVA
modelTrich<-gls(trichodea ~ spider.specialism + hosts_webs + hosts_activity + host_vert_stratification + head_width, data=waspHostData, correlation=corBrownian(1, phy=prunedTree), method="ML")
anova(modelTrich)
coef(modelTrich)
intervals(modelTrich)

#create a generalised linear model for the density of antennal sensilla basiconica, using fixed factors as per line 43, then run an ANOVA
modelBasiconica<-gls(basiconica ~ spider.specialism + hosts_webs + hosts_activity + host_vert_stratification + head_width, data=waspHostData, correlation=corBrownian(1, phy=prunedTree), method="ML")
anova(modelBasiconica)
coef(modelBasiconica)
intervals(modelBasiconica)

#create a generalised linear model for the density of antennal sensilla chaetica, using fixed factors as per line 43, then run an ANOVA
modelChaetica<-gls(chaetica ~ spider.specialism + hosts_webs + hosts_activity + host_vert_stratification + head_width, data=waspHostData, correlation=corBrownian(1, phy=prunedTree), method="ML")
anova(modelChaetica)
coef(modelChaetica)
intervals(modelChaetica)

#the effect of host activity time on wasp antennal sensilla chaetica density was p < 0.05 in the model, so create a boxplot to visualise that pattern
chaetica_activity_plot <- ggplot(waspHostData, aes(x = hosts_activity, 
                                                y = chaetica,
                                                fill = hosts_activity)) +
  geom_boxplot(width=0.6, position = position_dodge(0.7)) + 
  labs(y = "Sensilla chaetica density", x = "Host spider activity time") + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black", fill="black") +
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.1) +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line = element_line(size = 0.5, linetype = "solid"), 
        panel.grid.major = element_line(colour = NA), 
        axis.text = element_text(colour = "black", hjust = 0, face = "plain"), 
        axis.text.x = element_text(size = 10, colour = "black", face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10, colour = "black", face = "bold",  vjust = 1),
        panel.background = element_rect(fill = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10, face = "bold"), 
        legend.background = element_rect(fill = "white"), 
        legend.position = "none", legend.direction = "horizontal",
        plot.title = element_text(size = 10, face = "bold"))
chaetica_activity_plot

#the effect of host web building on wasp antennal sensilla chaetica density was p < 0.05 in the model, so create a boxplot to visualise that pattern
chaetica_web_plot <- ggplot(waspHostData, aes(x = hosts_webs, 
                                                   y = chaetica,
                                                   fill = hosts_webs)) +
  geom_boxplot(width=0.6, position = position_dodge(0.7)) + 
  labs(y = "Sensilla chaetica density", x = "Host spider web status") + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black", fill="black") +
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.1) +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line = element_line(size = 0.5, linetype = "solid"), 
        panel.grid.major = element_line(colour = NA), 
        axis.text = element_text(colour = "black", hjust = 0, face = "plain"), 
        axis.text.x = element_text(size = 10, colour = "black", face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10, colour = "black", face = "bold",  vjust = 1),
        panel.background = element_rect(fill = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10, face = "bold"), 
        legend.background = element_rect(fill = "white"), 
        legend.position = "none", legend.direction = "horizontal",
        plot.title = element_text(size = 10, face = "bold"))
chaetica_web_plot

#create a generalised linear model for the density of antennal sensilla placodea, using fixed factors as per line 43, then run an ANOVA
modelPlacodea<-gls(placodea ~ spider.specialism + hosts_webs + hosts_activity + host_vert_stratification + head_width, data=waspHostData, correlation=corBrownian(1, phy=prunedTree), method="ML")
anova(modelPlacodea)
coef(modelPlacodea)
intervals(modelPlacodea)

#the effect of host web building on wasp antennal sensilla placodea density was p < 0.05 in the model, so create a boxplot to visualise that pattern
placodea_web_plot <- ggplot(waspHostData, aes(x = hosts_webs, 
                                              y = placodea,
                                              fill = hosts_webs)) +
  geom_boxplot(width=0.6, position = position_dodge(0.7)) + 
  labs(y = "Sensilla placodea density", x = "Host spider web status") + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black", fill="black") +
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.1) +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line = element_line(size = 0.5, linetype = "solid"), 
        panel.grid.major = element_line(colour = NA), 
        axis.text = element_text(colour = "black", hjust = 0, face = "plain"), 
        axis.text.x = element_text(size = 10, colour = "black", face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10, colour = "black", face = "bold",  vjust = 1),
        panel.background = element_rect(fill = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10, face = "bold"), 
        legend.background = element_rect(fill = "white"), 
        legend.position = "none", legend.direction = "horizontal",
        plot.title = element_text(size = 10, face = "bold"))
placodea_web_plot

#create a generalised linear model for the density of antennal sensilla coeloconica, using fixed factors as per line 43, then run an ANOVA
modelCoeloconica<-gls(coeloconica ~ spider.specialism + hosts_webs + hosts_activity + host_vert_stratification + head_width, data=waspHostData, correlation=corBrownian(1, phy=prunedTree), method="ML")
anova(modelCoeloconica)
coef(modelCoeloconica)
intervals(modelCoeloconica)

#create a generalised linear model for the density of antennal sensilla ampullacea, using fixed factors as per line 43, then run an ANOVA
modelAmpullacea<-gls(ampullacea ~ spider.specialism + hosts_webs + hosts_activity + host_vert_stratification + head_width, data=waspHostData, correlation=corBrownian(1, phy=prunedTree), method="ML")
anova(modelAmpullacea)
coef(modelAmpullacea)
intervals(modelAmpullacea)

#the effect of host web building on wasp antennal sensilla ampullacea density was p < 0.05 in the model, so create a boxplot to visualise that pattern
ampullacea_web_plot <- ggplot(waspHostData, aes(x = hosts_webs, 
                                              y = ampullacea,
                                              fill = hosts_webs)) +
  geom_boxplot(width=0.6, position = position_dodge(0.7)) + 
  labs(y = "Sensilla ampullacea density", x = "Host spider web status") + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black", fill="black") +
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.1) +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line = element_line(size = 0.5, linetype = "solid"), 
        panel.grid.major = element_line(colour = NA), 
        axis.text = element_text(colour = "black", hjust = 0, face = "plain"), 
        axis.text.x = element_text(size = 10, colour = "black", face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10, colour = "black", face = "bold",  vjust = 1),
        panel.background = element_rect(fill = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10, face = "bold"), 
        legend.background = element_rect(fill = "white"), 
        legend.position = "none", legend.direction = "horizontal",
        plot.title = element_text(size = 10, face = "bold"))
ampullacea_web_plot
