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
library(cowplot)

#Open the data for analysing sensory organ morphology between kleptoparasitic and parasitoid spider wasps, and call this object "waspData"
waspData<-read.csv("pompilid_sensory_data.csv", header = TRUE)
#Read the phylogenetic tree for the speices of spider wasp in this dataset, and create an object from it called "waspTree"
waspTree<-read.tree("pompilinae_pruned.tre")
#Plot this phylogenetic tree
plot(waspTree, cex=0.7) # cex adjusts text size
d<-matrix(nrow=1, ncol=26)
colnames(d)<-waspData$species
#creates a list based on the species listed in datasheet with phenotypic traits
plot(waspTree)
#the first row of the dataset are headings for the columns
head(waspData)
#Check the names between the tree and dataset match
name.check(waspTree, waspData)
#create an object called "behaviour.simplenames" from column 3 of the waspHostData dataset, which will indicate if the wasp species is a kleptoparasite or a parasitoid
behaviour.simplenames<-waspData[,3]
names(behaviour.simplenames)<-rownames(waspData)
# assign colors to behaviour names (kleptoparasite, non-kleptoparasite)
cols<-setNames(palette()[1:length(unique(behaviour.simplenames))],sort(unique(behaviour.simplenames)))
# plot tree
plotTree(waspTree, fsize=0.8,ftype="i");
# add tip labels
tiplabels(pie=to.matrix(behaviour.simplenames,sort(unique(behaviour.simplenames))),piecol=cols,cex=0.4)
# add legend
add.simmap.legend(colors=cols,vertical=TRUE,prompt=TRUE,fsize=0.8)

#Statistical model for compound eye ommatidia diameter
#create a generalised linear model called "modelOmmat"; this model has ommatidia diamater ("ommatidia" column in the dataset) as the dependent variable with the variables (column names) behaviour.simple (kleptoparasite, non-kleptoparasite) and wasp head width (mm) as fixed factors
modelOmmat<-gls(ommatidia ~ behaviour.simple + head_width, data=waspData, correlation=corBrownian(1, phy=waspTree), method="ML")
#run an ANOVA on this linear model
anova(modelOmmat)
coef(modelOmmat)
intervals(modelOmmat)

#Statistical model for compound eye surface area with variables as per line 40
modelEyeSA<-gls(eye_sa ~ behaviour.simple + head_width, data=waspData, correlation=corBrownian(1, phy=waspTree), method="ML")
anova(modelEyeSA)
coef(modelEyeSA)
intervals(modelEyeSA)

#Plot for relationship between eye surface area and parasitism strategy, which was p < 0.05 in model
eyeSA_behav_plot <- ggplot(waspData, aes(x = behaviour.simple, 
                                         y = eye_sa,
                                         fill = behaviour.simple)) +
  geom_boxplot(width=0.6, position = position_dodge(0.7)) + 
  scale_fill_manual(values = c("#F0E442","#0072B2")) +
  labs(y = "Eye surface area", x = "Parasitism strategy") + 
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
eyeSA_behav_plot

#Statistical model for density of antennal sensilla trichodea
modelTrich<-gls(trichodea ~ behaviour.simple + head_width, data=waspData, correlation=corBrownian(1, phy=waspTree), method="ML")
anova(modelTrich)
coef(modelTrich)
intervals(modelTrich)

#Statistical model for density of antennal sensilla basiconica
modelBasiconica<-gls(basiconica ~ behaviour.simple + head_width, data=waspData, correlation=corBrownian(1, phy=waspTree), method="ML")
anova(modelBasiconica)
coef(modelBasiconica)
intervals(modelBasiconica)

#Statistical model for density of antennal sensilla chaetica
modelChaetica<-gls(chaetica ~ behaviour.simple + head_width, data=waspData, correlation=corBrownian(1, phy=waspTree), method="ML")
anova(modelChaetica)
coef(modelChaetica)
intervals(modelChaetica)

#Statistical model for density of antennal sensilla placodea
modelPlacodea<-gls(placodea ~ behaviour.simple + head_width, data=waspData, correlation=corBrownian(1, phy=waspTree), method="ML")
anova(modelPlacodea)
coef(modelPlacodea)
intervals(modelPlacodea)

#Statistical model for density of antennal sensilla coeloconica
modelCoeloconica<-gls(coeloconica ~ behaviour.simple + head_width, data=waspData, correlation=corBrownian(1, phy=waspTree), method="ML")
anova(modelCoeloconica)
coef(modelCoeloconica)
intervals(modelCoeloconica)

#Plot for relationship between sensilla coeloconica density and parasitism strategy, which was p < 0.05 in model
coeloconica_plot_behaviour <- ggplot(waspData, aes(x = behaviour.simple, 
                                            y = coeloconica,
                                            fill = behaviour.simple)) +
  geom_boxplot(width=0.6, position = position_dodge(0.7)) + 
  scale_fill_manual(values = c("#F0E442","#0072B2")) +
  labs(y = "Coeloconica density", x = "Parasitism strategy") + 
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
coeloconica_plot_behaviour

#Statistical model for density of antennal sensilla ampullacea
modelAmpullacea<-gls(ampullacea ~ behaviour.simple + head_width, data=waspData, correlation=corBrownian(1, phy=waspTree), method="ML")
anova(modelAmpullacea)
coef(modelAmpullacea)
intervals(modelAmpullacea)
