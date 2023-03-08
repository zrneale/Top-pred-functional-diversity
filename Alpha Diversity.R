
#Here I calculate the local (alpha) functional diversity of dragonfly communities in ponds with different top predators - invertebrates, salamanders, green sunfish, and largemouth bass. I'll try three 


#Load packages and data
library(tidyverse)
library(FD)
library(lme4)
library(MuMIn)
library(car)
library(emmeans)



Traitaverage <- read.csv("Data/Traitaverage.csv")
AbundancewideFD <- read.csv("Data/AbundancewideFD.csv")
Pondata <- read.csv("Data/Pondata.csv")%>%
  mutate_at(vars(Pondnum), as.factor)
Timekey <- read.csv("Data/Timekey.csv")
DomPredata <- read.csv("Data/PredType2.csv")



#Some quick data prep for functional diversity calculations. The function to calculate function diversity metrics requires a wide version of abundance data.
#Also, the spID in Traitaverage need to be assigned to row names

#ID variable needs to be moved to row names. Remove samples with fewer than 3 observations to reduce biased averages.
#Also removing samples with <3 observations to avoid biased averages.
AbundancewideFD <- AbundancewideFD%>%
  rowwise()%>%
  filter(sum(c_across(Ajun:Tramea))>=3)%>%
  column_to_rownames("ID")

#spID column in Traitaverage needs to be assigned to row names and removed
Traitaverage <- Traitaverage%>%
  column_to_rownames("spID")

#Calculate FD metrics
FD<-dbFD(Traitaverage, AbundancewideFD, w.abun=T)

#Compile the desired metrics into a data frame and format the predictor variables
FDistance <- FD%>%
  data.frame()%>%
  dplyr::select(FDis, FEve, FRic)%>% #Specify the metrics you want to analyze
  rownames_to_column("ID")%>% #Row names contain the united predictor variables
  separate(ID, c("dompred", "Pondnum", "year", "season"), convert=F)%>% 
  mutate_at(vars(Pondnum), as.factor)%>%
  mutate_if(is.character, as.factor)%>%
  left_join(dplyr::select(Pondata, -dompred), by = "Pondnum")


#Reorder dompred and season levels for graphing
FDistance$season <- factor(FDistance$season, levels = c("W","Sp","Su","F"))
FDistance$dompred <- factor(FDistance$dompred, levels = c("N","S","G","B"))


#Check for correlation between any of the predictor variables

FDistance%>%
  dplyr::select(dompred, area, season, year, perimeter, depth, veg, canopy)%>%
  mutate_all(as.numeric)%>%
  cor()

#None of them look overly correlated. Time to start analysis.


#Select appropriate model. I want to see if any of the environmental predictors improve the model fit.

model <- lmer(FDis ~ dompred + year + season + dompred:year + dompred:season  +  
                            area + depth + veg + canopy + (1|Pondnum), FDistance, na.action = na.fail)
dredge(model, beta = "none") #This function runs all possible models and prints an AIC table

#The model that includes dompred with the lowest AIC doesn't include any of the environmental predictors. That's the one I'll proceed with.

#Model results, including checking variance inflation


FDislmer <- lmer(FDis ~ dompred + year + season + dompred:year + dompred:season  +  (1|Pondnum), FDistance)


#look for multicollinearity (VIF>10 is too much)
vif(FDislmer)

Anova(FDislmer,type=3)
plot(FDislmer)



#Post-hoc to generate emmeans

FDis.dompred.posthoc <-emmeans(FDislmer,pairwise~dompred)
FDis.dompred.emmeans <- as.data.frame(FDis.dompred.posthoc$emmeans)%>%
  dplyr::mutate("LowerSE" = emmean-SE,"UpperSE"=emmean+SE)




FDis.season.posthoc <- emmeans(FDislmer,pairwise~season*dompred)
FDis.season.emmeans <- as.data.frame(FDis.season.posthoc$emmeans)%>%
  dplyr::mutate("LowerSE" = emmean-SE,"UpperSE"=emmean+SE)

FDis.all.posthoc <- emmeans(FDislmer, pairwise ~dompred * season * year)
FDis.all.emmeans <- as.data.frame(FDis.all.posthoc$emmeans)%>%
  dplyr::mutate("LowerSE" = emmean-SE,"UpperSE"=emmean+SE)%>%
  filter(!(year == 2008 & season == "Su" ),
         !(year == 2008 & season == "W"))


#Save emmeans for graphing in a separate workbook
FDis.dompred.emmeans%>%
  mutate(season = "Mean")%>%
  rbind(FDis.season.emmeans)%>%
  write.csv("Data/FDis.emmeans.csv", row.names = F)


#Colorblind friendly palette for graphing
cbPalette <- c("#CC79A7", "#E69F00", "#009E73","#56B4E9")


#FDis pred mean  graphs
Alpha.dompred.plot <- FDis.dompred.emmeans%>%
  ggplot(aes(x = dompred,y = emmean,color=dompred)) +
  geom_point(aes(size=16), show.legend = F)+
  geom_linerange(aes(ymin=lower.CL, ymax= upper.CL),size=0.5, show.legend = F) +
  scale_color_manual(labels = c("Invertebrate","Salamander","Sunfish","Bass"),
                     values=cbPalette) +
  labs(x = "Predator", y = expression(paste("FDis (",alpha," diversity) +/- SE"))) +
  scale_x_discrete(labels = c("Invertebrate", "Salamander", "Sunfish", "Bass")) +
  theme_classic()
Alpha.dompred.plot + 
  theme(axis.title= element_text(size=20),
        axis.text.x = element_text(size=16, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size = 14))
ggsave("Figures/Alpha.dompred.plot.tiff")


#Next is graph of predator means separated by season. Pred on x-axis.

season.labs <- c("Winter","Spring", "Summer", "Fall")
names(season.labs) <- c("W", "Sp", "Su", "F")


#dompred x season, dompred on x-axis
Alpha.plot <- FDis.season.emmeans%>%
  ggplot(aes(x = dompred, y = emmean, color = season))+
  geom_point(position=position_dodge(width=0.3),size=5)+
  geom_line(aes(group = season), position=position_dodge(width=0.3),size=0.5) +
  geom_linerange(aes(ymin=lower.CL, ymax= upper.CL),
                 position=position_dodge(width=0.3),size=0.5) +
  theme_classic() +
  labs(x = "Top predator",
       y = expression(paste("FDis (",alpha," diversity)")),
       color = "Season")+
  #scale_color_manual(labels = c("W","Sp","Su","F"),
  #                           values=cbPalette) +
  theme(legend.key.height = unit(0.9,"cm"), 
        legend.key.width = unit(0.9,"cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size =10),
        legend.position = c(0.95,1),
        legend.justification=c(.5,1.05),
        axis.title= element_text(size=14),
        axis.text = element_text(size=12))

ggsave("Figures/Alpha_dompred_plot.tiff")

#Dompred x season, season on x-axis

#Here's an alternate FDis line graph by season with season on x axis

Alpha.season.plot <- FDis.season.emmeans%>%
  ggplot(aes(x = season,y = emmean,color=dompred)) +
  geom_point(position=position_dodge(width=0.3),size=7)+
  geom_line(aes(group=dompred),position=position_dodge(width=0.3),size=0.5)+
  geom_linerange(aes(ymin=lower.CL, ymax= upper.CL),
                 position=position_dodge(width=0.3),size=0.5) +
  theme_classic() +
  labs(x = "Season", 
       y = expression(paste("FDis (",alpha," diversity)")),
       color = "Top predator") +
  scale_color_manual(labels = c("Invertebrate","Salamander","Sunfish","Bass"),
                     values=cbPalette) +
  scale_x_discrete(labels=c("Winter","Spring","Summer","Fall"))
Alpha.season.plot + theme(legend.key.height = unit(0.9,"cm"), 
                          legend.key.width = unit(0.9,"cm"),
                          legend.text = element_text(size = 8),
                          legend.title = element_text(size =10),
                          legend.position = c(0.7,1),
                          legend.justification=c(0,1.05),
                          axis.title= element_text(size=14),
                          axis.text = element_text(size=12))

#ggsave("Figures/Alpha.season.plot.tiff") #This will replace the above graph

#And here's a line graph of all dompred x year x season combos

#For some reason the below graph has two extra NA points tacked on to the end of the x axis. Trying to diagnose it

#FDis line graph pred x year x season
Alpha_all_plot <-FDis.all.emmeans%>%
  unite(season_year, c("season","year"))%>%
  na.omit()%>%
  ggplot(aes(x=season_year, y=emmean, color=dompred))+
  geom_point(position=position_dodge(width=0.3),size=5)+
  geom_line(aes(group=dompred),position=position_dodge(width=0.3),size=0.5)+
  geom_linerange(aes(ymin=lower.CL, ymax= upper.CL),
                 position=position_dodge(width=0.3),size=0.5) +
  theme_classic() +
  labs(x = element_blank(),
       y = expression(paste("FDis (",alpha," diversity)")),
       color = "Top predator") +
  scale_color_manual(labels = c("Invertebrate","Salamander","Sunfish","Bass"),
                     values=cbPalette) +
  scale_x_discrete(labels=Timekey$season_year)+
  theme(#legend.key.height = unit(0.9,"cm"), 
        #legend.key.width = unit(0.9,"cm"),
        #legend.text = element_text(size = 8),
        #legend.title = element_text(size =10),
        #legend.position = c(0.7,1.05),
        #legend.justification=c(0,1.05),
        legend.title = element_blank(),
        axis.title= element_text(size=14),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1))


ggsave("Figures/Alpha_all_plot.tiff")





#Run the analysis with resampled data to make sure sample sizes aren't influencing our results

#FDis resampling
n <-3 #number of ponds to draw from each predator type
numsim <- 1000 #number of simulations to run
FDis.rand <- tibble(Pondnum = factor(),
                    dompred = factor(levels = c("S","G","B")),
                    FDis = numeric(),
                    year = factor(),
                    season = factor(levels = c("W","Sp","Su","F")),
                    sim = NA) #dataframe to insert the mean distance values from randomizations

#Run the randomization. 
for(i in 1:numsim){
  samponds <- DomPredata%>%
    filter(dompred != "O")%>%
    group_by(dompred)%>%
    sample_n(n)
  
  FDis.rand <-FDistance%>%
    dplyr::select(Pondnum, dompred, FDis, year, season)%>%
    filter(Pondnum %in% samponds$Pondnum)%>%
    dplyr::mutate(sim = i)%>%
    rbind(FDis.rand)
}



#Run if need to save the randomized values
write.csv(FDis.rand, "Data/FDis.rand.csv", row.names = F)

#Run if file already saved and need to load it
#FDis.rand <- read.csv("Data/FDis.rand")

#Graph distribution of mean values by dompred
FDis.rand%>%
  group_by(dompred, sim)%>%
  dplyr::summarise(FDis = mean(FDis))%>%
  ggplot(aes(x = dompred, y = FDis, color = dompred)) +
  geom_violin(aes(fill = dompred)) +
  theme_classic() +
  scale_color_manual(values=cbPalette) +
  scale_fill_manual(values = cbPalette) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  ggtitle("Resampled Alpha Diversity") +
  labs(x = "Predator", 
       y = "Mean FDis") 

#Graph distribution of mean values by dompred*season
FDis.rand%>%
  group_by(dompred, season, sim)%>%
  dplyr::summarise(FDis = mean(FDis))%>%
  ggplot(aes(x = dompred, y = FDis, color = dompred)) +
  geom_violin(aes(fill = dompred)) +
  theme_classic() +
  scale_color_manual(values=cbPalette) +
  scale_fill_manual(values = cbPalette) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Resampled alpha Dissimilarity") +
  labs(x = "Predator", 
       y = "Mean FDis") +
  facet_grid(.~season)

ggsave("Figures/Alpha_rand.tiff")

