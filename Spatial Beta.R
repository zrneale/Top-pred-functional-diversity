#This file calculates the spatial beta diversity metrics

#Load data
library(tidyverse)
Pondtraitmeans <- read.csv("Data/Pondtraitmeans.csv")
Timekey <- read.csv("Data/Timekey.csv")%>%
  mutate_at(vars(year, season, time), as.factor)%>%
  mutate(season = factor(season, levels = c("W", "Sp", "Su", "F"))) #Reorder for graphing


#Calculate a distance matrix of the pond level average trait values and convert to a usable data frame
library(cluster)
Distmat <- Pondtraitmeans%>%
  column_to_rownames("ID")%>% #Daisy function requires nxp matrix.  Identity values moved to row names to be preserved
  daisy(metric="gower",stand=T)

#Convert the distance matrix to a data frame.  Need to distinguish between the time values of the two samples used for each pariwise combo
#Var1 and Var2 are the united predator, time ID, and pond ID for the two ponds compared. These will be separated to individual columns
library(reshape2)
Allpairdist <- Distmat%>%
  as.matrix()%>% 
  melt()%>% #Convert to a data frame
  separate(Var1, c("dompred1","time1","Pondnum1"))%>%
  separate(Var2, c("dompred2","time2","Pondnum2"))%>%
  filter(Pondnum1 >= Pondnum2)%>% #This gets rid of duplicates
  dplyr::rename("distance" = value)%>%
  mutate_at(vars(dompred1, dompred2, time1, time2, Pondnum1, Pondnum2), as.factor)%>%
  mutate(dompred2 = factor(dompred2, levels = c("N", "S", "G", "B")))%>% #reorder predator levels for graphing later
  mutate(dompred1 = factor(dompred1, levels = c("N", "S", "G", "B")))

#Uncomment to save
#write.csv(Allpairdist, "Data/Allpairdist.csv", row.names = F) 


#From the all-pairwise dataset, the comparisons between ponds of the same predator type, year, and season are selected here for the spatial beta diversity analysis

Spatialbeta <- Allpairdist %>%
  filter(time2 == time1, dompred1 == dompred2, Pondnum1 != Pondnum2)%>%
  dplyr::rename("dompred" = dompred1, "time" = time1)%>%
  left_join(Timekey, by =  "time")%>%
  mutate_at(vars(year), as.factor)%>%
  dplyr::select(-c("dompred2","time2"))%>%
  mutate(season = factor(season, levels = c("W", "Sp", "Su", "F"))) #Reorder season levels for graphing




#Gamma glm. Adding a small value to the response variable so there are no zeroes to fit the distribution.

library(lme4)
library(optimx)
Spatialbetagamma <- glmer(distance ~ dompred + season + year + dompred:year + dompred:season + 
                          (1|Pondnum1) + (1|Pondnum2), 
                          transform(Spatialbeta, distance = distance + 0.0001), family = "Gamma",
                          control = glmerControl(optimizer = "optimx", optCtrl = list(method= "bobyqa")))

library(car)
Anova(Spatialbetagamma, type = 3)


#Generate the estimated marginal means
library(emmeans)
Spatial.posthoc <- emmeans(Spatialbetagamma, pairwise ~ dompred|season,type="response")  
Spatial.emmeans <- Spatial.posthoc$emmeans%>%
  as.data.frame()

#Uncomment to save
#write.csv(Spatial.emmeans, "Data/Spatial.emmeans.csv", row.names = F)



#Randomization test to see if the results are influenced by the difference in sample sizes of the predator pond types.
#Import dataset of top predators in ponds
DomPredata<-read.csv("Data/PredType2.csv",header=T)%>%
  mutate_at(vars(dompred, Pondnum), as.factor)


n <-3 #number of ponds to draw from each predator type
numsim <- 1000 #number of simulations to run
spbeta.rand <- tibble(pondnum1 = factor(),
                      dompred = factor(levels = c("N", "S","G","B")),
                      pondnum2 = factor(),
                      distance = numeric(),
                      year = factor(),
                      season = factor(levels = c("W","Sp","Su","F")),
                      sim = NA) #dataframe to insert the mean distance values from randomizations

#Run the simulation. I'm removing all the invertebrate ponds because I'm drawing 5 ponds from each predator type and there are only 5 invert ponds. I'll substitue the observed values for inverts in the figure.

for(i in 1:numsim){
  samponds <- DomPredata%>%
    filter(dompred != "O")%>% #Don't include the ponds with "other" top predators
    group_by(dompred)%>%
    sample_n(n)
  
  spbeta.rand <-Spatialbeta%>%
    filter(Pondnum1 %in% samponds$Pondnum & Pondnum2 %in% samponds$Pondnum)%>%
    dplyr::mutate(sim = i)%>%
    rbind(spbeta.rand)
  
}

#Uncomment to save the data
#write.csv(spbeta.rand, "Data/spbeta.rand.csv", row.names = F)
