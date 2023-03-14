
#This file calculates and analyzes temporal beta diversity metrics

#Load data
# Data set of all pairwise gower's distances for all year x season x pond combinations. This was created in the Spatial Beta code
library(tidyverse)
Allpairdist <- read.csv("Data/Allpairdist.csv", header = T)%>%
  mutate_at(vars(time1, time2, distance), as.numeric)%>%
  mutate_at(vars(Pondnum1, dompred1, Pondnum2, dompred2), as.factor)
#Key of the numerical ID's given to each season x year combination
Timekey <- read.csv("Data/Timekey.csv", header = T)%>%
  mutate(season = factor(season, levels = c("W", "Sp", "Su", "F")))%>%
  mutate_at(vars(year), as.factor)


#All pairwise distances between pond level trait averages were calculated in the spatial beta code. Here I select just the pairwise distances I want for the temporal beta diversity - comparisons between ponds and themselves in the next time step.

Temporalbeta <- Allpairdist%>%
  filter(Pondnum1 == Pondnum2, time1 == time2 + 1)%>% #Restrict to distances between ponds and themselves one timestep in future
  dplyr::select(Pondnum1, dompred1, time1, time2, distance)%>% #pondnum's and dompred's 1 vs 2 are duplicates in each pairwise distance
  dplyr::rename("Pondnum" = Pondnum1, "dompred" = dompred1)%>%
  mutate(dompred = factor(dompred, levels = c("N", "S", "G", "B")))%>%
  left_join(Timekey,by = c("time1" = "time"))
  

#Recode the season levels to reflect the fact they are similarities between two different seasons
Temporalbeta$season <- recode_factor(Temporalbeta$season, "W" = "W_Sp", "Sp" = "Sp_Su", "Su" = "Su_F", "F" = "F_W")


#Gamma glm. Add a small value to response variable to fit in gamma distribution
library(lme4)
Temporalbetagamma <- glmer(distance ~ dompred + year + season + dompred*year + dompred*season + (1|Pondnum), 
                           transform(Temporalbeta, distance = distance + 0.0001), family = "Gamma",
                           control = glmerControl(optimizer = "optimx", optCtrl = list(method= "bobyqa")))

library(car)
Anova(Temporalbetagamma, type=3)

#Calculate emmeans for graphing
Temporal.posthoc <- emmeans(Temporalbetagamma, pairwise ~ dompred*season, type="response")

Temporal.emmeans <- Temporal.posthoc$emmeans%>%
  data.frame()

#Uncomment to save emmeans for graphing in separate file
#write.csv(Temporal.emmeans, "Data/Temporal.emmeans.csv")




#Randomization test to see if different sample sizes is biasing results

#Upload predator data for randomly selecting ponds
DomPredata<-read.csv("Data/PredType2.csv", header=T)

n <-3 #number of ponds to draw from each predator type
numsim <- 1000 #number of simulations to run

#Create dataframe to insert the mean distance values from randomizations
tempbeta.rand <- tibble(Pondnum = factor(),
                        dompred = factor(levels = c("N", "S","G","B")),
                        time1 = numeric(),
                        time2 = numeric(),
                        distance = numeric(),
                        year = factor(),
                        season = factor(levels = c("W","Sp","Su","F")),
                        sim = NA) 


#For loop to perform the randomizations
for(i in 1:numsim){
  samponds <- DomPredata%>%
    filter(dompred != "O")%>%
    group_by(dompred)%>%
    sample_n(n)
  
  tempbeta.rand <-Temporalbeta%>%
    dplyr::select(-c(year))%>%
    filter(Pondnum %in% samponds$Pondnum)%>%
    dplyr::mutate(sim = i)%>%
    rbind(tempbeta.rand)
}


#Uncomment to save for graphing
#write.csv(tempbeta.rand, "Data/tempbeta.rand.csv", row.names = F)
