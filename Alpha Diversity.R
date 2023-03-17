
#Here I calculate the local (alpha) functional diversity of dragonfly communities in ponds with different top predators - invertebrates, salamanders, green sunfish, and largemouth bass. I'll try three 

#Load data

library(tidyverse)
Traitaverage <- read.csv("Data/Traitaverage.csv")%>%
  column_to_rownames("spID")#Assign the spID column to row names for the FD function

AbundancewideFD <- read.csv("Data/AbundancewideFD.csv")
Pondata <- read.csv("Data/Pondata.csv")%>%
  mutate_at(vars(Pondnum, depthcon, veg, wood, canopy, substrate, litter), as.factor)
Timekey <- read.csv("Data/Timekey.csv")
DomPredata <- read.csv("Data/PredType2.csv")



#Some quick data prep for functional diversity calculations. The function to calculate function diversity metrics requires a wide version of abundance data.
#Also, the spID in Traitaverage df need to be assigned to row names

AbundancewideFD <- AbundancewideFD%>%
  rowwise()%>%
  filter(sum(c_across(Ajun:Tramea))>=3)%>% #Remove samples with fewer than 3 observations to reduce biased averages.
  column_to_rownames("ID")


#Calculate FD metrics
library(FD)
FD<-dbFD(Traitaverage, AbundancewideFD, w.abun=T)

#Create a data set with the FDis values and predictor variables and format
FDistance <- FD%>%
  data.frame()%>%
  dplyr::select(FDis)%>% #Specify the functional diversity metric(s) to be analyzed
  rownames_to_column("ID")%>% #Move the row names back to columns
  separate(ID, c("dompred", "Pondnum", "year", "season"), convert=F)%>% #ID values contain multiple predictors. Separate them
  #mutate_at(vars(Pondnum), as.factor)%>%
  mutate_if(is.character, as.factor)%>%
  left_join(dplyr::select(Pondata, -dompred), by = "Pondnum")%>%
  mutate(season = fct_relevel(season, c("W", "Sp", "Su", "F")), #Reorder seasons for visualization
         dompred = fct_relevel(dompred, c("N", "S", "G", "B"))) #Reorder dompred for visualization



#Select appropriate model. I want to see if any of the environmental predictors improve the model fit.
library(lme4)
model <- lmer(FDis ~ dompred + year + season + dompred:year + dompred:season  +  
                            area + depth + veg + canopy + (1|Pondnum), FDistance, na.action = na.fail)
#Compare AICc values of all possible models
library(MuMIn)
dredge(model, beta = "none")

#The model that includes dompred with the lowest AIC doesn't include any of the environmental predictors. That's the one we'll proceed with.

#Run the chosen model

FDislmer <- lmer(FDis ~ dompred + year + season + dompred:year + dompred:season  +  (1|Pondnum), FDistance)

library(car)
Anova(FDislmer, type = 3)
#Post-hoc to generate estimated marginal means (emmeans)

#Season-specific predator averages
FDis.posthoc <- emmeans(FDislmer,pairwise~season*dompred)
FDis.emmeans <- as.data.frame(FDis.season.posthoc$emmeans)%>%
  dplyr::mutate("LowerSE" = emmean-SE,"UpperSE"=emmean+SE)


#Uncomment to save data for graphing in another r file
#write.csv(FDis.emmeans,"Data/FDis.emmeans.csv", row.names = F)


#Run the analysis with resampled data to make sure sample sizes aren't influencing our results

#FDis resampling
n <-3 #number of ponds to draw from each predator type
numsim <- 1000 #number of simulations to run
#Make a data frame for the randomized values to be inserted into
FDis.rand <- tibble(Pondnum = factor(),
                    dompred = factor(levels = c("N", "S","G","B")),
                    FDis = numeric(),
                    year = factor(),
                    season = factor(levels = c("W","Sp","Su","F")),
                    sim = NA) 

#Run the randomization. 
for(i in 1:numsim){
  samponds <- DomPredata%>%
    filter(dompred != "O")%>% #Remove ponds with "other" dominant predators
    group_by(dompred)%>%
    sample_n(n)
  
  FDis.rand <-FDistance%>%
    dplyr::select(Pondnum, dompred, FDis, year, season)%>%
    filter(Pondnum %in% samponds$Pondnum)%>%
    dplyr::mutate(sim = i)%>%
    rbind(FDis.rand)
}


#Save the randomized data set
#write.csv(FDis.rand, "Data/FDis.rand.csv", row.names = F)

