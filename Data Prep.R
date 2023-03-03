#title: "Top pred func div - data prep"



#Required packages and input raw data
library(tidyverse)
library(fitdistrplus)
library(cluster)


#Load the data sets
Abundata<-read.csv("Data/Abundances.csv",header=T) #Dragonfly Abundances
Abundancewide<-read.csv("Data/Abundances_Wide.csv",header=T)#Abundances wide version
Traitdata<-read.csv("Data/Traits.csv",header=T) #Traits
Envdata<-read.csv("Data/env.csv",header=T) #Pond environmental characteristics
#Predata<-read.csv("Data/PredatorData.csv",header=T) #Presence/absence of all predators. I don't think I actually use this
DomPredata<-read.csv("Data/PredType2.csv",header=T) #Dominant predator
Permdata<-read.csv("Data/perm.csv",header=T) #Pond permanance data
Taxonomy<-read.csv("Data/taxonomy.csv",header=T) #List of taxnomic shorthand used
Libtraits<-read.csv("Data/Libellula\ measurements.csv", header=T) #Body length and area of all F-0 Libellula specimens



#Combine environmental data, dominant predator, and pondnum into one dataset for later use
Pondata <- left_join(DomPredata, Envdata, by = "pond")

#write.csv(Pondata, "Data/Pondata.csv", row.names = F)


#calculate average trait values by spp

Traitaverage<-Traitdata %>%
  filter(Stage== "F-0" | Stage== "F-0?")%>% #Include only the F-0 stages
  group_by(latin)%>%
  dplyr::summarise_at(vars("Body_A":"Eye_A","Gape_W":"Mentum_L"),mean,na.rm=T)%>%
  left_join(Taxonomy,by="latin") #Attach species abbreviations 


##Make decision on Libellula data, since there are 3 spp in trait dataset but abundance data only resolved to genus
#Histograms of body length and area of all F-0 samples
ggplot(Libtraits, aes(Body.Length)) + geom_histogram(binwidth=0.5) +
  theme_classic()
ggplot(Libtraits, aes(Body.Length)) + geom_histogram(binwidth=0.5) +
  facet_grid(Season~.) + theme_classic()

ggplot(Libtraits, aes(Body.Area)) + geom_histogram(binwidth=2) +
  facet_grid(Year~.) + theme_classic()
ggplot(Libtraits, aes(Body.Area)) + geom_histogram(binwidth=2) +
  facet_grid(Season~.) + theme_classic()


#Assign Libullela incesta value of "Libullela"  There were three libullela spp in trait dataset, but abundance data only resolved to genus for libullela
Traitaverage[22,14] <- "Libellula"
Traitaverage[22,15] <- "Anisoptera"
Traitaverage[22,16] <- "Libellulidae"

#Gomphus lividus only had 1 sample from which the size measurement was taken.  Here I check to see how many samples contained this species and how many other species were included in those samples
Gomph<-Traitdata%>%
  filter(Stage== "F-0" | Stage== "F-0?")%>% #Include only the F-0 stages
  group_by(latin)%>%
  dplyr::summarise(count = n())%>%
  dplyr::rename(spID = latin)%>%
  right_join(Abundata, by = "spID")%>%
  filter(count < 3)%>%view()
  dplyr::select(year, season, pond)%>%
  left_join(Abundata, by = c("year", "season", "pond"))
Gomph%>%
  filter(spID == "Gomphus lividus")%>%
  group_by(year, season, pond)%>%
  dplyr::summarise(count = n()) 

Abundata%>%
  filter()
#Remove the other two species of Libellula and Gomphus lividus
Traitaverage <- filter(Traitaverage, latin != "Libellula vibrans" 
                       & latin != "Libellula auripennis"
                       & latin != "Gomphus lividus")%>%
  dplyr::select(-c(latin, suborder, family)) # Don't need these

#Save Traitaverage file with rows arranged alphabetically. I'll need this file for the dbFD function
#Traitaverage%>%
  #arrange(species)%>%
  #write.csv("Data/Traitaverage.csv",row.names = F)

##Add the pond environmental data
#join dominant predator and environment datasets
Pondata<-left_join(DomPredata, Envdata,by="pond")

#change "na" values in ponddata to true "NA"'s
Pondata$substrate[Pondata$substrate == "na"] <- NA
Pondata$litter[Pondata$litter == "na"] <- NA

#Remove "O" predator ponds from ponddata and save
Pondata%>%
  filter(dompred != "O")%>%
  write.csv("Data/Pondata.csv")
#uncount abundance data
Abundanceuncount<-uncount(Abundata,abundance)

#Add trait values to abundance data
Abundancetraits<-left_join(Abundanceuncount,Traitaverage,by="species")%>%
  drop_na()

#add pond data to abundance/traits data. Dropping the "other" predator ponds.
Finaldata<-left_join(Abundancetraits,Pondata,by="pond")%>%
  group_by(dompred)%>%
  filter(dompred != "O")%>%
  droplevels()



#write.csv(Finaldata,"Data/Finaldata.csv", row.names = F)

#Some more prep for calculating FDis, which requires a wide version of abundance data

#Abundance data need to be wide for FD
AbundancewideFD <- dplyr::select(Abundancewide,-c(Gomphus.lividus,Arigomphus.lentulus,Arigomphus.maxwelli,C_igens,Celithemis_Unidentified,CelithemisA,Celithemis.elisa,
                                                  Didymops_transversa,Ebasidans,Epitheca_princeps,Iramburii,Lestes_inaequalis,Lestes_vigilax,N_pentacantha,
                                                  Sympetrum_ambiguum))%>%
  left_join(DomPredata, by="pond")%>%
  filter(dompred != "O")%>%
  unite("ID", c("dompred", "Pondnum", "year", "season"), sep = ".")%>%
  dplyr::select(c(ID, Ajun:Libellula, -X))


#Sort by column names and save as csv
#AbundancewideFD%>%
  #dplyr::select(ID, sort(colnames(.)))%>%

#For a couple of the analyses I'll need each season x year combination to be sequentially assigned a number.

#I'll need to assign each year x season combo a unique number for a few things.

Timekey <- Finaldata%>%
  ungroup()%>%
  mutate(season = factor(season, levels = c("W", "Sp", "Su", "F")))%>%
  arrange(year, season)%>%
  dplyr::select(year, season)%>%
  unique()%>%
  mutate(time=1:n())%>%
  unite(season_year, c(season, year), remove = F)%>%
  mutate(year = factor(year))

#write.csv(Timekey, "Data/Timekey.csv", row.names = F)

#Next some data prep for beta diversity analyses

#Make a data set of the average trait values for each pond

#Spatial beta diversity
#Calculate average trait values for each sample and combine necessary variables to attach as row names for dissimilarity matrix
Pondtraitmeans <- Finaldata%>%
  mutate_at(vars(year, season), as.factor)%>%
  left_join(Timekey, by = c("year", "season"))%>%
  group_by(dompred, time, Pondnum)%>%
  mutate(count = n())%>%
  filter(count >= 3)%>% #Filter out any ponds with fewer than 3 dragonflies to avoid spurious averages
  dplyr::summarize_at(vars(Body_A:Mentum_L), mean)%>%
  unite("ID",c(dompred, time, Pondnum),sep=".",remove=T)


#write.csv(Pondtraitmeans, "Data/Pondtraitmeans.csv", row.names = F)

