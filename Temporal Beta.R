
library(tidyverse)
library(lme4)
library(emmeans)

# Data set of all pairwise gower's distances for all year x season x pond combinations. This was created in the Spatial Beta code
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


#Exploratory data visualization

Temporalbeta %>%
  ggplot(aes(x=dompred, y=distance)) +
  geom_beeswarm()+
  facet_grid(~year) +
  theme_classic()
Temporalbeta %>%
  ggplot(aes(x=dompred, y=distance)) +
  geom_beeswarm()+
  facet_grid(~season) +
  theme_classic()


#GLM

Temporalbetalmer <- lmer(distance ~ dompred + year + season + dompred*year + dompred*season +
                          (1|Pondnum), Temporalbeta)
plot(Temporalbetalmer)
qqnorm(resid(Temporalbetalmer))
qqline(resid(Temporalbetalmer))

#Got some heteroscedasticity.  Checking variance for each  variable.

plot(Temporalbeta$dompred, resid(Temporalbetalmer))
plot(Temporalbeta$season, resid(Temporalbetalmer))
plot(Temporalbeta$year, resid(Temporalbetalmer))

#Looks ok.  Try some alternate distributions.


library(fitdistrplus)
descdist(Temporalbeta$distance)

#Looks like it's somewhere between lognormal and gamma.  I'll go with gamma first since that's what I used for the spatial beta

Temporalbetagamma <- glmer(distance ~ dompred + year + season + dompred*year + dompred*season + (1|Pondnum), 
                           transform(Temporalbeta, distance = distance + 0.0001), family = "Gamma",
                           control = glmerControl(optimizer = "optimx", optCtrl = list(method= "bobyqa")))

plot(Temporalbetagamma)
qqnorm(resid(Temporalbetagamma))
qqline(resid(Temporalbetagamma))
Anova(Temporalbetagamma, type=3)

#That converged, and the variance is better.  Pred main effect and pred*season interaction significant.  Next I'll do posthoc.


Temporal.all.posthoc <- emmeans(object=Temporalbetagamma, pairwise~dompred*year*season,type = "response")
Temporal.posthoc <- emmeans(object=Temporalbetagamma, pairwise~dompred,type = "response")
Temporal.season.posthoc <- emmeans(Temporalbetagamma, pairwise ~ dompred*season, type="response")

Temporal.emmeans <- Temporal.posthoc$emmeans%>%
  data.frame()%>%
  mutate(season = "Mean")%>%
  rbind(data.frame(Temporal.season.posthoc$emmeans))
write.csv(Temporal.emmeans, "Data/Temporal.emmeans.csv")



#Graphs


#Colorblind friendly palette for graphing
cbPalette <- c("#CC79A7", "#E69F00", "#009E73","#56B4E9")

#Graph of pred mean temporal beta
Temporal_dompred_plot <- as.data.frame(Temporal.posthoc$emmeans) %>%
  dplyr::rename(distance = response)%>%
  ggplot(aes(x = dompred, y= distance, color=dompred)) + 
  geom_point(size = 9, show.legend = F) +
  geom_linerange(aes(ymin=asymp.LCL, ymax = asymp.UCL), size = 0.5, show.legend=F) +
  scale_color_manual(values=cbPalette) +
  labs(x = "Predator", 
       y = expression(paste("Temporal Dissimilarity (",beta," diversity)"))) +
  scale_x_discrete(labels = c("Invertebrate", "Salamander", "Sunfish", "Bass")) +
  theme_classic()

#ggsave("Figures/Temporal_dompred_plot.tiff")

#Dompred on x axis, grouped by season
Temporal.dompred.plot <- Temporal.season.posthoc$emmeans%>%
  data.frame()%>%
  ggplot(aes(x = dompred, y = response, color = season)) +
  geom_point(position = position_dodge(width = 0.3), size = 7) +
  geom_line(aes(group = season), position = position_dodge(width = 0.3), size = 0.5) +
  geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                 position = position_dodge(width = 0.3), size = 0.5) +
  theme_classic() +
  labs(x = "Top Predator", 
              y = expression(paste("Temporal Dissimilarity (",beta," diversity)")),
              color = "Season") +
         #scale_color_manual(labels = c("W_Sp","Sp_Su","Su_F","F_W"),
                              #values=cbPalette) +
         scale_x_discrete(labels=c("N","S","G","B")) +
         theme(#legend.key.height = unit(0.9,"cm"), 
               #legend.key.width = unit(0.9,"cm"),
               #legend.text = element_blank(),
               legend.title = element_blank(),
               #legend.position = "top",,
               #legend.position = c(0.7,0.95),
               #legend.justification=c(0,0.9),
               axis.title= element_text(size=18),
               axis.text = element_text(size=16))
#ggsave("Figures/Temporal.dompred.plot.tiff")


#Season on x axis, grouped by dompred

Temporal_season_plot <-as.data.frame(Temporal.season.posthoc$emmeans) %>%
  ggplot(aes(x=season, y=response, color=dompred))+
         geom_point(position=position_dodge(width=0.3),size=7)+
         geom_line(aes(group=dompred),position=position_dodge(width=0.3),size=0.5)+
         geom_linerange(aes(ymin=asymp.LCL, ymax= asymp.UCL),
                 position=position_dodge(width=0.3),size=0.5) +
         theme_classic() +
         labs(x = "Season", 
              y = expression(paste("Temporal Dissimilarity (",beta," diversity)")),
              color = "Top predator") +
         scale_color_manual(labels = c("Invertebrate","Salamander","Sunfish","Bass"),
                              values=cbPalette) +
         scale_x_discrete(labels=c("Winter","Spring","Summer","Fall")) 
Temporal_season_plot +theme(legend.key.height = unit(0.9,"cm"), 
               legend.key.width = unit(0.9,"cm"),
               legend.text = element_text(size = 12),
               legend.title = element_text(size =14),
               legend.position = c(0.7,0.95),
               legend.justification=c(0,0.9),
               axis.title= element_text(size=18),
               axis.text = element_text(size=16))
#ggsave("Figures/Temporal_season_plot.tiff")

#Graph of all year x season combos

Temporal.all.plot <-as.data.frame(Temporal.all.posthoc$emmeans) %>%
  filter(!(year == 2008 & season == "W"))%>% #Remove emmeans for year x season that weren't sampled
filter(!(year == 2008 & season == "Sp"))%>%
  filter(!(year == 2011 & season == "F"))%>%
  left_join(Timekey, by=c("year","season"),na.rm=T)%>%
  ggplot(aes(x=season_year, y=response, color=dompred))+
  geom_point(position=position_dodge(width=0.3),size=7)+
  geom_line(aes(group=dompred),position=position_dodge(width=0.3),size=0.5)+
  geom_linerange(aes(ymin=asymp.LCL, ymax= asymp.UCL),
                 position=position_dodge(width=0.3),size=0.5) +
  theme_classic() +
  labs(x = "Season", 
       y = expression(paste("Temporal Dissimilarity (",beta," diversity)")),
       color = "Top predator") +
  scale_color_manual(labels = c("Invertebrate","Salamander","Sunfish","Bass"),
                     values=cbPalette) +
  scale_x_discrete(labels=Timekey$season_year) +
  theme(legend.key.height = unit(0.9,"cm"), 
        legend.key.width = unit(0.9,"cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size =10),
        legend.position = c(1,1.05),
        legend.justification=c(0,1.05),
        axis.title= element_text(size=14),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = unit(c(1,3,1,1),"cm")) 

ggsave("Figures/Temporal.all.plot.tiff", width=8) 





#Randomization test to see if different sample sizes is biasing results

#Upload 
DomPredata<-read.csv("Data/PredType2.csv", header=T) #Data set of the top predators in each pond

n <-3 #number of ponds to draw from each predator type
numsim <- 1000 #number of simulations to run
tempbeta.rand <- tibble(Pondnum = factor(),
                        dompred = factor(levels = c("N", "S","G","B")),
                        time1 = numeric(),
                        time2 = numeric(),
                        distance = numeric(),
                        year = factor(),
                        season = factor(levels = c("W","Sp","Su","F")),
                        sim = NA) #dataframe to insert the mean distance values from randomizations


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



#write.csv(tempbeta.rand, "Data/tempbeta.rand.csv", row.names = F)

#Dompred main effect
tempbeta.rand%>%
  group_by(dompred, sim)%>%
  dplyr::summarise(mean.dist = mean(distance))%>%
  ggplot(aes(x = dompred, y = mean.dist, color = dompred)) +
  geom_violin(aes(fill = dompred)) +
  theme_classic() +
  scale_color_manual(values=cbPalette) +
  scale_fill_manual(values = cbPalette) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  ggtitle("Randomized Temporal Dissimilarity") +
  labs(x = "Predator", 
       y = "Mean Temporal Dissimilarity") 

#ggsave("Figures/temp_rand.tiff")

#Dompred x season interaction

tempbeta.rand%>%
  group_by(dompred, season, sim)%>%
  dplyr::summarise(mean.dist = mean(distance))%>%
  ggplot(aes(x = dompred, y = mean.dist, color = dompred)) +
  geom_violin(aes(fill = dompred)) +
  theme_classic() +
  scale_color_manual(values=cbPalette) +
  scale_fill_manual(values = cbPalette) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Randomized Spatial Dissimilarity") +
  labs(x = "Predator", 
       y = "Mean Spatial Dissimilarity") +
  facet_grid(.~season)

#ggsave("Figures/temp_rand.tiff")


