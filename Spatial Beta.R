library(tidyverse)
library(cluster)
library(reshape2)
library(ggbeeswarm)
library(lme4)
library(fitdistrplus)
library(optimx)
Pondtraitmeans <- read.csv("Data/Pondtraitmeans.csv")
Timekey <- read.csv("Data/Timekey.csv")%>%
  mutate_at(vars(year, season, time), as.factor)%>%
  mutate(season = factor(season, levels = c("W", "Sp", "Su", "F"))) #Reorder for graphing


#Calculate a distance matrix of the pond level average trait values and convert to a usable data frame

Distmat <- Pondtraitmeans%>%
  column_to_rownames("ID")%>% #Daisy function requires nxp matrix.  Identity values moved to row names to be preserved
  daisy(metric="gower",stand=T)

#Convert the distance matrix to a data frame.  Need to distinguish between the time values of the two samples used for each pariwise combo
Allpairdist <- Distmat%>%
  as.matrix()%>%
  melt()%>%
  separate(Var1, c("dompred1","time1","Pondnum1"))%>%
  separate(Var2, c("dompred2","time2","Pondnum2"))%>%
  filter(Pondnum1 >= Pondnum2)%>% #This gets rid of one  of the triangle
  dplyr::rename("distance" = value)%>%
  mutate_at(vars(dompred1, dompred2, time1, time2, Pondnum1, Pondnum2), as.factor)%>%
  mutate(dompred2 = factor(dompred2, levels = c("N", "S", "G", "B")))%>% #reorder predator levels for graphing later
  mutate(dompred1 = factor(dompred1, levels = c("N", "S", "G", "B")))

#write.csv(Allpairdist, "Data/Allpairdist.csv", row.names = F) 

#From the all-pairwise dataset, the comparisons between ponds of the same predator type, year, and season were selected for the spatial beta diversity analysis


Spatialbeta <- Allpairdist %>%
  filter(time2 == time1, dompred1 == dompred2, Pondnum1 != Pondnum2)%>%
  dplyr::rename("dompred" = dompred1, "time" = time1)%>%
  left_join(Timekey, by =  "time")%>%
  mutate_at(vars(year), as.factor)%>%
  dplyr::select(-c("dompred2","time2"))%>%
  mutate(season = factor(season, levels = c("W", "Sp", "Su", "F"))) #Reorder season levels for graphing
#Some exploratory data visualization.  

Spatialbeta %>%
  ggplot(aes(x = dompred, y = distance)) +
  geom_beeswarm() +
  facet_wrap(~season) +
  theme_classic()
Spatialbeta %>%
  ggplot(aes(x = dompred, y = distance)) +
  geom_beeswarm() +
  facet_wrap(~year:season) +
  theme_classic()

#First I conducted a general linear model.  

Spatialbetalmer <- lmer(distance~ dompred + year + season + dompred:year + dompred:season + 
                          (1|Pondnum1) + (1|Pondnum2), Spatialbeta)

#Assumption checking demonstrated deviations from normality and heteroscedasicity.

plot(Spatialbetalmer)
qqnorm(resid(Spatialbetalmer))
qqline(resid(Spatialbetalmer))

#Look at the variance by factor

plot(Spatialbeta$dompred, resid(Spatialbetalmer))
plot(Spatialbeta$year, resid(Spatialbetalmer))
plot(Spatialbeta$season, resid(Spatialbetalmer))

#The variances by factor don't look too bad.



#I'll try fitdistrplus to see what the distribution is closest to.


descdist(Spatialbeta$distance)


#Looks like beta or gamma.  Let's try beta first.


SpatialbetaBetareg <- betareg(distance ~ dompred + year + season + dompred:year + dompred:season +
                                (1|pondnum1) + (1|pondnum2),
                              transform(Spatialbeta,distance = distance + 0.0001))
plot(SpatialbetaBetareg)



#Got an error that contrasts can only be applied to factors with 2 or more levels.  Here's a function I found online that can possibly shed light
#Couldn't get this to work
#DO NOT RUN

debug_contr_error <- function (dat, subset_vec = NULL) {
  if (!is.null(subset_vec)) {
    ## step 0
    if (mode(subset_vec) == "logical") {
      if (length(subset_vec) != nrow(dat)) {
        stop("'logical' `subset_vec` provided but length does not match `nrow(dat)`")
      }
      subset_log_vec <- subset_vec
    } else if (mode(subset_vec) == "numeric") {
      ## check range
      ran <- range(subset_vec)
      if (ran[1] < 1 || ran[2] > nrow(dat)) {
        stop("'numeric' `subset_vec` provided but values are out of bound")
      } else {
        subset_log_vec <- logical(nrow(dat))
        subset_log_vec[as.integer(subset_vec)] <- TRUE
      } 
    } else {
      stop("`subset_vec` must be either 'logical' or 'numeric'")
    }
    dat <- base::subset(dat, subset = subset_log_vec)
  } else {
    ## step 1
    dat <- stats::na.omit(dat)
  }
  if (nrow(dat) == 0L) warning("no complete cases")
  ## step 2
  var_mode <- sapply(dat, mode)
  if (any(var_mode %in% c("complex", "raw"))) stop("complex or raw not allowed!")
  var_class <- sapply(dat, class)
  if (any(var_mode[var_class == "AsIs"] %in% c("logical", "character"))) {
    stop("matrix variables with 'AsIs' class must be 'numeric'")
  }
  ind1 <- which(var_mode %in% c("logical", "character"))
  dat[ind1] <- lapply(dat[ind1], as.factor)
  ## step 3
  fctr <- which(sapply(dat, is.factor))
  if (length(fctr) == 0L) warning("no factor variables to summary")
  ind2 <- if (length(ind1) > 0L) fctr[-ind1] else fctr
  dat[ind2] <- lapply(dat[ind2], base::droplevels.factor)
  ## step 4
  lev <- lapply(dat[fctr], base::levels.default)
  nl <- lengths(lev)
  ## return
  list(nlevels = nl, levels = lev)
}
debug_contr_error(Spatialbeta)



#Everything has more than 2 levels.  Don't know what to do.


#Here's the glm with gamma function.  Since the response has zero's, I'll add a small value to them.

Spatialbetagamma <- glmer(distance ~ dompred + season + year + dompred:year + dompred:season + 
                          (1|Pondnum1) + (1|Pondnum2), 
                          transform(Spatialbeta, distance = distance + 0.0001), family = "Gamma",
                          control = glmerControl(optimizer = "optimx", optCtrl = list(method= "bobyqa")))

plot(Spatialbetagamma)
qqnorm(resid(Spatialbetagamma))
qqline(resid(Spatialbetagamma))
Anova(Spatialbetagamma, type = 3)



#Fit looks much better

#Generate the estimated marginal means
library(emmeans)
#Spatial.dompred.posthoc <- emmeans(Spatialbetagamma, pairwise ~ dompred, type = "response")
#Spatial.year.posthoc <- emmeans(Spatialbetagamma, pairwise ~ dompred|year, type="response")
Spatial.season.posthoc <- emmeans(Spatialbetagamma, pairwise ~ dompred|season,type="response")  

#Spatial.all.posthoc <- emmeans(Spatialbetagamma, pairwise ~ dompred|year|season, type = "response")

Spatial.emmeans <- Spatial.season.posthoc$emmeans%>%
  as.data.frame()

write.csv(as.data.frame(Spatial.emmeans), "Data/Spatial.emmeans.csv", row.names = F)


#Ok I got the emmeans and contrasts.  Now to graph it.  Not too sure how best to visualize it. I'll start with a figure with pred on x axis grouped by season

#colorblind friendly palet 
cbPalette <- c("#CC79A7", "#E69F00", "#009E73","#56B4E9")

Spatial.dompred.plot <- Spatial.emmeans%>%
  ggplot(aes(x = dompred, y = response, color = season)) +
  geom_point(position = position_dodge(width = 0.3), size = 7) +
  geom_line(aes(group = season), position = position_dodge(width = 0.3), size = 0.5) +
  geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                 position = position_dodge(width = 0.3), size = 0.5) +
  theme_classic() +
  labs(x = "Top Predator", 
              y = expression(paste("Spatial Dissimilarity (",beta," diversity)")),
              color = "Season") +
         #scale_color_manual(labels = c("W","Sp","Su","F"),
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

saveRDS(Spatial.dompred.plot, "Figures/Spatial.dompred.plot.rds")


#Group by season
Spatial_season_plot <-as.data.frame(Spatial.posthoc.season$emmeans) %>%
  ggplot(aes(x=season, y=response, color=dompred))+
  geom_point(position=position_dodge(width=0.3),size=7)+
  geom_line(aes(group=dompred),position=position_dodge(width=0.3),size=0.5)+
  geom_linerange(aes(ymin=asymp.LCL, ymax= asymp.UCL),
                 position=position_dodge(width=0.3),size=0.5) +
  theme_classic() +
  labs(x = "Season", 
       y = expression(paste("Spatial Dissimilarity (",beta," diversity)")),
       color = "Top predator") +
  scale_color_manual(labels = c("Invertebrate","Salamander","Sunfish","Bass"),
                     values=cbPalette) +
  scale_x_discrete(labels=c("Winter","Spring","Summer","Fall")) +
  theme(#legend.key.height = unit(0.9,"cm"), 
    #legend.key.width = unit(0.9,"cm"),
    #legend.text = element_blank(),
    legend.title = element_blank(),
    #legend.position = "top",,
    #legend.position = c(0.7,0.95),
    #legend.justification=c(0,0.9),
    axis.title= element_text(size=18),
    axis.text = element_text(size=16))
ggsave("Figures/Spatial_season_plot.tiff")
#Now year

#Group by year
Spatial_year_plot<-as.data.frame(Spatial.posthoc.year$emmeans) %>%
  ggplot(aes(x=year, y=response, color=dompred))+
  geom_point(position=position_dodge(width=0.3),size=7)+
  geom_line(aes(group=dompred),position=position_dodge(width=0.3),size=0.5)+
  geom_linerange(aes(ymin=asymp.LCL, ymax= asymp.UCL),
                 position=position_dodge(width=0.3),size=0.5) +
  theme_classic() +
  labs(x = "Year", 
       y = expression(paste("Spatial Dissimilarity (",beta," diversity)")),
       color = "Top predator") +
  scale_color_manual(labels = c("Invertebrate","Salamander","Sunfish","Bass"),
                     values=cbPalette) +
  theme(legend.key.height = unit(0.9,"cm"), 
        legend.key.width = unit(0.9,"cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size =14),
        legend.position = c(0.7,0.95),
        legend.justification=c(0,0.9),
        axis.title= element_text(size=18),
        axis.text = element_text(size=16))
ggsave("Figures/Spatial_year_plot.tiff")


#Plot of every year x season combo

Spatial_all_plot <-Spatial.posthoc.all$emmeans%>%
  as.data.frame()%>%
  left_join(Timekey, by=c("year","season"),na.rm=T)%>%
  na.omit()%>%
  ggplot(aes(x=as.factor(time), y=response, color=dompred))+
  geom_point(position=position_dodge(width=0.3),size=5)+
  geom_line(aes(group=dompred),position=position_dodge(width=0.3),size=0.5)+
  geom_linerange(aes(ymin=asymp.LCL, ymax= asymp.UCL),
                 position=position_dodge(width=0.3),size=0.5) +
  theme_classic() +
  labs(x = element_blank(),
       y = expression(paste("Spatial Dissimilarity (",beta," diversity)")),
       color = "Top predator") +
  scale_color_manual(labels = c("Invertebrate","Salamander","Sunfish","Bass"),
                     values=cbPalette) +
  scale_x_discrete(labels=Timekey$season_year)+
  theme(legend.key.height = unit(0.9,"cm"), 
        legend.key.width = unit(0.9,"cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size =10),
        legend.position = c(0.7,1.05),
        legend.justification=c(0,1.05),
        axis.title= element_text(size=14),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Figures/Spatial_all_plot.tiff")






#Randomization test to see if the results are influenced by the difference in sample sizes of the predator pond types.
#Import dataset of top predators in ponds
DomPredata<-read.csv("Data/PredType2.csv",header=T)%>%
  mutate_at(vars(dompred, Pondnum), as.factor)


n <-5 #number of ponds to draw from each predator type
numsim <- 500 #number of simulations to run
spbeta.rand <- tibble(pondnum1 = factor(),
                      dompred = factor(levels = c("S","G","B")),
                      pondnum2 = factor(),
                      distance = numeric(),
                      year = factor(),
                      season = factor(levels = c("W","Sp","Su","F")),
                      sim = NA) #dataframe to insert the mean distance values from randomizations

#Run the simulation. I'm removing all the invertebrate ponds because I'm drawing 5 ponds from each predator type and there are only 5 invert ponds. I'll substitue the observed values for inverts in the figure.

for(i in 1:numsim){
  samponds <- DomPredata%>%
    filter(dompred != "O", dompred != "N")%>% 
    group_by(dompred)%>%
    sample_n(n)
  
  spbeta.rand <-Spatialbeta%>%
    filter(Pondnum1 %in% samponds$Pondnum & Pondnum2 %in% samponds$Pondnum)%>%
    dplyr::mutate(sim = i)%>%
    rbind(spbeta.rand)
  
}

#Calculate average values and CI's for predator x season combinations
spbeta.rand.avg <- spbeta.rand%>%
  group_by(dompred, season)%>%
  summarise(avgdist = mean(distance), 
            asymp.LCL = quantile(distance, probs = 0.025),
            asymp.UCL = quantile(distance, probs = 0.975))





#Save the data
write.csv(spbeta.rand.avg, "Data/spbeta.rand.avg.csv", row.names = F)

#Graph distribution of mean values by dompred
spbeta.rand%>%
  group_by(dompred, sim)%>%
  dplyr::summarise(mean.dist = mean(distance))%>%
  ggplot(aes(x = dompred, y = mean.dist, color = dompred)) +
  geom_violin(aes(fill = dompred)) +
  theme_classic() +
  scale_color_manual(values=cbPalette) +
  scale_fill_manual(values = cbPalette) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  ggtitle("Randomized Spatial Dissimilarity") +
  labs(x = "Predator", 
       y = "Mean Spatial Dissimilarity") 

#Graph distribution of mean values by dompred*season
spbeta.rand%>%
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

ggsave("Figures/Spatial_rand.tiff")
#Graph distribution of mean values by dompred*year
spbeta.rand%>%
  group_by(dompred, year, sim)%>%
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
  facet_grid(.~year)






