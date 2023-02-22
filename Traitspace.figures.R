
library(tidyverse)
library(ape)
library(cluster)
library(ggrepel)
Finaldata <- read.csv("Data/Finaldata.csv")%>%
  mutate_at(vars(year), as.factor)%>%
  mutate(season = factor(season, levels = c("W", "Sp", "Su", "F")))%>% #Reorder for graphing
  mutate(dompred = factor(dompred, levels = c("N", "S", "G", "B")))
Traitaverage <- read.csv("Data/Traitaverage.csv")%>%
  column_to_rownames("species")

                     
#Conduct a PCoA on the species traits. PCoA was chosen because it's the method used in calculating FDis


##Need a distance matrix first
speciesdist <- daisy(Traitaverage, metric = "gower", stand = T)
PCOA <- pcoa(speciesdist, correction = "none")
biplot(PCOA, plot.axes = c(1,2))

##Attach the appropriate PCoA axis values to each of the observations in the final abundance data
Abundanceordi <- PCOA$vectors[,1:2]%>%
  as.data.frame()%>%
  rownames_to_column("species")%>%
  right_join(Finaldata, by = "species")


##Calculate the pond averages on the two axes
Averageordi <- Abundanceordi%>%
  group_by(dompred, year, season, Pondnum)%>%
  dplyr::summarise(count = n(), meanAxis1 = mean(Axis.1), meanAxis2 = mean(Axis.2),
                   SDAxis1 = sd(Axis.1), SDAxis2 = sd(Axis.2))%>%
  filter(count >= 3)%>%#Remove samples with fewer than 3 observations
  mutate_at(vars(Pondnum, year, dompred, season), as.factor)

##Now calculate the average of the pond averages within pred x year x season combos
Averageordi2 <- Averageordi%>%
  group_by(dompred, year, season)%>%
  dplyr::summarise(meanAxis1 = mean(meanAxis1), meanAxis2 = mean(meanAxis2), count = n())


##Extract and plot vector loadings using function from https://gist.github.com/Rekyt/ee15330639f8719d87aebdb8a5b095d4

compute_arrows = function(given_pcoa, trait_df) {
  
  # Keep only quantitative or ordinal variables
  # /!\ Change this line for different dataset
  #     or select only quantitative/ordinal var. /!\
  #trait_df = trait_df[, c(4:6, 20, 21)]
  
  n <- nrow(trait_df)
  points.stand <- scale(given_pcoa$vectors)
  
  # Compute covariance of variables with all axes
  S <- cov(trait_df, points.stand)
  
  
  # Select only positive eigenvalues
  pos_eigen = given_pcoa$values$Eigenvalues[seq(ncol(S))]
  
  # Standardize value of covariance (see Legendre & Legendre 1998)
  U <- S %*% diag((pos_eigen/(n - 1))^(-0.5))
  colnames(U) <- colnames(given_pcoa$vectors)
  
  # Add values of covariances inside object
  given_pcoa$U <- U
  
  return(given_pcoa)
}

###Run the function.
trait_loadings <- compute_arrows(PCOA, data.frame(apply(Traitaverage, 2, scale, center=TRUE, scale=TRUE)))

###Isolate the loadings in a data set and make the rownames (the traits) a column
arrows_df = as.data.frame(trait_loadings$U)%>%
  rownames_to_column("variable")

#Make plots First I'll create the PCOA loadings to be inserted as insets, then I'll make the traitspace figures

#colorblind friendly palette when faceting by season 
seasonPalette <- c('#EE6677', '#228833', '#4477AA', '#CCBB44', "grey35")
#Colorblind friendly palette wheb faceting by predator
predPalette <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00")

#Loading arrows plot
Loadings <- arrows_df%>%
  ggplot(aes(x=0, y =0, label = variable))+
  geom_segment(xend = arrows_df$Axis.1, yend = arrows_df$Axis.2, size = 0.1)+
  xlim(-5, 0) +
  ylim(-8, 10) +
  theme_classic() +
  geom_text_repel(data = arrows_df, inherit.aes = F, aes(x = Axis.1, y = Axis.2, label = variable),
                  size = 2.5)+
  labs(x = paste0("PCOA Axis 1 (", round(trait_loadings$values$Relative_eig[1] * 100, 2), "%)"),
       y = paste0("PCOA Axis 2 (", round(trait_loadings$values$Relative_eig[2] * 100, 2), "%)")) + 
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))





#Facet by dompred
Dompred.traitspace <- Averageordi2%>%
  ggplot(aes(x=meanAxis1,y=meanAxis2))+
  layer(data = filter(Averageordi2, dompred == "G"), stat = StatIdentity, position = PositionIdentity, geom = ggplot2::GeomCustomAnn,
        params = list(grob = ggplotGrob(Loadings), 
                      xmin = -0.14, xmax = -.025,
                      ymin = -.018, ymax = 0.04)) +
  geom_point(aes(color = season, shape=year), size = 2)+
  geom_polygon(aes(color = season, fill = season),alpha=0.1)+
  facet_wrap(~dompred, labeller =  as_labeller(c("N" = "Invertebrate", "S" = "Salamander", "G" = "Sunfish", "B" = "Bass"))) + 
  scale_color_manual(values=seasonPalette, labels = c("Winter","Spring","Summer","Fall")) +
  scale_fill_manual(values = seasonPalette, guide = NULL) +
  labs(color = "Season", shape = "Year") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 18)) +
  labs(x = paste0("PCOA Axis 1 (", round(trait_loadings$values$Relative_eig[1] * 100, 2), "%)"),
       y = paste0("PCOA Axis 2 (", round(trait_loadings$values$Relative_eig[2] * 100, 2), "%)")) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2))


#Save
ggsave("Figures/Traitspace.dompred.jpg", width = 15.32, height = 9.34)

#Facet by season
Season.traitspace <- Averageordi2%>%
  ggplot(aes(x=meanAxis1,y=meanAxis2))+
  layer(data = filter(Averageordi2, season == "Su"), stat = StatIdentity, position = PositionIdentity, geom = ggplot2::GeomCustomAnn,
        params = list(grob = ggplotGrob(Loadings), 
                      xmin = -0.14, xmax = -.025,
                      ymin = -.018, ymax = 0.04)) +
  geom_point(aes(color = dompred, shape=year))+
  geom_polygon(aes(color = dompred, fill = dompred),alpha=0.1)+
  facet_wrap(~season, labeller = as_labeller(c("W" = "Winter", "Sp" = "Spring", "Su" = "Summer", "F" = "Fall"))) + 
  scale_color_manual(values=predPalette, labels = c("Invertebrate", "Salamander", "Sunfish", "Bass")) +
  scale_fill_manual(values = predPalette, guide = NULL) +
  labs(color = "Top Predator", shape = "Year") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 18)) +
  labs(x = paste0("PCOA Axis 1 (", round(trait_loadings$values$Relative_eig[1] * 100, 2), "%)"),
       y = paste0("PCOA Axis 2 (", round(trait_loadings$values$Relative_eig[2] * 100, 2), "%)")) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2))

#Save
ggsave("Figures/Traitspace.season.jpg", width = 15.32, height = 9.34)


#Facet by year
Averageordi2%>%
  ggplot(aes(x=meanAxis1,y=meanAxis2))+
  geom_point(aes(color = dompred, shape=season))+
  geom_polygon(aes(color = dompred),alpha=0.1)+
  facet_wrap(~year) + 
  scale_color_manual(values=cbPalette) +
  scale_fill_manual(values = cbPalette, guide = NULL) +
  labs(color = "dompred") +
  theme_classic() 

#Un-comment to save
#ggsave("Figures/Traitspace.dompred.season.tiff")


#Plot averages of single traits to see if the vector loadings make sense. I'm going to try having one figure with separate panels for each trait. I'll need to convert the data to long form, group by pred and trait, and then summarize to calculate mean and SE. This will go in the supplement so I'm going to have separate figures for each season.

#Convert data to long form
Traitslong <- Finaldata%>%
  dplyr::select(pond, season, year, species, dompred, Body_A:Mentum_L)%>%
  gather(trait, value, Body_A:Mentum_L)%>%
  group_by(dompred, season, trait) 

library(ggbeeswarm)

PlotTraits <- function(seas){
  Traitslong%>%
  filter(season == seas)%>%
  summarize(traitmean = mean(value), count = n(), stdv = sd(value), se = stdv/sqrt(n))%>%
  ggplot(aes(x = dompred, y = traitmean, color = dompred)) +
    geom_point(size = 5) +
    geom_linerange(aes(ymin = traitmean - se, ymax = traitmean + se), size = 1) +
    #geom_beeswarm(data = filter(Traitslong, season == "W"), aes(x = dompred, y = value)) + #For some reason this won't run
facet_wrap(vars(trait), scales = "free_y") +
  scale_color_manual(values = predPalette) +
  labs(x = "Top Predator", y = "Trait mean (mm) +/- SE") +
  scale_x_discrete(labels = c("Inv", "Sal", "Sun", "Bass")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 16)) 


}

#Plot winter traits
PlotTraits("W") +
  ggtitle("Winter") +
  theme(plot.title = element_text(size = 20, face = "bold"))

ggsave("Figures/traits.winter.jpg", width = 10.19, height = 9.19)


#Plot spring traits
PlotTraits("Sp") +
  ggtitle("Spring") +
  theme(plot.title = element_text(size = 20, face = "bold"))

ggsave("Figures/traits.spring.jpg", width = 10.19, height = 9.19)


#Plot summer traits
PlotTraits("Su") +
  ggtitle("Summer") +
  theme(plot.title = element_text(size = 20, face = "bold"))

ggsave("Figures/traits.summer.jpg", width = 10.19, height = 9.19)

#Plot fall traits
PlotTraits("F") +
  ggtitle("Fall") +
  theme(plot.title = element_text(size = 20, face = "bold"))

ggsave("Figures/traits.fall.jpg", width = 10.19, height = 9.19)


