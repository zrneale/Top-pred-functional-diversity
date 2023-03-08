
library(tidyverse)
library(ggpubr)


#Data sets for the observed diversity results

FDis.emmeans <- read.csv("Data/FDis.emmeans.csv")%>%
  mutate(dompred = fct_relevel(dompred, c("N", "S", "G", "B")),
         season = fct_relevel(season, c("W", "Sp", "Su", "F", "Mean")))

Spatial.emmeans <- read.csv("Data/Spatial.emmeans.csv")%>%
  mutate(dompred = fct_relevel(dompred, c("N", "S", "G", "B")),
         season = fct_relevel(season, c("W", "Sp", "Su", "F")))

Temporal.emmeans <- read.csv("Data/Temporal.emmeans.csv")%>%
  mutate(dompred = fct_relevel(dompred, c("N", "S", "G", "B")),
         season = fct_relevel(season, c("W_Sp", "Sp_Su", "Su_F", "F_W")))

#Data sets for the randomized diversity results

FDis.rand <- read.csv("Data/FDis.rand.csv")%>%
  mutate(dompred = fct_relevel(dompred, c("N", "S", "G", "B")),
         season = fct_relevel(season, c("W", "Sp", "Su", "F")))

spbeta.rand <- read.csv("Data/spbeta.rand.csv")%>%
  mutate(dompred = fct_relevel(dompred, c("N", "S", "G", "B")),
         season = fct_relevel(season, c("W", "Sp", "Su", "F")))


tempbeta.rand <- read.csv("Data/tempbeta.rand.csv")%>%
  mutate(dompred = fct_relevel(dompred, c("N", "S", "G", "B")),
         season = fct_relevel(season, c("W_Sp", "Sp_Su", "Su_F", "F_W")))


#Specifications for plots

#colorblind friendly palette when faceting by season 
seasonPalette <- c('#EE6677', '#228833', '#4477AA', '#CCBB44', "grey35")
#Colorblind friendly palette wheb faceting by predator
predPalette <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00")
#Fonts. I'll only specify legend and y axis since only the bottom plot (temporal) will have x axis labels

Textsize <- theme(legend.text = element_text(size = 12),
                  legend.title = element_text(size = 14),
                  axis.title.y = element_text(size = 14),
                  axis.text.y = element_text(size = 12))



#Alpha diversity plot

#I'm going to represent predator means as horizontal lines rather than points. Here's a go at that


Alpha.plot <- FDis.emmeans%>%
  filter(season != "Mean")%>%
  ggplot(aes(x = dompred, y = emmean, color = season))+
  geom_point(position=position_dodge(width=0.5),size=5)+
  geom_line(aes(group = season), position=position_dodge(width=0.5), linewidth=0.5) +
  geom_linerange(aes(ymin=lower.CL, ymax= upper.CL),
                 position=position_dodge(width=0.5),size=0.5) +
  stat_summary(fun = "mean", geom = "point", color = "black", alpha = 0.6, shape = 17, size = 3) + 
  theme_classic() +
  labs(x = "Top predator",
       y = paste("FDis\n(α diversity)"),
       color = "Season")+
  scale_color_manual(values=seasonPalette, labels = c("Winter", "Spring", "Summer", "Fall")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  Textsize

#Spatial beta plot

Spatial.plot <- Spatial.emmeans%>%
  filter(season != "Mean") %>%
  ggplot(aes(x = dompred, y = response, color = season), data = .) +
  geom_point(position = position_dodge(width = 0.5), size = 5) +
  geom_line(aes(group = season), position = position_dodge(width = 0.5), linewidth = 0.5) +
  geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                 position = position_dodge(width = 0.5), size = 0.5) +
  stat_summary(fun = "mean", geom = "point", color = "black", alpha = 0.6, shape = 17, size = 3) + 
  theme_classic() +
  labs(x = "Top Predator", 
       y = paste("Spatial Dissimilarity\n(β diversity)"),
       color = "Season") +
  scale_color_manual(values=seasonPalette, labels = c("Winter", "Spring", "Summer", "Fall")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  Textsize


#Temporal beta plot

Temporal.plot <- Temporal.emmeans%>%
  filter(season != "Mean")%>%
  data.frame()%>%
  ggplot(aes(x = dompred, y = response, color = season)) +
  geom_point(position = position_dodge(width = 0.5), size = 5) +
  geom_line(aes(group = season), position = position_dodge(width = 0.5), linewidth = 0.5) +
  geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                 position = position_dodge(width = 0.5), size = 0.5) +
  stat_summary(fun = "mean", geom = "point", color = "black", alpha = 0.6, shape = 17, size = 3) + 
  theme_classic() +
  labs(x = "Top Predator", 
       y = paste("Temporal Dissimilarity\n(β diversity)"),
       color = "Season") +
  scale_color_manual(values=seasonPalette) +
  scale_x_discrete(labels=c("Invertebrate","Salamander","Sunfish","Bass")) +
  theme(axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12)) +
  Textsize


#I'm going to have one legend for the alpha and spatial beta plots and a separate one for the temporal beta. I think the best way to do this is to first make a combined plot for the first two sharing a legend, then combine that one with the temporal plot

library(ggpubr)
Alpha.spatial.plot <- ggarrange(Alpha.plot, Spatial.plot, common.legend = T, legend = "right", ncol = 1)

#Add the temporal plot
ggarrange(Alpha.spatial.plot, Temporal.plot, common.legend = F, ncol = 1, heights = c(2,1))
ggsave("Figures/Diversity.combined.figure.jpg", width = 7.33, height = 9.19)


####Randomization data plots


#Randomized alpha plot

library(scales) #for setting the number of digits in FDis values

Alpha.rand.plot <- FDis.rand%>%
  group_by(dompred, season, sim)%>%
  dplyr::summarise(FDis = mean(FDis))%>%
  ggplot(aes(x = dompred, y = FDis, color = dompred)) +
  geom_violin(aes(fill = dompred)) +
  theme_classic() +
  scale_color_manual(values=predPalette, labels = c("Invertebrate", "Salamander", "Sunfish", "Bass"), name = "Top Predator") +
  scale_fill_manual(values = predPalette, labels = c("Invertebrate", "Salamander", "Sunfish", "Bass"), name = "Top Predator") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Predator", 
       y = "Mean FDis") +
  facet_grid(.~season, labeller = labeller(season = c("W" = "Winter", "Sp" = "Spring", 
                                                      "Su" = "Summer", "F" = "Fall"))) +
  labs(x = "Top Predator", 
       y = paste("Mean FDis"),
       legend.title = "Top Predator") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  Textsize +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) 


#Randomized spatial plot

Spatial.rand.plot <- spbeta.rand%>%
  group_by(dompred, season, sim)%>%
  dplyr::summarise(mean.dist = mean(distance))%>%
  ggplot(aes(x = dompred, y = mean.dist, color = dompred)) +
  geom_violin(aes(fill = dompred)) +
  theme_classic() +
  scale_color_manual(values=predPalette, labels = c("Invertebrate", "Salamander", "Sunfish", "Bass"), name = "Top Predator") +
  scale_fill_manual(values = predPalette, labels = c("Invertebrate", "Salamander", "Sunfish", "Bass"), name = "Top Predator") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  labs(x = "Predator", 
       y = "Mean Spatial Dissimilarity") +
  facet_grid(.~season, labeller = labeller(season = c("W" = "Winter", "Sp" = "Spring", 
                                                      "Su" = "Summer", "F" = "Fall"))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  Textsize

#Randomized temporal plot

Temporal.rand.plot <-tempbeta.rand%>%
  group_by(dompred, season, sim)%>%
  dplyr::summarise(mean.dist = mean(distance))%>%
  ggplot(aes(x = dompred, y = mean.dist, color = dompred)) +
  geom_violin(aes(fill = dompred)) +
  theme_classic() +
  scale_color_manual(values=predPalette, labels = c("Invertebrate", "Salamander", "Sunfish", "Bass"), name = "Top Predator") +
  scale_fill_manual(values = predPalette, labels = c("Invertebrate", "Salamander", "Sunfish", "Bass"), name = "Top Predator") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Randomized Spatial Dissimilarity") +
  labs(x = "Predator", 
       y = "Mean Spatial Dissimilarity") +
  facet_grid(.~season, labeller = labeller(season = c("W_Sp" = "Winter_Spring", 
                                                      "Sp_Su" = "Spring_Summer", 
                                                      "Su_F" = "Summer_Fall", 
                                                      "F_W" = "Fall_Winter"))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  Textsize


#Combine the randomized plots

ggarrange(Alpha.rand.plot, Spatial.rand.plot, Temporal.rand.plot, common.legend = T, 
          legend = "right", ncol = 1)

ggsave("Figures/Diversity.rand.combined.jpg", width = 7.33, height = 9.19)
