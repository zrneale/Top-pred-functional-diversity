---
  title: "Functional diversity figures"
output: html_notebook
---
  
  Set WD
```{r, setup, include=FALSE}
knitr::opts_knit$set(root.dir = "/Users/Zoey/Dropbox/Dragonflies/Zoey/Top-pred-and-func-div")
```

```{r}
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

FDis.rand.avg <- read.csv("Data/FDis.rand.avg.csv")%>%
  mutate(dompred = factor(dompred))

spbeta.rand.avg <- read.csv("Data/spbeta.rand.avg.csv")%>%
  mutate(dompred = factor(dompred))

tempbeta.rand.avg <- read.csv("Data/tempbeta.rand.avg.csv")%>%
  mutate(dompred = factor(dompred),
         season = fct_relevel(season, c("W_Sp", "Sp_Su", "Su_F", "F_W")))
```

Specifications for plots

```{r}
#colorblind friendly palet 
cbPalette <- c("#CC79A7", "#E69F00", "#009E73","#56B4E9", "grey35")

#Fonts. I'll only specify legend and y axis since only the bottom plot (temporal) will have x axis labels

Textsize <- theme(legend.text = element_text(size = 12),
                  legend.title = element_text(size = 14),
                  axis.title.y = element_text(size = 14),
                  axis.text.y = element_text(size = 12))


```

Alpha plot

```{r}
#I'm going to represent predator means as horizontal lines rather than points. Here's a go at that


Alpha.plot <- FDis.emmeans%>%
  filter(season != "Mean")%>%
  ggplot(aes(x = dompred, y = emmean, color = season))+
  geom_point(position=position_dodge(width=0.5),size=5)+
  geom_line(aes(group = season), position=position_dodge(width=0.5),size=0.5) +
  geom_linerange(aes(ymin=lower.CL, ymax= upper.CL),
                 position=position_dodge(width=0.5),size=0.5) +
  stat_summary(fun = "mean", geom = "segment", 
               aes(xend = ..x.. - 0.25, yend = ..y..), color = "black") +
  stat_summary(fun = "mean", geom = "segment", 
               aes(xend = ..x.. + 0.25, yend = ..y..), color = "black") +
  theme_classic() +
  labs(x = "Top predator",
       y = paste("FDis\n(α diversity)"),
       color = "Season")+
  scale_color_manual(values=cbPalette, labels = c("Winter", "Spring", "Summer", "Fall")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  Textsize
```

Spatial beta plot
```{r}

Spatial.plot <- Spatial.emmeans%>%
  filter(season != "Mean") %>%
  ggplot(aes(x = dompred, y = response, color = season), data = .) +
  geom_point(position = position_dodge(width = 0.5), size = 5) +
  geom_line(aes(group = season), position = position_dodge(width = 0.5), size = 0.5) +
  geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                 position = position_dodge(width = 0.5), size = 0.5) +
  stat_summary(fun = "mean", geom = "segment", 
               aes(xend = ..x.. - 0.25, yend = ..y..), color = "black") +
  stat_summary(fun = "mean", geom = "segment", 
               aes(xend = ..x.. + 0.25, yend = ..y..), color = "black") +
  theme_classic() +
  labs(x = "Top Predator", 
       y = paste("Spatial Dissimilarity\n(β diversity)"),
       color = "Season") +
  scale_color_manual(values=cbPalette, labels = c("Winter", "Spring", "Summer", "Fall")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  Textsize


```

Temporal beta plot

```{r}
Temporal.plot <- Temporal.emmeans%>%
  filter(season != "Mean")%>%
  data.frame()%>%
  ggplot(aes(x = dompred, y = response, color = season)) +
  geom_point(position = position_dodge(width = 0.5), size = 5) +
  geom_line(aes(group = season), position = position_dodge(width = 0.5), size = 0.5) +
  geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                 position = position_dodge(width = 0.5), size = 0.5) +
  stat_summary(fun = "mean", geom = "segment", 
               aes(xend = ..x.. - 0.25, yend = ..y..), color = "black") +
  stat_summary(fun = "mean", geom = "segment", 
               aes(xend = ..x.. + 0.25, yend = ..y..), color = "black") +
  theme_classic() +
  labs(x = "Top Predator", 
       y = paste("Temporal Dissimilarity\n(β diversity)"),
       color = "Season") +
  scale_color_manual(values=cbPalette) +
  scale_x_discrete(labels=c("Invertebrate","Salamander","Sunfish","Bass")) +
  theme(axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12)) +
  Textsize
```

```{r, fig.width = 7, fig.height = 10}
#I'm going to have one legend for the alpha and spatial beta plots and a separate one for the temporal beta. I think the best way to do this is to first make a combined plot for the first two sharing a legend, then combine that one with the temporal plot

Alpha.spatial.plot <- ggarrange(Alpha.plot, Spatial.plot, common.legend = T, legend = "right", ncol = 1)

#Add the temporal plot
ggarrange(Alpha.spatial.plot, Temporal.plot, common.legend = F, ncol = 1, heights = c(2,1))
ggsave("Figures/Diversity.combined.figure.tiff", width = 7.33, height = 9.19)
```

Randomized alpha plot

```{r}
#The randomized only includes salamander, sunfish, and bass ponds because we drew five ponds in each randomization to match the invertebrate sample size which means each randomization would give the same invertebrate values. For the figure I'll graph the observed invert mean + SE beside the resampled mean + SE for the other three pond types.


#Combine the randomized mean values with the observed values for invert ponds and make the plot. 

Alpha.rand.plot <- FDis.emmeans%>%
  filter(dompred == "N", season != "Mean")%>%
  rename(avgFDis = emmean)%>%
  dplyr::select(dompred, season, avgFDis, lower.CL, upper.CL)%>%
  rbind(FDis.rand.avg)%>%
  ggplot(aes(x = dompred, y = avgFDis, color = season), data = .) +
  geom_point(position = position_dodge(width = 0.5), size = 5) +
  #geom_line(aes(group = season), position = position_dodge(width = 0.5), size = 0.5) +
  geom_linerange(aes(ymin = lower.CL, ymax = upper.CL),
                 position = position_dodge(width = 0.5), size = 0.5) +
  theme_classic() +
  labs(x = "Top Predator", 
       y = paste("FDis\n(α diversity) +/- CI"),
       color = "Season") +
  scale_color_manual(values=cbPalette, labels = c("Winter", "Spring", "Summer", "Fall")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  Textsize

```

Randomized spatial plot
```{r}

#The randomized only includes salamander, sunfish, and bass ponds because we drew five ponds in each randomization to match the invertebrate sample size which means each randomization would give the same invertebrate values. For the figure I'll graph the observed invert mean + SE beside the resampled mean + SE for the other three pond types.


#Combine the randomized mean values with the observed values for invert ponds and make the plot. 

Spatial.rand.plot <- Spatial.emmeans%>%
  filter(dompred == "N")%>%
  rename(avgdist = response)%>%
  dplyr::select(dompred, season, avgdist, asymp.LCL, asymp.UCL)%>%
  rbind(spbeta.rand.avg)%>%
  ggplot(aes(x = dompred, y = avgdist, color = season), data = .) +
  geom_point(position = position_dodge(width = 0.5), size = 5) +
  #geom_line(aes(group = season), position = position_dodge(width = 0.5), size = 0.5) +
  geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                 position = position_dodge(width = 0.5), size = 0.5) +
  theme_classic() +
  labs(x = "Top Predator", 
       y = paste("Spatial Dissimilarity\n(β diversity) +/- CI"),
       color = "Season") +
  scale_color_manual(values=cbPalette, labels = c("Winter", "Spring", "Summer", "Fall")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  Textsize

```

Randomized temporal plot
```{r}

#The randomized only includes salamander, sunfish, and bass ponds because we drew five ponds in each randomization to match the invertebrate sample size which means each randomization would give the same invertebrate values. For the figure I'll graph the observed invert mean + SE beside the resampled mean + SE for the other three pond types.


#Combine the randomized mean values with the observed values for invert ponds and make the plot. 

Temporal.rand.plot <- Temporal.emmeans%>%
  filter(dompred == "N", season != "Mean")%>%
  rename(avgdist = response)%>%
  dplyr::select(dompred, season, avgdist, asymp.LCL, asymp.UCL)%>%
  rbind(tempbeta.rand.avg)%>%
  ggplot(aes(x = dompred, y = avgdist, color = season), data = .) +
  geom_point(position = position_dodge(width = 0.5), size = 5) +
  #geom_line(aes(group = season), position = position_dodge(width = 0.5), size = 0.5) +
  geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                 position = position_dodge(width = 0.5), size = 0.5) +
  theme_classic() +
  labs(x = "Top Predator", 
       y = paste("Temporal Dissimilarity\n(β diversity) +/- CI"),
       color = "Season") +
  scale_color_manual(values=cbPalette) +
  #theme(axis.title.x = element_blank(),
  #axis.text.x = element_blank(),
  #axis.ticks.x = element_blank()) +
  scale_x_discrete(labels=c("Invertebrate","Salamander","Sunfish","Bass")) +
  theme(axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12)) +
  Textsize

```

Combine the randomized plots
```{r, fig.width = 7, fig.height = 10}
#I'm going to have one legend for the alpha and spatial beta plots and a separate one for the temporal beta. I think the best way to do this is to first make a combined plot for the first two sharing a legend, then combine that one with the temporal plot

Alpha.spatial.rand.plot <- ggarrange(Alpha.rand.plot, Spatial.rand.plot, common.legend = T, legend = "right", ncol = 1)

#Add the temporal plot
ggarrange(Alpha.spatial.rand.plot, Temporal.rand.plot, common.legend = F, ncol = 1, heights = c(2,1))
ggsave("Figures/Diversity.rand.combined.figure.tiff", width = 7.33, height = 9.19)
```