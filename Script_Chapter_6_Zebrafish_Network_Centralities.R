# Chapter 6: Zebrafish Network metrics - Giant component ####
# Author: R. Liscovsky (rliscovs@ed.ac.uk)

# 1. Install & Load Packages ####

if (!require("ggpubr")) install.packages("ggpubr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")

library(ggpubr)
library(ggplot2)
library(dplyr)

# 2. Read Data In ####

networ_centralities_bt <- read.csv("Data/degree_distribution2_betweenness.csv",
                          header=TRUE,
                          stringsAsFactors=FALSE, 
                          check.names = FALSE, 
                          sep = ",")

networ_centralities_wd <- read.csv("Data/degree_distribution2_.csv",
                                header=TRUE,
                                stringsAsFactors=FALSE, 
                                check.names = FALSE, 
                                sep = ",")


head(networ_centralities_wd)


# 3. Run visual ####

ggplot(data = networ_centralities_wd, aes(continent_iso_2, Weighted_Degree)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue") +
  labs(title = "Count of papers per year in each volume",
       subtitle = "Volumes 1 to 49",
       y = "Count of papers", x = "Year of publication") + 
  facet_wrap(~ Period)

## Final visual with statistics ####

### Weighted_Degree  ####

net.summary <- networ_centralities_wd %>%
  group_by(Period, Community) %>%
  summarise(
    sd = sd(Weighted_Degree, na.rm = TRUE),
    Weighted_Degree = mean(Weighted_Degree)
  )

net.summary

net.summary2 <- net.summary
net.summary2$Weighted_Degree <- as.numeric(net.summary2$Weighted_Degree )


ggplot(data = networ_centralities_wd, aes(Community, Weighted_Degree)) +
  geom_jitter(aes(color = "steelblue")) +
  geom_line(aes(group = "Period"), data = net.summary) +
  geom_errorbar(
    aes(ymin = Weighted_Degree-sd, ymax = Weighted_Degree+sd),
    data = net.summary, width = 0.2) +
  geom_point(data = net.summary, size = 2) +
  theme(legend.position = "none") +
  labs(subtitle="Weighted degree distribution", 
       y="Weighted degree", 
       x="Community", 
       title="Score distribution in the 6 regional communities", 
       caption = "Source: own preparation based on Scopus data")+
  scale_color_manual(guide=FALSE, values="darkgrey") + #turn off the legend, define the colors
  theme_bw()+
  facet_wrap(~ Period)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

### betweenesscentrality####

net.summary <- networ_centralities_bt %>%
  group_by(Period, Community) %>%
  summarise(
    sd = sd(betweenesscentrality, na.rm = TRUE),
    betweenesscentrality = mean(betweenesscentrality)
  )
net.summary

ggplot(data = networ_centralities_bt, aes(Community, betweenesscentrality)) +
  geom_jitter(aes(color = "steelblue")) +
  geom_line(aes(group = "Period"), data = net.summary) +
  geom_errorbar(
    aes(ymin = betweenesscentrality-sd, ymax = betweenesscentrality+sd),
    data = net.summary, width = 0.2) +
  geom_point(data = net.summary, size = 2) +
  theme(legend.position = "none") +
  labs(subtitle="Betweenness centrality distribution", 
       y="Betweenness (normalized)", 
       x="Community", 
       title="Score distribution in the 6 regional communities", 
       caption = "Source: own preparation based on Scopus data")+
  scale_color_manual(guide=FALSE, values="darkgrey") + #turn off the legend, define the colors
  theme_bw()+
  facet_wrap(~ Period)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
