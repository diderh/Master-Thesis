#To install biomonitoR package, please install devtools first and then use the following code.
#library(devtools)
#install_github("alexology/biomonitoR", ref = "develop", build_vignettes = TRUE)

library(biomonitoR)
library(vegan)
library(ade4)
library(StatMatch)
library(ggplot2)
library(tidyverse)

########################################
# Data table indexing 
########################################

## Insert our data table into r environment

library(readxl)
benthos_data <- read_excel("Benthos_final_updatedfile_Tanypodinae.xlsx")


## Reorganize the data table into a more elongated format, where individual rows represent distinct taxa names. 
## This restructuring will facilitate the comparison between our taxonomic data and 
## the reference data required for trait data generation.

longer_format_taxa <- benthos_data %>%
  pivot_longer(cols = starts_with(c("Gammarus", "Chironomidae", "Tanypodinae", "Baetis", "Ephemerellidae", 
                                    "Tipulidae", "Oligochaeta")),
               names_to = "Taxa",
               values_to = "Abundance")
glimpse(longer_format_taxa)

#### grouping

most.abun.taxa <- longer_format_taxa %>%
  select(Treatment, Time_point, Position, Taxa, Abundance)

total.abundance <- longer_format_taxa %>%
  select(Treatment, Taxa, Abundance) %>%
  group_by(Treatment) %>%
  summarise(Total_Abundance = sum(Abundance), .groups = "drop") 
  
data.with.total <- most.abun.taxa %>%
  left_join(total.abundance, by = "Treatment") %>%
  mutate(Rel.abun = Abundance / Total_Abundance)

## Rename time point factor levels

data.with.total$Time_point <- recode(data.with.total$Time_point, "week -1" = "Before Treatment", "week_6" = "After Treatment")

### Converting character variable to factor

data.with.total$Position <- as.factor(data.with.total$Position)
data.with.total$Treatment <- as.factor(data.with.total$Treatment)
data.with.total$Time_point <- as.factor(data.with.total$Time_point)

str(data.with.total)


## Reorders the levels of factor variables

data.with.total$Treatment <- relevel(data.with.total$Treatment, ref="Control")
levels(data.with.total$Treatment)

data.with.total$Time_point <- relevel(data.with.total$Time_point, ref="Before Treatment")
levels(data.with.total$Time_point)

data.with.total$Position <- relevel(data.with.total$Position, ref="Upstream")
levels(data.with.total$Position)

## Final data

data.with.total <- data.with.total %>%
  mutate(interaction = interaction(Treatment, Time_point))


## Reorders the levels of factor variables
data.with.total$interaction <- relevel(data.with.total$interaction, ref="Crayfish + ALAN.After Treatment")
data.with.total$interaction <- relevel(data.with.total$interaction, ref="Crayfish + ALAN.Before Treatment")
data.with.total$interaction <- relevel(data.with.total$interaction, ref="Crayfish.After Treatment")
data.with.total$interaction <- relevel(data.with.total$interaction, ref="Crayfish.Before Treatment")
data.with.total$interaction <- relevel(data.with.total$interaction, ref="ALAN.After Treatment")
data.with.total$interaction <- relevel(data.with.total$interaction, ref="ALAN.Before Treatment")
data.with.total$interaction <- relevel(data.with.total$interaction, ref="Control.After Treatment")
data.with.total$interaction <- relevel(data.with.total$interaction, ref="Control.Before Treatment")
levels(data.with.total$interaction)


## Plotting relative abundance of taxa

rel.abu.plot <- ggplot(data.with.total, aes(x = interaction, y = (Rel.abun)*100, fill = Taxa, alpha = Time_point)) +
  geom_col() + 
  facet_wrap(~Position, scales = "free_y") +
  scale_alpha_manual(values=c(0.7,1)) +
  guides(alpha=guide_legend(override.aes=list(fill=hcl(c(15,195),200,0,alpha=c(0.2,1))))) +
  theme_classic()+
  scale_fill_manual(values=c("green", "darkgreen", "#E69F00", "red", "darkred", "#56B4E9", "darkblue"))+
  labs(x = "Treatment", 
       y = "Relative Abundance (%)",
       alpha = "Time point",
       tag = "a)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,family = "Arial", size = 12, colour = c("black", "black", "darkgreen", "darkgreen", "blue", "blue", "darkred", "darkred"), margin = margin(b = 15)),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, family = "Arial", face = "bold"),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(t = 12),family = "Arial", size = 14, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 16),family = "Arial", size = 14, face = "bold"),
        legend.title = element_text(color = "darkblue",margin = margin(r = 15),family = "Arial", size = 14, face = "bold"),
        legend.text = element_text(size = 13, family = "Arial"),
        legend.position = "top",
        legend.background = element_rect(fill = "transparent"),
        axis.ticks.length = unit(.2, "cm"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold", color = "white")) 

rel.abu.plot

ggsave("RelativeAbundancePlot.png", rel.abu.plot, width = 14, height = 8, dpi = 1600, device = "png")
  

