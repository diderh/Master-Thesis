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
                                    "Heptageniidae", "Leptophlebiidae", "Hydropsychidae", "Rhyacophilidae", 
                                    "Simuliidae", "Tipulidae", "Limoniidae", "Dixidae", "Culicidae", 
                                    "Ceratopogonidae", "Psychodidae", "Empididae", "Ephydridae", 
                                    "Pediciidae", "Elmis", "Limnius", "Hydraenidae", "Chrysomelidae", 
                                    "Oligochaeta")),
               names_to = "Taxa",
               values_to = "Abundance")
glimpse(longer_format_taxa)

## Generate a column of samples corresponding to the names specified in our lab samples.

longer_format_taxa$Sample_name <- paste(longer_format_taxa$Channel,
                                        longer_format_taxa$Treatment,
                                        longer_format_taxa$Position,
                                        longer_format_taxa$Time_point, sep = "_")

## Creating final data table
Final_data <- longer_format_taxa%>%
  select(Sample_name, Taxa, Abundance) %>%
  pivot_wider(names_from = Sample_name, values_from = Abundance, values_fill = 0)

## Select a subset of our data table that was utilized in the process of generating trait data.

datafortraitgen <- Final_data
datafortraitgen[1:6, 1:6]



##############################
# Assign trait data 
##############################

## Applying as_biomonitor function to our species abundance data,  
## this function will merge our data with the reference taxonomic data table.
## I selected group = "mi" which means that "as_biomonitor" function will merge our data
## with the reference macro-invertebrate taxonomic data table

biomonitor_data <- as_biomonitor(datafortraitgen, group = "mi", traceB = FALSE)


## Aggregate taxa

taxa_aggregation <- aggregate_taxa(biomonitor_data)
print(taxa_aggregation)

tree <- taxa_aggregation$Tree
## Structure of macro-invertebrate community. For each observation (n = 64) the relative abundance is reported. 
## Observations were ordered using an agglomeration hierarchical clustering with 
## the Ward.D2 method to enhance pattern detection based on community similarity.

plot(taxa_aggregation)

## Assign trait

traitscailing_data <- assign_traits( taxa_aggregation, tax_lev = "Taxa")
colnames(traitscailing_data)

## Calculating average traits for selected traits

All_trait_av <- average_traits(traitscailing_data)

colnames(All_trait_av)

## Selecting our desired traits

selected_trait_av <- All_trait_av[, c(1, 52:114)]

## Ommittig NA values in our dataset
selected_trait_av[is.na(selected_trait_av)] <- 0

colnames(selected_trait_av)


###################################
# Calculating Functional index
###################################

library(ade4)

## we are using fuzzy coded traits. For this reason, we set the number of categories for each trait.
Trait_category_nr <- c(9,8,8,5,7,5,4,4,2,3,8)


## Functional dispersion

Functional_dispersion <- f_disp(taxa_aggregation, trait_db = selected_trait_av, nbdim = 2, type = "F", col_blocks = Trait_category_nr, correction = "cailliez")
Functional_dispersion <- as.data.frame(Functional_dispersion)
Functional_dispersion$Sample_name <- colnames(datafortraitgen[-1])
glimpse(Functional_dispersion)


## Functional diversity

Functional_diversity <- f_divs(taxa_aggregation, trait_db = selected_trait_av, type = "F", col_blocks = Trait_category_nr, correction = "cailliez")
Functional_diversity <- as.data.frame(Functional_diversity)
Functional_diversity$Sample_name <- colnames(datafortraitgen[-1])
glimpse(Functional_diversity)


## Functional redundancy

Functional_redundancy <- f_red(taxa_aggregation, trait_db = selected_trait_av, type = "F", col_blocks = Trait_category_nr, correction = "cailliez")
Functional_redundancy <- as.data.frame(Functional_redundancy)
Functional_redundancy$Sample_name <- colnames(datafortraitgen[-1])
glimpse(Functional_redundancy)


## Functional evenness

Functional_evenness <- f_eve(taxa_aggregation, trait_db = selected_trait_av, type = "F", nbdim = 2, col_blocks = Trait_category_nr, correction = "cailliez")
Functional_evenness <- as.data.frame(Functional_evenness)
Functional_evenness$Sample_name <- colnames(datafortraitgen[-1])
glimpse(Functional_evenness)


## Functional richness

Functional_richness <- f_rich(taxa_aggregation, trait_db = selected_trait_av, type = "F", nbdim = 2, col_blocks = Trait_category_nr, correction = "cailliez")
Functional_richness <- as.data.frame(Functional_richness)
Functional_richness$Sample_name <- colnames(datafortraitgen[-1])
glimpse(Functional_richness)

## Community weighted mean trait values

CWM <- cwm(taxa_aggregation, trait_db = selected_trait_av, tax_lev = "Taxa", trans = log1p, traceB = FALSE)
CWM <- as.data.frame(CWM)
CWM$Sample_name <- colnames(datafortraitgen[-1])
glimpse(CWM)



###################################################
# Final data frame for functional diversity indices
###################################################

## Gather functional indices for each traits

FD_indices <- merge(Functional_richness, Functional_dispersion, by = "Sample_name")
FD_indices <- merge(FD_indices, Functional_evenness, by = "Sample_name")
FD_indices <- merge(FD_indices, Functional_diversity, by = "Sample_name")
FD_indices <- merge(FD_indices, Functional_redundancy, by = "Sample_name")

## Removing NA values in our data table
FD_indices[is.na(FD_indices)] <- 0

glimpse(FD_indices)

################################################################################################################

#####################################################
## Multivariate analysis for macro invertebrate functional traits
#####################################################

###############
###############################
## NMDS 
###############################
library(FD)
library(ade4)
library(factoextra)
library(FactoMineR)
library(corrr)
library(ggcorrplot)
library(vegan)
colnames(CWM)
NMDS.data <- CWM%>%
  select(-Sample) %>%
  separate(Sample_name, into = c("Channel", "Treatment", "Position", "Time_point"), sep = "_", remove = FALSE)

##  Rename factor levels
NMDS.data$Time_point <- recode(NMDS.data$Time_point, "week -1" = "Before Treatment", "week" = "After Treatment")

## Final PCA data
NMDS.data <- NMDS.data %>%
  mutate(time.pos =paste(Time_point, Position, sep = "_"))
colnames(CWM)
rownames(NMDS.data) <- NMDS.data$Sample_name
colnames(NMDS.data)
colnames(NMDS.data)[1:63] <- c("Aquatic_passive", "Aquatic_active", "Aerial_passive", "Aerial_active",
                              "Absorber", "Deposit_feeder", "Shredder", "Scraper", "Filter_feeder", "Piercer", "Predator", "Parasite",
                              "Microorganism", "Detritus <1 mm", "Dead plant ≥1 mm", "Living_microphytes", "Living_macrophytes", "Dead animal ≥1 mm", "Living_microinvertebrates", "Living_macroinvertebrates", "Vertebrates",
                              "≤1 year", ">1 year",
                              "Flier", "Surface_swimmer", "F_Water_swimmer", "Crawler", "Burrower", "Interstitial", "Temp_attached", "Perm_attached",
                              "Ovoviviparity", "I_free_eggs", "I_cemented_eggs", "C_cemented", "C_free", "C_in_vegetation", "C_terrestrial", "Asexual",
                              "Eggs/statoblasts", "Cocoons", "Housings_against_desiccation", "Diapause/dormancy", "None",
                              "Tegument", "Gill", "Plastron", "Spiracle", "Hydrostatic_vesicle",
                              "≤0.25 cm", ">0.25–0.5 cm", ">0.5–1 cm", ">1–2 cm", ">2–4 cm", ">4–8 cm", ">8 cm",
                              "Egg", "Larva", "Nymph", "Adult",
                              "<1", "1", ">1")

## Converting character variables to factor variables
NMDS.data$Time_point <- factor(NMDS.data$Time_point)
NMDS.data$Treatment <- factor(NMDS.data$Treatment)
NMDS.data$Position <- factor(NMDS.data$Position)
NMDS.data$time.pos  <- factor(NMDS.data$time.pos )
str(NMDS.data)

## Reorders the levels of factor variable

NMDS.data$Treatment <- relevel(NMDS.data$Treatment, ref="Control")
levels(NMDS.data$Treatment)

NMDS.data$Time_point <- relevel(NMDS.data$Time_point, ref="Before Treatment")
levels(NMDS.data$Time_point)

NMDS.data$Position <- relevel(NMDS.data$Position, ref="Upstream")
levels(NMDS.data$Position)

# Create the distance matrix with ade4 functions.
trait.matrix <- NMDS.data[, 1:63]
colnames(trait.matrix)
trait.modalities <- c(4,8,9,2,8,8,5,5,7,4,3)

traits_prep <- prep.fuzzy(trait.matrix, col.blocks = trait.modalities)
traits_dist <- ktab.list.df(list(traits_prep))
traits_dist <- dist.ktab(traits_dist, type = "F")
str(traits_dist)

## NMDS
set.seed(123)
nmds <- metaMDS(traits_dist, k = 2, trymax = 1000)
NMDS_result <- metaMDS(traits_dist, k = 2, trymax = 1000)%>%
  scores() %>%
  as_tibble(rownames = "Sample_name")

NMDS.data.result <- inner_join(NMDS.data, NMDS_result, by = "Sample_name")


######################################
## Plotting NMDS for factor variables
######################################

### Before treatment at upstream position
centroid_bu <- NMDS.data.result %>%
  filter(Time_point == "Before Treatment", Position == "Upstream") %>%
  group_by(Treatment) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))


plot_bu <- NMDS.data.result %>%
  filter(Time_point == "Before Treatment", Position == "Upstream")

bu <- ggplot(plot_bu, aes(x = NMDS1, y = NMDS2, color = Treatment)) +
  geom_point() +
  geom_point(data = centroid_bu, size = 5, shape = 21, color = "black",
             aes(fill = Treatment), show.legend = FALSE) +
  geom_text(aes(x = -0.24, y = 0.15,
                label = "Stress = 0.146"),
            stat = "unique", family ="Arial", size = 4.5, color = "black") +
  coord_fixed(ratio = 1) +
  theme_bw() + 
  theme(axis.text=element_text(size=14,family="Arial"),
        axis.title = element_text(size=16,family="Arial", face = "bold"),
        legend.text = element_text(size=14,family="Arial"),
        legend.title = element_text(size=16,family="Arial", face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 16, family = "Arial", face = "bold", color = "darkblue"),
        plot.tag = element_text(size = 18, family = "Arial", face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        aspect.ratio = 0.56) +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0.15))),
               data = plot_bu[plot_bu$Treatment != "versicolor",],
               show.legend = FALSE)+
  scale_color_discrete(name = "Treatment") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title = "Before treatment at upstream position", tag = "a)")

bu

## After treatment at upstream position
centroid_au <- NMDS.data.result %>%
  filter(Time_point == "After Treatment", Position == "Upstream") %>%
  group_by(Treatment) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))

plot_au <- NMDS.data.result %>%
  filter(Time_point == "After Treatment", Position == "Upstream")

au <- ggplot(plot_au, aes(x = NMDS1, y = NMDS2, color = Treatment)) +
  geom_point() +
  geom_point(data = centroid_au, size = 5, shape = 21, color = "black",
             aes(fill = Treatment), show.legend = FALSE) +
  coord_fixed(ratio = 1) +
  geom_text(aes(x = -0.18, y = -0.14,
                label = "Stress = 0.146"),
            stat = "unique", family ="Arial", size = 4.5, color = "black") +
  theme_bw() + 
  theme(axis.text=element_text(size=14,family="Arial"),
        axis.title = element_text(size=16,family="Arial", face = "bold"),
        strip.text = element_text(size = 15, family = "Arial"),
        legend.title = element_text(size=16,family="Arial", face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 16, family = "Arial", face = "bold", color = "darkblue"),
        plot.tag = element_text(size = 18, family = "Arial", face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        aspect.ratio = 0.56) +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0.15))),
               data = plot_au[plot_au$Treatment != "versicolor",],
               show.legend = FALSE)+
  scale_color_discrete(name = "Treatment") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title = "After treatment at upstream position", tag = "b)")

au

## Before treatment at downstream position
centroid_bd <- NMDS.data.result %>%
  filter(Time_point == "Before Treatment", Position == "Downstream") %>%
  group_by(Treatment) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))

plot_bd <- NMDS.data.result %>%
  filter(Time_point == "Before Treatment", Position == "Downstream")

bd <- ggplot(plot_bd, aes(x = NMDS1, y = NMDS2, color = Treatment)) +
  geom_point() +
  geom_point(data = centroid_bd, size = 5, shape = 21, color = "black",
             aes(fill = Treatment), show.legend = FALSE) +
  coord_fixed(ratio = 1) +
  geom_text(aes(x = 0.32, y = 0.35,
                label = "Stress = 0.146"),
            stat = "unique", family ="Arial", size = 4.5, color = "black") +
  theme_bw() + 
  theme(axis.text=element_text(size=14,family="Arial"),
        axis.title = element_text(size=16,family="Arial", face = "bold"),
        legend.text = element_text(size=14,family="Arial"),
        legend.title = element_text(size=16,family="Arial", face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 16, family = "Arial", face = "bold", color = "darkred"),
        plot.tag = element_text(size = 18, family = "Arial", face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        aspect.ratio = 0.56) +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0.15))),
               data = plot_bd[plot_bd$Treatment != "versicolor",],
               show.legend = FALSE)+
  scale_color_discrete(name = "Treatment") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title = "Before treatment at downstream position", tag = "c)")
bd

## After treatment at downstream position
centroid_ad <- NMDS.data.result %>%
  filter(Time_point == "After Treatment", Position == "Downstream") %>%
  group_by(Treatment) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))

plot_ad <- NMDS.data.result %>%
  filter(Time_point == "After Treatment", Position == "Downstream")

ad <- ggplot(plot_ad, aes(x = NMDS1, y = NMDS2, color = Treatment)) +
  geom_point() +
  geom_point(data = centroid_ad, size = 5, shape = 21, color = "black",
             aes(fill = Treatment), show.legend = FALSE) +
  geom_text(aes(x = 0.21, y = -0.22,
                label = "Stress = 0.146"),
            stat = "unique", family ="Arial", size = 4.5, color = "black") +
  coord_fixed(ratio = 1) +
  theme_bw() + 
  theme(axis.text=element_text(size=14,family="Arial"),
        axis.title = element_text(size=16,family="Arial", face = "bold"),
        legend.text = element_text(size=14,family="Arial"),
        legend.title = element_text(size=16,family="Arial", face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 16, family = "Arial", face = "bold", color = "darkred"),
        plot.tag = element_text(size = 18, family = "Arial", face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        aspect.ratio = 0.56) +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0.15))),
               data = plot_ad[plot_ad$Treatment != "versicolor",],
               show.legend = FALSE)+
  scale_color_discrete(name = "Treatment") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title = "After treatment at downstream position", tag = "d)")
ad

library(ggpubr)
nmds_plot <- ggarrange(bu, bd, au, ad, nrow = 2, ncol = 2, common.legend = TRUE, legend = "top")
ggsave("NMDS_plot.png", nmds_plot, width = 16, height = 8, dpi = 1600, device = "png")



###################################
## Overall PERMANOVA for functional traits
###################################
## Performing PERMANOVA
set.seed(42)

perm <- adonis2(traits_dist ~ Treatment*Time_point*Position,
                 data = NMDS.data, permutations = 999)
perm

#############################################################################
## PERMDISP Treatment
#############################################################################

##  Rename factor levels
NMDS.data$Time_point <- recode(NMDS.data$Time_point, "Before Treatment" = "Week(-1)",  "After Treatment" = "Week6")
NMDS.data$time.pos <- recode(NMDS.data$time.pos, "Before Treatment_Upstream" = "Week(-1)_Upstream", "Before Treatment_Downstream" = "Week(-1)_Downstream",
                             "After Treatment_Upstream" = "Week6_Upstream", "After Treatment_Downstream" = "Week6_Downstream")

str(NMDS.data)

## Reorders the levels of factor variable
NMDS.data$time.pos <- factor(NMDS.data$time.pos)
NMDS.data$time.pos <- relevel(NMDS.data$time.pos, ref="Week(-1)_Downstream")
levels(NMDS.data$time.pos)

NMDS.data$time.pos <- relevel(NMDS.data$time.pos, ref="Week(-1)_Upstream")
levels(NMDS.data$time.pos)

## PERMDISP
library(vegan)
Trt <- betadisper(traits_dist, NMDS.data$Treatment)
anova(Trt)
permutest(Trt, pairwise = TRUE)
TukeyHSD(Trt)

## Boxplot for PERMDISP
tr.permdisp.data <- data.frame(Distance_to_centroid = Trt$distances, Group = NMDS.data$Treatment)

perm1 <- ggplot(data=tr.permdisp.data,aes(x=Group,y=Distance_to_centroid, fill = Group))+
  scale_fill_manual(values = c("grey", "#009E73", "#D55E00", "lightblue"))+
  geom_boxplot()+
  theme_bw() +
  theme(axis.text.x=element_text(margin = margin(b = 10),colour = "black", size = 12, family = "Arial"),
        axis.text.y=element_text(margin = margin(l = 10),colour = "black",size=12,family="Arial"),
        axis.text=element_text(size=14,family="Arial"),
        axis.title = element_text(size=14,family="Arial", face = "bold"),
        panel.grid = element_blank(),
        plot.title=element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold", color = "black"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold"),
        legend.position = "none",
        axis.ticks.length = unit(.2, "cm")) +
  labs(x = "Treatment", y = "Distance to centroid", tag = "c)")

#############################################################################
## PERMDISP for Time point
#############################################################################

## PERMDISP
tp <- betadisper(traits_dist, NMDS.data$Time_point)
anova(tp)
permutest(tp, pairwise = TRUE)
TukeyHSD(tp)



## Boxplot for PERMDISP
tp.permdisp.data <- data.frame(Distance_to_centroid = tp$distances, Group = NMDS.data$Time_point)

perm2 <- ggplot(data=tp.permdisp.data,aes(x=Group,y=Distance_to_centroid, fill = Group))+
  scale_fill_manual(values = c("grey", "lightblue"))+
  geom_boxplot()+
  theme_bw() +
  theme(axis.text.x=element_text(margin = margin(b = 10),colour = "black", size = 12, family = "Arial"),
        axis.text.y=element_text(margin = margin(l = 10),colour = "black",size=12,family="Arial"),
        axis.text=element_text(size=14,family="Arial"),
        axis.title = element_text(size=14,family="Arial", face = "bold"),
        panel.grid = element_blank(),
        plot.title=element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold", color = "black"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold"),
        legend.position = "none",
        axis.ticks.length = unit(.2, "cm")) +
  labs(x = "Time point", y = "Distance to centroid", tag = "e)")


####################################################################
## PERMDISP for Position
####################################################################
## PERMDISP
pos <- betadisper(traits_dist, NMDS.data$Position)
anova(pos)
permutest(pos, pairwise = TRUE)
TukeyHSD(pos)

## Boxplot for PERMDISP
pos.permdisp.data <- data.frame(Distance_to_centroid = pos$distances, Group = NMDS.data$Position)

perm3 <- ggplot(data=pos.permdisp.data,aes(x=Group,y=Distance_to_centroid, fill = Group))+
  scale_fill_manual(values = c("grey", "lightblue"))+
  geom_boxplot()+
  theme_bw() +
  theme(axis.text.x=element_text(margin = margin(b = 10),colour = "black", size = 12, family = "Arial"),
        axis.text.y=element_text(margin = margin(l = 10),colour = "black",size=12,family="Arial"),
        axis.text=element_text(size=14,family="Arial"),
        axis.title = element_text(size=14,family="Arial", face = "bold"),
        panel.grid = element_blank(),
        plot.title=element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold", color = "black"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold"),
        legend.position = "none",
        axis.ticks.length = unit(.2, "cm")) +
  labs(x = "Position", y = "Distance to centroid", tag = "d)")



#############################################################################
## PERMDISP for time point and position interaction
#############################################################################

## PERMDISP
tp.pos <- betadisper(traits_dist, NMDS.data$time.pos)
anova(tp.pos)
permutest(tp.pos, pairwise = TRUE)
TukeyHSD(tp.pos)

## Boxplot for PERMDISP
tp.pos.permdisp.data <- data.frame(Distance_to_centroid = tp.pos$distances, Group = NMDS.data$time.pos)

perm4 <- ggplot(data=tp.pos.permdisp.data,aes(x=Group,y=Distance_to_centroid, fill = Group))+
  scale_fill_manual(values = c("lightgrey", "lightgrey", "lightblue", "lightblue"))+
  geom_boxplot()+
  theme_bw() +
  theme(axis.text.x=element_text(margin = margin(b = 10),colour = "black", size = 12, family = "Arial"),
        axis.text.y=element_text(margin = margin(l = 10),colour = "black",size=12,family="Arial"),
        axis.text=element_text(size=14,family="Arial"),
        axis.title = element_text(size=14,family="Arial", face = "bold"),
        panel.grid = element_blank(),
        plot.title=element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold", color = "black"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold"),
        legend.position = "none",
        axis.ticks.length = unit(.2, "cm")) +
  labs(x = "Time point x Position", y = "Distance to centroid", tag = "b)")


NMDS.data <- NMDS.data %>%
  mutate(interaction =paste(Treatment, Time_point, Position, sep = "_"))

NMDS.data$interaction <- factor(NMDS.data$interaction)

NMDS.data$interaction <- relevel(NMDS.data$interaction, ref="Control_Week6_Downstream")
NMDS.data$interaction <- relevel(NMDS.data$interaction, ref="Control_Week6_Upstream")
NMDS.data$interaction <- relevel(NMDS.data$interaction, ref="Control_Week(-1)_Downstream")
NMDS.data$interaction <- relevel(NMDS.data$interaction, ref="Control_Week(-1)_Upstream")
levels(NMDS.data$interaction)

###################################################################################
## PERMDISP for interaction between treatment, time point and position
###################################################################################

## PERMDISP
interaction <- betadisper(traits_dist, NMDS.data$interaction)
anova(interaction)
permutest(interaction, pairwise = TRUE)
TukeyHSD(interaction)

## Boxplot for PERMDISP
interaction.permdisp.data <- data.frame(Distance_to_centroid = interaction$distances, Group = NMDS.data$interaction)

## Convert character variable to factor variable

interaction.permdisp.data$Group <- factor(interaction.permdisp.data$Group)

## Reorder the levels of factor variable

interaction.permdisp.data$Group <- relevel(interaction.permdisp.data$Group, ref="Crayfish + ALAN_Week6_Downstream")
interaction.permdisp.data$Group <- relevel(interaction.permdisp.data$Group, ref="Crayfish + ALAN_Week6_Upstream")
interaction.permdisp.data$Group <- relevel(interaction.permdisp.data$Group, ref="Crayfish + ALAN_Week(-1)_Downstream")
interaction.permdisp.data$Group <- relevel(interaction.permdisp.data$Group, ref="Crayfish + ALAN_Week(-1)_Upstream")
interaction.permdisp.data$Group <- relevel(interaction.permdisp.data$Group, ref="Crayfish_Week6_Downstream")
interaction.permdisp.data$Group <- relevel(interaction.permdisp.data$Group, ref="Crayfish_Week6_Upstream")
interaction.permdisp.data$Group <- relevel(interaction.permdisp.data$Group, ref="Crayfish_Week(-1)_Downstream")
interaction.permdisp.data$Group <- relevel(interaction.permdisp.data$Group, ref="Crayfish_Week(-1)_Upstream")
interaction.permdisp.data$Group <- relevel(interaction.permdisp.data$Group, ref="ALAN_Week6_Downstream")
interaction.permdisp.data$Group <- relevel(interaction.permdisp.data$Group, ref="ALAN_Week6_Upstream")
interaction.permdisp.data$Group <- relevel(interaction.permdisp.data$Group, ref="ALAN_Week(-1)_Downstream")
interaction.permdisp.data$Group <- relevel(interaction.permdisp.data$Group, ref="ALAN_Week(-1)_Upstream")
interaction.permdisp.data$Group <- relevel(interaction.permdisp.data$Group, ref="Control_Week6_Downstream")
interaction.permdisp.data$Group <- relevel(interaction.permdisp.data$Group, ref="Control_Week6_Upstream")
interaction.permdisp.data$Group <- relevel(interaction.permdisp.data$Group, ref="Control_Week(-1)_Downstream")
interaction.permdisp.data$Group <- relevel(interaction.permdisp.data$Group, ref="Control_Week(-1)_Upstream")
levels(interaction.permdisp.data$Group)


perm5 <- ggplot(data=interaction.permdisp.data,aes(x=Group,y=Distance_to_centroid, fill = Group))+
  scale_fill_manual(values = c("grey", "grey", "grey", "grey",
                               "#009E73", "#009E73", "#009E73", "#009E73",
                               "lightblue", "lightblue", "lightblue", "lightblue",
                               "#D55E00", "#D55E00", "#D55E00", "#D55E00"))+
  geom_boxplot()+
  theme_bw() +
  theme(axis.text.x=element_text(angle = 50, vjust = 1, hjust = 1, margin = margin(b = 10),
                                 colour = c("darkgrey","darkgrey","black", "black",
                                            "#009E73", "#009E73","darkgreen", "darkgreen",
                                            "#56B4E9","#56B4E9","darkblue", "darkblue",
                                            "red", "red","darkred", "darkred"),
                                 size = 12, family = "Arial"),
        axis.text.y=element_text(margin = margin(l = 10),colour = "black",size=12,family="Arial"),
        axis.text=element_text(size=14,family="Arial"),
        axis.title = element_text(size=14,family="Arial", face = "bold"),
        panel.grid = element_blank(),
        plot.title=element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold", color = "black"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold"),
        legend.position = "none") +
  labs(x = "Treatment x Time point x Position", y = "Distance to centroid", tag = "a)")


library(cowplot)
permdisp_all <- plot_grid(perm5, plot_grid(perm4, perm1, perm3, perm2), nrow = 2)
ggsave("permall.png", permdisp_all, width = 16, height = 13, dpi = 1200, device = "png")

#############################################################################################

####################################
## SIMPER
####################################

## SIMPER for Treatment
simper.treatment <- simper(trait.matrix, NMDS.data$Treatment)

simp.tr.trait1 <- simper.treatment$ALAN_Control %>%
  as.data.frame() %>%
  rownames_to_column(var = "Trait") %>%
  filter(cusum < 0.7, p < 0.05, ratio > 0)

simp.tr.trait2 <- simper.treatment$Crayfish_Control %>%
  as.data.frame() %>%
  rownames_to_column(var = "Trait") %>%
  filter(cusum < 0.7, p < 0.05, ratio > 0)

simp.tr.trait3 <- simper.treatment$`Crayfish + ALAN_Control` %>%
  as.data.frame() %>%
  rownames_to_column(var = "Trait") %>%
  filter(cusum < 0.7, p < 0.05, ratio > 0)

## SIMPER for Time point
simper.time <- simper(trait.matrix, NMDS.data$Time_point)
simper.time
simp.time.trait <- simper.time$`Week(-1)_Week6` %>%
  as.data.frame() %>%
  rownames_to_column(var = "Trait") %>%
  filter(cusum < 0.7, p < 0.05, ratio > 0)

## SIMPER for Position
simper.position <- simper(trait.matrix, NMDS.data$Position)
simper.position
simp.position.trait <- simper.position$`Upstream_Downstream` %>%
  as.data.frame() %>%
  rownames_to_column(var = "Trait") %>%
  filter(cusum < 0.7, p < 0.05, ratio > 0)


#############################################################################################
### Uni-variate analysis
##############################################################
## Linear mixed effect model for functional diversity indices 
##############################################################
library(brms)
library(bayesplot)
library(rstanarm)
library(lme4)
library(lmerTest)
library(tidyverse)
library(ggridges)
library(shinystan)
library(tidybayes)
library(ggmcmc)
library(sjstats)
library(tidybayes)
library(showtext)
library(ggthemes)

FD_indices <- FD_indices %>%
  separate(Sample_name, into = c("Channel", "Treatment", "Position", "Time_point"), sep = "_", remove = FALSE)

##  Rename factor levels
FD_indices$Time_point <- recode(FD_indices$Time_point, "week -1" = "Before Treatment", "week" = "After Treatment")

## Converting character variables to factor variables
FD_indices$Time_point <- factor(FD_indices$Time_point)
FD_indices$Treatment <- factor(FD_indices$Treatment)
FD_indices$Position <- factor(FD_indices$Position)

## Reorders the levels of factor variable
FD_indices$Treatment <- relevel(FD_indices$Treatment, ref="Control")
levels(FD_indices$Treatment)

FD_indices$Time_point <- relevel(FD_indices$Time_point, ref="Before Treatment")
levels(FD_indices$Time_point)

FD_indices$Position <- relevel(FD_indices$Position, ref="Upstream")
levels(FD_indices$Position)


############################
## Functional richness
############################

# Checking the range of functional richness
range(FD_indices$Functional_richness)  # range is in between 0 to 0.8 and zero inflated

# Checking the mean value of functional richness value
mean(FD_indices$Functional_richness) # mean value = 0.34

# Calculating standard deviation of functional richness
sd(FD_indices$Functional_richness) # sd = 0.19

 # Checking the number of zero values in functional richness
any(FD_indices$Functional_richness == 0)

## Setting prior

get_prior(sqrt(Functional_richness) ~ Treatment*Time_point + (1|Channel), data = FD_indices, family = zero_inflated_beta())

## Model fitting
rich_model <- brms::brm(Functional_richness ~ Treatment*Time_point*Position + (1|Channel),
                        data = FD_indices, family = zero_inflated_beta(),
                        prior = set_prior("normal(0.34, 0.19)", class = "b"),
                        chains = 3,iter = 3000, warmup = 1500, control = list(adapt_delta = 0.95, max_treedepth = 15), seed = 123)

## Model summary
summary(rich_model)
conditional_effects(rich_model, ask = FALSE)

## Look at the chains and the convergence of the model
plot(rich_model, ask = FALSE)

color_scheme_set("mix-brightblue-teal")
color_scheme_view()

rich.trace <- mcmc_trace(rich_model, pars = c("b_Intercept",
                                "b_TreatmentALAN",
                                "b_TreatmentCrayfish",
                                "b_TreatmentCrayfishPALAN",
                                "b_TreatmentCrayfishPALAN:PositionDownstream",
                                "b_TreatmentCrayfishPALAN:Time_pointAfterTreatment")) +
  theme_bw() + 
  labs(title = "Trace plot (FRich model)", tag = "a)") +
  theme(axis.text.x = element_text(family = "Arial", size = 12, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5, family = "Arial", size = 16, color = "black", face = "bold"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold"))

ggsave("rich_trace.png", rich.trace, width = 12, height = 8, dpi = 800, device = "png")



## Check how well the model fits the data

color_scheme_set("brightblue")
rich_post_check <- pp_check(rich_model, type = "dens_overlay", ndraws = 100) +
  theme_classic() +
  labs(title = "Posterior predictive checks (FRich model)", tag = "a)")+
  theme(axis.text.x = element_text(family = "Arial", size = 12, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5, family = "Arial", size = 16, color = "black", face = "bold"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold"))

rich_post_check

ggsave("rich_post_check.png", rich_post_check, width = 12, height = 8, dpi = 800, device = "png")


#########################################
## Plotting the model conditional effects
#########################################
predicted.data.rich <- FD_indices %>%
  add_predicted_draws(rich_model) 

alpha_max <- 1
alpha_min <- 0.3
alpha_vals <- c(
  seq(alpha_max, alpha_min, length.out = 2), 
  seq(alpha_min, alpha_max, length.out = 2)[-1]
)
alpha_vals

frich <- qplot(interaction(Time_point, Treatment), .prediction, data=predicted.data.rich, geom= "boxplot", fill=Treatment, alpha = Time_point) +
  scale_fill_manual(values = c("grey", "#009E73","#56B4E9", "#E69F00"))+
  facet_wrap(~ Position, nrow = 1) +
  scale_alpha_manual(values=c(0.2,0.7)) +
  guides(alpha=guide_legend(override.aes=list(fill=hcl(c(15,195),200,0,alpha=c(0.2,1))))) +
  theme_bw() +
  labs(x = "Interaction(Treatment x Time point)", 
       y = "Functional richness (predicted)",
       fill = "Before-After Treatment",
       alpha = "Time point",
       tag = "a)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, family = "Arial", size = 12, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, family = "Arial", face = "bold"),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(t = 12),family = "Arial", size = 14, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 16),family = "Arial", size = 14, face = "bold"),
        legend.title = element_text(color = "darkblue",margin = margin(b = 0),family = "Arial", size = 14, face = "bold"),
        legend.text = element_text(size = 13, family = "Arial"),
        legend.position = "top",
        legend.background = element_rect(fill = "transparent"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold", color = "white"))    

ggsave("frich.png", frich, width = 14, height = 8, dpi = 1600, device = "png")
##############################################
## Posterior distribution of model parameters
##############################################
library(bayesplot)
posterior_rich <- as.array(rich_model)
posterior_rich
dim(posterior_rich)

## Central posterior uncertainty intervals

color_scheme_set("viridisA")
rich_pd <- mcmc_areas(posterior_rich, pars = c("b_Intercept",
                                    "b_TreatmentCrayfishPALAN",
                                    "b_TreatmentCrayfish:Time_pointAfterTreatment:PositionDownstream",
                                    "b_TreatmentALAN:Time_pointAfterTreatment:PositionDownstream",
                                    "b_TreatmentCrayfishPALAN:Time_pointAfterTreatment:PositionDownstream"),
           prob = 0.95, # 95% intervals
           prob_outer = 0.99, # 99%
           point_est = "median") + 
  scale_y_discrete(
    labels = c("b_Intercept" = "Intercept",
               "b_TreatmentCrayfishPALAN" = "Crayfish+ALAN",
               "b_TreatmentCrayfish:Time_pointAfterTreatment:PositionDownstream" = "Crayfish x Week6 x Downstream",
               "b_TreatmentALAN:Time_pointAfterTreatment:PositionDownstream" = "ALAN x Week6 x Downstream",
               "b_TreatmentCrayfishPALAN:Time_pointAfterTreatment:PositionDownstream" = "(Crayfish+ALAN) x Week6 x Downstream")) +
  labs(title = "Posterior distribution (FRich model)", tag = "a)") 
rich.post <-rich_pd + theme_classic() +
  theme(axis.text.x = element_text(family = "Arial", size = 12, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black"),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(family = "Arial", size = 16, color = "black", face = "bold"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold")) 


rich.post
### Groupwise comparison
library(emmeans)
emmeans_rich_model <- emmeans(rich_model, ~ Treatment*Time_point*Position)
pairs(emmeans_rich_model, simple = "each")

##############################
## Functional dispersion model
##############################

## Range of trait values
range(FD_indices$Functional_dispersion) # in between 0.001 to 0.17 and zero inflated

## Calculating standard deviation of functional dispersion
sd(FD_indices$Functional_dispersion) # sd value = 0.049

## Mean value of functional dispersion
mean(FD_indices$Functional_dispersion) ## mean value = 0.05 

## Checking the number of zero values in functional dispersion
any(FD_indices$Functional_dispersion == 0)

## Density plot for functional dispersion
ggplot(FD_indices, aes(x = Functional_dispersion)) +
  geom_density(fill = "lightblue", alpha = 0.6) +
  theme_bw() +
  theme(axis.text=element_text(size=12,family="Arial"),
        axis.title = element_text(size=12,family="Arial", face = "bold"),
        legend.text = element_text(size=12,family="Arial", face = "italic"),
        legend.title = element_text(size=12,family="Arial", face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold", color = "darkred"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold"),
        panel.grid = element_blank(),
        legend.position = "none") +
  labs(x = "Functional dispersion", y = "Density", title = "Density plot for functional dispersion")

## Model fitting
disp_model <- brms::brm(Functional_dispersion ~ Treatment*Time_point*Position + (1|Channel),
                        data = FD_indices, family = skew_normal(), # Functional dispersion is left skewed
                        prior = set_prior("normal(0.05, 0.049)", class = "b"),
                        chains = 3,iter = 3000, warmup = 1500, control = list(adapt_delta = 0.95, max_treedepth = 15), seed = 123)

## Model summary
summary(disp_model)

## Look at the chains and the convergence of the model
color_scheme_set("mix-brightblue-teal")
color_scheme_view()
plot(disp_model, ask = FALSE)

disp.trace <- mcmc_trace(disp_model, pars = c("b_Intercept",
                                              "b_PositionDownstream",
                                              "b_TreatmentALAN:PositionDownstream",
                                              "b_TreatmentCrayfish:PositionDownstream",
                                              "b_TreatmentCrayfishPALAN:PositionDownstream",
                                              "b_Time_pointAfterTreatment:PositionDownstream")) +
  theme_bw() +
  labs(title = "Trace plot (FDisp model)", tag = "b)") +
  theme(axis.text.x = element_text(family = "Arial", size = 12, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5, family = "Arial", size = 16, color = "black", face = "bold"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold"))


ggsave("disp_trace.png", disp.trace, width = 12, height = 8, dpi = 800, device = "png")


## Check how well the model fits the data
color_scheme_set("brightblue")
disp_post_check <- pp_check(disp_model, type = "dens_overlay", ndraws = 100) +
  theme_classic() +
  labs(title = "Posterior predictive checks (FDisp model)", tag = "b)")+
  theme(axis.text.x = element_text(family = "Arial", size = 12, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5, family = "Arial", size = 16, color = "black", face = "bold"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold"))


disp_post_check

ggsave("disp_post_check.png", disp_post_check, width = 12, height = 8, dpi = 800, device = "png")

## Plotting the model conditional effects

predicted.data.disp <- FD_indices %>%
  add_predicted_draws(disp_model)

fdisp <- qplot(interaction(Time_point, Treatment), .prediction, data=predicted.data.disp, geom= "boxplot", fill=Treatment, alpha = Time_point) +
  scale_fill_manual(values = c("grey", "#009E73","#56B4E9", "#E69F00"))+
  facet_wrap(~ Position, nrow = 1) +
  scale_alpha_manual(values=c(0.2,0.7)) +
  guides(alpha=guide_legend(override.aes=list(fill=hcl(c(15,195),200,0,alpha=c(0.2,1))))) +
  theme_bw() +
  labs(x = "Interaction(Treatment x Time point)", 
       y = "Functional disperion (predicted)",
       fill = "Before-After Treatment",
       alpha = "Time point",
       tag = "a)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, family = "Arial", size = 12, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, family = "Arial", face = "bold"),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(t = 12),family = "Arial", size = 14, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 16),family = "Arial", size = 14, face = "bold"),
        legend.title = element_text(color = "darkblue",margin = margin(b = 0),family = "Arial", size = 14, face = "bold"),
        legend.text = element_text(size = 13, family = "Arial"),
        legend.position = "top",
        legend.background = element_rect(fill = "transparent"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold", color = "white"))  


ggsave("fdisp.png", fdisp, width = 14, height = 8, dpi = 1600, device = "png")

## Posterior distribution of model parameters
posterior_d <- as.array(disp_model)
dim(posterior_d)

## Central posterior uncertainty intervals

color_scheme_set("viridisA")
disp_pd <- mcmc_areas(posterior_d, pars = c("b_Intercept",
                                               "b_PositionDownstream",
                                               "b_TreatmentCrayfish:Time_pointAfterTreatment:PositionDownstream",
                                               "b_TreatmentALAN:Time_pointAfterTreatment:PositionDownstream",
                                               "b_TreatmentCrayfishPALAN:Time_pointAfterTreatment:PositionDownstream"),
                      prob = 0.95, # 95% intervals
                      prob_outer = 0.99, # 99%
                      point_est = "median") + 
  scale_y_discrete(
    labels = c("b_Intercept" = "Intercept",
               "b_PositionDownstream" = "Downstream",
               "b_TreatmentCrayfish:Time_pointAfterTreatment:PositionDownstream" = "Crayfish x Week6 x Downstream",
               "b_TreatmentALAN:Time_pointAfterTreatment:PositionDownstream" = "ALAN x Week6 x Downstream",
               "b_TreatmentCrayfishPALAN:Time_pointAfterTreatment:PositionDownstream" = "(Crayfish+ALAN) x Week6 x Downstream")) +
  labs(tag = "b)", title = "Posterior distribution (FDisp model)")
disp.post <- disp_pd + theme_classic() +
  theme(axis.text.x = element_text(family = "Arial", size = 12, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black"),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(family = "Arial", size = 16, color = "black", face = "bold"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold")) 


disp.post

## Groupwise comparison for brms model
emmeans_disp_model <- emmeans(disp_model, ~ Treatment*Time_point*Position)
pairs(emmeans_disp_model, simple = "each")

###########################################################################################################################

############################
## Functional evenness model
############################

## Range of trait values
range(FD_indices$Functional_evenness) # in between 0 to 0.94 and zero inflated

## Calculating standard deviation of functional evenness
sd(FD_indices$Functional_evenness) # sd value = 0.21

## Mean value of functional evenness
mean(FD_indices$Functional_evenness) ## mean value = 0.29

## Checking the number of zero values in functional evenness
any(FD_indices$Functional_evenness == 0)

## Model fitting
even_model <- brms::brm(Functional_evenness ~ Treatment*Time_point*Position + (1|Channel),
                        data = FD_indices, family = zero_inflated_beta(),
                        prior = set_prior("normal(0.29, 0.21)", class = "b"),
                        chains = 3,iter = 3000, warmup = 1500, control = list(adapt_delta = 0.95, max_treedepth = 15), seed = 123)

## Model summary
summary(even_model)

## Look at the chains and the convergence of the model
plot(even_model, ask = FALSE)

color_scheme_set("mix-brightblue-teal")
even.trace <- mcmc_trace(even_model, pars = c("b_Intercept",
                                              "b_TreatmentALAN",
                                              "b_TreatmentCrayfish",
                                              "b_TreatmentCrayfishPALAN",
                                              "b_TreatmentCrayfishPALAN:PositionDownstream",
                                              "b_TreatmentCrayfishPALAN:Time_pointAfterTreatment")) +
  theme_bw() + 
  labs(title = "Trace plot (FEve model)", tag = "c)") +
  theme(axis.text.x = element_text(family = "Arial", size = 12, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5, family = "Arial", size = 16, color = "black", face = "bold"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold"))


even.trace

ggsave("even_trace.png", even.trace, width = 12, height = 8, dpi = 800, device = "png")

## Check how well the model fits the data
color_scheme_set("brightblue")
eve_post_check <- pp_check(even_model, type = "dens_overlay", ndraws = 100) +
  theme_classic() +
  labs(title = "Posterior predictive checks (FEven model)", tag = "c)")+
  theme(axis.text.x = element_text(family = "Arial", size = 12, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5, family = "Arial", size = 16, color = "black", face = "bold"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold"))


eve_post_check

ggsave("eve_post_check.png", eve_post_check, width = 12, height = 8, dpi = 800, device = "png")

## Plotting the model conditional effects

predicted.data.even <- FD_indices %>%
  add_predicted_draws(even_model)

feve <- qplot(interaction(Time_point, Treatment), .prediction, data=predicted.data.even, geom= "boxplot", fill=Treatment, alpha = Time_point) +
  scale_fill_manual(values = c("grey", "#009E73","#56B4E9", "#E69F00"))+
  facet_wrap(~ Position, nrow = 1) +
  scale_alpha_manual(values=c(0.2,0.7)) +
  guides(alpha=guide_legend(override.aes=list(fill=hcl(c(15,195),200,0,alpha=c(0.2,1))))) +
  theme_bw() +
  labs(x = "Interaction(Treatment x Time point)", 
       y = "Functional evenness (predicted)",
       fill = "Before-After Treatment",
       alpha = "Time point",
       tag = "a)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, family = "Arial", size = 12, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, family = "Arial", face = "bold"),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(t = 12),family = "Arial", size = 14, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 16),family = "Arial", size = 14, face = "bold"),
        legend.title = element_text(color = "darkblue",margin = margin(b = 0),family = "Arial", size = 14, face = "bold"),
        legend.text = element_text(size = 13, family = "Arial"),
        legend.position = "top",
        legend.background = element_rect(fill = "transparent"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold", color = "white"))    

ggsave("feve.png", feve, width = 14, height = 8, dpi = 1600, device = "png")

## Posterior distribution of model parameters
posterior_e <- as.array(even_model)
dim(posterior_e)

## Central posterior uncertainty intervals

color_scheme_set("viridisA")
even_pd <- mcmc_areas(posterior_e, pars = c("b_Intercept",
                                               "b_TreatmentCrayfish:Time_pointAfterTreatment:PositionDownstream",
                                               "b_TreatmentALAN:Time_pointAfterTreatment:PositionDownstream",
                                               "b_TreatmentCrayfishPALAN:Time_pointAfterTreatment:PositionDownstream"),
                      prob = 0.95, # 95% intervals
                      prob_outer = 0.99, # 99%
                      point_est = "median") + 
  scale_y_discrete(
    labels = c("b_Intercept" = "Intercept",
               "b_TreatmentCrayfish:Time_pointAfterTreatment:PositionDownstream" = "Crayfish x Week6 x Downstream",
               "b_TreatmentALAN:Time_pointAfterTreatment:PositionDownstream" = "ALAN x Week6 x Downstream",
               "b_TreatmentCrayfishPALAN:Time_pointAfterTreatment:PositionDownstream" = "(Crayfish+ALAN) x Week6 x Downstream")) +
  labs(tag = "c)", title = "Posterior distribution (FEve model)")
eve.post <- even_pd + theme_classic() +
  theme(axis.text.x = element_text(family = "Arial", size = 12, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black"),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(family = "Arial", size = 16, color = "black", face = "bold"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold")) 

eve.post

## Groupwise comparison for brms model
emmeans_even_model <- emmeans(even_model, ~ Treatment*Time_point*Position)
pairs(emmeans_even_model, simple = "each")

###########################################################################################################################

############################
## raoQ model
############################
## density plot for raoQ
ggplot(FD_indices, aes(x = raoQ)) +
  geom_density(fill = "lightblue", color = "red") +
  theme_bw() +
  labs(x = "Rao's Q", y = "Density", title = "Density plot for Rao's Q")

## Range of trait values
range(FD_indices$raoQ) # in between 0.003 to 0.55 and left skewed

## Calculating standard deviation of raoQ
sd(FD_indices$raoQ) # sd value = 0.149

## Mean value of raoQ
mean(FD_indices$raoQ) ## mean value = 0.148

## Checking the number of zero values in raoQ
any(FD_indices$raoQ == 0)

## Model fitting
raoQ_model <- brms::brm(raoQ ~ Treatment*Time_point*Position + (1|Channel),
                        data = FD_indices, family = skew_normal(), # raoQ is left skewed
                        prior = set_prior("normal(0.148, 0.149)", class = "b"),
                        chains = 3,iter = 3000, warmup = 1500, control = list(adapt_delta = 0.95, max_treedepth = 15), seed = 123)

## Model summary
summary(raoQ_model)

## Look at the chains and the convergence of the model
plot(raoQ_model, ask = FALSE) 

color_scheme_set("mix-brightblue-teal")

raoQ.trace <- mcmc_trace(raoQ_model, pars = c("b_Intercept",
                                              "b_PositionDownstream",
                                              "b_TreatmentALAN:PositionDownstream",
                                              "b_TreatmentCrayfish:PositionDownstream",
                                              "b_TreatmentCrayfishPALAN:PositionDownstream",
                                              "b_Time_pointAfterTreatment:PositionDownstream")) +
  theme_bw() + 
  labs(title = "Trace plot (RaoQ model)", tag = "d)") +
  theme(axis.text.x = element_text(family = "Arial", size = 12, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5, family = "Arial", size = 16, color = "black", face = "bold"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold"))


raoQ.trace

ggsave("raoQ_trace.png", raoQ.trace, width = 12, height = 8, dpi = 800, device = "png")

## Check how well the model fits the data
color_scheme_set("brightblue")
raoq_post_check <- pp_check(raoQ_model, type = "dens_overlay", ndraws = 100) +
  theme_classic() +
  labs(title = "Posterior predictive checks (RaoQ model)", tag = "d)")+
  theme(axis.text.x = element_text(family = "Arial", size = 12, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5, family = "Arial", size = 16, color = "black", face = "bold"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold"))


raoq_post_check

library(gridExtra)
predict.plot <- grid.arrange(rich_post_check, disp_post_check, eve_post_check, raoq_post_check, ncol = 2)

ggsave("predict.plot.png", predict.plot, width = 17, height = 13, dpi = 1200, device = "png")

## Plotting the model conditional effects

predicted.data.raoQ <- FD_indices %>%
  add_predicted_draws(raoQ_model)

raoq <- qplot(interaction(Time_point, Treatment), .prediction, data=predicted.data.raoQ, geom= "boxplot", fill=Treatment, alpha = Time_point) +
  scale_fill_manual(values = c("grey", "#009E73","#56B4E9", "#E69F00"))+
  facet_wrap(~ Position, nrow = 1) +
  scale_alpha_manual(values=c(0.2,0.7)) +
  guides(alpha=guide_legend(override.aes=list(fill=hcl(c(15,195),200,0,alpha=c(0.2,1))))) +
  theme_bw() +
  labs(x = "Interaction(Treatment x Time point)", 
       y = "Rao's Q (predicted)",
       fill = "Before-After Treatment",
       alpha = "Time point",
       tag = "a)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, family = "Arial", size = 12, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, family = "Arial", face = "bold"),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(t = 12),family = "Arial", size = 14, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 16),family = "Arial", size = 14, face = "bold"),
        legend.title = element_text(color = "darkblue",margin = margin(b = 0),family = "Arial", size = 14, face = "bold"),
        legend.text = element_text(size = 13, family = "Arial"),
        legend.position = "top",
        legend.background = element_rect(fill = "transparent"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold", color = "white"))    

ggsave("raoQ.png", raoq, width = 14, height = 8, dpi = 1600, device = "png")

## Posterior distribution of model parameters
posterior_raoQ <- as.array(raoQ_model)
dim(posterior_raoQ)

## Central posterior uncertainty intervals

color_scheme_set("viridisA")
raoQ_pd <- mcmc_areas(posterior_raoQ, pars = c("b_Intercept",
                                               "b_PositionDownstream",
                                               "b_TreatmentCrayfish:Time_pointAfterTreatment:PositionDownstream",
                                               "b_TreatmentALAN:Time_pointAfterTreatment:PositionDownstream",
                                               "b_TreatmentCrayfishPALAN:Time_pointAfterTreatment:PositionDownstream"),
                      prob = 0.95, # 95% intervals
                      prob_outer = 0.99, # 99%
                      point_est = "median") + 
  scale_y_discrete(
    labels = c("b_Intercept" = "Intercept",
               "b_PositionDownstream" = "Downstream",
               "b_TreatmentCrayfish:Time_pointAfterTreatment:PositionDownstream" = "Crayfish x Week6 x Downstream",
               "b_TreatmentALAN:Time_pointAfterTreatment:PositionDownstream" = "ALAN x Week6 x Downstream",
               "b_TreatmentCrayfishPALAN:Time_pointAfterTreatment:PositionDownstream" = "(Crayfish+ALAN) x Week6 x Downstream")) +
  labs(tag = "d)", title = "Posterior distribution (RaoQ model)")
raoq.post <- raoQ_pd + theme_classic() +
  theme(axis.text.x = element_text(family = "Arial", size = 12, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black"),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(family = "Arial", size = 16, color = "black", face = "bold"),
        plot.tag = element_text(size = 16, family = "Arial", face = "bold")) 


raoq.post


ggsave("rich.post.png", rich.post, width = 12, height = 8, dpi = 800, device = "png")
ggsave("disp.post.png", disp.post, width = 12, height = 8, dpi = 800, device = "png")
ggsave("even.post.png", eve.post, width = 12, height = 8, dpi = 800, device = "png")
ggsave("raoq.post.png", raoq.post, width = 12, height = 8, dpi = 800, device = "png")


## Groupwise comparison for brms model
emmeans_raoQ_model <- emmeans(raoQ_model, ~ Treatment*Time_point*Position)
pairs(emmeans_raoQ_model, simple = "each")


##########################################################################################################################
## Thank You
##########################################################################################################################
