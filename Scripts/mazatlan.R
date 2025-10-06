##################################################################
######### Analisis de Core Microbiome - Lagunas Costeras #########
##################################################################



############ Elaboró JORGE ALEXANDER ROJAS-VARGAS ######################


## Caso de estudio: Mazatlan

# Location of work directory
setwd("/Users/jorgerv/Documents/lagunas_costeras/mazatlan")




###########################
##### ABUNDANCE PLOTS #####
###########################

## Run Pavian to visualize the Kraken reports and download count and abundance tables
library(pavian)
pavian::runApp(port=5000)

## Read the count tables
count_class <- read.csv("A_arc_bac_count_class.csv", row.names = 1)
count_genus <- read.csv("A_arc_bac_count_genus.csv", row.names = 1)
count_species <- read.csv("A_arc_bac_count_species.csv", row.names = 1)

## Convierto NAs en ceros
count_class[is.na(count_class)] <- 0
count_genus[is.na(count_genus)] <- 0
count_species[is.na(count_species)] <- 0

save(count_class, file = "count_arc_bac_class.RData")
save(count_genus, file = "count_arc_bac_genus.RData")
save(count_species, file = "count_arc_bac_species.RData")

## Function to process relative abundance (RA) tables
process_abundance_table <- function(count_table) {
  abundance_table <- apply(count_table, 2, function(i) i / sum(i))
  abundance_table <- as.data.frame(abundance_table)
  abundance_table <- abundance_table[rowSums(abundance_table) > 0, , drop = FALSE]
  return(abundance_table)
}

## Apply the abundance function to count tables
abund_arc_bac_class <- process_abundance_table(count_class)*100
abund_arc_bac_genus <- process_abundance_table(count_genus)*100
abund_arc_bac_species <- process_abundance_table(count_species)*100

## Eliminate taxa with RA < 0.00001%
abund_arc_bac_class <- abund_arc_bac_class[apply(abund_arc_bac_class, 1, function(x) any(x > 0.00001)), , drop = FALSE]
abund_arc_bac_genus <- abund_arc_bac_genus[apply(abund_arc_bac_genus, 1, function(x) any(x > 0.00001)), , drop = FALSE]
abund_arc_bac_species <- abund_arc_bac_species[apply(abund_arc_bac_species, 1, function(x) any(x > 0.00001)), , drop = FALSE]

## Re-calculate RAs
abund_arc_bac_class <- process_abundance_table(abund_arc_bac_class)*100
abund_arc_bac_genus <- process_abundance_table(abund_arc_bac_genus)*100
abund_arc_bac_species <- process_abundance_table(abund_arc_bac_species)*100

save(abund_arc_bac_class, file = "abund_arc_bac_class.RData")
save(abund_arc_bac_genus, file = "abund_arc_bac_genus.RData")
save(abund_arc_bac_species, file = "abund_arc_bac_species.RData")



## Read metadata
metadata <- read.csv("metadata_mazatlan.csv")
save(metadata, file = "metadata_ZMO.RData")

require(ggplot2)
require(reshape2)


###### CLASS LEVEL ######

load("abund_arc_bac_class.RData")
load("abund_arc_bac_genus.RData")
load("metadata_ZMO.RData")

## Se eliminarán las abundancias menores a 1% para el barstack plot
abund_arc_bac_class_stack <- abund_arc_bac_class[apply(abund_arc_bac_class, 1, function(x) any(x >= 1, na.rm = TRUE)), ]
abund_arc_bac_genus_stack <- abund_arc_bac_genus[apply(abund_arc_bac_genus, 1, function(x) any(x >= 1, na.rm = TRUE)), ]


## Añado fila de "Other" para recuperar la taxa filtrada
column_sums <- colSums(abund_arc_bac_class_stack)
diff_values <- 100 - column_sums
other_row <- as.data.frame(t(diff_values))
rownames(other_row) <- "Other"
abund_arc_bac_class_stack <- rbind(other_row, abund_arc_bac_class_stack)
new_column_names <- metadata$ID_2
colnames(abund_arc_bac_class_stack) <- new_column_names

## Organizo las clases segun su abundancia media de menor a mayor
abund_arc_bac_class_stack_filt <- apply(abund_arc_bac_class_stack, 1, median)
top_class <- names(sort(abund_arc_bac_class_stack_filt, decreasing = FALSE)[1:20])
top_class <- c("Other", setdiff(top_class, "Other")) # Reorganizar para que "Other" quede de primero
top_class

## Renombrar las columnas de abund_arc_bac_class_stack
#colnames(abund_arc_bac_class_stack) <- metadata$ID_2[match(colnames(abund_arc_bac_class_stack), metadata$Sample)]

## Ordenar las columnas alfabéticamente
abund_arc_bac_class_stack <- abund_arc_bac_class_stack[, order(colnames(abund_arc_bac_class_stack))]


## Transformo tabla para graficar
dat_m <- melt(data.matrix(abund_arc_bac_class_stack))
colnames(dat_m)<-c("Class", "Sample", "Abundance")

## Agrego metadatos
Sector <- metadata$Sector[match(dat_m$Sample, metadata$ID_2)]
Season <- metadata$Season[match(dat_m$Sample, metadata$ID_2)]
dat_m$Sector <- Sector
dat_m$Season <- Season

## Ordeno las Taxa segun el orden de top_class
dat_m$Class <- factor(dat_m$Class, levels = top_class)

## Define color palette for taxa
my_colors <- c(
  "Other" = "#808080",
  "Thermoleophilia" = "#2d1650",
  "Acidimicrobiia" = "#5a2ca0",
  "Ignavibacteria" = "#8d5fd3",
  "Dehalococcoidia" = "#c6afe9",
  "Cytophagia" = "#162d50",
  "Spirochaetia" = "#214478",
  "Bacteroidia" = "#3771c8",
  "Flavobacteriia" = "#5f8dd3",
  "Gemmatimonadetes" = "#afc6e9",
  "Bacilli" = "#806600",
  "Planctomycetia" = "#d4aa00",
  "Nitrospira" = "#ffd42a",
  "Clostridia" = "#ffe680",
  "Anaerolineae" = "#fff6d5",       
  "Betaproteobacteria" = "#445016",
  "Actinobacteria" = "#677821",
  "Alphaproteobacteria" = "#89a02c",
  "Gammaproteobacteria" = "#abc837",
  "Deltaproteobacteria" = "#cdde87"
)

svglite("Abundance_bacteria_class_ZMO.svg", width=9, height=8)
ggplot(dat_m, aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") + xlab("Samples") + ylab("Relative Abundance (%)") + theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 18), axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20, face = "bold"), 
    axis.title.y = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16), 
    legend.title = element_text(size = 18, face = "bold"),
    strip.text = element_text(size = 14, face = "bold")
  ) +  
  #theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) + 
  #ggtitle("Prokaryotic class level") +
  scale_fill_manual(
    values = my_colors,
    guide = guide_legend(nrow = 29)
  ) +
  #facet_grid(~ Sector, scales = "free", space = "free")
  facet_grid(~ Season, scales = "free", space = "free")
dev.off()





#######################################
##### ANALYSIS OF ALPHA DIVERSITY #####
#######################################


###### GENUS LEVEL ######


count_genus_arc_bac <- count_genus

## Phyloseq object
taxa_table <- row.names(count_genus_arc_bac)
setdiff(rownames(count_genus_arc_bac),taxa_table)
taxa_table =  as.matrix(taxa_table)
colnames(taxa_table) <- "Genus"
rownames(taxa_table) <- taxa_table
OTU_taxa <- otu_table(count_genus_arc_bac, taxa_are_rows = TRUE)
TAX_taxa <- tax_table(taxa_table)
sampledata <- sample_data(metadata)
rownames(sampledata) <- metadata$Sample
ps_count_genus = phyloseq(OTU_taxa, TAX_taxa, sampledata)
# Keep only taxa with positive sums
ps_count_genus <- prune_taxa(taxa_sums(ps_count_genus)>0, ps_count_genus)
print(ps_count_genus)
save(ps_count_genus, file = "ps_count_genus_ZMO.RData")


## Start the analysis

load("ps_count_genus_ZMO.RData")

plot_richness(ps_count_genus)

a_my_comparisons <- list(c("La_Nina", "El_Nino"), 
                         c("A_Open_ocean","B_Ocean_slope","C_Coastal_slope", "D_Coast"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","ns"))

ad_sector <- plot_richness(ps_count_genus, x="Sector", measures = c("Observed","Shannon"), color = "Sector") + theme_bw() + 
  geom_boxplot(lwd = 0.7, alpha = 0.6) + scale_color_manual(values = c("#007ea7","#6a994e","#c8ab37ff")) 

ad_season <- plot_richness(ps_count_genus, x="Season", measures = c("Observed","Shannon"), color = "Season") + theme_bw() + 
  geom_boxplot(lwd = 0.7, alpha = 0.6)+ scale_color_manual(values = c("#c83737","#2c5aa0"))

svglite("Alpha_Diversity_genus_ZMO.svg", width=5, height=6)
ad_sector / ad_season
dev.off()


library(phyloseq)
library(dplyr)
library(ggpubr)
library(rstatix)

# Extraer la matriz de riqueza
alpha_df <- estimate_richness(ps_count_genus, measures = c("Observed", "Shannon"))

# Agregar metadata de los grupos
metadata <- data.frame(sample_data(ps_count_genus))

# Unir riqueza con metadata
alpha_df <- cbind(metadata, alpha_df)

kruskal_sector_observed <- kruskal_test(Observed ~ Sector, data = alpha_df)
kruskal_sector_shannon <- kruskal_test(Shannon ~ Sector, data = alpha_df)

# Mostrar resultados
kruskal_sector_observed
kruskal_sector_shannon

kruskal_season_observed <- kruskal_test(Observed ~ Season, data = alpha_df)
kruskal_season_shannon <- kruskal_test(Shannon ~ Season, data = alpha_df)

# Mostrar resultados
kruskal_season_observed
kruskal_season_shannon



###### SPECIES LEVEL ######


count_species_arc_bac <- count_species

## Phyloseq object
taxa_table <- row.names(count_species_arc_bac)
setdiff(rownames(count_species_arc_bac),taxa_table)
taxa_table =  as.matrix(taxa_table)
colnames(taxa_table) <- "species"
rownames(taxa_table) <- taxa_table
OTU_taxa <- otu_table(count_species_arc_bac, taxa_are_rows = TRUE)
TAX_taxa <- tax_table(taxa_table)
sampledata <- sample_data(metadata)
rownames(sampledata) <- metadata$Sample
ps_count_species = phyloseq(OTU_taxa, TAX_taxa, sampledata)
# Keep only taxa with positive sums
ps_count_species <- prune_taxa(taxa_sums(ps_count_species)>0, ps_count_species)
print(ps_count_species)
save(ps_count_species, file = "ps_count_species_ZMO.RData")



## Start the analysis

load("ps_count_species_ZMO.RData")

plot_richness(ps_count_species)

a_my_comparisons <- list(c("La_Nina", "El_Nino"), 
                         c("A_Open_ocean","B_Ocean_slope","C_Coastal_slope", "D_Coast"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","ns"))

ad_sector <- plot_richness(ps_count_species, x="Sector", measures = c("Observed","Shannon"), color = "Sector") + theme_bw() + 
  geom_boxplot(lwd = 0.7, alpha = 0.6) + scale_color_manual(values = c("#007ea7","#6a994e","#c8ab37ff")) 

ad_season <- plot_richness(ps_count_species, x="Season", measures = c("Observed","Shannon"), color = "Season") + theme_bw() + 
  geom_boxplot(lwd = 0.7, alpha = 0.6)+ scale_color_manual(values = c("#c83737","#2c5aa0"))

svglite("Alpha_Diversity_species_ZMO.svg", width=5, height=6)
ad_sector / ad_season
dev.off()


library(phyloseq)
library(dplyr)
library(ggpubr)
library(rstatix)

# Extraer la matriz de riqueza
alpha_df <- estimate_richness(ps_count_species, measures = c("Observed", "Shannon"))

# Agregar metadata de los grupos
metadata <- data.frame(sample_data(ps_count_species))

# Unir riqueza con metadata
alpha_df <- cbind(metadata, alpha_df)

kruskal_sector_observed <- kruskal_test(Observed ~ Sector, data = alpha_df)
kruskal_sector_shannon <- kruskal_test(Shannon ~ Sector, data = alpha_df)

# Mostrar resultados
kruskal_sector_observed
kruskal_sector_shannon

kruskal_season_observed <- kruskal_test(Observed ~ Season, data = alpha_df)
kruskal_season_shannon <- kruskal_test(Shannon ~ Season, data = alpha_df)

# Mostrar resultados
kruskal_season_observed
kruskal_season_shannon




######################################
##### ANALYSIS OF BETA DIVERSITY #####
######################################


## Load data
load("abund_arc_bac_class.RData")
load("abund_arc_bac_genus.RData")
load("abund_arc_bac_species.RData")
load("metadata_ZMO.RData")



#### AT GENUS LEVEL ####


## ------------------------------------------
##
## USING BRAY-CURTIS DISSIMILARITY
##
## ------------------------------------------

# Se pretende capturar las diferencias en la abundancia relativa y la presencia/ausencia de especies,
# y no en los efectos de la naturaleza composicional de los datos (normalización de las proporciones)
# ni se tienen en cuenta los efectos de la escala.
# La fórmula básicamente mide cuánto se diferencian las composiciones de las muestras en términos de 
# las proporciones relativas de cada taxón, lo que significa que las diferencias en las abundancias 
# relativas de los taxones afectan directamente el valor de la distancia de Bray-Curtis.


## Calculate Bray-Curtis distance
bray_curtis_dist_genus <- vegan::vegdist(t(abund_arc_bac_genus), method = "bray")

bray_curtis_pcoa_genus <- ecodist::pco(bray_curtis_dist_genus)
bray_curtis_pcoa_df_genus <- data.frame(pcoa1 = bray_curtis_pcoa_genus$vectors[,1], 
                                        pcoa2 = bray_curtis_pcoa_genus$vectors[,2])
dim(bray_curtis_pcoa_df_genus)
dim(metadata)

bray_curtis_pcoa_df_genus <- cbind(bray_curtis_pcoa_df_genus, metadata)

# Extraer varianza explicada por cada eje
eig_vals <- bray_curtis_pcoa_genus$values / sum(bray_curtis_pcoa_genus$values) * 100
pc1_var <- round(eig_vals[1], 2)  # % explicado por PC1
pc2_var <- round(eig_vals[2], 2)  # % explicado por PC2

# Definir colores
season_colors <- c("El_Nino" = "#c83737", "La_Nina" = "#2c5aa0")
sector_colors <- c("A_Open" = "#007ea7", "B_Intermedia" = "#6a994e", "C_Coast" = "#c8ab37ff")

#"A_Open_ocean","B_Ocean_slope","C_Coastal_slope", "D_Coast"

# Plot con ajustes
p1 <- ggplot(bray_curtis_pcoa_df_genus, aes(x = pcoa1, y = pcoa2, color = Season, shape = Sector)) +
  theme_bw() +
  geom_point(size = 4) +  # Ajustar tamaño de los puntos
  stat_ellipse(aes(x = pcoa1, y = pcoa2, group = Season, color = Season), 
               geom = "path", linewidth = 0.5) +  # Agrupar solo por Season
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"), 
    y = paste0("PC2 (", pc2_var, "%)"),
    title = "Bray-Curtis PCoA Genus - Season"
  ) + 
  theme(
    title = element_text(size = 14, face = "bold", hjust = 0.5), 
    axis.title.x = element_text(size = 12, face = "bold"), 
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12), 
    axis.text.y = element_text(size = 12), 
    legend.text = element_text(size = 12)
  )

print(p1)


p2 <- ggplot(bray_curtis_pcoa_df_genus, aes(x = pcoa1, y = pcoa2, color = Season)) +
  theme_bw() +
  geom_point(size = 4) +  # Ajustar tamaño de los puntos
  stat_ellipse(aes(x = pcoa1, y = pcoa2, group = Season, color = Season), 
               geom = "path", linewidth = 0.5
               #, alpha = 0.2
               ) +
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"), 
    y = paste0("PC2 (", pc2_var, "%)"),
    title = "Bray-Curtis PCoA Genus - Season"
  ) + 
  theme(
    title = element_text(size = 14, face = "bold", hjust = 0.5), 
    axis.title.x = element_text(size = 12, face = "bold"), 
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12), 
    axis.text.y = element_text(size = 12), 
    legend.text = element_text(size = 12)
  )

print(p2)


## Statistics

# Asegurar que los metadatos tienen las variables de interés
metadata <- data.frame(mtdt_filt) # Convertir a data.frame si es necesario

# PERMANOVA para Cohort
set.seed(123)
permanova_cohort <- adonis2(bray_curtis_dist_genus ~ Cohort, data = metadata, permutations = 9999)
print(permanova_cohort)
# global p-value 0.022



svglite("Beta_Diversity_genus_ZMO_opc1.svg", width=6, height=4)
p1
dev.off()



## ------------------------------------------
##
## USING AITCHISON DISTANCE
##
## ------------------------------------------


# Para comparar las proporciones relativas de especies en lugar de las cantidades absolutas, es decir,
# la composición relativa de las comunidades (cómo cambian las proporciones de taxones en lugar de sus abundancias)

## ------------------------------------------
## Build the phyloseq object
## ------------------------------------------

#Set a seed for the purposes of reproducibility in this document.
set.seed(1)

#Perform a CLR transformation - #We are imputing zeroes using the 'const' method 
#Essentially, we replace zeroes with 65% of the next lowest value - see Lubbe et al 2021.
count_genus.clr <- clr_c(count_genus)

#ps object of count table
#Used for Aitchison distance
taxa_table <- as.matrix(row.names(count_genus.clr))
colnames(taxa_table) <- "Genus"
rownames(taxa_table) <- taxa_table
OTU_taxa <- otu_table(count_genus.clr, taxa_are_rows = TRUE)
TAX_taxa <- tax_table(taxa_table)
sampledata <- sample_data(metadata)
rownames(sampledata) <- metadata$Sample

ps_count_genus.clr = phyloseq(OTU_taxa, TAX_taxa, sampledata)
save(ps_count_genus.clr, file = "ps_count_genus.clr.RData")

## ------------------------------------------
## Calculate Aitchison distance
## ------------------------------------------

load("count_PEDS_genus.RData")

#ps_c_clr <- microbiome::transform(ps_count_class, "clr")
ps_c_clr <- microbiome::transform(ps_count_genus, "clr")
#ps_c_clr <- microbiome::transform(ps_count_species, "clr")

#Generate Aitchison distance matrix
clr_c_matrix <- phyloseq::distance(ps_c_clr, method = "euclidean")

### We choose the Aitchison distance using CLR transformed data at GENUS level
#Generate distance matrix using Aitchison distances
dist_matrix <- phyloseq::distance(ps_count_genus.clr, method = "euclidean")

## https://microbiome.github.io/course_2021_radboud/beta-diversity.html
# Does principal coordinate analysis
euclidean_genus_pcoa <- ecodist::pco(dist_matrix)
# Creates a data frame from principal coordinates
euclidean_genus_pcoa_df <- data.frame(
  pcoa1 = euclidean_genus_pcoa$vectors[,1], 
  pcoa2 = euclidean_genus_pcoa$vectors[,2])
# Creates a plot
euclidean_genus_plot <- ggplot(data = euclidean_genus_pcoa_df,
                                aes(x=pcoa1, y=pcoa2)) +
  geom_point() +
  labs(x = "PC1",
       y = "PC2",
       title = "Aitchison distances at Genus level") +  
  theme(title = element_text(size = 12)) # makes titles smaller

euclidean_genus_plot
#####


## ------------------------------------------
## Dispersion test with Season
## ------------------------------------------

#Calculate the beta dispersion
dispr <- vegan::betadisper(dist_matrix, phyloseq::sample_data(ps_count_genus.clr)$Season)

#Check the results of the dispersion test
dispr

#PERMANOVA test to dispersion
#Set a seed for the purposes of reproducibility
set.seed(123)
vegan::adonis2(dist(dispr$distances) ~dispr$group) # Global p-value 0.179

#Pairwise PERMANOVA test with fdr correction of dispersion
set.seed(123)
pairwise_test <- vegan::permutest(dispr, pairwise = TRUE)
pairwise_pvalues <- pairwise_test$pairwise$permuted 
pairwise_pvalues
pairwise_pvalues_fdr <- p.adjust(pairwise_pvalues, method = "fdr")
pairwise_pvalues_fdr 
# q = 0.893

# PERMANOVA test to microbiota beta diversity (communities structure)
set.seed(123)
vegan::adonis2(dist_matrix ~ phyloseq::sample_data(ps_count_genus.clr)$Season, data = metadata, permutations = 999)
# global p-value 0.079

# Extraer varianza explicada por cada eje
eig_vals <- dispr$eig / sum(dispr$eig) * 100
pc1_var <- round(eig_vals[1], 2)  # % explicado por PC1
pc2_var <- round(eig_vals[2], 2)  # % explicado por PC2

#PCA and dispersion plot
#svglite("disp_pca_Season_genus_ZMO.svg", width=4, height=4)
plot(dispr, main = "", sub = "",
     xlab = paste0("PC1 (", pc1_var, "%)"), 
     ylab = paste0("PC2 (", pc2_var, "%)"),
     col = c("El_Nino" = "#c83737", "La_Nina" = "#2c5aa0"), 
     cex = 1.5)
#dev.off()

## ------------------------------------------
## Dispersion test with Sector
## ------------------------------------------

#Calculate the beta dispersion
dispr <- vegan::betadisper(dist_matrix, phyloseq::sample_data(ps_count_genus.clr)$Sector)

#Check the results of the dispersion test
dispr

#PERMANOVA test to dispersion
#Set a seed for the purposes of reproducibility
set.seed(123)
vegan::adonis2(dist(dispr$distances) ~dispr$group) # Global p-value 0.191

#Pairwise PERMANOVA test with fdr correction of dispersion
set.seed(123)
pairwise_test <- vegan::permutest(dispr, pairwise = TRUE)
pairwise_pvalues <- pairwise_test$pairwise$permuted 
pairwise_pvalues
pairwise_pvalues_fdr <- p.adjust(pairwise_pvalues, method = "fdr")
pairwise_pvalues_fdr 
# A_Open-C_Coast q = 0.012

# PERMANOVA test to microbiota beta diversity (communities structure)
set.seed(123)
vegan::adonis2(dist_matrix ~ phyloseq::sample_data(ps_count_genus.clr)$Sector, data = metadata, permutations = 999)
# global p-value 0.079

# Extraer varianza explicada por cada eje
eig_vals <- dispr$eig / sum(dispr$eig) * 100
pc1_var <- round(eig_vals[1], 2)  # % explicado por PC1
pc2_var <- round(eig_vals[2], 2)  # % explicado por PC2

#PCA and dispersion plot
svglite("disp_pca_Sector_genus_ZMO.svg", width=4, height=4)
plot(dispr, main = "", sub = "",
     xlab = paste0("PC1 (", pc1_var, "%)"), 
     ylab = paste0("PC2 (", pc2_var, "%)"),
     col = c("A_Open" = "#007ea7", "B_Intermedia" = "#6a994e", "C_Coast" = "#c8ab37ff"), 
     cex = 1.5)
dev.off()

# "A_Open_ocean","B_Ocean_slope","C_Coastal_slope", "D_Coast"


#### AT SPECIES LEVEL ####


## ------------------------------------------
##
## USING BRAY-CURTIS DISSIMILARITY
##
## ------------------------------------------

# Se pretende capturar las diferencias en la abundancia relativa y la presencia/ausencia de especies,
# y no en los efectos de la naturaleza composicional de los datos (normalización de las proporciones)
# ni se tienen en cuenta los efectos de la escala.
# La fórmula básicamente mide cuánto se diferencian las composiciones de las muestras en términos de 
# las proporciones relativas de cada taxón, lo que significa que las diferencias en las abundancias 
# relativas de los taxones afectan directamente el valor de la distancia de Bray-Curtis.


## Calculate Bray-Curtis distance
bray_curtis_dist_species <- vegan::vegdist(t(abund_arc_bac_species), method = "bray")

bray_curtis_pcoa_species <- ecodist::pco(bray_curtis_dist_species)
bray_curtis_pcoa_df_species <- data.frame(pcoa1 = bray_curtis_pcoa_species$vectors[,1], 
                                          pcoa2 = bray_curtis_pcoa_species$vectors[,2])
dim(bray_curtis_pcoa_df_species)
dim(metadata)

bray_curtis_pcoa_df_species <- cbind(bray_curtis_pcoa_df_species, metadata)

# Extraer varianza explicada por cada eje
eig_vals <- bray_curtis_pcoa_species$values / sum(bray_curtis_pcoa_species$values) * 100
pc1_var <- round(eig_vals[1], 2)  # % explicado por PC1
pc2_var <- round(eig_vals[2], 2)  # % explicado por PC2

# Definir colores
season_colors <- c("El_Nino" = "#c83737", "La_Nina" = "#2c5aa0")
sector_colors <- c("A_Open" = "#007ea7", "B_Intermedia" = "#6a994e", "C_Coast" = "#c8ab37ff")

#"A_Open_ocean","B_Ocean_slope","C_Coastal_slope", "D_Coast"

# Plot con ajustes
p1 <- ggplot(bray_curtis_pcoa_df_species, aes(x = pcoa1, y = pcoa2, color = Season, shape = Sector)) +
  theme_bw() +
  geom_point(size = 4) +  # Ajustar tamaño de los puntos
  stat_ellipse(aes(x = pcoa1, y = pcoa2, group = Season, color = Season), 
               geom = "path", linewidth = 0.5) +  # Agrupar solo por Season
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"), 
    y = paste0("PC2 (", pc2_var, "%)"),
    title = "Bray-Curtis PCoA species - Season"
  ) + 
  theme(
    title = element_text(size = 14, face = "bold", hjust = 0.5), 
    axis.title.x = element_text(size = 12, face = "bold"), 
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12), 
    axis.text.y = element_text(size = 12), 
    legend.text = element_text(size = 12)
  )

print(p1)


p2 <- ggplot(bray_curtis_pcoa_df_species, aes(x = pcoa1, y = pcoa2, color = Season)) +
  theme_bw() +
  geom_point(size = 4) +  # Ajustar tamaño de los puntos
  stat_ellipse(aes(x = pcoa1, y = pcoa2, group = Season, color = Season), 
               geom = "path", linewidth = 0.5
               #, alpha = 0.2
  ) +
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"), 
    y = paste0("PC2 (", pc2_var, "%)"),
    title = "Bray-Curtis PCoA species - Season"
  ) + 
  theme(
    title = element_text(size = 14, face = "bold", hjust = 0.5), 
    axis.title.x = element_text(size = 12, face = "bold"), 
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12), 
    axis.text.y = element_text(size = 12), 
    legend.text = element_text(size = 12)
  )

print(p2)


## Statistics

# Asegurar que los metadatos tienen las variables de interés
metadata <- data.frame(metadata) # Convertir a data.frame si es necesario

# PERMANOVA para Season
set.seed(123)
permanova_season <- adonis2(bray_curtis_dist_species ~ Season, data = metadata, permutations = 999)
print(permanova_season)
# global p-value 0.079


# PERMANOVA para Sector
set.seed(123)
permanova_sector <- adonis2(bray_curtis_dist_species ~ Sector, data = metadata, permutations = 999)
print(permanova_sector)
# global p-value 0.562


svglite("Beta_Diversity_species_ZMO_opc1.svg", width=6, height=4)
p1
dev.off()

svglite("Beta_Diversity_species_ZMO_opc2.svg", width=6, height=4)
p2
dev.off()



## ------------------------------------------
##
## USING AITCHISON DISTANCE
##
## ------------------------------------------


# Para comparar las proporciones relativas de especies en lugar de las cantidades absolutas, es decir,
# la composición relativa de las comunidades (cómo cambian las proporciones de taxones en lugar de sus abundancias)

## ------------------------------------------
## Build the phyloseq object
## ------------------------------------------

#Set a seed for the purposes of reproducibility in this document.
set.seed(1)

#Perform a CLR transformation - #We are imputing zeroes using the 'const' method 
#Essentially, we replace zeroes with 65% of the next lowest value - see Lubbe et al 2021.
count_species.clr <- clr_c(count_species)


#ps object of count table
#Used for Aitchison distance
taxa_table <- as.matrix(row.names(count_species.clr))
colnames(taxa_table) <- "species"
rownames(taxa_table) <- taxa_table
OTU_taxa <- otu_table(count_species.clr, taxa_are_rows = TRUE)
TAX_taxa <- tax_table(taxa_table)
sampledata <- sample_data(metadata)
rownames(sampledata) <- metadata$Sample

ps_count_species.clr = phyloseq(OTU_taxa, TAX_taxa, sampledata)
save(ps_count_species.clr, file = "ps_count_species.clr.RData")

## ------------------------------------------
## Calculate Aitchison distance
## ------------------------------------------

# Load phyloseq file
load("ps_count_species.clr.RData")

### We choose the Aitchison distance using CLR transformed data at species level
#Generate distance matrix using Aitchison distances
dist_matrix <- phyloseq::distance(ps_count_species.clr, method = "euclidean")

## ------------------------------------------
## Dispersion test with Season
## ------------------------------------------

#Calculate the beta dispersion
dispr <- vegan::betadisper(dist_matrix, phyloseq::sample_data(ps_count_species.clr)$Season)

#Check the results of the dispersion test
dispr

#PERMANOVA test to dispersion
#Set a seed for the purposes of reproducibility
set.seed(123)
vegan::adonis2(dist(dispr$distances) ~dispr$group) # Global p-value 0.179

#Pairwise PERMANOVA test with fdr correction of dispersion
set.seed(123)
pairwise_test <- vegan::permutest(dispr, pairwise = TRUE)
pairwise_pvalues <- pairwise_test$pairwise$permuted 
pairwise_pvalues
pairwise_pvalues_fdr <- p.adjust(pairwise_pvalues, method = "fdr")
pairwise_pvalues_fdr 
# A_Open-C_Coast q = 0.923

# PERMANOVA test to microbiota beta diversity (communities structure)
set.seed(123)
vegan::adonis2(dist_matrix ~ phyloseq::sample_data(ps_count_species.clr)$Season, data = metadata, permutations = 9999)
# global p-value 0.079

# Extraer varianza explicada por cada eje
eig_vals <- dispr$eig / sum(dispr$eig) * 100
pc1_var <- round(eig_vals[1], 2)  # % explicado por PC1
pc2_var <- round(eig_vals[2], 2)  # % explicado por PC2

#PCA and dispersion plot
svglite("disp_pca_Season_species_ZMO.svg", width=4, height=4)
plot(dispr, main = "", sub = "",
     xlab = paste0("PC1 (", pc1_var, "%)"), 
     ylab = paste0("PC2 (", pc2_var, "%)"),
     col = c("El_Nino" = "#c83737", "La_Nina" = "#2c5aa0"), 
     cex = 1.5)
dev.off()

## ------------------------------------------
## Dispersion test with Sector
## ------------------------------------------

#Calculate the beta dispersion
dispr <- vegan::betadisper(dist_matrix, phyloseq::sample_data(ps_count_species.clr)$Sector)

#Check the results of the dispersion test
dispr

#PERMANOVA test to dispersion
#Set a seed for the purposes of reproducibility
set.seed(123)
vegan::adonis2(dist(dispr$distances) ~dispr$group) # Global p-value 0.191

#Pairwise PERMANOVA test with fdr correction of dispersion
set.seed(123)
pairwise_test <- vegan::permutest(dispr, pairwise = TRUE)
pairwise_pvalues <- pairwise_test$pairwise$permuted 
pairwise_pvalues
pairwise_pvalues_fdr <- p.adjust(pairwise_pvalues, method = "fdr")
pairwise_pvalues_fdr 
# A_Open-C_Coast Q = 0.012

# PERMANOVA test to microbiota beta diversity (communities structure)
set.seed(123)
vegan::adonis2(dist_matrix ~ phyloseq::sample_data(ps_count_species.clr)$Sector, data = metadata, permutations = 999)
# global p-value 0.375

# Extraer varianza explicada por cada eje
eig_vals <- dispr$eig / sum(dispr$eig) * 100
pc1_var <- round(eig_vals[1], 2)  # % explicado por PC1
pc2_var <- round(eig_vals[2], 2)  # % explicado por PC2

#PCA and dispersion plot
#svglite("disp_pca_Sector_species_ZMO.svg", width=4, height=4)
plot(dispr, main = "", sub = "",
     xlab = paste0("PC1 (", pc1_var, "%)"), 
     ylab = paste0("PC2 (", pc2_var, "%)"),
     col = c("A_Open" = "#007ea7", "B_Intermedia" = "#6a994e", "C_Coast" = "#c8ab37ff"), 
     cex = 1.5)
#dev.off()

# "A_Open_ocean","B_Ocean_slope","C_Coastal_slope", "D_Coast"







################################
#### FUNCTIONAL COMPARISONS ####
################################



# Definir el directorio
base_path <- "/Users/jorgerv/Documents/lagunas_costeras/Results_databases/Second_version/1-8/"

# Leer los archivos usando row.names para que la primera columna sean los nombres de las filas
total_cds_quintin <- read.csv(paste0(base_path, "1-8_count_cds.tsv"), sep = "\t", row.names = 1)
valid_samples <- rownames(total_cds_quintin)

## -------------------------
##
## Análisis con genes
##
## -------------------------

metal_genes <- read.csv(paste0(base_path, "count_BACMET_Metal_genes.csv"), sep = ",", row.names = 1)
drug_genes <- read.csv(paste0(base_path, "count_CARD_Drug_Class_genes.csv"), sep = ",", row.names = 1)
drug_class <- read.csv(paste0(base_path, "count_CARD_Drug_Class.csv"), sep = ",", row.names = 1)
nitrogen_genes <- read.csv(paste0(base_path, "count_NITROGEN_pathway_genes.csv"), sep = ",", row.names = 1)
sulfur_genes <- read.csv(paste0(base_path, "count_SULFUR_pathway_genes.csv"), sep = ",", row.names = 1)
virulence_genes <- read.csv(paste0(base_path, "count_VFDB_type_genes.csv"), sep = ",", row.names = 1)

# Función para filtrar las tablas de genes por las muestras en valid_samples
# y eliminar las columnas numéricas que sumen cero, preservando las filas por nombres
filtrar_y_eliminar_columnas_cero <- function(genes, valid_samples) {
  # Filtrar por muestras en valid_samples
  genes_filtered <- genes[rownames(genes) %in% valid_samples, ]
  # Eliminar las columnas cuyo conteo suma cero
  genes_filtered <- genes_filtered[, colSums(genes_filtered) != 0]
  
  return(genes_filtered)
}

# Aplicar la función a cada tabla de genes
metal_genes_filtered <- filtrar_y_eliminar_columnas_cero(metal_genes, valid_samples)
drug_genes_filtered <- filtrar_y_eliminar_columnas_cero(drug_genes, valid_samples)
drug_class_filtered <- filtrar_y_eliminar_columnas_cero(drug_class, valid_samples)
nitrogen_genes_filtered <- filtrar_y_eliminar_columnas_cero(nitrogen_genes, valid_samples)
sulfur_genes_filtered <- filtrar_y_eliminar_columnas_cero(sulfur_genes, valid_samples)
virulence_genes_filtered <- filtrar_y_eliminar_columnas_cero(virulence_genes, valid_samples)

# Renombrar la columna de 'total_cds_quintin' si es necesario
# colnames(total_cds_quintin)[1] <- "sample"  # No es necesario porque ya se usaron rownames

# Función para calcular abundancias con transformación logarítmica
calcular_abundancias_log <- function(genes_filtered, total_cds) {
  # Unir las tablas por las rownames
  genes_abund <- cbind(genes_filtered, Conteo_CDS = total_cds[rownames(genes_filtered), ])
  
  # Dividir los conteos por el total de CDS y aplicar la transformación logarítmica
  genes_abund <- genes_abund %>%
    dplyr::mutate(across(-Conteo_CDS, ~ ((.x / Conteo_CDS)*100))) %>%
    dplyr::select(-Conteo_CDS)
  
  # Eliminar las columnas cuya suma es cero
  genes_abund <- genes_abund %>%
    dplyr::select(where(~ sum(.) != 0))
  
  return(genes_abund)
}

# Aplicar la función a cada tabla de genes
metal_genes_abund <- calcular_abundancias_log(metal_genes_filtered, total_cds_quintin)
drug_genes_abund <- calcular_abundancias_log(drug_genes_filtered, total_cds_quintin)
nitrogen_genes_abund <- calcular_abundancias_log(nitrogen_genes_filtered, total_cds_quintin)
sulfur_genes_abund <- calcular_abundancias_log(sulfur_genes_filtered, total_cds_quintin)
virulence_genes_abund <- calcular_abundancias_log(virulence_genes_filtered, total_cds_quintin)

# Transpose of abundance tables
metal_genes_abund <- as.data.frame(t(metal_genes_abund))
drug_genes_abund <- as.data.frame(t(drug_genes_abund))
nitrogen_genes_abund <- as.data.frame(nitrogen_genes_abund)
sulfur_genes_abund <- as.data.frame(sulfur_genes_abund)
virulence_genes_abund <- as.data.frame(t(virulence_genes_abund))

save(metal_genes_abund, file = "genes_abund_metal.RData")
save(drug_genes_abund, file = "genes_abund_drug.RData")
save(nitrogen_genes_abund, file = "genes_abund_nitrogen.RData")
save(sulfur_genes_abund, file = "genes_abund_sulfur.RData")
save(virulence_genes_abund, file = "genes_abund_virulence.RData")

# Most abundant genes
# Calcular el porcentaje promedio de cada columna
metal_genes_abund$mean <- rowMeans(metal_genes_abund)
drug_genes_abund$mean <- rowMeans(drug_genes_abund)
nitrogen_genes_abund$mean <- rowMeans(nitrogen_genes_abund)
sulfur_genes_abund$mean <- rowMeans(sulfur_genes_abund)
virulence_genes_abund$mean <- rowMeans(virulence_genes_abund)


## --------------------------
##
## Análisis con pathways
##
## --------------------------

metal_pathways <- read.csv(paste0(base_path, "count_BACMET_Metal.csv"), sep = ",", row.names = 1)
drug_class <- read.csv(paste0(base_path, "count_CARD_Drug_Class.csv"), sep = ",", row.names = 1)
nitrogen_pathways <- read.csv(paste0(base_path, "count_NITROGEN_pathway.csv"), sep = ",", row.names = 1)
sulfur_pathways <- read.csv(paste0(base_path, "count_SULFUR_pathway.csv"), sep = ",", row.names = 1)
virulence_pathways <- read.csv(paste0(base_path, "count_VFDB_type.csv"), sep = ",", row.names = 1)


# Función para filtrar las tablas de pathways por las muestras en valid_samples
# y eliminar las columnas numéricas que sumen cero, preservando las filas por nombres
filtrar_y_eliminar_columnas_cero <- function(pathways, valid_samples) {
  # Filtrar por muestras en valid_samples
  pathways_filtered <- pathways[rownames(pathways) %in% valid_samples, ]
  # Eliminar las columnas cuyo conteo suma cero
  pathways_filtered <- pathways_filtered[, colSums(pathways_filtered) != 0]
  
  return(pathways_filtered)
}

# Aplicar la función a cada tabla de pathways
metal_pathways_filtered <- filtrar_y_eliminar_columnas_cero(metal_pathways, valid_samples)
drug_pathways_filtered <- filtrar_y_eliminar_columnas_cero(drug_class, valid_samples)
nitrogen_pathways_filtered <- filtrar_y_eliminar_columnas_cero(nitrogen_pathways, valid_samples)
sulfur_pathways_filtered <- filtrar_y_eliminar_columnas_cero(sulfur_pathways, valid_samples)
virulence_pathways_filtered <- filtrar_y_eliminar_columnas_cero(virulence_pathways, valid_samples)

# Renombrar la columna de 'total_cds_quintin' si es necesario
# colnames(total_cds_quintin)[1] <- "sample"  # No es necesario porque ya se usaron rownames

# Función para calcular abundancias con transformación logarítmica FUNCIONO
calcular_abundancias_log <- function(pathways_filtered, total_cds) {
  total_cds <- as.data.frame(total_cds)  # Asegurar que es un dataframe
  pathways_abund <- cbind(pathways_filtered, Conteo_CDS = total_cds[rownames(pathways_filtered), , drop = FALSE])
  
  # Comprobar si 'Conteo_CDS' está en el dataframe
  if (!"Conteo_CDS" %in% colnames(pathways_abund)) {
    stop("Error: 'Conteo_CDS' no está en pathways_abund. Revisa los rownames.")
  }
  
  print("Columnas antes de select:")  
  print(colnames(pathways_abund))  # Debugging
  
  pathways_abund <- pathways_abund %>%
    mutate(across(-Conteo_CDS, ~ ((.x / Conteo_CDS) * 100))) %>%
    dplyr::select(-Conteo_CDS)  # Asegurar que se usa el select de dplyr
  
  pathways_abund <- pathways_abund %>%
    dplyr::select(where(~ sum(.) != 0))
  
  return(pathways_abund)
}


# Aplicar la función a cada tabla de pathways
metal_pathways_abund <- calcular_abundancias_log(metal_pathways_filtered, total_cds_quintin)
drug_pathways_abund <- calcular_abundancias_log(drug_pathways_filtered, total_cds_quintin)
nitrogen_pathways_abund <- calcular_abundancias_log(nitrogen_pathways_filtered, total_cds_quintin)
sulfur_pathways_abund <- calcular_abundancias_log(sulfur_pathways_filtered, total_cds_quintin)
virulence_pathways_abund <- calcular_abundancias_log(virulence_pathways_filtered, total_cds_quintin)

# Transform to data frames the abundance tables
metal_pathways_abund <- as.data.frame(metal_pathways_abund)
drug_pathways_abund <- as.data.frame(drug_pathways_abund)
nitrogen_pathways_abund <- as.data.frame(nitrogen_pathways_abund)
sulfur_pathways_abund <- as.data.frame(sulfur_pathways_abund)
virulence_pathways_abund <- as.data.frame(virulence_pathways_abund)

summary(metal_pathways_abund)
summary(drug_pathways_abund)
summary(nitrogen_pathways_abund)
summary(sulfur_pathways_abund)
summary(virulence_pathways_abund)

# Save data
save(drug_pathways_abund, file = "pathways_abund_drug.RData")
save(metal_pathways_abund, file = "pathways_abund_metal.RData")
save(nitrogen_pathways_abund, file = "pathways_abund_nitrogen.RData")
save(sulfur_pathways_abund, file = "pathways_abund_sulfur.RData")
save(virulence_pathways_abund, file = "pathways_abund_virulence.RData")


load("pathways_abund_metal.RData")
load("pathways_abund_virulence.RData")

summary(metal_pathways_abund)
mean(as.numeric(as.matrix(metal_pathways_abund)))
max(as.numeric(as.matrix(metal_pathways_abund))) - mean(as.numeric(as.matrix(metal_pathways_abund)))
summary(virulence_pathways_abund)
mean(as.numeric(as.matrix(virulence_pathways_abund)))
max(as.numeric(as.matrix(virulence_pathways_abund))) - mean(as.numeric(as.matrix(virulence_pathways_abund)))

## -----------------------------
## 
## BOX PLOTS ABUNDANCE PATHWAYS
##
## -----------------------------


library(ggplot2)
library(dplyr)
library(tidyr)

# Load data
load("pathways_abund_drug.RData")
load("pathways_abund_metal.RData")
load("pathways_abund_nitrogen.RData")
load("pathways_abund_sulfur.RData")
load("pathways_abund_virulence.RData")

load("metadata_ZMO.RData")

# Transponer la tabla de abundancia
drug_pathways_abund$Sample <- rownames(drug_pathways_abund)
metal_pathways_abund$Sample <- rownames(metal_pathways_abund)
nitrogen_pathways_abund$Sample <- rownames(nitrogen_pathways_abund)
sulfur_pathways_abund$Sample <- rownames(sulfur_pathways_abund)
virulence_pathways_abund$Sample <- rownames(virulence_pathways_abund)

# Unir las tablas de abundancia y metadatos
combined_data <- merge(drug_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(metal_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(nitrogen_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(sulfur_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(virulence_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")

# Definir colores
my_colors <- c("A_La_Nina" = "#2c5aa0",
               "B_El_Nino" = "#c83737",
               "A_Open_ocean" = "#007ea7", 
               "B_Ocean_slope" = "#6a994e", 
               "C_Coastal_slope" = "#a05a2c",
               "D_Coast" = "#c8ab37ff")


# Crear una lista de columnas de interés
#Nitrogen
columnas_de_interes <- c(2, 3, 4, 5, 6, 7, 8)
# Metals
columnas_de_interes <- c(2,3,4,8,9,10,11,13,14,15,17,18,19,20,24)
# Sulfur
columnas_de_interes <- c(2, 3, 4, 5, 6, 7, 8)
# Drug
columnas_de_interes <- c(2:28)
# Virulence
columnas_de_interes <- c(2:17)

# Crear las comparaciones adecuadas
# Para ver resultados
a_my_comparisons <- list(
  c("A_La_Nina", "B_El_Nino"),
  c("A_Open_ocean", "B_Ocean_slope"),
  c("A_Open", "C_Coastal_slope"),
  c("A_Open", "D_Coast"),
  c("B_Ocean_slope", "C_Coastal_slope"),
  c("B_Ocean_slope", "D_Coast"),
  c("C_Coastal_slope", "D_Coast")
)

# Para graficar
a_my_comparisons <- list(
  c("B_Ocean_slope", "C_Coastal_slope"),
  c("B_Ocean_slope", "D_Coast")
)

# Bucle para generar y almacenar las figuras en variables a1, a2, a3, ...
for (i in 1:length(columnas_de_interes)) {
  columna_index <- columnas_de_interes[i]
  columna_de_interes <- colnames(combined_data)[columna_index]
  # Transformación a formato largo
  df_long <- combined_data %>% 
    dplyr::select(.data[[columna_de_interes]], Sector, Season) %>%           
    pivot_longer(cols = c(Sector, Season),
                 names_to  = "GroupType",
                 values_to = "Group")
  
  # Fijamos el orden de las siete categorías en el eje x
  level_order <- c("A_La_Nina","B_El_Nino", 
                   "A_Open_ocean","B_Ocean_slope","C_Coastal_slope","D_Coast")
  df_long$Group <- factor(df_long$Group, levels = level_order)

  # Crear el boxplot con el t test usando la nueva columna combinada
  figura <- ggplot(df_long, aes(x = Group, y = .data[[columna_de_interes]], fill = Group)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.75)) +
    geom_jitter(width = 0.0, size = 0.8, alpha = 0.5, colour = "black") +
    scale_fill_manual(values = my_colors) +
    stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.format", 
                       position = position_dodge(0.75)) +
    labs(x = "", y = "", fill = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 14),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste(columna_de_interes)) #+
  #ylim(0, 7.5)  # Fijar los límites del eje y
  
  # Almacenar la figura en una variable dinámica (a1, a2, a3, ...)
  assign(paste0("a", i), figura)
}


svglite("annot_nitrogen_metabolism.svg", width=29, height=4)
((a1 | a2 | a3 | a4 | a5 | a6 | a7))
dev.off()

svglite("annot_metal_metabolism.svg", width=29, height=8)
((a1 | a2 | a3 | a4 | a5 | a6 | a7) / (a8 | a9 | a10| a11 | a12 | a13 | a14))
dev.off()

svglite("annot_sulfur_metabolism.svg", width=29, height=4)
((a1 | a2 | a3 | a4 | a5 | a6 | a7))
dev.off()

svglite("annot_drug_metabolism.svg", width=29, height=16)
((a1 | a2 | a3 | a4 | a5 | a6 | a7) / (a8 | a9 | a10| a11 | a12 | a13 | a14) / (a15 | a16 | a17| a18 | a19 | a20 | a21) / (a22 | a23 | a24 | a25 | a26 | a1 | a2))
dev.off()

svglite("annot_virulence_metabolism.svg", width=24.7, height=12)
((a1 | a2 | a3 | a4 | a5 | a6) / (a7 | a8 | a9 | a10 | a11 | a12) / (a13 | a14 | a15 | a16 | a1 | a12))
dev.off()




### ---------------------------
## Diferencias entre pathways
## ----------------------------


# Cargar paquetes necesarios
library(dplyr)
library(tidyr)
library(FSA)        # para dunnTest()
library(rstatix)    # para ajustar p-values si se prefiere

# Unir las tablas de abundancia y metadatos
combined_data <- merge(drug_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(metal_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(nitrogen_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(sulfur_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(virulence_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")

# Ajusta el modelo con los tres factores y todas sus interacciones
combined_data$Sector <- as.factor(combined_data$Sector)
combined_data$Season <- as.factor(combined_data$Season)

# 1) Columnas de abundancia
pathways_cols <- combined_data %>%
  dplyr::select(-Sample, -Sector, -Season, -ID_1, -ID_2) %>%
  names()

# 2) Calcular KW p-values para cada pathway
kw_results <- lapply(pathways_cols, function(var){
  kw <- kruskal.test(combined_data[[var]] ~ combined_data$Sector)
  data.frame(
    Pathway = var,
    P.value = kw$p.value
  )
}) %>%
  bind_rows() %>%
  arrange(P.value)

# 3) Kruskal–Wallis + Dunn para Sector
kw_dunn_results <- lapply(pathways_cols, function(var){
  # 1) Kruskal–Wallis
  kw <- kruskal.test(combined_data[[var]] ~ combined_data$Sector)
  
  # 2) Si KW sale p ≤ 0.05, lanza Dunn
  if(kw$p.value <= 0.1){
    d <- dunnTest(
      x      = combined_data[[var]],
      g      = combined_data$Sector,
      method = "bh"
    )$res
    
    d %>%
      dplyr::transmute(
        Pathway    = var,
        Comparison = Comparison,
        P.unadj    = P.unadj,
        P.adj      = P.adj
      )
  } else {
    NULL
  }
}) %>%
  bind_rows()

# 4) Wilcoxon para Season
wilcox_results <- lapply(c("Season"), function(var){
  lapply(pathways_cols, function(col){
    wt <- wilcox.test(
      combined_data[[col]] ~ combined_data[[var]],
      data = combined_data,
      exact = FALSE
    )
    data.frame(
      Pathway  = col,
      Variable = var,
      P.value  = wt$p.value
    )
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  # corrección BH por grupo de Variable
  group_by(Variable) %>%
  mutate(P.adj = p.adjust(P.value, method = "BH")) %>%
  ungroup()

# 5) Mostrar todos los resultados
kw_results
print(wilcox_results, n = 100)

# 6) Filtrar resultados “tendencia” 0.05 < p ≤ 0.1
kw_dunn_sig   <- kw_dunn_results   %>% dplyr::filter(P.adj <= 0.1)
wilcox_sig    <- wilcox_results    %>% dplyr::filter(P.adj <= 0.1)

# 7) Listado final
list(
  Dunn_Sector     = kw_dunn_sig,
  Wilcox_SeasonHabitat = wilcox_sig
)

# Nitrogen: no significant results
# Sulfur: no significant results
# Drugs: no significant results
# Metals: no significant results
# Virulence: no significant results






## ----------------------------
## 
## HEATMAPS GENES
##
## ----------------------------


library(ggplot2)
library(dplyr)
library(tidyr)

# Load data
load("genes_abund_metal.RData")
load("genes_abund_drug.RData")
load("genes_abund_nitrogen.RData")
load("genes_abund_sulfur.RData")
load("genes_abund_virulence.RData")

load("metadata_ZMO.RData")

# Función para convertir matriz en formato largo, unir con metadata, y calcular medias
procesar_genes <- function(genes_abund, metadata) {
  # Convertir la matriz en data.frame
  genes_abund_df <- as.data.frame(t(genes_abund))
  
  # Convertir al formato largo
  genes_long <- genes_abund_df %>%
    rownames_to_column(var = "Gene") %>%
    pivot_longer(-Gene, names_to = "Sample", values_to = "Abundance")
  
  # Unir con metadata
  genes_metadata <- genes_long %>%
    left_join(metadata, by = c("Sample" = "Sample"))
  
  # Calcular medias por Sector
  sector_means <- genes_metadata %>%
    group_by(Gene, Sector) %>%
    dplyr::summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    pivot_wider(names_from = Sector, values_from = mean_abundance)
  
  # Calcular medias por Season
  season_means <- genes_metadata %>%
    group_by(Gene, Season) %>%
    dplyr::summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    pivot_wider(names_from = Season, values_from = mean_abundance)
  
  # Unir todas las medias en un único dataframe
  genes_means <- genes_abund_df %>%
    rownames_to_column(var = "Gene") %>%
    left_join(sector_means, by = "Gene") %>%
    left_join(season_means, by = "Gene")
  
  return(genes_means)
}

# Aplicar la función a cada una de las matrices
metal_means <- procesar_genes(metal_genes_abund, metadata)
drug_means <- procesar_genes(drug_genes_abund, metadata)
nitrogen_means <- procesar_genes(nitrogen_genes_abund, metadata)
sulfur_means <- procesar_genes(sulfur_genes_abund, metadata)
virulence_means <- procesar_genes(virulence_genes_abund, metadata)

# Calculo los valores máximos
#max_value_drug <- max(drug_means[ , (ncol(drug_means) - 6):ncol(drug_means)], na.rm = TRUE)
max_value_nitrogen <- max(nitrogen_means[ , (ncol(nitrogen_means) - 6):ncol(nitrogen_means)], na.rm = TRUE)
max_value_sulfur <- max(sulfur_means[ , (ncol(sulfur_means) - 6):ncol(sulfur_means)], na.rm = TRUE)

# Función para generar los heatmaps de un archivo *_means
crear_heatmaps <- function(means_data, max_val) {
  
  # Convertir a formato largo para Sector and Season
  sector_data <- means_data %>%
    dplyr::select(Gene, starts_with("A_Open_ocean"), starts_with("B_"), 
                  starts_with("C_"), starts_with("D_")) %>%
    pivot_longer(-Gene, names_to = "Sector", values_to = "Abundance")
  
  season_data <- means_data %>%
    dplyr::select(Gene, starts_with("A_La_Nina"), starts_with("B_El_Nino")) %>%
    pivot_longer(-Gene, names_to = "Season", values_to = "Abundance")
  
  # Crear heatmap para Sector
  heatmap_sector <- ggplot(sector_data, aes(x = Sector, y = Gene, fill = Abundance)) +
    geom_tile(color = "white", linewidth = 0.5) +
    #scale_fill_gradient2(low = "white", high = "darkred", 
    scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),
                         na.value = "grey50", limit = c(0, max_val),
                         name = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "", fill = "Abundance") +
    coord_fixed(ratio = 0.2) 
  
  # Crear heatmap para Season
  heatmap_season <- ggplot(season_data, aes(x = Season, y = Gene, fill = Abundance)) +
    geom_tile(color = "white", linewidth = 0.5) +
    #scale_fill_gradient2(low = "white", high = "darkred", 
    scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),
                         na.value = "grey50", limit = c(0, max_val),
                         name = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),  # Ocultar nombres de los genes en el eje y
          axis.ticks.y = element_blank()) +
    labs(title = "", fill = "Abundance") +
    coord_fixed(ratio = 0.2) 
  
  # Unir los dos heatmaps en un solo gráfico utilizando patchwork
  library(patchwork)
  
  # Combinar los tres heatmaps en una figura
  combined_plot <- (heatmap_sector | heatmap_season) + 
    plot_layout(guides = "collect") & theme(legend.position = "right")
  
  return(combined_plot)
}

# Generar heatmaps para cada archivo *_means
#heatmap_drug <- crear_heatmaps(drug_means, max_value_drug)
heatmap_nitrogen <- crear_heatmaps(nitrogen_means, max_value_nitrogen)
heatmap_sulfur <- crear_heatmaps(sulfur_means, max_value_sulfur)

# Mostrar el plot de cada uno
#print(heatmap_drug)
print(heatmap_nitrogen)
print(heatmap_sulfur)

# Definir el número de genes en cada dataset
num_genes_drug <- nrow(drug_means)
num_genes_nitrogen <- nrow(nitrogen_means)
num_genes_sulfur <- nrow(sulfur_means)

# Ajustar la altura dinámicamente
height_factor <- 0.3  # Ajusta este valor según el espaciado que necesites

# Guardar las figuras con dimensiones ajustadas
svglite("heatmap_CARD_genes.svg", width = 8, height = max(10, num_genes_drug * height_factor))
heatmap_drug
dev.off()

svglite("heatmap_NITROGEN_genes.svg", width = 8, height = max(8, num_genes_nitrogen * height_factor))
heatmap_nitrogen
dev.off()

svglite("heatmap_SULFUR_genes.svg", width = 8, height = max(8, num_genes_sulfur * height_factor))
heatmap_sulfur
dev.off()


## ----------------------------
## Diferencias entre genes
## ----------------------------


# Cargar paquetes necesarios
library(FSA)      # Para Dunn test
library(dplyr)    # Para manipulación de datos
library(tidyr)    # Para transformar datos
library(rstatix)  # Para múltiples comparaciones ajustadas

# Transponer la tabla de abundancia
drug_genes_abund$Sample <- rownames(drug_genes_abund)
metal_genes_abund$Sample <- rownames(metal_genes_abund)
nitrogen_genes_abund$Sample <- rownames(nitrogen_genes_abund)
sulfur_genes_abund$Sample <- rownames(sulfur_genes_abund)
virulence_genes_abund$Sample <- rownames(virulence_genes_abund)

# Unir las tablas de abundancia y metadatos
combined_data <- merge(drug_genes_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(metal_genes_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(nitrogen_genes_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(sulfur_genes_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(virulence_genes_abund, metadata, by.x = "Sample", by.y = "Sample")

# Ajusta el modelo con los tres factores y todas sus interacciones
combined_data$Sector <- as.factor(combined_data$Sector)
combined_data$Season <- as.factor(combined_data$Season)

# Seleccionar solo las columnas numéricas (excluyendo "Sample", "Sector", "Season", "Habitat")
genes_cols <- combined_data %>%
  dplyr::select(-Sample, -Sector, -Season, -ID_1, -ID_2)

# Función para aplicar Dunn test (Sector: 3 niveles)
dunn_results <- lapply(names(genes_cols), function(gene) {
  test <- dunnTest(combined_data[[gene]] ~ combined_data$Sector, method = "bh")
  result <- data.frame(Gene = gene, test$res)
  return(result)
}) %>%
  bind_rows() %>%
  filter(P.adj <= 0.1)  # Filtrar p-values ajustados < 0.05

# Función para aplicar Wilcoxon test (Season: 2 niveles)
wilcox_results <- lapply(c("Season"), function(variable) {
  lapply(names(genes_cols), function(gene) {
    test <- wilcox.test(combined_data[[gene]] ~ combined_data[[variable]])
    data.frame(Gene = gene, Variable = variable, p.value = test$p.value)
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  filter(p.value <= 0.1)  # Filtrar p-values <= 0.05


# Ajuste de p-values global (sobre todos los tests)
wilcox_results <- wilcox_results %>%
  dplyr::mutate(
    p.adj = p.adjust(p.value, method = "BH")
  ) %>%
  arrange(p.adj)

# Opcional: quedarte solo con los tests significativos tras ajuste
sig_wilcox <- wilcox_results %>%
  filter(p.adj <= 0.1)

# Mostrar los resultados significativos
list(Dunn_Test_Significativo = dunn_results, Wilcoxon_Significativo = sig_wilcox )

## NITROGEN

# $Dunn_Test_Significativo
# Gene            Comparison       Z    P.unadj      P.adj
# <0 rows> (or 0-length row.names)

# $Wilcoxon_Significativo
# Gene Variable    p.value      p.adj
# 1 amoA_A   Season 0.05714286 0.05714286
# 2   nirS   Season 0.05714286 0.05714286
# 3   nosZ   Season 0.05714286 0.05714286
# 4   narZ   Season 0.05714286 0.05714286
# 5   narB   Season 0.05714286 0.05714286
# 6   narC   Season 0.05714286 0.05714286
# 7   nifH   Season 0.05714286 0.05714286
# 8   ureB   Season 0.05714286 0.05714286
# 9   glsA   Season 0.05714286 0.05714286


## SULFUR

# $Dunn_Test_Significativo
# Gene             Comparison         Z    P.unadj      P.adj
# 1 ddhA  A_Open - B_Intermedia  2.333333 0.01963066 0.05889197
# 2 metY  A_Open - B_Intermedia  2.138090 0.03250944 0.09752833
# 3 rdlA       A_Open - C_Coast -2.160247 0.03075356 0.04613034
# 4 rdlA B_Intermedia - C_Coast -2.366432 0.01796048 0.05388143

# $Wilcoxon_Significativo
# Gene Variable    p.value      p.adj
# 1    sir   Season 0.05714286 0.05714286
# 2   dsrL   Season 0.05714286 0.05714286
# 3   dsrM   Season 0.05714286 0.05714286
# 4   dsrN   Season 0.05714286 0.05714286
# 5   hydD   Season 0.05714286 0.05714286
# 6   acuI   Season 0.05714286 0.05714286
# 7   dddC   Season 0.05714286 0.05714286
# 8   dmsA   Season 0.05714286 0.05714286
# 9   dmsC   Season 0.05714286 0.05714286
# 10  dsyB   Season 0.05714286 0.05714286
# 11   gdh   Season 0.05714286 0.05714286
# 12   mdh   Season 0.05714286 0.05714286
# 13  slcC   Season 0.05714286 0.05714286
# 14  cysE   Season 0.05714286 0.05714286
# 15 hdrA1   Season 0.05714286 0.05714286
# 16 hdrA2   Season 0.05714286 0.05714286
# 17  iseJ   Season 0.03191120 0.05714286



## ----------------------------
## 
## COMPLETENESS PATHWAYS
##
## ----------------------------

# Load data
load("genes_abund_nitrogen.RData")
load("genes_abund_sulfur.RData")

load("metadata_ZMO.RData")

# Transponer la tabla de abundancia
nitrogen_genes_abund$Sample <- rownames(nitrogen_genes_abund)
sulfur_genes_abund$Sample <- rownames(sulfur_genes_abund)

# Unir las tablas de abundancia y metadatos
combined_data <- merge(nitrogen_genes_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(sulfur_genes_abund, metadata, by.x = "Sample", by.y = "Sample")

# Ajusta el modelo con los tres factores y todas sus interacciones
combined_data$Sector <- as.factor(combined_data$Sector)
combined_data$Season <- as.factor(combined_data$Season)

## Para NITROGEN

# 1) Definir la lista de genes por pathway, incluyendo pasajes alternativos
genes_list <- list(
  Anammox = c("hzo", "hzsA", "hzsB", "hzsC", "hdh"),
  
  Nitrogen_fixation = c("anfG", "nifD", "nifH", "nifK", "nifW"),
  
  Nitrification = list(
    groupA = c("amoA_A","amoB_A","amoC_A"),
    groupB = c("amoA_B","amoB_B","amoC_B"),
    common = c("hao","nxrA","nxrB")
  ),
  
  Denitrification = list(
    step1 = list(
      c("napA","napB","napC"),
      c("narH","narJ","narI","narG"),
      c("narV","narW","narY","narZ")
    ),
    step2 = c("nirK","nirS"),
    step3 = c("norB","norC","norZ"),
    step4 = "nosZ"
  ),
  
  Assimilatory_nitrate_reduction = list(
    step1 = list(
      c("narB","narC"),
      "NR",
      c("nasA","nasB")
    ),
    step2 = "nirA"
  ),
  
  Dissimilatory_nitrate_reduction = list(
    step1 = list(
      c("napA","napB","napC"),
      c("narH","narJ","narI","narG"),
      c("narV","narW","narY","narZ")
    ),
    step2 = list(
      c("nirB","nirD"),
      c("nrfA","nrfB","nrfC","nrfD")
    )
  ),
  
  Organic_degradation_and_synthesis = list(
    route1 = c("glsA","glnA","asnB","ansB",
               "gs_K00264","gs_K00265","gs_K00266","gs_K00284"),
    route2 = c("ureA","ureB","ureC","nao","nmo",
               "gdh_K00260","gdh_K00261","gdh_K00262","gdh_K15371")
  )
)

# Transponer y unir
nitrogen_genes_abund$Sample <- rownames(nitrogen_genes_abund)
combined_data <- merge(nitrogen_genes_abund, metadata, by = "Sample")

# 3) Crear data.frame de presencia/ausencia (0/1) para **todos** los genes 
#    Si un gen no existe en combined_data, se considera 0
desired_genes <- unlist(
  lapply(genes_list, function(x) if(is.list(x)) unlist(x) else x)
) %>% unique()

present_genes <- intersect(desired_genes, names(combined_data))
missing_genes <- setdiff(desired_genes, names(combined_data))

# Construir presence_df con una columna Sample + todos los genes
presence_df <- combined_data %>%
  dplyr::select(Sample) %>%
  bind_cols(
    combined_data %>%
      dplyr::select(all_of(present_genes)) %>%
      mutate(across(everything(), ~ as.integer(. > 0)))
  )

# Añadir los genes faltantes como columnas de ceros
for(g in missing_genes){
  presence_df[[g]] <- 0L
}

# Reordenar columnas: Sample + desired_genes
presence_df <- presence_df %>%
  dplyr::select(Sample, all_of(desired_genes))

# 4) Función auxiliar para pasos OR en pipelines alternativas
step_present <- function(df_row, groups){
  # df_row: vector numérico de 0/1
  # groups: lista de vectores de nombres
  any(sapply(groups, function(g) any(df_row[g] == 1))) * 1
}

# 5) Calcular completitud de cada pathway (%) por muestra
completeness_df <- presence_df %>%
  rowwise() %>%
  mutate(
    Anammox = sum(c_across(all_of(genes_list$Anammox))) / 5 * 100,
    
    Nitrogen_fixation = if_else(
      c_across("anfG") == 1,
      100,
      sum(c_across(setdiff(genes_list$Nitrogen_fixation, "anfG"))) / 4 * 100
    ),
    
    Nitrification = {
      cntA <- sum(c_across(genes_list$Nitrification$groupA))
      cntB <- sum(c_across(genes_list$Nitrification$groupB))
      cntC <- sum(c_across(genes_list$Nitrification$common))
      (max(cntA, cntB) + cntC) / 6 * 100
    },
    
    Denitrification = {
      s1 <- step_present(cur_data(), genes_list$Denitrification$step1)
      s2 <- any(c_across(genes_list$Denitrification$step2)) * 1
      s3 <- any(c_across(genes_list$Denitrification$step3)) * 1
      s4 <- c_across(genes_list$Denitrification$step4)
      (s1 + s2 + s3 + s4) / 4 * 100
    },
    
    Assimilatory_nitrate_reduction = {
      s1 <- step_present(cur_data(), genes_list$Assimilatory_nitrate_reduction$step1)
      s2 <- c_across(genes_list$Assimilatory_nitrate_reduction$step2)
      (s1 + s2) / 2 * 100
    },
    
    Dissimilatory_nitrate_reduction = {
      s1 <- step_present(cur_data(), genes_list$Dissimilatory_nitrate_reduction$step1)
      s2 <- step_present(cur_data(), genes_list$Dissimilatory_nitrate_reduction$step2)
      (s1 + s2) / 2 * 100
    },
    
    Organic_degradation_and_synthesis = {
      c1 <- sum(c_across(genes_list$Organic_degradation_and_synthesis$route1)) /
        length(genes_list$Organic_degradation_and_synthesis$route1) * 100
      c2 <- sum(c_across(genes_list$Organic_degradation_and_synthesis$route2)) /
        length(genes_list$Organic_degradation_and_synthesis$route2) * 100
      max(c1, c2)
    }
  ) %>%
  ungroup() %>%
  dplyr::select(
    Sample,
    Anammox,
    Nitrogen_fixation,
    Nitrification,
    Denitrification,
    Assimilatory_nitrate_reduction,
    Dissimilatory_nitrate_reduction,
    Organic_degradation_and_synthesis
  )

# 6) Ver resultados
head(completeness_df)
# (Opcional) unirlo a combined_data:
combined_data <- combined_data %>% left_join(completeness_df, by = "Sample")

## Calculo las medias para el completeness

# 2) Función para convertir tu data.frame de completeness en medias por factor
procesar_completeness <- function(completeness_df, metadata) {
  # pivotar a “largo”
  long <- completeness_df %>%
    pivot_longer(
      cols      = -Sample,
      names_to  = "Pathway",
      values_to = "Completeness"
    ) %>%
    left_join(metadata, by = "Sample")
  
  # medias por Sector
  sector_means <- long %>%
    group_by(Pathway, Sector) %>%
    summarise(mean_completeness = mean(Completeness, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = Sector, values_from = mean_completeness)
  
  # medias por Season
  season_means <- long %>%
    group_by(Pathway, Season) %>%
    summarise(mean_completeness = mean(Completeness, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = Season, values_from = mean_completeness)
  
  # unir todo
  full_means <- sector_means %>%
    left_join(season_means,  by = "Pathway")
  
  return(full_means)
}

# 3) Calcular medias de nitrogen completeness
nitrogen_means <- procesar_completeness(completeness_df, metadata)

# 4) Función para graficar heatmaps de medias de completeness
crear_heatmaps <- function(means_df) {
  # long formato para cada factor
  sec <- means_df %>%
    pivot_longer(-Pathway, names_to="Sector",    values_to="Completeness") %>%
    filter(Sector %in% c("A_Open_ocean","B_Ocean_slope","C_Coastal_slope", "D_Coast"))
  
  seas <- means_df %>%
    pivot_longer(-Pathway, names_to="Season",    values_to="Completeness") %>%
    filter(Season %in% c("A_La_Nina","B_El_Nino"))
  
  
  # helper
  mk <- function(df, xvar){
    ggplot(df, aes_string(x = xvar, y = "Pathway", fill = "Completeness")) +
      geom_tile(color = "white", linewidth = 0.5) +
      scale_fill_gradientn(
        colors = c("#ffd166","#619b8a","#3e5c76","#6a4c93","#a4161a"),
        limits = c(0,100),
        na.value = "grey80"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = if(xvar=="Sector") element_text() else element_blank(),
        axis.ticks.y = if(xvar=="Sector") element_line() else element_blank()
      ) +
      labs(x = NULL, y = NULL) +
      coord_fixed(ratio = 0.2)
  }
  
  p1 <- mk(sec,  "Sector")
  p2 <- mk(seas,  "Season")
  
  (p1 | p2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
}

# Generar el heatmap de nitrogen completeness
heatmap_nitrogen <- crear_heatmaps(nitrogen_means)
print(heatmap_nitrogen)

# Ajustar la altura dinámicamente
height_factor <- 0.5  # Ajusta este valor según el espaciado que necesites
num_pathways_nitrogen <- 7

# Guardar las figuras con dimensiones ajustadas
svglite("heatmap_NITROGEN_completeness.svg", width = 8, height = max(8, num_pathways_nitrogen * height_factor))
heatmap_nitrogen
dev.off()



## Para SULFUR

# 0) Cargar librerías
library(dplyr)

# 1) Definir la lista de genes por pathway de azufre
genes_list_sulfur <- list(
  Assimilatory_sulfate_reduction = list(
    step1 = list(c("sat"),
                 c("cysN_cysC","cysN","cysD")),
    step2 = c("cysN_cysC","cysC"),
    step3 = list(c("cysH"), c("nrnA")),
    step4 = list(c("cysJ","cysI"), c("sir"))
  ),
  
  Dissimilatory_sulfur_reduction_and_oxidation = list(
    step1 = c("sat"),
    step2 = list(c("aprA","aprB"),
                 c("qmoA","qmoB","qmoC")),
    step3 = c("dsrA","dsrB","dsrL"),
    step4 = c("dsrC"),
    step5 = c("dsrJ","dsrK","dsrM","dsrN","dsrO","dsrP","dsrT"),
    step6 = c("dsrE","dsrF","dsrH"),
    step7 = c("dsrD","rdsr")
  ),
  
  Sulfur_reduction = list(
    step1 = list(c("asrA","asrB","asrC"),
                 c("fsr"),
                 c("mccA")),
    step2 = list(c("hydB","hydD","hydG"),
                 c("shyA","shyB","shyC","shyD")),
    step3 = c("sudA","sudB","rdlA"),
    step4 = list(c("otr"),
                 c("ttrA","ttrB","ttrC")),
    step5 = list(c("psrA","psrB","psrC"),
                 c("sreA","sreB","sreC"))
  ),
  
  SOX_systems = list(
    step1 = list(c("soxA","soxX"),
                 c("soxY","soxZ")),
    step2 = c("soxB"),
    step3 = c("soxC","soxD")
  ),
  
  Sulfur_oxidation = list(
    step1 = list(c("doxA","doxD"),
                 c("tsda","tsdB")),
    step2 = c("fccA","fccB","sqr"),
    step3 = c("glpE","sseA"),
    step4 = list(c("sorA","sorB"),
                 c("soeA","soeB","soeC"))
  ),
  
  Sulfur_disproportionation = list(
    step1 = c("phsA","phsB","phsC"),
    step2 = c("tetH"),
    step3 = c("sor")
  ),
  
  Organic_sulfur_transformation = list(
    route1 = c("dsyB","dddA","dddC","dddD"),
    route2 = c("dsyB","dddK","dddL","dddP","dddQ","dddT","dddW","dddY","prpE"),
    route3 = c("dsyB","dmdA","dmdB","dmdC","dmdD"),
    route4 = list(c("dsyB"),
                  c("dddD"),
                  c("dddK","dddL","dddP","dddQ","dddT","dddW","dddY"),
                  c("dmoA"),
                  c("mddA")),
    route5 = list(c("dsyB"),
                  c("dddD"),
                  c("dddK","dddL","dddP","dddQ","dddT","dddW","dddY"),
                  c("dmsA","dmsB","dmsC"),
                  c("ddhA","ddhB","ddhC")),
    route6 = c("gah","hpsN","hpsO","hpsP","iseJ","isfD","mdh","mtsA","mtsB",
               "pta","sfnG","slcC","slcD","sqdB","sqdD","sqdX","tauX","tauY",
               "tmm","toa","tpa","yihQ")
  ),
  
  Link_between_inorganic_and_organic_sulfur_transformation = c(
    "cuyA","cysE","cysK","cysM","cysO",
    "hdrA1","hdrA2","hdrB1","hdrB2","hdrC1","hdrC2","hdrD","hdrE",
    "mccB","metA","metB","metC","metX","metY","metZ",
    "msmA","msmB","mtoX","ssuD","ssuE","suyA","suyB","tauD",
    "tbuB","tbuC","tmoC","tmoF","touC","touF","xsc"
  )
)

# 2) Preparar presence_df para estos genes
desired_genes <- unique(unlist(genes_list_sulfur))
present_genes <- intersect(desired_genes, names(sulfur_genes_abund))
missing_genes <- setdiff(desired_genes, names(sulfur_genes_abund))

presence_sulfur <- combined_data %>%
  dplyr::select(Sample) %>%
  bind_cols(
    combined_data %>%
      dplyr::select(all_of(present_genes)) %>%
      mutate(across(everything(), ~ as.integer(. > 0)))
  )

for(g in missing_genes) {
  presence_sulfur[[g]] <- 0L
}

presence_sulfur <- presence_sulfur %>%
  dplyr::select(Sample, all_of(desired_genes))

# 3) Función auxiliar para pasos OR
step_present <- function(df_row, groups) {
  any(sapply(groups, function(g) any(df_row[g] == 1))) * 1
}

# 4) Calcular completeness para cada pathway de azufre
completeness_sulfur <- presence_sulfur %>%
  rowwise() %>%
  mutate(
    Assimilatory_sulfate_reduction = {
      s1 <- step_present(cur_data(), genes_list_sulfur$Assimilatory_sulfate_reduction$step1)
      s2 <- all(c_across(genes_list_sulfur$Assimilatory_sulfate_reduction$step2) == 1) * 1
      s3 <- step_present(cur_data(), genes_list_sulfur$Assimilatory_sulfate_reduction$step3)
      s4 <- step_present(cur_data(), genes_list_sulfur$Assimilatory_sulfate_reduction$step4)
      (s1 + s2 + s3 + s4) / 4 * 100
    },
    Dissimilatory_sulfur_reduction_and_oxidation = {
      s1 <- any(c_across(genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step1)) * 1
      s2 <- step_present(cur_data(), genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step2)
      s3 <- all(c_across(genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step3) == 1) * 1
      s4 <- all(c_across(genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step4) == 1) * 1
      s5 <- all(c_across(genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step5) == 1) * 1
      s6 <- all(c_across(genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step6) == 1) * 1
      s7 <- all(c_across(genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step7) == 1) * 1
      (s1 + s2 + s3 + s4 + s5 + s6 + s7) / 7 * 100
    },
    Sulfur_reduction = {
      sum(sapply(genes_list_sulfur$Sulfur_reduction, function(gr) step_present(cur_data(), gr))) / 5 * 100
    },
    SOX_systems = {
      sum(sapply(genes_list_sulfur$SOX_systems, function(gr) step_present(cur_data(), gr))) / 3 * 100
    },
    Sulfur_oxidation = {
      sum(sapply(genes_list_sulfur$Sulfur_oxidation, function(gr) step_present(cur_data(), gr))) / 4 * 100
    },
    Sulfur_disproportionation = {
      # Paso 1: phsA, phsB, phsC todos presentes
      s1 <- all(c_across(genes_list_sulfur$Sulfur_disproportionation$step1) == 1) * 1
      # Paso 2: tetH presente
      s2 <- c_across(genes_list_sulfur$Sulfur_disproportionation$step2)
      # Paso 3: sor presente
      s3 <- c_across(genes_list_sulfur$Sulfur_disproportionation$step3)
      # Completitud
      (s1 + s2 + s3) / 3 * 100
    },
    Organic_sulfur_transformation = {
      # Calcular cada ruta y tomar la máxima
      r1 <- sum(c_across(genes_list_sulfur$Organic_sulfur_transformation$route1)) /
        length(genes_list_sulfur$Organic_sulfur_transformation$route1) * 100
      r2 <- sum(c_across(genes_list_sulfur$Organic_sulfur_transformation$route2)) /
        length(genes_list_sulfur$Organic_sulfur_transformation$route2) * 100
      r3 <- sum(c_across(genes_list_sulfur$Organic_sulfur_transformation$route3)) /
        length(genes_list_sulfur$Organic_sulfur_transformation$route3) * 100
      r4 <- step_present(cur_data(), genes_list_sulfur$Organic_sulfur_transformation$route4) / 5 * 100
      r5 <- step_present(cur_data(), genes_list_sulfur$Organic_sulfur_transformation$route5) / 5 * 100
      r6 <- sum(c_across(genes_list_sulfur$Organic_sulfur_transformation$route6)) /
        length(genes_list_sulfur$Organic_sulfur_transformation$route6) * 100
      max(r1,r2,r3,r4,r5,r6)
    },
    Link_between_inorganic_and_organic_sulfur_transformation = {
      sum(c_across(genes_list_sulfur$Link_between_inorganic_and_organic_sulfur_transformation)) /
        length(genes_list_sulfur$Link_between_inorganic_and_organic_sulfur_transformation) * 100
    }
  ) %>%
  ungroup()

# 5) Ver resultados
pathway_cols <- c("Assimilatory_sulfate_reduction", "Dissimilatory_sulfur_reduction_and_oxidation",
                  "Sulfur_reduction", "SOX_systems", "Sulfur_oxidation", "Sulfur_disproportionation",
                  "Organic_sulfur_transformation", "Organic_sulfur_transformation",
                  "Link_between_inorganic_and_organic_sulfur_transformation")
nonzero_pathways <- completeness_sulfur %>%
  dplyr::select(all_of(pathway_cols)) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
  unlist() %>%
  `>`(0) %>%                # vector lógico
  which() %>%               # índices de TRUE
  names()          

completeness_df <- completeness_sulfur %>%
  dplyr::select(Sample, all_of(nonzero_pathways))


# (Opcional) Unir a combined_data:
combined_data <- combined_data %>% left_join(completeness_sulfur, by="Sample")



# 3) Calcular medias de sulfur completeness
sulfur_means <- procesar_completeness(completeness_df, metadata)

# 4) Función para graficar heatmaps de medias de completeness
crear_heatmaps <- function(means_df) {
  # long formato para cada factor
  sec <- means_df %>%
    pivot_longer(-Pathway, names_to="Sector",    values_to="Completeness") %>%
    filter(Sector %in% c("A_Open_ocean","B_Ocean_slope","C_Coastal_slope", "D_Coast"))
  
  seas <- means_df %>%
    pivot_longer(-Pathway, names_to="Season",    values_to="Completeness") %>%
    filter(Season %in% c("B_El_Nino", "A_La_Nina"))
  
  # helper
  mk <- function(df, xvar){
    ggplot(df, aes_string(x = xvar, y = "Pathway", fill = "Completeness")) +
      geom_tile(color = "white", linewidth = 0.5) +
      scale_fill_gradientn(
        colors = c("#ffd166","#619b8a","#3e5c76","#6a4c93","#a4161a"),
        limits = c(0,100),
        na.value = "grey80"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = if(xvar=="Sector") element_text() else element_blank(),
        axis.ticks.y = if(xvar=="Sector") element_line() else element_blank()
      ) +
      labs(x = NULL, y = NULL) +
      coord_fixed(ratio = 0.2)
  }
  
  p1 <- mk(sec,  "Sector")
  p2 <- mk(seas,  "Season")
  
  (p1 | p2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
}

# Generar el heatmap de nitrogen completeness
heatmap_sulfur <- crear_heatmaps(sulfur_means)
print(heatmap_sulfur)

# Ajustar la altura dinámicamente
height_factor <- 0.5  # Ajusta este valor según el espaciado que necesites
num_pathways_sulfur <- 7

# Guardar las figuras con dimensiones ajustadas
svglite("heatmap_SULFUR_completeness.svg", width = 8, height = max(8, num_pathways_sulfur * height_factor))
heatmap_sulfur
dev.off()







### PHYSICOCHEMICAL

library(ARTool)

# Unir las tablas de abundancia y metadatos
combined_data <- merge(drug_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(metal_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(nitrogen_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(sulfur_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(virulence_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")

# Ajusta el modelo con los tres factores y todas sus interacciones
combined_data$Sector <- as.factor(combined_data$Sector)
combined_data$Season <- as.factor(combined_data$Season)

# Ajusta el modelo con los tres factores y todas sus interacciones
art_mod <- art(Sand ~ Sector * Season, data = combined_data)
anova(art_mod)
as.data.frame(anova(art_mod)) %>%
  tibble::rownames_to_column("Term") %>%
  mutate(
    label = paste0(
      Term, ": F(", Df, ",", Df.res, ")=", round(`F value`,2),
      ", p=", signif(`Pr(>F)`,2)
    )
  )

# Luego en ggplot, añadir con annotate() o patchwork::plot_annotation()
plot_annotation(
  title = "Sand abundance: ART( Sector×Season×Habitat )",
  subtitle = paste(anova_tab$label, collapse = "  |  ")
)


art_mod <- art(Silt ~ Sector * Season * Habitat, data = combined_data)
anova(art_mod)
art_mod <- art(pH ~ Sector * Season * Habitat, data = combined_data)
anova(art_mod)
art_mod <- art(TOC_umol ~ Sector * Season * Habitat, data = combined_data)
anova(art_mod)
art_mod <- art(TN_umol ~ Sector * Season * Habitat, data = combined_data)
anova(art_mod)
art_mod <- art(NH4_g ~ Sector * Season * Habitat, data = combined_data)
anova(art_mod)
art_mod <- art(NO3_ug ~ Sector * Season * Habitat, data = combined_data)
anova(art_mod)
art_mod <- art(NO2_ug ~ Sector * Season * Habitat, data = combined_data)
anova(art_mod)
art_mod <- art(Fe_III_mg ~ Sector * Season * Habitat, data = combined_data)
anova(art_mod)
art_mod <- art(Fe_II_mg ~ Sector * Season * Habitat, data = combined_data)
anova(art_mod)









