#tentukan folder kerja 
setwd("D:/DMBP_Undip/Pelatihan_Lab_Data_COI/phyloseq") # silakan ganti sesuai address file masing-masing PC

# install package yang dibutuhkan namun jika sudah terinstall bisa lewatkan tahapan ini
install.packages("tidyverse")

if (!requireNamespace("devtools", quietly = TRUE))
  install.packages('devtools')

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

BiocManager::install("phyloseq")

install.packages('vegan')

devtools::install_github("gauravsk/ranacapa")

# memanggil package yang telah terinstal
library(tidyverse)
library(phyloseq)
library(vegan)
library(ranacapa)

# Read in OTU table
otu_table=read.csv("otu_table_galaxy.csv",sep=",",row.names=1)
otu_table=as.matrix(otu_table)

#   Read in taxonomy
# Seperate by kingdom phylum, class, order, family, genus, species
taxonomy <- read.csv("taxonomy_galaxy.csv",header=T, sep=",", row.names=1)
taxonomy=as.matrix(taxonomy)

#   Read in metadata
metadata=read.csv("map.csv",sep=",",row.names=1)

#   Import as phyloseq objects
OTU = otu_table(otu_table, taxa_are_rows=TRUE)
TAX = tax_table(taxonomy)
META = sample_data(metadata)

workshop = phyloseq(OTU, TAX, META)
workshop

#Remove Singletons and Zeros
workshop_clean = filter_taxa(workshop, function(x) sum(x) > 1, TRUE)
workshop_clean

ssum <- sum(sample_sums(workshop_clean)); ssum

#subset taxa
workshop_clean_1 <- workshop_clean%>% subset_taxa(
  Kingdom   == "Eukaryota"); workshop_clean_1

workshop_clean_2 <- workshop_clean_1 %>% subset_taxa(
  Phylum   == "Chordata" |
    Phylum   == "Mollusca" |
    Phylum   == "Echinodermata"); workshop_clean_2

#hapus taxa yang tidak diinginkan
workshop_clean_3 <- workshop_clean_2 %>% subset_taxa(
  Order   != "Primates" &
    Order   != "Rodentia" &
    Class   != "Aves" &
    Class   != "Amphibia" &
    Family   != "Bovidae"); workshop_clean_3

ssum <- sum(sample_sums(workshop_clean_3)); ssum
sample_sums(workshop_clean_3)

# rarefraction curve
plot(ranacapa::ggrare(workshop_clean_3, color = "InputFileName", se = FALSE, parallel=TRUE)) + geom_line(size= 1.5)

#----------------------
### 1. Stacked Barplot (Class)

workshop_phylum <- workshop_clean_3 %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Class)                                      # Sort data frame alphabetically by phylum

phylum_colors <- c("#FFFF00","#CBD588","#5F7FC7","orange","#DA5724","#CD9BCD","#AD6F3B","#990000","#673770","#D14285","#C84248","#8569D5","#5E738F","#599861")

ggplot(workshop_phylum, aes(x = Sample_Origin, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = phylum_colors) +
  theme(axis.text.x = element_text(size=7, angle=0, hjust=0.5, vjust=0)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  theme(legend.position="right") +
  ylab("Taxa Abundance")

#### Alpha Diversity
p1 <- plot_richness(workshop_clean_3, x= "Sample_Origin", measures= c("Observed", "Shannon", "Simpson")) +
  theme(axis.text.x = element_text(angle = 0)) +
  geom_boxplot() +
  theme(legend.position="bottom") ; p1

#### Beta Diversity
## Bray-curtis INDEX
set.seed(1)

bray <- phyloseq::distance(workshop_clean_3, "bray")
bray.ord <-ordinate(workshop_clean_3, method = "NMDS", bray)
sampledf <- data.frame(sample_data(workshop_clean_3))

## Bray-sample origin
Adonis2 <- adonis2(bray ~ Sample_Origin, data = sampledf, permutation=999); Adonis2

beta <- betadisper(bray, sampledf$Sample_Origin)
permutest(beta, pairwise = T)

plot_bray <- plot_ordination(workshop_clean_3, bray.ord, shape = "Sample_Origin")
plot_bray + geom_point(size = 5)