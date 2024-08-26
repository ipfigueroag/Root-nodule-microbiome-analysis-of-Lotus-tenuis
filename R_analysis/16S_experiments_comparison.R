### Plantas de campo - Analysis 
# Siguiendo tutoria MICROECO disponible en: https://chiliubio.github.io/microeco_tutorial/basic-class.html
### Instalaci?n y carga de librerias
#install.packages("microeco") # la principal
library(microeco)

#### CRAN packages ####
# allow more waiting time to download each package
options(timeout = 1000)
# If a package is not installed, it will be installed from CRAN
# First select the packages of interest
tmp <- c("microeco", "MASS", "GUniFrac", "ggpubr", "randomForest", "ggdendro", "ggrepel", "agricolae", "igraph", "picante", "pheatmap", "rgexf", 
         "ggalluvial", "ggh4x", "rcompanion", "FSA", "gridExtra", "aplot", "NST", "GGally", "ggraph", "networkD3", "poweRlaw", "ggtern")
# Now check or install
for(x in tmp){
  if(!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
  }
}

## Bioconductor packages ##
install.packages("BiocManager")
install.packages("file2meco", repos = BiocManager::repositories())
install.packages("MicrobiomeStat", repos = BiocManager::repositories())
install.packages("WGCNA", repos = BiocManager::repositories())
BiocManager::install("ggtree")
BiocManager::install("metagenomeSeq")
BiocManager::install("ALDEx2")
BiocManager::install("ANCOMBC")
library(BiocManager); library(file2meco); library(metagenom)

## Github packages ##
# downloading link of the compressed archive
url <- "https://github.com/ChiLiubio/microeco_dependence/releases/download/v0.17.0/microeco_dependence.zip"
# allow more time to download the zip file in R
options(timeout = 2000)
# If failed or too slow, please open the upper url in browser or other tool to download the zip file and move it to the current R working directory
download.file(url = url, destfile = "microeco_dependence.zip")
# uncompress the file in R
tmp <- "microeco_dependence"
unzip(paste0(tmp, ".zip"))
# install devtools
if(!require("devtools", character.only = TRUE)){install.packages("devtools", dependencies = TRUE)}
# run these one by one
devtools::install_local(paste0(tmp, "/", "SpiecEasi-master.zip"), dependencies = TRUE)
devtools::install_local(paste0(tmp, "/", "mixedCCA-master.zip"), dependencies = TRUE)
devtools::install_local(paste0(tmp, "/", "SPRING-master.zip"), dependencies = TRUE)
devtools::install_local(paste0(tmp, "/", "NetCoMi-main.zip"), repos = BiocManager::repositories())
devtools::install_local(paste0(tmp, "/", "beem-static-master.zip"), dependencies = TRUE)
devtools::install_local(paste0(tmp, "/", "chorddiag-master.zip"), dependencies = TRUE)
devtools::install_local(paste0(tmp, "/", "ggradar-master.zip"), dependencies = TRUE)
devtools::install_local(paste0(tmp, "/", "ggnested-main.zip"), dependencies = TRUE)
devtools::install_local(paste0(tmp, "/", "ggcor-1-master.zip"), dependencies = TRUE)

## Tax4Fun ##
install.packages("RJSONIO")
install.packages(system.file("extdata", "biom_0.3.12.tar.gz", package="microeco"), repos = NULL, type = "source")
install.packages(system.file("extdata", "qiimer_0.9.4.tar.gz", package="microeco"), repos = NULL, type = "source")
install.packages(system.file("extdata", "Tax4Fun_0.3.1.tar.gz", package="microeco"), repos = NULL, type = "source")

## Tax4Fun2 ##
# Either seqinr or Biostrings package should be installed for reading and writing fasta file
install.packages("seqinr", dependencies = TRUE)
# or install Biostrings from bioconductor https://bioconductor.org/packages/release/bioc/html/Biostrings.html
# Now we show how to read the fasta file
# see https://github.com/ChiLiubio/file2meco to install file2meco
rep_fasta_path <- system.file("extdata", "rep.fna", package="file2meco")
rep_fasta <- seqinr::read.fasta(rep_fasta_path)
# or use Biostrings package
rep_fasta <- Biostrings::readDNAStringSet(rep_fasta_path)
# try to create a microtable object with rep_fasta
data("otu_table_16S")
# In microtable class, all the taxa names should be necessarily included in rep_fasta
otu_table_16S <- otu_table_16S[rownames(otu_table_16S) %in% names(rep_fasta), ]
test <- microtable$new(otu_table = otu_table_16S, rep_fasta = rep_fasta)
test

library(ggplot2)  # para los graficos

# load the example data; 16S rRNA gene amplicon sequencing dataset
data(sample_info_16S)
data(otu_table_16S)
data(taxonomy_table_16S)

###### Cargando mis propios datos #####

#lugar de trabajo actual
getwd()

#setwd("G:/Mi unidad/Doctorado UNSAM/16S_sequencing/R_analysis")
#### cargando datos ####

metadata <- read.delim2("metadata_expcomprs.txt"); head(metadata) # metadata solo para comparar experimentos
tax_data <- read.delim2("taxonomy_table_expcomprs.txt"); head(tax_data)
otu_data <- read.delim2("otu_count_table_expcomprs.txt"); head(otu_data)


library(dplyr)
row.names(tax_data)<-tax_data[,1]
row.names(otu_data)<-otu_data[,1]
head(otu_data)
head(tax_data)
row.names(metadata)<-metadata$sample.id
head(metadata)

otu_data<-otu_data[-c(1)]
tax_data<-tax_data[-c(1)]

# make the taxonomic information unified, very important
library(magrittr) # Define la funci?n %<>%
tax_data %<>% tidy_taxonomy
head(tax_data)

#Before creating microtable object, make sure that the rownames of sample information table are sample names.
metadata[1:5, ]; class(metadata)
tax_data[1:5, 1:3]; class(tax_data)
otu_data[1:5, 1:5]; class(otu_data)

#To make the OTU and sample information consistent across all files in the dataset object, we use function tidy_dataset to trim the dataset.

dataset<- microtable$new(otu_table = otu_data, sample_table = metadata, tax_table = tax_data)
class(dataset)

dataset$tidy_dataset()
print(dataset)
#Then we use sample_sums() to check the sequence numbers in each sample.
dataset$sample_sums() %>% range # Resultado: 127910 - 130082 seqs per sample


#Sometimes, in order to reduce the effects of sequencing depth on the diversity measurements, we need to perform the resampling to make the sequence number equal for each sample. The function rarefy_samples can invoke the function tidy_dataset automatically before and after the rarefying.
# As an example, we use 50000 sequences in each sample
#meco_dataset$rarefy_samples(sample.size = 50000)
#19639 OTUs were removed because they are no longer present in any sample after random subsampling ...
#57324 eran los otus de salida de MOTHUR, despu?s del RESAMPLING SE RESTAN 19639 y quedan 37685


#Then, we calculate the taxa abundance at each taxonomic rank using cal_abund().
# use default parameters
dataset$cal_abund()  ## The result is stored in dataset$taxa_abund
# return meco_dataset$taxa_abund
head(dataset$taxa_abund)
class(dataset$taxa_abund)

#The function save_abund() can be used to save the taxa abundance file to a local place easily.
dataset$save_abund(dirpath = "taxa_abund")

#Then, let's calculate the alpha diversity.
# If you want to add Faith's phylogenetic diversity, use PD = TRUE, this will be a little slow
dataset$cal_alphadiv(PD = FALSE)
print(dataset$alpha_diversity)
write.table(dataset$alpha_diversity, file = "alpha_diversity/indices_alpha_div.txt", quote = FALSE, sep = "\t")
# return dataset$alpha_diversity
class(dataset$alpha_diversity)

# save dataset$alpha_diversity to a directory
dataset$save_alphadiv(dirpath = "alpha_diversity")

#Let's go on to beta diversity with function cal_betadiv(). We provide four most frequently used indexes: Bray-curtis, Jaccard, weighted Unifrac and unweighted unifrac.
# If you do not want to calculate unifrac metrics, use unifrac = FALSE
# require GUniFrac package installed
dataset$cal_betadiv(unifrac = FALSE)
# return dataset$beta_diversity
class(dataset$beta_diversity)
# save dataset$beta_diversity to a directory
dataset$save_betadiv(dirpath = "beta_diversity")

###### Chapter 4 Composition-based class  #####
#Organizar la etiqueta Description en Sample_table para que queden en el orden l?gico que queremos
#dataset$sample_table$Type_soil.Site <- factor(dataset$sample_table$Type_soil.Site, levels = c())
#dataset$sample_table$Type_soil
###### Barplot #####
# create trans_abund object
# use 10 Phyla with the highest abundance in the dataset.

#CARGAR LAS TABLAS TAXONOMY Y ABUNDANCIA DE OTUS PARA GRAFICAR A LAS BACTERIAS EXCLUSIVAS

totalphyla<-trans_abund$new(dataset = dataset, taxrank = "Phylum")

t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10)
# t1 object now include the transformed abundance data t1$abund_data and other elements for the following plotting
library(ggplot2)

barp<-t1$plot_bar(others_color = "grey70", facet = "exp_type", xtext_keep = FALSE, legend_text_italic = FALSE, strip_text = "15", x_axis_name = "sample.id") +
  theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.title = element_text(size = 0), legend.margin=margin(0,0,10,0), legend.box.margin=margin(-10,-10,-10,-10))

#590 x 424 JPJ, PNG, TIFF
#pdf 7 x 4 landscape
#Genus
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 10)
# t1 object now include the transformed abundance data t1$abund_data and other elements for the following plotting
library(ggplot2)

bgenus<-t1$plot_bar(color_values = RColorBrewer::brewer.pal(9, "Set1"), others_color = "grey70", facet = "exp_type", xtext_keep = FALSE, legend_text_italic = FALSE, strip_text = "15", x_axis_name = "sample.id") +
  theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.title = element_text(size = 0), legend.margin=margin(0,0,10,0), legend.box.margin=margin(-10,-10,-10,-10))

#plot tern
install("ggtern")
library(ggtern)
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 10, groupmean = "site")
t1$plot_tern()

######## AHORA FILTRO, PARA ELIMINAR A PROTEOBACTERIA Y VER LAS ABUNDANCIAS DE LOS DEM?S FILOS ########
# use R subset function to filter taxa in tax_table
dataset_sin_p_proteobacteria <- clone(dataset)

dataset_sin_p_proteobacteria$tax_table %<>% base::subset(Phylum != "p__Proteobacteria")

dataset_sin_p_proteobacteria$tidy_dataset()
print(dataset_sin_p_proteobacteria)

#Then, we calculate the taxa abundance at each taxonomic rank using cal_abund().
# use default parameters
dataset_sin_p_proteobacteria$cal_abund()  ## The result is stored in dataset$taxa_abund
# return meco_dataset$taxa_abund
head(dataset_sin_p_proteobacteria$taxa_abund)
class(dataset_sin_p_proteobacteria$taxa_abund)
## phylum
t1_sin_proteob <- trans_abund$new(dataset_sin_p_proteobacteria, taxrank = "Phylum", ntaxa = 7)

#t1_sin_proteob$plot_donut(label = FALSE)
#devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE); library(ggradar)
#t1_sin_proteob$plot_radar(values.radar = c("0%", "25%", "50%"))
#install.packages("ggbreak")
library(ggbreak)

t1_sin_proteob$plot_box(group = "exp_type") +
  scale_y_continuous(breaks = round(seq(0, 12, by=2),12)) +
  scale_y_break(breaks = c(12, 80), scales = 0.5, ticklabels = c(80, 90), space = 0.1, expand = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 15), axis.title.y = element_text(angle = 90, size = 20), legend.title =  element_text(size = 0), legend.text = element_text(size = 14), axis.text = element_text(size = 12))
# 740 x 550

## Genus
t1_sin_proteob <- trans_abund$new(dataset_sin_p_proteobacteria, taxrank = "Genus", ntaxa = 8, groupmean = "exp_type")

t1_sin_proteob$plot_donut(label = FALSE)
# another way with grepl function

####### VUELVO A CARGAR LOS DATOS, LUEGO FILTRO A MESORHIZOBIUM PARA VER LAS ABUNDANCIAS DE LOS DEM?S G?NEROS ########

dataset_sin_g_Mesorhizobium <- clone(dataset)

dataset_sin_g_Mesorhizobium$tax_table %<>% base::subset(Genus != "g__Mesorhizobium")

dataset_sin_g_Mesorhizobium$tidy_dataset()
print(dataset_sin_g_Mesorhizobium)

#Then, we calculate the taxa abundance at each taxonomic rank using cal_abund().
# use default parameters
dataset_sin_g_Mesorhizobium$cal_abund()  ## The result is stored in dataset$taxa_abund
# return meco_dataset$taxa_abund
head(dataset_sin_g_Mesorhizobium$taxa_abund)
class(dataset_sin_g_Mesorhizobium$taxa_abund)
## Genus
t1_sin_meso <- trans_abund$new(dataset_sin_g_Mesorhizobium, taxrank = "Genus", ntaxa = 7)

library(ggpubr)
#a<-t1_sin_meso$plot_donut(label = FALSE, facet_nrow = 2)

t1_sin_meso$plot_box(group = "exp_type") +
  scale_y_continuous(breaks = round(seq(0, 30, by=10),30)) +
  scale_y_break(breaks = c(3.5, 4, 30, 35), scales = c(0.4), ticklabels = c(40, 60, 80), space = 0.1, expand = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 15), axis.title.y = element_text(angle = 90, size = 20), legend.title =  element_text(size = 0), legend.text = element_text(size = 14), axis.text = element_text(size = 12))
# 740 x 550

t1_sin_meso$plot_pie(facet_nrow = 3)
## ternary plot
t1_sin_meso <- trans_abund$new(dataset = dataset_sin_g_Mesorhizobium, taxrank = "Genus", ntaxa = 10, groupmean = "site")
t1_sin_meso$plot_tern()

#HASTA AQU? PEGU? DEL SCRIPT DE SANTAMARIA

#t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 15)
#box<-t1$plot_box(group = "Type_soil.Site", xtext_size = 11) +
#  theme(legend.position = c(0.8, 0.7), legend.title =  element_text(size = 0), axis.text.x = element_text(angle = 45, hjust=1), axis.title.y = element_text(size = 14), panel.background = element_rect(fill = "white"), panel.border = element_rect(linetype = "solid", fill = NA), plot.margin = unit(c(0.3,0.3,0.1,1.0), "cm"))
#590 x 424

#t1 <- trans_abund$new(dataset = meco_dataset, taxrank = "Genus", ntaxa = 40)
#t1$plot_heatmap(facet = "Salinity.level", xtext_keep = FALSE, withmargin = FALSE)

#write.table(dataset$otu_table, file = "dataset_OTUtableabundance.txt", quote = FALSE, sep = "\t")
#write.table(dataset$tax_table, file = "dataset_TAXtableabundance.txt", quote = FALSE, sep = "\t")
#Then, we calculate the taxa abundance at each taxonomic rank using cal_abund().
# use default parameters


# Escoger los bacterias nitrificantes
####### AOB ### Gr?ficos de abund relativa descriptivos #####
A <- t1$data_abund
A %<>% filter(grepl("\\bNitrosomonas\\b|\\bNitrosospira\\b|\\bNitrosococcus\\b|\\bNitrosolobus\\b|\\bNitrosovibrio\\b", A$Taxonomy, ignore.case = TRUE))

A$Sample <- factor(A$Sample, levels = c("Low", "Medium", "High"))

a<-ggplot(data=A, aes(x=Sample, y=Abundance, fill=Taxonomy)) +
  geom_bar(stat="identity", position = "fill") +
  #Change the y axis labels as percent
  scale_y_continuous(labels = scales::percent_format(suffix = ""), expand = expand_scale(mult = c(0, 0))) +
  ylab("\nRelative abundance of\n AOB phylotypes (%)") + xlab("") +
  theme(legend.position = "bottom", axis.text = element_text(size = 14, colour = "black"), axis.title.y = element_text(size = 18), legend.text = element_text(size = 10)) +
  scale_fill_discrete(name=NULL) + 
  theme(panel.border = element_rect(linetype = "solid", fill = NA), plot.margin = unit(c(1.5,0.5,0.3,0.3), "cm")) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

#400x380

NOB <- t1$data_abund
NOB %<>% filter(grepl("\\bNitrospira\\b|\\bNitrococcus\\b|\\bNitrotoga\\b|\\bNitrospina\\b|\\Nitrobacter\\b|\\Nitrolancetus\\b|\\bNitrolancea\\b", NOB$Taxonomy, ignore.case = TRUE))

NOB$Sample <- factor(NOB$Sample, levels = c("Low", "Medium", "High"))

b<-ggplot(data=NOB, aes(x=Sample, y=Abundance, fill=Taxonomy)) +
  geom_bar(stat="identity", position = "fill") +
  #Change the y axis labels as percent
  scale_y_continuous(labels = scales::percent_format(suffix = ""), expand = expand_scale(mult = c(0, 0))) +
  ylab("\nRelative abundance of\n NOB phylotypes (%)") + xlab("") +
  theme(legend.position = "bottom", axis.text = element_text(size = 14, colour = "black"), axis.title.y = element_text(size = 18), legend.text = element_text(size = 10)) +
  scale_fill_discrete(name=NULL) + 
  theme(panel.border = element_rect(linetype = "solid", fill = NA),  plot.margin = unit(c(1.5,0.5,0.3,0.3), "cm")) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

ab<-ggarrange(a, b, labels = c("A", "B"), ncol = 2, nrow = 1)
annotate_figure(ab)
# 750 x 400

####### sin calcular ninguna abundancia, hacer lo siguiente:
aob_com <- clone(meco_dataset)

aob_com$tax_table %<>% .[grepl("Nitrosomonas|Nitrosospira|Nitrosococcus|Nitrosolobus|Nitrosovibrio", .$Genus), ]
aob_com$tax_table$Genus <- "g__AOB"
#a <- aob_com$tax_table
#a %<>% filter(grepl("\\bNitrosomonas\\b|\\bNitrosospira\\b|\\bNitrosococcus\\b|\\bNitrosolobus\\b|\\bNitrosovibrio\\b", a$Taxonomy, ignore.case = TRUE))

aob_t1 <- trans_abund$new(dataset = aob_com, taxrank = "Genus")

aob_t1$plot_box(group = "Salinity.level", xtext_angle = 30)

aob_t1 <- trans_diff$new(dataset = aob_com, method = "lefse", group = "Salinity.level", metastat_taxa_level = "Genus", p_adjust_method = "none")


library(RColorBrewer)
ggplot(data=b, aes(x=Group, y=freq, fill=Group)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_brewer(palette="Dark2") +
  geom_errorbar(aes(ymin=freq, ymax=freq+se), width=.2) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.2))) + 
  theme(panel.border = element_rect(linetype = "solid", fill = NA)) +
  ylab("Relative frequency of AOB\n in total reads (%)") + xlab("") +
  theme(legend.position = "none", axis.text = element_text(size = 14, colour = "black"), axis.title.y = element_text(size = 18))


c <- t5$data_abund
c %<>% filter(grepl("\\bNitrospira\\b|\\bNitrococcus\\b|\\bNitrotoga\\b|\\bNitrospina\\b|\\Nitrobacter\\b|\\Nitrolancetus\\b|\\bNitrolancea\\b", c$Taxonomy, ignore.case = TRUE))

abunmedc <- c(c$Abundance[1:3], c$Abundance[10:12])
abunlowc <- c(c$Abundance[4:6], c$Abundance[13:15])
abunhighc <- c(c$Abundance[7:9], c$Abundance[16:18])
semc <- sd(abunmedc) / sqrt(length(abunmedc))
selc <- sd(abunlowc) / sqrt(length(abunlowc))
sehc <- sd(abunhighc) / sqrt(length(abunhighc))

d <- cbind(c(c$Group[1], c$Group[4], c$Group[7]), rbind(mean(c(c$Abundance[1:3], c$Abundance[10:12], c$Abundance[19:21], c$Abundance[28:30])), mean(c(c$Abundance[4:6], c$Abundance[13:15], c$Abundance[22:24], c$Abundance[31:33])), mean(c(c$Abundance[7:9], c$Abundance[16:18], c$Abundance[25:27], c$Abundance[34:36]))), c(semc, selc, sehc))
colnames(d) <- c("Group", "freq", "se"); d <- as.data.frame(d)
d$Group <- c("Medium", "Low", "High")
d$Group<- factor(d$Group, levels = c("Low", "Medium", "High"))

ggplot(data=d, aes(x=Group, y=freq, fill=Group)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_brewer(palette="Dark2") +
  geom_errorbar(aes(ymin=freq, ymax=freq+se), width=.2) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.2))) + 
  theme(panel.border = element_rect(linetype = "solid", fill = NA)) +
  ylab("Relative frequency of NOB\n in total reads (%)") + xlab("") +
  theme(legend.position = "none", axis.text = element_text(size = 14, colour = "black"), axis.title.y = element_text(size = 18))


######## an?lisis diferencial AOB: g?nero por g?nero ######
aob <- clone(meco_dataset)

# use R subset function to filter taxa in tax_table
#aob$tax_table %<>% base::subset(Genus == "g__Nitrosomonas" | Genus == "g__Nitrosospira" | Genus == "g__Nitrosococcus" | Genus == "g__Nitrosolobus" | Genus == "g__Nitrosovibrio")
# another way with grepl function
aob$tax_table %<>% .[grepl("Nitrosomonas|Nitrosospira|Nitrosococcus|Nitrosolobus|Nitrosovibrio", .$Genus), ]
aob

#aob$tidy_dataset()
print(aob)

aob_t1 <- trans_abund$new(dataset = aob, taxrank = "Genus", groupmean = "Salinity.level")

aob_t1 <- trans_diff$new(dataset = aob, method = "wilcox", group = "Salinity.level", metastat_taxa_level = "Genus", p_adjust_method = "none")

aob_t1$plot_diff_abund(group_order = c("Low", "Medium", "High"), color_values = c("#1B9E77", "#D95F02", "#7570B3"), coord_flip = F, add_sig = TRUE, select_taxa = c("g__Nitrosomonas", "g__Nitrosococcus")) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black"), axis.title = element_text(size = 18), legend.text = element_text(size = 13), legend.position = c(0.78, 0.15), legend.title = element_text(size = 0), legend.background = element_rect(fill = 0 ))  +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)))

#### NOB

nob <- clone(meco_dataset)

# use R subset function to filter taxa in tax_table
#aob$tax_table %<>% base::subset(Genus == "g__Nitrosomonas" | Genus == "g__Nitrosospira" | Genus == "g__Nitrosococcus" | Genus == "g__Nitrosolobus" | Genus == "g__Nitrosovibrio")
# another way with grepl function
nob$tax_table %<>% .[grepl("Nitrospira|Nitrococcus|Nitrotoga|Nitrospina|Nitrobacter|Nitrolancetus|Nitrolancea", .$Genus), ]

nob

print(nob)
install.packages("randomForest"); library("randomForest")
nob_t1 <- trans_abund$new(dataset = nob, taxrank = "Genus")
nob_t1$plot_box(group = "Salinity.level")

nob_t1 <- trans_diff$new(dataset = nob, method = "lefse", group = "Salinity.level", alpha = 1, taxa_level = "Genus") #nob_t1 en realidad contiene todo, NOB y AOB community

#nob_t1$plot_diff_abund(add_sig = "TRUE", select_taxa = "g__Nitrospira", add_sig_label = "P.adj", coord_flip = F)
nob_t1$plot_diff_bar()

### Creando los boxplot ###

res_diff<- nob_t1$res_diff
res_diff %<>% filter(grepl("Nitrospira|Nitrococcus|Nitrosococcus|Nitrospina|Nitrosomonas|Nitrolancea", res_diff$Taxa, ignore.case = TRUE))

res_abund<- nob_t1$res_abund
res_abund %<>% filter(grepl("Nitrospira|Nitrococcus|Nitrosococcus|Nitrospina|Nitrosomonas|Nitrolancea", res_abund$Taxa, ignore.case = TRUE))

write.table(res_diff, file = "Differential_abundance/res_diff_nobaob.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
write.table(res_abund, file = "Differential_abundance/res_abund_nobaob.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

### boxplots ###
nob_t1 <- trans_abund$new(dataset = nob, taxrank = "Genus")
nob_t1$plot_box(group = "Salinity.level", show_point = TRUE)
data_abund<- nob_t1$data_abund #
data_abund %<>% filter(grepl("Nitrospira|Nitrococcus|Nitrosococcus|Nitrospina|Nitrosomonas|Nitrolancea", data_abund$Taxonomy, ignore.case = TRUE))
write.table(data_abund, file = "Differential_abundance/data_abund.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
nob_t1$data_abund<-data_abund
b<-nob_t1$plot_box(group = "Salinity.level")

nitrifiers <- read.delim("Differential_abundance/box_plot_nitrifiers.txt"); head(nitrifiers)
#nitrospira
nitrospira <- nitrifiers
nitrospira %<>% filter(grepl("Nitrospira", nitrospira$Taxonomy, ignore.case = TRUE))


install.packages("devtools")
library(devtools)
install_github("kassambara/easyGgplot2")
library(easyGgplot2)

nitrospira$Site<-factor(nitrospira$Site, levels = c("Low", "Medium", "High"))
library(ggpubr)

nob1<-ggplot(nitrospira, aes(x=Taxonomy, y=Abundance, fill=Site))+
  stat_boxplot(geom="errorbar", width = 0.3, position = position_dodge(width = 0.75)) +
  geom_boxplot(show.legend = "none") +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) + 
  stat_compare_means(label = "p.format", label.x = 1.5, label.y = 1.35) +
  theme(axis.text.x = element_text(size = 15, color = "black"), axis.title.x = element_text(size = 0), axis.title.y = element_text(size = 12), panel.background = element_rect(fill = 'white', color = 'black'), panel.grid = element_blank())

# 340 x 300 jpg, tiff, png
# pdf 3.97 x 2.81 lanscape

#nitrospina
nitrospina <- nitrifiers
nitrospina %<>% filter(grepl("Nitrospina", nitrospina$Taxonomy, ignore.case = TRUE))

nitrospina$Site<-factor(nitrospina$Site, levels = c("Low", "Medium", "High"))

nob2<-ggplot(nitrospina, aes(x=Taxonomy, y=Abundance, fill=Site))+
  stat_boxplot(geom="errorbar", width = 0.3, position = position_dodge(width = 0.75)) +
  geom_boxplot(show.legend = "none") +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) + 
  stat_compare_means(label = "p.format", label.x = 1.5, label.y = 0.008) +
  theme(axis.text.x = element_text(size = 15, color = "black"), axis.title.x = element_text(size = 0), axis.title.y = element_text(size = 12), panel.background = element_rect(fill = 'white', color = 'black'), panel.grid = element_blank())

# 350 x 300 jpg, tiff, png
# pdf 3.97 x 2.81 lanscape

#nitrosococcus
nitrococcus <- nitrifiers
nitrococcus %<>% filter(grepl("Nitrococcus", nitrococcus$Taxonomy, ignore.case = TRUE))

nitrococcus$Site<-factor(nitrococcus$Site, levels = c("Low", "Medium", "High"))

nob3<-ggplot(nitrococcus, aes(x=Taxonomy, y=Abundance, fill=Site))+
  stat_boxplot(geom="errorbar", width = 0.3, position = position_dodge(width = 0.75)) +
  geom_boxplot(show.legend = "none") +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) + 
  stat_compare_means(label = "p.format", label.x = 1.5, label.y = 0.035) +
  theme(axis.text.x = element_text(size = 15, color = "black"), axis.title.x = element_text(size = 0), axis.title.y = element_text(size = 12), panel.background = element_rect(fill = 'white', color = 'black'), panel.grid = element_blank())

# 350 x 300 jpg, tiff, png
# pdf 3.97 x 2.81 lanscape

#nitrolancea
nitrolancea <- nitrifiers
nitrolancea %<>% filter(grepl("Nitrolancea", nitrolancea$Taxonomy, ignore.case = TRUE))

nitrolancea$Site<-factor(nitrolancea$Site, levels = c("Low", "Medium", "High"))

nob4<-ggplot(nitrolancea, aes(x=Taxonomy, y=Abundance, fill=Site))+
  stat_boxplot(geom="errorbar", width = 0.3, position = position_dodge(width = 0.75)) +
  geom_boxplot(show.legend = "none") +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) + 
  stat_compare_means(label = "p.format", label.x = 1.5, label.y = 0.0032) +
  theme(axis.text.x = element_text(size = 15, color = "black"), axis.title.x = element_text(size = 0), axis.title.y = element_text(size = 12), panel.background = element_rect(fill = 'white', color = 'black'), panel.grid = element_blank())
# 350 x 300 jpg, tiff, png
# pdf 3.97 x 2.81 lanscape
nob1234<-ggarrange(nob1, nob2, nob3, nob4, 
          labels = c("E", "F", "G", "H"),
          ncol = 2, nrow = 2)
annotate_figure(nob1234, top = text_grob("NOB community", color = "black", face = "bold", size = 15))

# 500 x 440 jpg, tiff, png
# pdf 3.97 x 2.81 lanscape

#### AOB community 
nitrosomonas <- nitrifiers
nitrosomonas %<>% filter(grepl("Nitrosomonas", nitrosomonas$Taxonomy, ignore.case = TRUE))

nitrosomonas$Site<-factor(nitrosomonas$Site, levels = c("Low", "Medium", "High"))

aob2<-ggplot(nitrosomonas, aes(x=Taxonomy, y=Abundance, fill=Site))+
  stat_boxplot(geom="errorbar", width = 0.3, position = position_dodge(width = 0.75)) +
  geom_boxplot(show.legend = "none") +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) + 
  stat_compare_means(label = "p.format", label.x = 1.5, label.y = 0.021) +
  theme(axis.text.x = element_text(size = 15, color = "black"), axis.title.x = element_text(size = 0), axis.title.y = element_text(size = 12), panel.background = element_rect(fill = 'white', color = 'black'), panel.grid = element_blank())

# 350 x 300 jpg, tiff, png
# pdf 3.97 x 2.81 lanscape

#nitrosococcus
nitrosococcus <- nitrifiers
nitrosococcus %<>% filter(grepl("Nitrosococcus", nitrosococcus$Taxonomy, ignore.case = TRUE))

nitrosococcus$Site<-factor(nitrosococcus$Site, levels = c("Low", "Medium", "High"))

aob1<-ggplot(nitrosococcus, aes(x=Taxonomy, y=Abundance, fill=Site))+
  stat_boxplot(geom="errorbar", width = 0.3, position = position_dodge(width = 0.75)) +
  geom_boxplot(show.legend = "none") +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) + 
  stat_compare_means(label = "p.format", label.x = 1.5, label.y = 0.036) +
  theme(axis.text.x = element_text(size = 15, color = "black"), axis.title.x = element_text(size = 0), axis.title.y = element_text(size = 12), panel.background = element_rect(fill = 'white', color = 'black'), panel.grid = element_blank())
# 350 x 300 jpg, tiff, png
# pdf 3.97 x 2.81 lanscape
aob12<-ggarrange(aob1, aob2, 
                   labels = c("C", "D"),
                   ncol = 1, nrow = 2)
annotate_figure(aob12, top = text_grob("AOB community", color = "black", face = "bold", size = 15))

# 500 x 440 jpg, tiff, png
# pdf 5.00 x 4.00 lanscape


ggplot(data=nitrifiers, aes(x=Taxa, y=Mean, fill=Group)) +
  geom_boxplot(stat="identity", color="black") +
  scale_fill_brewer(palette="Dark2") +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.2))) + 
  theme(panel.border = element_rect(linetype = "solid", fill = NA)) +
  ylab("Relative frequency of AOB\n in total reads (%)") + xlab("") +
  theme(legend.position = "none", axis.text = element_text(size = 14, colour = "black"), axis.title.y = element_text(size = 18))


nob_t1 <- trans_diff$new(dataset = nob, method = "metastat", group = "Salinity.level", p_adjust_method = "none", "Genus")

nob_t1$plot_diff_abund(group_order = c("Low", "Medium", "High"), color_values = c("#1B9E77", "#D95F02", "#7570B3"), coord_flip = F, add_sig = TRUE) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black"), axis.title = element_text(size = 18), legend.text = element_text(size = 13), legend.position = c(0.78, 0.15), legend.title = element_text(size = 0), legend.background = element_rect(fill = 0 ))  +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)))


#Then, let's calculate the alpha diversity.
# If you want to add Faith's phylogenetic diversity, use PD = TRUE, this will be a little slow
meco_dataset$cal_alphadiv(PD = FALSE)
print(meco_dataset$alpha_diversity)
write.table(meco_dataset$alpha_diversity, file = "alpha_diversity/indices_alpha_div.txt", quote = FALSE, sep = "\t")
# return dataset$alpha_diversity
class(meco_dataset$alpha_diversity)

# save dataset$alpha_diversity to a directory
dataset$save_alphadiv(dirpath = "alpha_diversity")

#Let's go on to beta diversity with function cal_betadiv(). We provide four most frequently used indexes: Bray-curtis, Jaccard, weighted Unifrac and unweighted unifrac.
# If you do not want to calculate unifrac metrics, use unifrac = FALSE
# require GUniFrac package installed
meco_dataset$cal_betadiv(unifrac = FALSE)
# return dataset$beta_diversity
class(meco_dataset$beta_diversity)
# save dataset$beta_diversity to a directory
meco_dataset$save_betadiv(dirpath = "Ciclo_nitrogeno_beta_diversity")


##### Correlation plot Beta-diversity and nitrifiers community ######
meco_dataset_nitrifiers <- clone(meco_dataset)
nitrifiers <- meco_dataset_nitrifiers$tax_table
nitrifiers %<>% filter(grepl("\\bg__Nitrospira\\b|\\bg__Nitrococcus\\b|\\bg__Nitrotoga\\b|\\bg__Nitrospina\\b|\\Nitrobacter\\b|\\g__Nitrolancetus\\b|\\bg__Nitrolancea\\b|\\bNitrosomonas\\b|\\bNitrosospira\\b|\\bNitrosococcus\\b|\\bNitrosolobus\\b|\\bNitrosovibrio\\b", nitrifiers$Genus, ignore.case = TRUE))

meco_dataset_nitrifiers$tax_table <- nitrifiers

meco_dataset_nitrifiers$cal_betadiv(unifrac = FALSE)
# return dataset$beta_diversity
class(meco_dataset$beta_diversity)

meco_dataset_nitrifiers$beta_diversity$bray

bray_sampletable <-cbind(meco_dataset_nitrifiers$beta_diversity$bray, meco_dataset_nitrifiers$sample_table)

#corbrayn <- (meco_dataset$beta_diversity %<>% subset (Measure == "bray"))

ggplot(data = bray_sampletable, aes(x = CE, y = '2A')) + 
  geom_point(size = 3) + 
  geom_smooth(method = "glm") +
  stat_cor(size = 5, color = "black") +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black"), axis.title = element_text(size = 15), legend.position = "none") +
  labs(y= "Nitrifiers community \n(Bray-Curtis similarity)", x = "Salinity")



##### Box plots #####
# show 15 taxa at Class level
t1 <- trans_abund$new(dataset = meco_dataset, taxrank = "Class", ntaxa = 15)
t1$plot_box(group = "Salinity.level")
# show 15 taxa at Order level
t1 <- trans_abund$new(dataset = meco_dataset, taxrank = "Order", ntaxa = 15)
t1$plot_box(group = "Salinity.level")

t1 <- trans_abund$new(dataset = meco_dataset, taxrank = "Family", ntaxa = 15)
t1$plot_box(group = "Salinity.level")

# show 40 taxa at Genus level
t1 <- trans_abund$new(dataset = meco_dataset, taxrank = "Genus", ntaxa = 15)
box<-t1$plot_box(group = "Salinity.level", xtext_size = 11) +
  theme(legend.position = c(0.8, 0.7), legend.title =  element_text(size = 0), axis.text.x = element_text(angle = 45, hjust=1), axis.title.y = element_text(size = 14), panel.background = element_rect(fill = "white"), panel.border = element_rect(linetype = "solid", fill = NA), plot.margin = unit(c(0.3,0.3,0.1,1.0), "cm"))

t1 <- trans_abund$new(dataset = meco_dataset, taxrank = "Genus", ntaxa = 40)
t1$plot_heatmap(facet = "Salinity.level", xtext_keep = FALSE, withmargin = FALSE)

#####Then, we show the pie chart with the group mean values #####
t1 <- trans_abund$new(dataset = meco_dataset, taxrank = "Order", ntaxa = 6, groupmean = "Salinity.level")
# all pie chart in one row
t1$plot_pie(facet_nrow = 1)


##### 4.2 trans_venn class #####
# merge samples as one community for each group
#dataset1 <- meco_dataset$merge_samples(use_group = "Salinity.level")
# dataset1 is a new microtable object
# create trans_venn object, see other available options in ratio

##### VENN DIAGRAM  #####

dataset1 <- dataset$merge_samples(group = "experiment")
t1 <- trans_venn$new(dataset1, ratio = "numratio")
v<-t1$plot_venn(text_size = 5, text_name_size = 5, petal_plot = FALSE) 
  
##### transform venn results to the sample-species table, here do not consider abundance, only use presence/absence.#####
t2 <- t1$trans_comm(use_frequency = FALSE)
# t2 is a new microtable class, each part is considered a sample
t2$otu_table <- t2$otu_table[,1:9]
class(t2)
# calculate taxa relative abundance, that is, the relative frequency
t2$cal_abund() #Aqu? la abundancia relativa ya est? en %
# transform and plot
t2$sample_table$Group <- c("Medium", "Medium", "Medium", "Low", "Low", "Low", "High", "High", "High")
t2$sample_table$Group<- factor(t2$sample_table$Group, levels = c("Low", "Medium", "High"))

#write.table(t2$taxa_abund$Genus, file = "taxa_abund/taxa_abund_t2.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

t4 <- clone(t2)

t5 <- trans_abund$new(dataset = t4, taxrank = "Genus"); head(t5$data_abund)

a <- t5$data_abund
a %<>% filter(grepl("\\bNitrosomonas\\b|\\bNitrosospira\\b|\\bNitrosococcus\\b|\\bNitrosolobus\\b|\\bNitrosovibrio\\b", a$Taxonomy, ignore.case = TRUE))

abunmed <- c(a$Abundance[1:3], a$Abundance[10:12])
abunlow <- c(a$Abundance[4:6], a$Abundance[13:15])
abunhigh <- c(a$Abundance[7:9], a$Abundance[16:18])
sem <- sd(abunmed) / sqrt(length(abunmed))
sel <- sd(abunlow) / sqrt(length(abunlow))
seh <- sd(abunhigh) / sqrt(length(abunhigh))

b <- cbind(c(a$Group[1], a$Group[4], a$Group[7]), rbind(mean(c(a$Abundance[1:3], a$Abundance[10:12])), mean(c(a$Abundance[4:6], a$Abundance[13:15])), mean(c(a$Abundance[7:9], a$Abundance[16:18]))), c(sem, sel, seh))
colnames(b) <- c("Group", "freq", "se"); b <- as.data.frame(b)
b$Group <- c("Medium", "Low", "High")
b$Group<- factor(b$Group, levels = c("Low", "Medium", "High"))

library(RColorBrewer)
ggplot(data=b, aes(x=Group, y=freq, fill=Group)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_brewer(palette="Dark2") +
  geom_errorbar(aes(ymin=freq, ymax=freq+se), width=.2) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.2))) + 
  theme(panel.border = element_rect(linetype = "solid", fill = NA)) +
  ylab("Relative frequency of AOB\n in total reads (%)") + xlab("") +
  theme(legend.position = "none", axis.text = element_text(size = 14, colour = "black"), axis.title.y = element_text(size = 18))
  

c <- t5$data_abund
c %<>% filter(grepl("\\bNitrospira\\b|\\bNitrococcus\\b|\\bNitrotoga\\b|\\bNitrospina\\b|\\Nitrobacter\\b|\\Nitrolancetus\\b|\\bNitrolancea\\b", c$Taxonomy, ignore.case = TRUE))

abunmedc <- c(c$Abundance[1:3], c$Abundance[10:12])
abunlowc <- c(c$Abundance[4:6], c$Abundance[13:15])
abunhighc <- c(c$Abundance[7:9], c$Abundance[16:18])
semc <- sd(abunmedc) / sqrt(length(abunmedc))
selc <- sd(abunlowc) / sqrt(length(abunlowc))
sehc <- sd(abunhighc) / sqrt(length(abunhighc))

d <- cbind(c(c$Group[1], c$Group[4], c$Group[7]), rbind(mean(c(c$Abundance[1:3], c$Abundance[10:12], c$Abundance[19:21], c$Abundance[28:30])), mean(c(c$Abundance[4:6], c$Abundance[13:15], c$Abundance[22:24], c$Abundance[31:33])), mean(c(c$Abundance[7:9], c$Abundance[16:18], c$Abundance[25:27], c$Abundance[34:36]))), c(semc, selc, sehc))
colnames(d) <- c("Group", "freq", "se"); d <- as.data.frame(d)
d$Group <- c("Medium", "Low", "High")
d$Group<- factor(d$Group, levels = c("Low", "Medium", "High"))

ggplot(data=d, aes(x=Group, y=freq, fill=Group)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_brewer(palette="Dark2") +
  geom_errorbar(aes(ymin=freq, ymax=freq+se), width=.2) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.2))) + 
  theme(panel.border = element_rect(linetype = "solid", fill = NA)) +
  ylab("Relative frequency of NOB\n in total reads (%)") + xlab("") +
  theme(legend.position = "none", axis.text = element_text(size = 14, colour = "black"), axis.title.y = element_text(size = 18))




t3 <- trans_diff$new(dataset = t2)

outAOVAOB <-aov(Abundance ~ Group  + Taxonomy, data = AOB)

#+
  geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
                position=position_dodge(.9))
# The integer is OTU number
# The percentage data is the sequence number/total sequence number

#####Sometimes, we are interested in the composition of unique and shared species.
t1 <- trans_venn$new(dataset1, ratio = "seqratio")
# transform venn results to the sample-species table, here do not consider abundance, only use presence/absence.
t1 <- t1$trans_venn_com(use_OTUs_frequency = TRUE)
# t2 is a new microtable class, each part is considered a sample
class(t1)

########### Chapter 5 Diversity-based class #########
##### 5.1 trans_alpha class #####

t1 <- trans_alpha$new(dataset = dataset, group = "experiment")
# return t1$alpha_stat
t1$data_stat
t1$data_alpha
write.table(t1$data_stat, file = "indices_alpha_div_Estadisticos.txt", quote = FALSE, sep = "\t")
#Then, we test the differences among groups using the KW (Kruscal Wallis) rank sum test and anova with multiple comparisons.
t1$cal_diff(method = "KW")
# return t1$res_alpha_diff
t1$res_diff
write.table(t1$res_diff, file = "indices_alpha_div_KW.txt", quote = FALSE, sep = "\t")

#install.packages("agricolae")
library(agricolae)
t1$cal_diff(method = "anova")
# return t1$res_alpha_diff
t1$res_diff
write.table(t1$res_alpha_diff, file = "indices_alpha_div_ANOVA.txt", quote = FALSE, sep = "\t")

#We can also use the boxplot to show the paired comparisons directly.
#install.packages("ggpubr")
library(ggpubr)
chao<-
  t1$plot_alpha(measure = "Chao1", color_values = RColorBrewer::brewer.pal(8, "Set2"), boxplot_add = FALSE) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA), axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = c(0.15, 0.1), legend.direction = "vertical") +
  geom_boxplot(aes(fill=experiment)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_boxplot(geom="errorbar", width=0.2) +
  guides(colour = "none") +
  stat_compare_means(label = "p.format", hjust = -1, vjust = 0.0) #+
  #scale_x_discrete(labels= c("I", "I", "I2", "I2", "M", "M"))

shanon<-t1$plot_alpha(measure = "Shannon", color_values = RColorBrewer::brewer.pal(8, "Set2"), boxplot_add = FALSE) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA), legend.position = c(0.15, 0.1), legend.direction = "horizontal", legend.title = element_text(size = 0), axis.text.x = element_text(angle = 0, hjust = 0)) +
  geom_boxplot(aes(fill=experiment)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_boxplot(geom="errorbar", width=0.2) +
  guides(colour = "none") +
  stat_compare_means(label = "p.format", hjust = 0.25, vjust = 0.6) #+ scale_x_discrete(labels= c("I", "I", "I2", "I2", "M", "M"))

ace<-t1$plot_alpha(measure = "ACE", color_values = RColorBrewer::brewer.pal(8, "Set2"), boxplot_add = FALSE) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA), legend.position = c(0.15, 0.1), legend.direction = "horizontal", legend.title = element_text(size = 0)) +
  geom_boxplot(aes(fill=experiment)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_boxplot(geom="errorbar", width=0.2) +
  guides(colour = "none") +
  stat_compare_means(label = "p.format", hjust = 0.25, vjust = 0.6) #+
  #scale_x_discrete(labels= c("I", "I", "I2", "I2", "M", "M"))

simp <- t1$plot_alpha(measure = "Simpson", color_values = RColorBrewer::brewer.pal(8, "Set2"), boxplot_add = FALSE) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA), legend.position = c(0.35, 0.1), legend.direction = "horizontal", legend.title = element_text(size = 0), axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_boxplot(aes(fill=experiment)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_boxplot(geom="errorbar", width=0.2) +
  guides(colour = "none") +
  stat_compare_means(label = "p.format", hjust = 0.25, vjust = 0.6) #+
  #scale_x_discrete(labels= c("I", "I", "I2", "I2", "M", "M"))

t1$plot_alpha(measure = "Fisher", color_values = RColorBrewer::brewer.pal(8, "Set2"), boxplot_add = FALSE) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA), legend.position = c(0.15, 0.1), legend.direction = "horizontal", legend.title = element_text(size = 0)) +
  geom_boxplot(aes(fill=experiment)) +
  geom_jitter(width = 0.1) +
  stat_boxplot(geom="errorbar", width=0.2) +
  guides(colour = "none") +
  stat_compare_means(label = "p.format", hjust = 0.5, vjust = 0.6)# +
#  scale_x_discrete(labels= c("I", "I", "I2", "I2", "M", "M"))

t1$plot_alpha(measure = "Observed", color_values = RColorBrewer::brewer.pal(8, "Set2"), boxplot_add = FALSE) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA), legend.position = c(0.15, 0.1), legend.direction = "horizontal", legend.title = element_text(size = 0)) +
  geom_boxplot(aes(fill=experiment)) +
  geom_jitter(width = 0.1) +
  stat_boxplot(geom="errorbar", width=0.2) +
  guides(colour = "none") +
  stat_compare_means(label = "p.format", hjust = 0.5, vjust = 0.6)# +
  #scale_x_discrete(labels= c("I", "I", "I2", "I2", "M", "M"))

t1$plot_alpha(measure = "InvSimpson", color_values = RColorBrewer::brewer.pal(8, "Set2"), boxplot_add = FALSE) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA), legend.position = c(0.35, 0.1), legend.direction = "horizontal", legend.title = element_text(size = 0)) +
  geom_boxplot(aes(fill=experiment)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_boxplot(geom="errorbar", width=0.2) +
  guides(colour = "none") #+
#  scale_x_discrete(labels= c("I", "I", "I2", "I2", "M", "M"))

alphad <-ggarrange(simp, shanon, labels = c("A", "B"), ncol = 2, nrow = 1, common.legend = TRUE)
annotate_figure(alphad)
#692 x 390


#The plot_alpha function add the significance label by searching the results in object$res_diff instead of calculating the significance again


t1$plot_alpha(measure = "InvSimpson", add_sig_text_size = 6, order_x_mean = FALSE, add_letter = TRUE, xtext_size =15, boxplot_add = "dotplot")

t1$plot_alpha(measure = "Fisher", add_sig_text_size = 6, order_x_mean = FALSE, add_letter = TRUE, xtext_size =15, boxplot_add = "dotplot")


###### Correlation plot #####
library(ggplot2)

t1 <- trans_alpha$new(dataset = dataset, group = "experiment")

corshannon <- (t1$data_alpha %<>% subset (Measure == "Shannon"))

s<-ggplot(data = corshannon, aes(x = EC, y = Value, color = Measure)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "glm") +
  stat_cor(size = 4, color = "black") +
  theme(axis.text.x = element_text(size = 10, color = "black"), axis.text.y = element_text(size = 10, color = "black"), axis.title = element_text(size = 14), legend.position = "none", panel.border = element_rect(linetype = "solid", fill = NA),  plot.margin = unit(c(1.5,0.5,0.3,0.3), "cm"), panel.background = element_rect(fill = 'white', color = 'black')) +
  labs(y= "Shannon index", x = "Salinity")
#remover leyenda legend.position = "none"

t1 <- trans_alpha$new(dataset = dataset, group = "Site")
simpson <- (t1$data_alpha %<>% subset (Measure == "Simpson"))

s<-ggplot(data = simpson, aes(x = EC, y = Value, color = Measure)) + 
  geom_point(size = 3, color = "3") + 
  geom_smooth(method = "glm", color = "3") +
  stat_cor(size = 4, color = "black") +
  theme(axis.text.x = element_text(size = 10, color = "black"), axis.text.y = element_text(size = 10, color = "black"), legend.position = "none", axis.title = element_text(size = 14), panel.border = element_rect(linetype = "solid", fill = NA),  plot.margin = unit(c(1.5,0.5,0.3,0.3), "cm"), panel.background = element_rect(fill = 'white', color = 'black')) +
  labs(y= "Simpson index", x = "Salinity")

t1 <- trans_alpha$new(dataset = dataset, group = "Site")
observed <- (t1$data_alpha %<>% subset (Measure == "Observed"))

o<-ggplot(data = observed, aes(x = EC, y = Value, color = Measure)) + 
  geom_point(size = 3, color = "3") + 
  geom_smooth(method = "glm", color = "3") +
  stat_cor(size = 4, color = "black") +
  theme(axis.text.x = element_text(size = 10, color = "black"), axis.text.y = element_text(size = 10, color = "black"), legend.position = "none", axis.title = element_text(size = 14), panel.border = element_rect(linetype = "solid", fill = NA),  plot.margin = unit(c(1.5,0.5,0.3,0.3), "cm"), panel.background = element_rect(fill = 'white', color = 'black')) +
  labs(y= "Observed OTUS", x = "Salinity")


t1 <- trans_alpha$new(dataset = dataset, group = "Site")
fisher <- (t1$data_alpha %<>% subset (Measure == "Fisher"))

f<-ggplot(data = fisher, aes(x = EC, y = Value, color = Measure)) + 
  geom_point(size = 3, color = "3") + 
  geom_smooth(method = "glm", color = "3") +
  stat_cor(size = 4, color = "black") +
  theme(axis.text.x = element_text(size = 10, color = "black"), axis.text.y = element_text(size = 10, color = "black"), legend.position = "none", axis.title = element_text(size = 14), panel.border = element_rect(linetype = "solid", fill = NA),  plot.margin = unit(c(1.5,0.5,0.3,0.3), "cm"), panel.background = element_rect(fill = 'white', color = 'black')) +
  labs(y= "Fisher", x = "Salinity")

osf<-ggarrange(o, s, f, p, tr, v, labels = c("A", "B", "C", "D", "E", "F"), ncol = 3, nrow = 2)
annotate_figure(osf)
#850x500 jpg
# 10 x 7 landscape
##### 5.2 trans_beta class #####
# we first create an trans_beta object
# measure parameter can invoke the distance matrix of bray in dataset$beta_diversity

t1 <- trans_beta$new(dataset = dataset, group = "experiment", measure = "bray")
# use PCoA as an example, PCA or NMDS is also available
t1$cal_ordination(method = "NMDS")
# t1$res_ordination is the ordination result list
class(t1$res_ordination)

# plot the PCoA result
p<-t1$plot_ordination(plot_color = "experiment", plot_type = "point", plot_shape = "site") +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black"), axis.title = element_text(size = 14), legend.text = element_text(size = 10), legend.title = element_text(size = 15), legend.background = element_rect(fill = 0), panel.border = element_rect(linetype = "solid", fill = NA),  plot.margin = unit(c(0.5,0.5,0.3,0.3), "cm"), panel.background = element_rect(fill = 'white', color = 'black')) #+
 # stat_ellipse(geom = "polygon", aes(group = experiment, color = experiment, fill = experiment), alpha = 0.3) +
  annotate("text", label = paste0("stress: ", format(t1$res_ordination$model$stress, digits = 2)))

#nmds <- t1$plot_ordination(plot_color = "Salinity.level", plot_group_ellipse = TRUE)
#nmds + 
#  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black"), axis.title = element_text(size = 15), legend.text = element_text(size = 13), legend.position = c(0.1, 0.15), legend.title = element_text(size = 0), legend.background = element_rect(fill = 0 ))

# Clustering plot is also a frequently used method.
# use replace_name to set the label name, group parameter used to set the color
#install.packages("ggdendro")
library(ggdendro)
install("ggtree"); install("treeio")
library(ggtree)
library(treeio)

tr<-t1$plot_clustering(group = "experiment", replace_name = "site_exp_type") +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.title = element_text(size = 14), line = element_line(linewidth = 1), panel.background = element_rect(fill = 'white'), legend.position = "bottom", legend.title = element_text(size = 0), legend.text = element_text(size = 12))


#448 x 327
# calculate and plot sample distances within groups
t1$cal_group_distance(within_group = TRUE); t1$res_group_distance
# return t1$res_group_distance
# perform Wilcoxon Rank Sum and Signed Rank Tests
t1$cal_group_distance_diff(method = "wilcox"); t1$res_group_distance_diff
# plot_group_order parameter can be used to adjust orders in x axis
t1$plot_group_distance(boxplot_add = "mean") +
  theme(panel.border = element_rect(linetype = "solid", fill = NA), legend.position = "none") +
  geom_boxplot(aes(fill=experiment)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_boxplot(geom="errorbar", width=0.2) +
  guides(colour = "none")

#perMANOVA(Anderson 2001) is often used in the differential test of distances among groups.
# manova for all groups
t1$cal_manova(cal_manova_all = TRUE)
t1$res_manova

# manova for each paired groups
t1$cal_manova(cal_manova_paired = TRUE)
t1$res_manova

# manova for specified group set: here "Group + Type"
t1$cal_manova(cal_manova_set = "Site")
t1$res_manova

# PERMDISP(Anderson et al. 2011) is also implemented to check multivariate homogeneity of groups dispersions (variances).
# PERMDISP for the whole comparison and for each paired groups
t1$cal_betadisper()
t1$res_betadisper

###### trans_diff class #####

### Metastat #####
# metastat analysis at Genus level
t1 <- trans_diff$new(dataset = meco_dataset, method = "metastat", group = "experiment", metastat_taxa_level = "Genus")
# t1$res_metastat is the result
# t1$res_metastat_group_matrix is the group comparisons order for plotting
# plot the first paired groups, choose_group = 1
t1$plot_metastat(use_number = 1:10, qvalue = 0.05, choose_group = 1)

########### LEFSE #############
#LEfSe combines the non-parametric test and linear discriminant analysis (Segata et al. 2011).
#dataset$sample_table$Type_soil <- factor(meco_dataset$sample_table$Salinity.level, levels = c("Low", "Medium", "High"))
dataset$sample_table$Type_soil
#install.packages("metagenomeSeq"); library(metagenomeSeq)
#BiocManager::install("metagenomeSeq", force=TRUE)
#install.packages("agricolae"); library(agricolae)
# LEFSE
t1 <-trans_diff$new(dataset = dataset, method = "lefse", group = "Type_soil", p_adjust_method = "none")
# t1$res_lefse is the result
# t1$res_abund is the abundance information
library(ggbreak)

t1$plot_diff_abund(color_values = c("#1B9E77", "#D95F02"), group_order = c("ML", "SLL"), coord_flip = T, add_sig = TRUE, use_number = 1:7) +
  theme(axis.text.x = element_text(size = 18, color = "black"), axis.text.y = element_text(size = 20, color = "black"), axis.title = element_text(size = 18), legend.text = element_text(size = 15), legend.title = element_text(size = 0), legend.background = element_rect(fill = 0 ))  +
  scale_y_continuous(breaks = round(seq(0, 0.01, by=0.001), 0.01)) +
  scale_y_break(breaks = c(0.01, 0.1), scales = 0.5, ticklabels = c(0.01, 0.25, 0.50, 0.75, 1), space = 0.1, expand = TRUE) +
scale_x_discrete(expand = c(0, 0))

#ticklabels = c(0.50, 0.75, 1.00), 
#, breaks = round(seq(0, 0.125, by=0.05),0.125)

head(t1$res_diff); head(t1$res_abund$Taxa)

#LDA BAR
t1$plot_diff_bar(threshold = 2, color_values = c("#1B9E77", "#D95F02"), group_order = c("ML", "SLL")) +
  theme(axis.text.x = element_text(size = 18, color = "black"), axis.text.y = element_text(size = 20, color = "black"), axis.title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 0), legend.background = element_rect(fill = 0 ))


####### selecciono los grupos de bacterias a graficar ####
library(dplyr)
n1 <-t1$res_abund
n1 %<>% filter(grepl("Nitrosomonas", Taxa))
n1rf <-t1$res_diff
n1rf %<>% filter(grepl("Nitrosomonas", Taxa))


n3 <-t1$res_abund
n3 %<>% filter(grepl("Nitrosococcus", Taxa))
n3rf <-t1$res_diff
n3rf %<>% filter(grepl("Nitrosococcus", Taxa))

n13 <-merge(n1, n3, all = TRUE); n13rf <-merge(n1rf, n3rf, all = TRUE)
t3$res_abund <-n13
t3$res_diff <-n13rf

t2 <- t1
t3 <- t2

t3$res_abund <- t1$res_abund

t1$plot_diff_abund()
t3$res_diff
t3$res_abund
a<-t3$plot_diff_abund(coord_flip = F) 

t3$res_abund %>%
  ggplot( aes(x=Taxa, y=Mean, fill=Taxa)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9)

#lefse
t1 <-trans_diff$new(dataset = meco_dataset, method = 'lefse', group = "Salinity.level")
# t1$res_lefse is the LEfSe result
# t1$res_abund is the abundance information
t1$plot_diff_abund(color_values = RColorBrewer::brewer.pal(8, "Dark2"), coord_flip = F, add_sig = TRUE)
head(t1$res_diff); head(t1$res_abund$Taxa)
# selecciono los grupos de bacterias a graficar
library(dplyr)
n1 <-t1$res_abund
n1 %<>% filter(grepl("Nitrosomonas", Taxa))
n1rf <-t1$res_diff
n1rf %<>% filter(grepl("Nitrosomonas", Taxa))


n3 <-t1$res_abund
n3 %<>% filter(grepl("Nitrosococcus", Taxa))
n3rf <-t1$res_diff
n3rf %<>% filter(grepl("Nitrosococcus", Taxa))

n13 <-merge(n1, n3, all = TRUE); n13rf <-merge(n1rf, n3rf, all = TRUE)
t3$res_abund <-n13
t3$res_diff <-n13rf

t2 <- t1
t3 <- t2

t3$res_abund <- t1$res_abund

t1$plot_diff_abund()
t3$res_diff
t3$res_abund
a<-t3$plot_diff_abund(coord_flip = F) 

t3$res_abund %>%
  ggplot( aes(x=Taxa, y=Mean, fill=Taxa)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9)
#add_significance = T, group_order = c("Low", "Medium", "High"), color_values = c("#7570B3", "#D95F02", "#1B9E77"))

#### otra cosa
t1$res_lefse
t1$plot_lefse_bar(LDA_score = 3.1, group_order = c("Low", "Medium", "High"))

t1$res_lefse[1:5, ]
dim(t1$res_lefse)


# clade_label_level 5 represent phylum level in this analysis
# require ggtree package
#install.packages("tidytree")
library(tidytree)
#install.packages("ggtree"); library(ggtree)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ggtree")

t1$plot_lefse_cladogram(use_taxa_num = 200, use_feature_num = 50, clade_label_level = 5)

# choose some taxa according to the positions in the previous picture; those taxa labels have minimum overlap
use_labels <- c("p__Actinobacteria", "p__Ignavibacteriae", "c__Alphaproteobacteria", "c__Gammaproteobacteria", "c__Bacteroidia", "o__Chromatiales", "o__Euzebyales", "o__Myxococcales", "c__Ignavibacteria", "g__Desulfopila", "g__Salegentibacter", "g__Salisaeta")

# then use parameter select_show_labels to show
t1$plot_lefse_cladogram(use_taxa_num = 150, use_feature_num = 45, select_show_labels = use_labels, group_order = c("Low", "Medium", "High"))
# Now we can see that more taxa names appear in the tree

######### trans_network class ###########
## We first introduce the correlation-based network
# When the OTU number is large, using R WGCNA package to replace "base" will be faster
# require WGCNA package, see https://chiliubio.github.io/microeco_tutorial/intro.html#wgcna for the installation
library("impute")
library("WGCNA")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("impute")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("preprocessCore")
library(preprocessCore)

#t1 <- trans_network$new(dataset = meco_dataset, cal_cor = "WGCNA", taxa_level = "OTU", filter_thres = 0.0001, cor_method = "spearman")

library(devtools)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)
# SparCC method, require SpiecEasi package, see https://chiliubio.github.io/microeco_tutorial/intro.html#spieceasi for the installation
# SparCC is very slow, so consider filtering more species with low abundance
###### Required packages #####
install.packages("devtools")
install.packages("BiocManager")

# Install NetCoMi
devtools::install_github("stefpeschel/NetCoMi", 
                         dependencies = c("Depends", "Imports", "LinkingTo"),
                         repos = c("https://cloud.r-project.org/",
                                   BiocManager::repositories()))

devtools::install_github("zdk123/SpiecEasi")
devtools::install_github("GraceYoon/SPRING")


########### Creando la red ##############

dataset

#red <- trans_network$new(dataset = dataset,cor_method = "spearman", taxa_level = "OTU")

red <- trans_network$new(dataset = dataset, cor_method = "pearson", filter_thres = 0.00001) #se crea la red

red$res_cor_p
head(red$res_cor_p$cor)
head(red$res_cor_p$p)

# construct network; require igraph package
red$cal_network(COR_p_thres = 0.01, add_taxa_name = "Genus", COR_cut = 0.2)
#Then, let's partition modules for the network.
red$cal_module(method = "cluster_fast_greedy")

# require rgexf package to be installed
red$save_network(filepath = "network_genus_corcut0.2.gexf")

library(igraph)

# construct network; require igraph package
#red$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE, add_taxa_name = "Phylum")
#The optimized COR threshold: 0.24...
#The parameter COR_cut can be used to select the correlation threshold. Furthermore, COR_optimization = TRUE can be used to find the optimized coefficient threshold (potential transition point of network eigenvalues) instead of the COR_cut based on the RMT theory (Deng et al. 2012).
# return 
red$res_network

# add modules in the network
# invoke igraph cluster_fast_greedy function for this undirected network 
#red$cal_module(method = "cluster_fast_greedy")
# save network
# open the gexf file using Gephi(https://gephi.org/)
# require rgexf package
#install.packages("rgexf"); 

library(rgexf)

write.table(red$res_cor_p, file="Network/resumen_network_genuscorcut0.2.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

red$save_network(filepath = "Network/network_sparcc_genuscorcut0.2.gexf")

# Now, we show the node colors with the Phylum information and the edges colors with the positive and negative correlations. 

# calculate network attributes
#red$cal_network_attr()
# return red$res_network_attr
######### otros calculos de la red ########
red$cal_node_type()
# return red$res_node_type

red$cal_eigen()
#red$res_eigen_expla

# plot node roles with phylum information
red$plot_taxa_roles(use_type = 1)

# plot node roles with phylum information
red$plot_taxa_roles(use_type =2)

#we show the eigengene analysis of modules
#The eigengene of a module represents the main variance of the abundance in the species of the module.
red$cal_eigen()
# return red$res_eigen

# create trans_env object like the above operation
t2 <- trans_env$new(dataset = meco_dataset, add_data =  meco_dataset$sample_table[, 2:7])
# calculate correlations
t2$cal_cor(add_abund_table = t1$res_eigen)
# plot the correlation heatmap: Then we perform correlation heatmap to show the relationships between eigengenes and environmental factors.
t2$plot_cor()


##### Network ####
#t1 <- trans_network$new(dataset = dataset, cal_cor = "spearman")
#t1$res_cor_p
#t1$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE)

library(igraph)

# construct network; require igraph package
#red$cal_network(p_thres = 0.05, COR_optimization = TRUE, add_taxa_name = "Phylum")
# return 
#red$res_network
# add modules in the network
#red$cal_module()
# save network
# open the gexf file using Gephi(https://gephi.org/)
# require rgexf package
#install.packages("rgexf"); 

#library(rgexf)
#red$save_network(filepath = "network_Otusphylum.gexf")

# Now, we show the node colors with the Phylum information and the edges colors with the positive and negative correlations. 
# calculate network attributes
#red$cal_network_attr()
# return red$res_network_attr

### Muestras SLL ###
sub_sll <- red$subset_network(node = dataset$otu_table %>% .[.[, c("I_BSA1", "I_BSA2", "I_BSA3", "I_BSA4", "I2_BSA1", "I2_BSA2", "I2_BSA3", "I2_BSA4", "M_BSA1", "M_BSA2", "M_BSA3", "M_BSA4")] != 0, ] %>% rownames, rm_single = TRUE)

red_sll <- clone(red)
red_sll$res_network <- sub_sll
# then t2 have a network for 'S1' and can be used for further analysis
red_sll$cal_module()
red_sll$save_network("Network/red_sll_genuscorcut0.2.gexf")
red_sll$cal_network_attr()
red_sll$res_network_attr

### Muestras ML ###
sub_ml <- red$subset_network(node = dataset$otu_table %>% .[.[, c("I_ML1", "I_ML2", "I_ML3", "I_ML4", "I2_ML1", "I2_ML2", "I2_ML3", "I2_ML4", "M_ML1", "M_ML2", "M_ML3", "M_ML4")] != 0, ] %>% rownames, rm_single = TRUE)

red_ml <- clone(red)
red_ml$res_network <- sub_ml
# then t2 have a network for 'S1' and can be used for further analysis
red_ml$cal_module()
red_ml$save_network("Network/red_ml_genuscorcut0.2.gexf")
red_ml$cal_network_attr()
red_ml$res_network_attr

### Redes por sitios ###
#redes sll de M


#redes ml de M


#redes sll de I


#redes ml de I


#redes sll de I2


#redes ml de I2


######### Creando una nueva red pero con el m?todo Pearson

redpearson <- trans_network$new(dataset = dataset, cor_method = "pearson", filter_thres = 0) #se crea la red
# construct network; require igraph package
redpearson$cal_network(COR_p_thres = 0.01, add_taxa_name = "Genus", COR_optimization = TRUE) #The optimized COR threshold: 0.14
#Then, let's partition modules for the network.
redpearson$cal_module(method = "cluster_fast_greedy")
sub_sll_p <- redpearson$subset_network(node = dataset$otu_table %>% .[.[, c("I_BSA1", "I_BSA2", "I_BSA3", "I_BSA4", "I2_BSA1", "I2_BSA2", "I2_BSA3", "I2_BSA4", "M_BSA1", "M_BSA2", "M_BSA3", "M_BSA4")] != 0, ] %>% rownames, rm_single = TRUE)

red_sll_p <- clone(redpearson)
red_sll_p$res_network <- sub_sll_p
# then t2 have a network for 'S1' and can be used for further analysis
red_sll_p$cal_module()
red_sll_p$save_network("Network/red_sll_genuscorcut0.2_p.gexf")
red_sll_p$cal_network_attr()
red_sll_p$res_network_attr

### Muestras ML ###
sub_ml_p <- redpearson$subset_network(node = dataset$otu_table %>% .[.[, c("I_ML1", "I_ML2", "I_ML3", "I_ML4", "I2_ML1", "I2_ML2", "I2_ML3", "I2_ML4", "M_ML1", "M_ML2", "M_ML3", "M_ML4")] != 0, ] %>% rownames, rm_single = TRUE)

red_ml_p <- clone(redpearson)
red_ml$res_network <- sub_ml_p
# then t2 have a network for 'S1' and can be used for further analysis
red_ml_p$cal_module()
red_ml_p$save_network("Network/red_ml_genuscorcut0.2_p.gexf")
red_ml_p$cal_network_attr()
red_ml_p$res_network_attr

######### Creando una nueva red pero con el m?todo Pearson y nivel PHYLUM

redpearson <- trans_network$new(dataset = dataset, cor_method = "pearson", filter_thres = 0) #se crea la red
# construct network; require igraph package
redpearson$cal_network(COR_p_thres = 0.01, add_taxa_name = "Phylum", COR_cut = 0.6) #The optimized COR threshold: default 0.6
#Then, let's partition modules for the network.
redpearson$cal_module(method = "cluster_fast_greedy")
sub_sll_p <- redpearson$subset_network(node = dataset$otu_table %>% .[.[, c("I_BSA1", "I_BSA2", "I_BSA3", "I_BSA4", "I2_BSA1", "I2_BSA2", "I2_BSA3", "I2_BSA4", "M_BSA1", "M_BSA2", "M_BSA3", "M_BSA4")] != 0, ] %>% rownames, rm_single = TRUE)

red_sll_p <- clone(redpearson)
red_sll_p$res_network <- sub_sll_p
# then t2 have a network for 'S1' and can be used for further analysis
red_sll_p$cal_module()
red_sll_p$save_network("Network/red_sll_phylumcorcut0.6_p.gexf")
red_sll_p$cal_network_attr()
red_sll_p$res_network_attr

### Muestras ML ###
sub_ml_p <- redpearson$subset_network(node = dataset$otu_table %>% .[.[, c("I_ML1", "I_ML2", "I_ML3", "I_ML4", "I2_ML1", "I2_ML2", "I2_ML3", "I2_ML4", "M_ML1", "M_ML2", "M_ML3", "M_ML4")] != 0, ] %>% rownames, rm_single = TRUE)

red_ml_p <- clone(redpearson)
red_ml$res_network <- sub_ml_p
# then t2 have a network for 'S1' and can be used for further analysis
red_ml_p$cal_module()
red_ml_p$save_network("Network/red_ml_phylumcorcut0.6_p.gexf")
red_ml_p$cal_network_attr()
red_ml_p$res_network_attr

######### Creando una nueva red pero con el m?todo Spearman y nivel PHYLUM

redspearman <- trans_network$new(dataset = dataset, cor_method = "spearman", filter_thres = 0.000001) #se crea la red
# construct network; require igraph package
redspearman$cal_network(COR_p_thres = 0.01, add_taxa_name = "Phylum", COR_optimization = TRUE) #The optimized COR threshold: default 0.6
#Then, let's partition modules for the network.
redspearman$cal_module(method = "cluster_fast_greedy")
sub_sll_p <- redspearman$subset_network(node = dataset$otu_table %>% .[.[, c("I_BSA1", "I_BSA2", "I_BSA3", "I_BSA4", "I2_BSA1", "I2_BSA2", "I2_BSA3", "I2_BSA4", "M_BSA1", "M_BSA2", "M_BSA3", "M_BSA4")] != 0, ] %>% rownames, rm_single = TRUE)

red_sll_p <- clone(redspearman)
red_sll_p$res_network <- sub_sll_p
# then t2 have a network for 'S1' and can be used for further analysis
red_sll_p$cal_module()
red_sll_p$save_network("Network/red_sll_phylumcorcut0.6_spearman.gexf")
red_sll_p$cal_network_attr()
red_sll_p$res_network_attr

### Muestras ML ###
sub_ml_p <- redspearman$subset_network(node = dataset$otu_table %>% .[.[, c("I_ML1", "I_ML2", "I_ML3", "I_ML4", "I2_ML1", "I2_ML2", "I2_ML3", "I2_ML4", "M_ML1", "M_ML2", "M_ML3", "M_ML4")] != 0, ] %>% rownames, rm_single = TRUE)

red_ml_p <- clone(redspearman)
red_ml$res_network <- sub_ml_p
# then t2 have a network for 'S1' and can be used for further analysis
red_ml_p$cal_module()
red_ml_p$save_network("Network/red_ml_phylumcorcut0.6_spearman.gexf")
red_ml_p$cal_network_attr()
red_ml_p$res_network_attr

#######otra cosa
t2$res_node_type
# plot node roles with phylum information
t2$plot_taxa_roles(use_type = 1)

#we show the eigengene analysis of modules
t2$cal_eigen()
# return t1$res_eigen

# create trans_env object like the above operation
#t2 <- trans_env$new(dataset = meco_dataset, add_data =  meco_dataset$sample_table[, 4:11])
# calculate correlations
#t2$cal_cor(add_abund_table = t1$res_eigen)
# plot the correlation heatmap
#t2$plot_cor()

# extract a sub network that contains all nodes in module M1
#t1$subset_network(node = t1$res_node_type %>% .[.$module == "M1", ] %>% rownames, rm_single = TRUE)
# return a new network with igraph class
# extract sub network in which all edge labels are "+", i.e. positive edges
#t1$subset_network(edge = "+")

####### Network para ciertas muestras extract the sub-network of sample 'S1' #######

red_sll$save_network("Network/low_3ABC_26may.gexf")
red_sll$cal_network_attr()
# return red_sll$res_network_attr

### Muestras Medium 2ABC ###
sub_medABC <- red$subset_network(node = meco_dataset$otu_table %>% .[.[, c("2A", "2B", "2C")] != 0, ] %>% rownames, rm_single = TRUE)

# see https://chiliubio.github.io/microeco_tutorial/notes.html#clone for the 'clone' function explanation

red_med_sal <- clone(red)
red_med_sal$res_network <- sub_medABC
# then t2 have a network for 'S1' and can be used for further analysis
red_med_sal$cal_module()

red_med_sal$save_network("Network/med_2ABC_26may.gexf")
red_med_sal$cal_network_attr()
# return red_med_sal$res_network_attr

### Muestras High 4ABC ###
sub_highABC <- red$subset_network(node = meco_dataset$otu_table %>% .[.[, c("4A", "4B", "4C")] != 0, ] %>% rownames, rm_single = TRUE)

# see https://chiliubio.github.io/microeco_tutorial/notes.html#clone for the 'clone' function explanation

red_high_sal <- clone(red)
red_high_sal$res_network <- sub_highABC
# then t2 have a network for 'S1' and can be used for further analysis
red_high_sal$cal_module()

red_high_sal$save_network("Network/high_4ABC_26may.gexf")
red_high_sal$cal_network_attr()
# return red_high_sal$res_network_attr

## Graficar en R
library("ggraph")
# default parameter represents using igraph plot.igraph function
red_low_sal$plot_network()
plot_network(red_low_sal)
# use ggraph method; require ggraph package
# If ggraph is not installed; first install it with command: install.packages("ggraph")
red_low_sal$plot_network(method = "ggraph", node_color = "Phylum")

# please use a loop for more samples

# extract the sub-network of sample 'S1'
sub_high4A <- red$subset_network(node = meco_dataset$otu_table %>% .[.[, "4A"] != 0, ] %>% rownames, rm_single = TRUE)

# see https://chiliubio.github.io/microeco_tutorial/notes.html#clone for the 'clone' function explanation

red_high_sal <- clone(red)
red_high_sal$res_network <- sub_high4A
# then t2 have a network for 'S1' and can be used for further analysis
red_high_sal$cal_module()
red_high_sal$save_network("high_4A.gexf")
# please use a loop for more samples

########### 7.1 trans_env class #############
###### RDA analysis (db-RDA and RDA) ####

# add_data is used to add the environmental data

#t3$data_env$TN <- as.numeric(t3$data_env$TN)
#t3$data_env$C.N.ratio <- as.numeric(t3$data_env$C.N.ratio)
#t3$data_env$pH<- as.numeric(t3$data_env$pH)
#t3$data_env$EC<- as.numeric(t3$data_env$EC)
#t3$data_env$RC<- as.numeric(t3$data_env$RC)
#t3$data_env$OM<- as.numeric(t3$data_env$OM)
#t3$data_env$NO3<- as.numeric(t3$data_env$NO3)
#t3$data_env$N.NO3<- as.numeric(t3$data_env$N.NO3)
#t3$data_env$P<- as.numeric(t3$data_env$P)
#t3$data_env$Ca<- as.numeric(t3$data_env$Ca)
#t3$data_env$Mg<- as.numeric(t3$data_env$Mg)
#t3$data_env$Na<- as.numeric(t3$data_env$Na)
#t3$data_env$K<- as.numeric(t3$data_env$K)

env_table<-cbind(dataset$sample_table$EC, dataset$sample_table$RC, dataset$sample_table$OM, dataset$sample_table$Ca, dataset$sample_table$Mg, dataset$sample_table$Na, dataset$sample_table$K); colnames(env_table)<-c("EC", "RC", "OM", "Ca", "Mg", "Na", "K"); row.names(env_table)<-row.names(dataset$sample_table)
env_table<-as.data.frame(env_table)
env_table$EC<- as.numeric(env_table$EC)
env_table$RC<- as.numeric(env_table$RC)
env_table$OM<- as.numeric(env_table$OM)
env_table$Ca<- as.numeric(env_table$Ca)
env_table$Mg<- as.numeric(env_table$Mg)
env_table$Na<- as.numeric(env_table$Na)
env_table$K<- as.numeric(env_table$K)

#t3 <- trans_env$new(dataset = dataset, env_cols = 2:14) ya no voy a usar todas estas variables, sino las que escog? en la 1337

t3<- trans_env$new(dataset = dataset, add_data = env_table)
t3$data_env

# Generally, it is beneficial to analyze environmental variables in order to better use more methods. The cal_diff function is used to test the significance of variables across groups like we have shown in trans_alpha and trans_diff class parts.
t3$cal_diff(group = "Type_soil"); t3$res_diff

t3$cal_diff(method = "anova", group = "Type_soil")
# place all the plots into a list
tmp <- list()
for(i in colnames(t3$data_env)){
  tmp[[i]] <- t3$plot_diff(measure = i, add_sig_text_size = 5, xtext_size = 12) + theme(plot.margin = unit(c(0.1, 0, 0, 1), "cm"))
}
plot(gridExtra::arrangeGrob(grobs = tmp, ncol = 3))

# 800 x 600


# use bray-curtis distance for dbRDA
t3$cal_ordination(method = "dbRDA", use_measure = "bray")
# t1$res_rda is the result list stored in the object
t3$res_ordination
t3$res_ordination_R2

t3$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1.5)
t3$trans_ordination()
t3$res_ordination_trans


t3$cal_ordination_anova()
t3$res_ordination_terms
t3$res_ordination_axis

t3$cal_ordination_envfit()
t3$res_ordination_envfit

t3$cal_ordination(method = "dbRDA", use_measure = "bray", feature_sel = TRUE)

c("TN", "pH", "EC", "RC", "OM", "Ca", "Mg", "Na", "K")

# t1$res_rda_trans is the transformed result for plotting
rda<- t3$plot_ordination(plot_color = "Type_soil", plot_type = c("point", "ellipse"), env_text_size = 6)  +
  theme(axis.text.x = element_text(size = 20, color = "black"), axis.text.y = element_text(size = 20, color = "black"), axis.title = element_text(size = 25), legend.text = element_text(size = 20), legend.position = c(0.15, 0.78), legend.title = element_text(size = 0), legend.background = element_rect(fill = 0 ), plot.margin = unit(c(1.0,0.5,0.2,0.5), "cm"))
# 750x500
# contribuci?n de cada variable al ordination plot
t3$cal_ordination_envfit()
t3$res_ordination_envfit
####
t3$res_rda_envsquare

# use Genus
t3$cal_rda(use_dbrda = TRUE, taxa_level = "Order")
# As the main results of RDA are related with the projection and angles between different arrows
#The rda total result is stored in object$res_rda ...
#The R2 is stored in object$res_rda_R2 ...
# terms anova result is stored in object$res_rda_terms ...
#The axis anova result is stored in object$res_rda_axis ...
#t3$res_rda
#t3$res_rda_terms

# we adjust the length of the arrow to show them clearly using several parameters.
t3$trans_rda(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1500, max_perc_tax = 3000, min_perc_env = 200, min_perc_tax = 300)
# t1$res_rda_trans is the transformed result for plotting
t3$plot_rda(plot_color = "Salinity.level")
t3$cal_rda_envsquare()
# Result is stored in object$res_rda_envsquare ...
t3$res_rda_envsquare

#Mantel test can be used to check whether there is significant correlations between environmental variables and distance matrix.
t3$cal_mantel(use_measure = "bray")
# return t1$res_mantel
t3$res_mantel

write.table(t3$res_mantel, file = "Env/mantel_test.txt", quote = FALSE, sep = "\t")


# correlation heatmap using Genus level data.
t3$cal_cor(use_data = "Genus", p_adjust_method = "fdr")
# return t3$res_cor

# default ggplot2 method with clustering
t3$plot_cor()
# filter genera that do not have at least one ***
t3$plot_cor(filter_feature = c(""))


# correlation heatmap using Phylum level data.
t3$cal_cor(use_data = "Phylum", p_adjust_method = "fdr")
# return t3$res_cor

# default ggplot2 method with clustering
t3$plot_cor()
# filter genera that do not have at least one ***
t3$plot_cor(filter_feature = c(""))


########### METAGENOMA #############
library(microeco)
library(file2meco)
library(magrittr)

sample_file_path <- system.file("extdata", "example_metagenome_sample_info.tsv", package="file2meco")
match_file_path <- system.file("extdata", "example_metagenome_match_table.tsv", package="file2meco")
# use KEGG pathway based HUMAnN result
abund_file_path <- system.file("extdata", "example_HUMAnN_KEGG_abund.tsv", package="file2meco")
# match_table parameter can be used to replace sample names
test <- humann2meco(abund_table = abund_file_path, db = "KEGG", sample_data = sample_file_path, match_table = match_file_path)
# remove the unclassified pathway in the top level
test$tax_table %<>% subset(level1 != "unclassified")
test$tidy_dataset()
# rel = FALSE donot use relative abundance, use the raw RPK
test$cal_abund(select_cols = 1:3, rel = FALSE)
test1 <- trans_abund$new(test, taxrank = "level2", ntaxa = 10)
test1$plot_bar(facet = "Group", ylab_title = "Abundance (RPK)", use_colors = RColorBrewer::brewer.pal(12, "Set3"))


