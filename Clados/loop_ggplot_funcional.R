library(rlang)
library(tidyverse)
library(ggplot2)
library(ggtreeExtra)
library(ggtree)
library(treeio)
library(ggstar)
library(ggnewscale)
library(TDbook)

# CARGAMOS DATOS
#-------------------------------------
setwd('/home/lalibelulalo/R_2/FIGURES_complete_chr/')
tabla = read.table(file="refseq_no_atypical_complete_chr_269_filtrado_Octanuc.txt", header=TRUE, sep="\t")
#***
tabla <- tabla%>%
  filter(spp=='Calothrix_sp_336-3_336-3'|
           spp=='Calothrix_sp_NIES-3974_NIES-3974'|
           spp=='Calothrix_sp_PCC_6303_PCC_6303'|
           spp=='Calothrix_sp_PCC_7716_PCC_7716'|
           spp=='Calothrix_sp_NIES-4071_NIES-4071'|
           spp=='Calothrix_sp_NIES-4105_NIES-4105')
write.table(tabla,file='Count_table_calothrix_clade_filtered.txt',sep="\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
#***
tabla = read.table(file="refseq_no_atypical_complete_chr_269_filtrado_Octanuc.txt", header=TRUE, sep="\t")
tabla <- tabla%>%
  filter(spp=='Candidatus_Atelocyanobacterium_thalassa_isolate_ALOHA'|
           spp=='cyanobacterium_endosymbiont_of_Braarudosphaera_bigelowii'|
           spp=='Crocosphaera_subtropica_ATCC_51142_ATCC_51142'|
           spp=='Rippkaea_orientalis_PCC_8801_PCC_8801'|
           spp=='cyanobacterium_endosymbiont_of_Epithemia_turgida_isolate_EtSB_Lake_Yunoko'|
           spp=='cyanobacterium_endosymbiont_of_Rhopalodia_gibberula')

write.table(tabla,file='Count_table_cyanobacterium_clade_filtered.txt',sep="\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
#***
tabla = read.table(file="refseq_no_atypical_complete_chr_269_filtrado_Octanuc.txt", header=TRUE, sep="\t")
tabla <- tabla%>%
  filter(spp=='Geminocystis_sp_NIES-3708_NIES-3708'|#NIES-3708
           spp=='Geminocystis_sp_NIES-3709_NIES-3709'|#NIES-3709
           spp=='Stanieria_sp_NIES-3757_NIES-3757'|#NIES-3757
           spp=='Chondrocystis_sp_NIES-4102_NIES-4102'|#NIES-4102
           spp=='Geminocystis_herdmanii_PCC_6308_PCC_6308'|#PCC_6308
           spp=='Chamaesiphon_minutus_PCC_6605_PCC_6605'|#PCC_6605
           spp=='Halothece_sp_PCC_7418_PCC_7418'|#PCC_7418
           spp=='Stanieria_cyanosphaera_PCC_7437_PCC_7437'|#PCC_7437
           spp=='Dactylococcopsis_salina_PCC_8305_PCC_8305'|#PCC_8305
           spp=='Cyanobacterium_aponinum_PCC_10605_PCC_10605'|#PCC_10605
           spp=='Euhalothece_natronophila_Z-M001_Z-M001')#Z-M001

write.table(tabla,file='Count_table_geminocystis_clade_filtered.txt',sep="\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
#***
tabla = read.table(file="refseq_no_atypical_complete_chr_269_filtrado_Octanuc.txt", header=TRUE, sep="\t")
tabla <- tabla%>%
  filter(spp=='Pseudanabaena_sp_ABRG5-3_ABRG5-3'|#ABRG5-3
           spp=='Synechococcus_sp_JA-2-3Ba2-13_JA-2-3Ba2-13'|#JA-2-3Ba
           spp=='Synechococcus_sp_JA-3-3Ab_JA-3-3Ab'|#JA-3-3Ab
           spp=='Synechococcus_sp_PCC_7336_PCC_7336'|#PCC_7336
           spp=='Pseudanabaena_sp_PCC_7367_PCC_7367'|#PCC_7367
           spp=='Synechococcus_sp_PCC_7502_PCC_7502')#PCC_7502

write.table(tabla,file='Count_table_pseudoanabaena_clade_filtered.txt',sep="\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
#***
tabla = read.table(file="refseq_no_atypical_complete_chr_269_filtrado_Octanuc.txt", header=TRUE, sep="\t")
tabla <- tabla%>%
  filter(spp=='Synechococcus_sp_CBW1002_CBW1002'|#CBW1002
           spp=='Synechococcus_sp_CBW1004_CBW1004'|#CBW1004
           spp=='Synechococcus_sp_CBW1006_CBW1006'|#CBW1006
           spp=='Synechococcus_sp_CBW1108_CBW1108'|#CBW1108
           spp=='Cyanobium_sp_M30B3'|#M30B3
           spp=='Cyanobium_sp_NIES-981_NIES-981'|#NIES-981
           spp=='Cyanobium_sp_NS01_NS01'|#NS01
           spp=='Cyanobium_gracile_PCC_6307_PCC_6307')#PCC_6307

write.table(tabla,file='Count_table_SynechococcusCyanobium_clade_filtered.txt',sep="\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
#***
tabla = read.table(file="refseq_no_atypical_complete_chr_269_filtrado_Octanuc.txt", header=TRUE, sep="\t")
tabla <- tabla%>%
  filter(spp=='Thermosynechococcus_vestitus_BP-1_BP-1'|#BP-1
           spp=='Thermosynechococcus_sp_CL-1_CL-1'|#CL-1
           spp=='Gloeomargarita_lithophora_Alchichica-D10_D10'|#D10
           spp=='Thermosynechococcus_sp_HN-54_HN-54'|#HN-54
           spp=='Acaryochloris_marina_MBIC11017_MBIC11017'|#MBIC11017
           spp=='Acaryochloris_sp_Moss_Beach_Moss_Beach'|#Moss_Beach
           spp=='Thermosynechococcus_sp_NK55a_NK55'|#NK55a
           spp=='Synechococcus_sp_PCC_6312_PCC_6312'|#PCC_6312
           spp=='Thermostichus_lividus_PCC_6715_PCC_6715'|#PCC_6715
           spp=='Acaryochloris_marina_S15_S15'|#S15
           spp=='Thermosynechococcus_elongatus_PKUAC-SCTE542_PKUAC-SCTE542'|#SCTE542
           spp=='Thermosynechococcus_sp_TA-1_TA-1')#TA-1

write.table(tabla,file='Count_table_thermosynechococcus_clade_filtered.txt',sep="\t",row.names = FALSE,col.names = TRUE,quote = FALSE)


#-------------------------------------

setwd('/home/lalibelulalo/R_2/FIGURES_pico/')
tabla = read.table(file="Markov_count_pico_2022_gbff_2022-10-28_10hrs29mins_Octanuc_M3_.txt", header=TRUE, sep="\t")
tabla <- tabla%>%
  filter(spp=='Parasynechococcus_marenigrum_WH_8102'|
           spp=='Synechococcus_sp_WH_8103'|
           spp=='Synechococcus_sp_A18-46_1'|
           spp=='Synechococcus_sp_BOUM118'|
           spp=='Synechococcus_sp_RS9915'|
           spp=='Synechococcus_sp_A18-40'|
           spp=='Synechococcus_sp_A15-24'|
           spp=='Synechococcus_sp_A15-28')
write.table(tabla,file='Count_table_A18-40_clade_filtered.txt',sep="\t",row.names = FALSE,col.names = TRUE,quote = FALSE)

##  EDITO EL NOMBRE DE LAS ESPECIES
setwd('/home/lalibelulalo/R_2/FIGURES_complete_chr/')

tabla = read.table(file="Count_table_calothrix_clade_filtered.txt", header=TRUE, sep="\t")
tabla[tabla=='Calothrix_sp_336-3_336-3'] <- '336-3'
tabla[tabla=='Calothrix_sp_NIES-3974_NIES-3974'] <- 'NIES-3974'
tabla[tabla=='Calothrix_sp_PCC_6303_PCC_6303'] <- 'PCC_6303'
tabla[tabla=='Calothrix_sp_PCC_7716_PCC_7716'] <- 'PCC_7716'
tabla[tabla=='Calothrix_sp_NIES-4071_NIES-4071'] <- 'NIES-4071'
tabla[tabla=='Calothrix_sp_NIES-4105_NIES-4105'] <- 'NIES-4105'
tree <- read.tree("/home/lalibelulalo/PIPELINES_2023/ASR_GENOMES/Callothrix_clade/SpeciesTree_rooted.txt")
#***
tabla = read.table(file="Count_table_A18-40_clade_filtered.txt", header=TRUE, sep="\t")
tree <- read.tree("A18-40_clade_tree_rooted.txt")
#***
tabla = read.table(file="Count_table_cyanobacterium_clade_filtered.txt", header=TRUE, sep="\t")
tabla[tabla=='Candidatus_Atelocyanobacterium_thalassa_isolate_ALOHA'] <- 'ALOHA'
tabla[tabla=='cyanobacterium_endosymbiont_of_Braarudosphaera_bigelowii'] <- 'bigelowii'
tabla[tabla=='Crocosphaera_subtropica_ATCC_51142_ATCC_51142'] <- 'ATCC_51142'
tabla[tabla=='Rippkaea_orientalis_PCC_8801_PCC_8801'] <- 'PCC_8801'
tabla[tabla=='cyanobacterium_endosymbiont_of_Epithemia_turgida_isolate_EtSB_Lake_Yunoko'] <- 'Yunoko'
tabla[tabla=='cyanobacterium_endosymbiont_of_Rhopalodia_gibberula'] <- 'gibberula'
tree <- read.tree("/home/lalibelulalo/PIPELINES_2023/ASR_GENOMES/Cyanobacterium_clade/SpeciesTree_rooted.txt")
#***
tabla = read.table(file="Count_table_geminocystis_clade_filtered.txt", header=TRUE, sep="\t")
tabla[tabla=='Geminocystis_sp_NIES-3708_NIES-3708'] <- 'NIES-3708'
tabla[tabla=='Geminocystis_sp_NIES-3709_NIES-3709'] <- 'NIES-3709'
tabla[tabla=='Stanieria_sp_NIES-3757_NIES-3757'] <- 'NIES-3757'
tabla[tabla=='Chondrocystis_sp_NIES-4102_NIES-4102'] <- 'NIES-4102'
tabla[tabla=='Geminocystis_herdmanii_PCC_6308_PCC_6308'] <- 'PCC_6308'
tabla[tabla=='Chamaesiphon_minutus_PCC_6605_PCC_6605'] <- 'PCC_6605'
tabla[tabla=='Halothece_sp_PCC_7418_PCC_7418'] <- 'PCC_7418'
tabla[tabla=='Stanieria_cyanosphaera_PCC_7437_PCC_7437'] <- 'PCC_7437'
tabla[tabla=='Dactylococcopsis_salina_PCC_8305_PCC_8305'] <- 'PCC_8305'
tabla[tabla=='Cyanobacterium_aponinum_PCC_10605_PCC_10605'] <- 'PCC_10605'
tabla[tabla=='Euhalothece_natronophila_Z-M001_Z-M001'] <- 'Z-M001'
tree <- read.tree("/home/lalibelulalo/PIPELINES_2023/ASR_GENOMES/Geminocystis_clade/SpeciesTree_rooted.txt")
#***
tabla = read.table(file="Count_table_pseudoanabaena_clade_filtered.txt", header=TRUE, sep="\t")
tabla[tabla=='Pseudanabaena_sp_ABRG5-3_ABRG5-3'] <- 'ABRG5-3'
tabla[tabla=='Synechococcus_sp_JA-2-3Ba2-13_JA-2-3Ba2-13'] <- 'JA-2-3Ba'
tabla[tabla=='Synechococcus_sp_JA-3-3Ab_JA-3-3Ab'] <- 'JA-3-3Ab'
tabla[tabla=='Synechococcus_sp_PCC_7336_PCC_7336'] <- 'PCC_7336'
tabla[tabla=='Pseudanabaena_sp_PCC_7367_PCC_7367'] <- 'PCC_7367'
tabla[tabla=='Synechococcus_sp_PCC_7502_PCC_7502'] <- 'PCC_7502'
tree <- read.tree("/home/lalibelulalo/PIPELINES_2023/ASR_GENOMES/Pseudoanabaena_clade/SpeciesTree_rooted.txt")
#***
tabla = read.table(file="Count_table_SynechococcusCyanobium_clade_filtered.txt", header=TRUE, sep="\t")
tabla[tabla=='Synechococcus_sp_CBW1002_CBW1002'] <- 'CBW1002'
tabla[tabla=='Synechococcus_sp_CBW1004_CBW1004'] <- 'CBW1004'
tabla[tabla=='Synechococcus_sp_CBW1006_CBW1006'] <- 'CBW1006'
tabla[tabla=='Synechococcus_sp_CBW1108_CBW1108'] <- 'CBW1108'
tabla[tabla=='Cyanobium_sp_M30B3'] <- 'M30B3'
tabla[tabla=='Cyanobium_sp_NIES-981_NIES-981'] <- 'NIES-981'
tabla[tabla=='Cyanobium_sp_NS01_NS01'] <- 'NS01'
tabla[tabla=='Cyanobium_gracile_PCC_6307_PCC_6307'] <- 'PCC_6307'
tree <- read.tree("/home/lalibelulalo/PIPELINES_2023/ASR_GENOMES/Synechococcus_cyanobium_clade/SpeciesTree_rooted.txt")
#***
tabla = read.table(file="Count_table_thermosynechococcus_clade_filtered.txt", header=TRUE, sep="\t")
tabla[tabla=='Thermosynechococcus_vestitus_BP-1_BP-1'] <- 'BP-1'
tabla[tabla=='Thermosynechococcus_sp_CL-1_CL-1'] <- 'CL-1'
tabla[tabla=='Gloeomargarita_lithophora_Alchichica-D10_D10'] <- 'D10'
tabla[tabla=='Thermosynechococcus_sp_HN-54_HN-54'] <- 'HN-54'
tabla[tabla=='Acaryochloris_marina_MBIC11017_MBIC11017'] <- 'MBIC11017'
tabla[tabla=='Acaryochloris_sp_Moss_Beach_Moss_Beach'] <- 'Moss_Beach'
tabla[tabla=='Thermosynechococcus_sp_NK55a_NK55'] <- 'NK55a'
tabla[tabla=='Synechococcus_sp_PCC_6312_PCC_6312'] <- 'PCC_6312'
tabla[tabla=='Thermostichus_lividus_PCC_6715_PCC_6715'] <- 'PCC_6715'
tabla[tabla=='Acaryochloris_marina_S15_S15'] <- 'S15'
tabla[tabla=='Thermosynechococcus_elongatus_PKUAC-SCTE542_PKUAC-SCTE542'] <- 'SCTE542'
tabla[tabla=='Thermosynechococcus_sp_TA-1_TA-1'] <- 'TA-1'
tree <- read.tree("/home/lalibelulalo/PIPELINES_2023/ASR_GENOMES/Thermosynechococcus_clade/SpeciesTree_rooted.txt")

#------------------------------------------------
#tabla = read.table(file="Markov_count_pico_2022_gbff_2022-10-28_10hrs29mins_Octanuc_M3_.txt", header=TRUE, sep="\t")
#tree <- read.tree("pico_tree_165_rooted.txt")

tabla = read.table(file="Count_table_A18-40_clade_filtered.txt", header=TRUE, sep="\t")
tree <- read.tree("A18-40_clade_tree_rooted.txt")

#tabla = read.table(file="Count_table_calothrix_clade_filtered.txt", header=TRUE, sep="\t")
#tree <- read.tree("/home/lalibelulalo/PIPELINES_2023/ASR_GENOMES/Callothrix_clade/SpeciesTree_rooted.txt")
pwidths <- c(4.55,2.50,1.30,0.80)*0.17
angles <- c(70,70,70,70)
TextSize <- c(1.527,2.20,2.85,2.88)*3.0
HeatOffset <- c(0.320,0.230,0.190,0.180)*1.1#0.600,0.397,0.310,0.278
BarOffset <- c(0.15,0.13,0.13,0.13)*1.75#
BarPlotPwidths<- c(0.35,0.35,0.35,0.35)*1.5#0.360,0.335,0.285,0.250
PLOT_OFFSET = 1.7
title.Alt = 0.12 # 0.01
sppLab = 4*1.4 #2
HeatSize = 0.09#0.09
contador = 1
significatives=1
# tabla[14] => FrecObs
# tabla[15] => O/E
CountType = 14
set = "cyanobacterium"
kmerI = 8
GuideW = 1.0
GuideH = 1.0

##  ------------------------------------------------------------------------------------------------------------------------------
##                                  pwidths         angles          TextSize              HeatOffset              BarOffset
##  OCTANUC-216_refseq_all  4.00,2.40,1.30,0.68    0,0,0,0   0.600,1.20,1.80,2.20   0.850,0.565,0.390,0.310    0.05,0.05,0.05,0.05
##  HEXANUC-216_refseq_all  2.70,1.40,0.70,0.40    0,0,0,0   1.000,1.60,2.60,3.00   0.600,0.400,0.315,0.275    0.05,0.05,0.05,0.05
##  PENTANUC-216_refseq_all 1.80,0.90,0.30,0.10    0,0,0,0   2.200,2.20,2.20,2.20   0.450,0.335,0.258,0.235    0.05,0.05,0.05,0.05
##  TETRANUC-216_refseq_all 1.10,0.90,0.50,0.30    0,0,0,0   2.200,2.20,2.20,2.20   0.360,0.335,0.285,0.263    0.05,0.05,0.05,0.05
##  HEXANUC-269_refseq_chr                                                                                     0.370,0.370,0.370,0.370
##  TETRANUC-269_refseq_chr                                                                                    0.365,0.365,0.285,0.250
##  ------------------------------------------------------------------------------------------------------------------------------
##  OCTANUC-269_refseq_chr  4.55,2.50,1.30,0.80    0,0,0,0   0.527,1.20,1.85,1.88   0.927,0.577,0.390,0.328    0.05,0.05,0.05,0.05
##  HEXANUC-269_refseq_chr  2.70,1.40,0.70,0.40    0,0,0,0   1.000,1.60,2.60,3.00   0.600,0.400,0.315,0.275    0.05,0.05,0.05,0.05
##  PENTANUC-269_refseq_chr 1.95,0.90,0.40,0.10    0,0,0,0   2.150,2.20,2.20,2.20   0.485,0.335,0.275,0.235    0.05,0.05,0.05,0.05
##  TETRANUC-269_refseq_chr 1.10,0.90,0.50,0.30    0,0,0,0   2.200,2.20,2.20,2.20   0.360,0.335,0.285,0.260    0.05,0.05,0.05,0.05
##  HEXANUCNOHIP-269_refseq_chr                                                                                0.370,0.370,0.370,0.370
##  TETRANUCNOHIP-269_refseq_chr                                                                               0.365,0.365,0.285,0.250
##  CALOTHIRX CLADE         4.55,2.50,1.30,0.80    70,70,70,70  1.527,2.20,2.85,2.88)*3   0.320,0.230,0.190,0.180   0.15,0.13,0.13,0.13   0.35,0.35,0.35,0.35   0.2   0.09
##  c(4.55,2.50,1.30,0.80)*4.5 c(70,70,70,70) c(1.527,2.20,2.85,2.88)*2.7 c(0.320,0.230,0.190,0.180)*4.2  c(0.15,0.13,0.13,0.13)*3.2  c(0.35,0.35,0.35,0.35)*5
##  GEMINOCYSTIS|PSEUDOANABAENA  (4.55,2.50,1.30,0.80)*0.26  (0,0,0,0) (1.527,2.20,2.85,2.88)*2.0  (0.320,0.230,0.190,0.180)*0.43  (0.15,0.13,0.13,0.13)*0.6 (0.35,0.35,0.35,0.35) 0.11 
##  SYNECHOCOCCUS/CYANOBIUM  (4.55,2.50,1.30,0.80)*0.26  (0,0,0,0)  (1.527,2.20,2.85,2.88)*4.0 (0.320,0.230,0.190,0.180)*0.8 (0.15,0.13,0.13,0.13)*2.0 (0.35,0.35,0.35,0.35) 0.11 
##  THERMOSYNECHOCOCCUS (4.55,2.50,1.30,0.80)*0.05  (0,0,0,0) c(1.527,2.20,2.85,2.88)*4.0 (0.320,0.230,0.190,0.180)*0.4 (0.15,0.13,0.13,0.13)*0.8 (0.35,0.35,0.35,0.35)
##  ------------------------------------------------------------------------------------------------------------------------------
##  OCTANUC-297_elhai       4.00,2.00,0.90,0.68    0,0,0,0   0.600,1.45,1.80,2.20   0.433,0.270,0.190,0.175    0.05,0.05,0.05,0.05
##  HEXANUC-297_elhai       2.70,1.20,0.70,0.60    0,0,0,0   1.000,1.60,2.60,3.00   0.313,0.205,0.175,0.169    0.05,0.05,0.05,0.05
##  PENTANUC-297_elhai      1.843,0.510,0.3,0.1    0,0,0,0   2.200,2.45,3.00,2.20   0.250,0.165,0.148,0.137    0.05,0.05,0.05,0.05
##  TETRANUC-297_elhai      0.80,0.60,0.50,0.40    0,0,0,0   2.200,2.20,2.20,2.20   0.185,0.175,0.161,0.155    0.05,0.05,0.05,0.05
##  HEXANUCNOHIP-297_elhai                                                                                     0.195,0.195,0.195,0.195
##  TETRANUCNOHIP-297_elhai                                                                                    0.195,0.195,0.195,0.195
##  ------------------------------------------------------------------------------------------------------------------------------
##  OCTANUC-165_pico        0.69,0.30,0.10,0.10    0,0,0,0   2.200,2.20,2.20,2.20   0.222,0.188,0.171,0.171    0.05,0.05,0.05,0.05
##  HEXANUC-165_pico
##  PENTANUC-165_pico       0.69,0.30,0.10,0.10    0,0,0,0   2.200,2.20,2.20,2.20   0.222,0.188,0.171,0.171    0.05,0.05,0.05,0.05
##  TETRANUC-165_pico       0.98,0.50,0.10,0.10    0,0,0,0   2.200,2.20,3.20,3.20   0.249,0.210,0.586,0.571    0.05,0.05,0.255,0.255
##  HEXANUCNOHIP-165_pico                                                                                      0.260,0.260,0.260,0.260
##  TETRANUCNOHIP-165_pico                                                                                     0.260,0.260,0.260,0.260
##
###

########
{
if (CountType==14){
  CountType2 = "FrecObs"
  CountType2H = "HigherFrecObs"
}else if(CountType==15){
  CountType2 = "OE"
  CountType2H = "HigherOE"
}

if(kmerI==4){
  kmerS = "Tetranuc"
  PAL = "GATC"
}else if(kmerI==5){
  kmerS = "Pentanuc"
  PAL = "CAATG"#GCGGC
}else if(kmerI==6){
  kmerS = "Hexanuc"
  PAL = "CGATCG"
}else if(kmerI==8){
  kmerS = "Octanuc"
  PAL = "GCGATCGC"
}else if(kmerI==41){
  kmerS = "TetranucNoHIP"
  kmerI = 4
  PAL = "GATC"
}else if(kmerI==61){
  kmerS = "HexanucNoHIP"
  kmerI = 6
  PAL = "CGATCG"
}else if(kmerI==7){
  kmerS = "Heptanuc"
  kmerI = 7
  PAL = "AGGGCCT"
}else if(kmerI==73){
  kmerS = "Heptanuc1ntNoPal"
  kmerI = 7
  PAL = "GCGATCG"
}else if(kmerI==82){
  kmerS = "OctanucHIP2ntPal"
  kmerI = 8
  PAL = "GCGATCGC"
}else if(kmerI==83){
  kmerS = "OctanucHIP1ntNoPal"
  kmerI = 8
  PAL = "GCGATCGC"
}
###
}
# CALCULAMOS PVALUES
if (kmerI==4 | kmerI==41){
  tabla$pval = pbinom((tabla$obs-1),
                      (tabla$genomesize-kmerI+1),
                      (tabla$markov2/(tabla$genomesize-kmerI+1)),lower.tail = FALSE) # log.p = "TRUE" ## EDITAR markov3 o markov2
  # CALCULAMOS FDRS, FRECUENCIAS OBSERVADAS, OE y CONTENIDO GC
  tabla$fdrs <- p.adjust(tabla$pval, method="fdr")
  tabla$frecObs <- (1000*(tabla$obs)/(tabla$genomesize))
  tabla$OE <- ((tabla$obs)/(tabla$markov2)) ## EDITAR markov3 o markov2
  tabla$GC <-(((tabla$G)+(tabla$C))/((tabla$A)+(tabla$Th)+(tabla$G)+(tabla$C)))
}else{
  tabla$pval = pbinom((tabla$obs-1),
                      (tabla$genomesize-kmerI+1),
                      (tabla$markov3/(tabla$genomesize-kmerI+1)),lower.tail = FALSE) # log.p = "TRUE" ## EDITAR markov3 o markov2
  # CALCULAMOS FDRS, FRECUENCIAS OBSERVADAS, OE y CONTENIDO GC
  tabla$fdrs <- p.adjust(tabla$pval, method="fdr")
  tabla$frecObs <- (1000*(tabla$obs)/(tabla$genomesize))
  tabla$OE <- ((tabla$obs)/(tabla$markov3)) ## EDITAR markov3 o markov2
  tabla$GC <-(((tabla$G)+(tabla$C))/((tabla$A)+(tabla$Th)+(tabla$G)+(tabla$C)))
}


## IDENTIFICAMOS LOS CONJUNTOS DE PALINDROMOS SIGNIFICATIVOS
## Graficamos
## Y CREAMOS UNA TABLA PARA CADA PVALUE 
{prefijo = paste0("Markov_",set,"_",kmerS,"_",CountType2,"_")#"Markov_elhai_97_Tetranuc_"
pvalues <- c(1e-32,1e-64,1e-128,1e-256)
dataframes <- c("sel32","sel64","sel128","sel256")
for (pvalue in pvalues) {
  significatives.pvalue <- c()
  for(i in 1:nrow(tabla)) {
    if(tabla[i,12] < pvalue){
      significatives.pvalue <- rbind(significatives.pvalue,c(tabla[i,]))
      write.table(significatives.pvalue, file=paste(prefijo,pvalue,".txt",sep=""), sep="\t", row.names = F)
    }
  }
}

# CARGAMOS TABLAS PARA CADA PVALUE DE CORTE
if (significatives==2){
  sel32=read.table(file=paste0(prefijo,"1e-32.txt"), header=TRUE, sep="\t")
  sel64=read.table(file=paste0(prefijo,"1e-64.txt"), header=TRUE, sep="\t")
  df2 <- list(sel32,sel64)
  dataframes <- c("sel32","sel64")
}else if(significatives==3){
  sel32=read.table(file=paste0(prefijo,"1e-32.txt"), header=TRUE, sep="\t")
  sel64=read.table(file=paste0(prefijo,"1e-64.txt"), header=TRUE, sep="\t")
  sel128=read.table(file=paste0(prefijo,"1e-128.txt"), header=TRUE, sep="\t")
  df2 <- list(sel32,sel64,sel128)
  dataframes <- c("sel32","sel64","sel128")
}else if(significatives==4){
  sel32=read.table(file=paste0(prefijo,"1e-32.txt"), header=TRUE, sep="\t")
  sel64=read.table(file=paste0(prefijo,"1e-64.txt"), header=TRUE, sep="\t")
  sel128=read.table(file=paste0(prefijo,"1e-128.txt"), header=TRUE, sep="\t")
  sel256=read.table(file=paste0(prefijo,"1e-256.txt"), header=TRUE, sep="\t")
  df2 <- list(sel32,sel64,sel128,sel256)
  dataframes <- c("sel32","sel64","sel128","sel256")  
}else{
  sel32=read.table(file=paste0(prefijo,"1e-32.txt"), header=TRUE, sep="\t")
  df2 <- list(sel32)
  dataframes <- c("sel32")  
  }


# GRAFICAMOS CADA TABLA
#library(tidyverse)

#Create dataframes(In this example n = 3)
SigNames <- function(x) {
  d <- substitute(x)
  n <- sapply(d[-1],deparse)
  return(n)
}

# GUARDAMOS DATAFRAMES EN UNA LISTA
if (significatives==2){
  df_1 <- sel32  
  df_2 <- sel64
  all_data <- list(sel32,sel64)
  df_names <- SigNames(c(sel32,sel64))
}else if(significatives==3){
  df_1 <- sel32  
  df_2 <- sel64
  df_3 <- sel128
  all_data <- list(sel32,sel64,sel128)
  df_names <- SigNames(c(sel32,sel64,sel128))
}else if(significatives==4){
  df_1 <- sel32  
  df_2 <- sel64
  df_3 <- sel128
  df_4 <- sel256
  all_data <- list(sel32,sel64,sel128,sel256)
  df_names <- SigNames(c(sel32,sel64,sel128,sel256))
}else{
  df_1 <- sel32
  all_data <- list(sel32)
  df_names <- SigNames(c(sel32))
}

names(all_data)<-lapply(df_names, function(x) paste0("df_",set,"_",kmerS,"_",CountType2,"_",x))

#Graph and save for each dataframe

if(kmerI==4){
  for (i in 1:length(all_data)){
    df_i <- all_data[[i]]
    benp <-  
      #df_i %>%
      ggplot(data = tabla,
             aes(x = log10(frecObs),
                 y = log10(obs/markov2))) + ## EDITAR markov3 o markov2
      geom_point(color = "gray", size = 0.5, alpha = 1/5) +
      geom_point(data=df_i, aes(x = log10(frecObs), y = log10(obs/markov2), col = palindrome)) + ## EDITAR markov3 o markov2
      labs(title=paste0("Palindromos ", names(all_data)[i]), subtitle="") 
    
    ggsave(benp, file=paste0(names(all_data)[i],"_significative-palindromes.png"),width=7, height=6, units="in", scale=1.5)
  }
}else{
    for (i in 1:length(all_data)){
      df_i <- all_data[[i]]
      benp <-  
        #df_i %>%
        ggplot(data = tabla,
               aes(x = log10(frecObs),
                   y = log10(obs/markov3))) + ## EDITAR markov3 o markov2
        geom_point(color = "gray", size = 0.5, alpha = 1/5) +
        geom_point(data=df_i, aes(x = log10(frecObs), y = log10(obs/markov3), col = palindrome)) + ## EDITAR markov3 o markov2
        labs(title=paste0("Palindromos ", names(all_data)[i]), subtitle="") 
      
      ggsave(benp, file=paste0(names(all_data)[i],"_significative-palindromes.png"),width=7, height=6, units="in", scale=1.5)
    }
  }
}

###
###
## Obtenemos los Archivos
###
###
{## Funcion para obtener nombres
SigNames <- function(x) {
  d <- substitute(x)
  n <- sapply(d[-1],deparse)
  return(n)
}

if (significatives==2){
  df_sel <- list(sel32,sel64)
  df_names <- SigNames(c(sel32,sel64))
}else if(significatives==3){
  df_sel <- list(sel32,sel64,sel128)
  df_names <- SigNames(c(sel32,sel64,sel128))
}else if (significatives==4){
  df_sel <- list(sel32,sel64,sel128,sel256)
  df_names <- SigNames(c(sel32,sel64,sel128,sel256))  
}else{
  df_sel <- list(sel32)
  df_names <- SigNames(c(sel32))  
  }

names(df_sel)<-lapply(df_names, function(x) paste0("df_",x))

for (l in 1:length(df_sel)){
  df_l <- df_sel[[l]]
  ### *** ###
  
  # HACEMOS UNA LISTA DE ESPECIES 
  spp_sel <- c()
  for(i in 1:nrow(df_l)) {
    spp_sel <- append(spp_sel,df_l[i,1])
  }
  # QUITAMOS LAS REPETICIONES
  spp_sel <-unique(spp_sel)
  spp_sel
  
  # CREAMOS UNA LISTA VACIA QUE CONTENDRÁ AL PALINDROMO MAS ABUNDANTE POR ESPECIE
  highest_counts <- c()
  
  # BUSCAMOS EL PALINDROMO MAS ABUNDANTE:
  
  ##
  ##  tabla[14] => FrecObs
  ##  tabla[15] => O/E
  
  for(j in 1:length(spp_sel)) {         # PARA CADA ESPECIE:
    spp_pals_counts <- c()              # CREAMOS UNA LISTA DE LOS PALINDROMOS EXISTENTES
    print(spp_sel[j])
    for(i in 1:length(df_l[,1])) {     # PARA CADA RENGLON DE LA TABLA DE SIGNIFICATIVOS
      if (isTRUE(spp_sel[j] == df_l[i,1])==TRUE){ # SI LA ESPECIE COINCIDE CON LA QUE ESTAMOS CHECANDO
        spp_pals_counts <- append(spp_pals_counts, c(df_l[i,CountType])) # ENTONCES GUARDAMOS EL OE DEL PALINDROMO
      }                                                            # EN LA LISTA DE PALINDROMOS EXISTENTES
    }
    for (k in 1:length(df_l[,1])) {    # PARA CADA REGLON DE LA TABLA DE SIGNIFICATIVOS
      if (isTRUE((df_l[k,CountType])==spp_pals_counts[which.max(abs(spp_pals_counts))])==TRUE){ # SI EL OE EN CUESTION ES EL OE MAS GRANDE
        highest_counts <- rbind(highest_counts, c(df_l[k,])) # ENTONCES GUARDAMOS LA LINEA EN  HIGHEST COUNTS
        highest_counts <- unique(highest_counts)  # QUITAMOS LOS REPETIDOS
      }
    }
  }
  
  # CREAMOS LA TABLA CON LOS PALINDROMOS MAS SOBREREPRESENTADOS
  write.table(highest_counts, file=paste0("Highest_counts_",set,"_",kmerS,"_",CountType2,"_",names(df_sel)[l],".txt"), sep="\t", row.names = F)
  
  # ABRIMOS LA TABLA PARA CREAR EL ARCHIVO
  tabla_l <- read.csv(file=paste0("Highest_counts_",set,"_",kmerS,"_",CountType2,"_",names(df_sel)[l],".txt"), header=TRUE, sep="\t")
  tabla_l <- unique(tabla_l)
  
  # TABLA DE ESPECIES COMPLETAS CON 0's
  # En elcaso de las especies sin palindromos les agrego HIP con valor de 0
  # este no se cuenta en el grafico
  df_barplot_attr_zeros <- c()
  for (i in 1:length(tabla[,1])){
    vec_zeros <- c(tabla[i,1], PAL, 0) ## EDITAR ESTE CERO DEPENDIENDO EL ANALISIS  ## GCGGC
    df_barplot_attr_zeros <- rbind(df_barplot_attr_zeros, vec_zeros)
  }
  df_barplot_attr_zeros <- unique(df_barplot_attr_zeros)
  
  pico_df_barplot_attr <- c()
  # CREAMOS df_barplot_attr.txt
  for (i in 1:length(tabla_l[,1])){
    vec <- c(tabla_l[i,1], tabla_l[i,3], tabla_l[i,CountType])
    pico_df_barplot_attr <- rbind(pico_df_barplot_attr, vec)
  }
  # Combino las tablas con especies que tienen palindromos y las que no
  for (i in 1:length(df_barplot_attr_zeros[,1])){
    for (j in 1:length(pico_df_barplot_attr[,1])){
      if ((isTRUE(df_barplot_attr_zeros[i,1]==pico_df_barplot_attr[j,1])==TRUE)){
        df_barplot_attr_zeros[i,2] = pico_df_barplot_attr[j,2]
        df_barplot_attr_zeros[i,3] = pico_df_barplot_attr[j,3]
      }
    }
  }
  colnames(df_barplot_attr_zeros) = c("Spp", "Palindromes", CountType2H)
  write.table(df_barplot_attr_zeros, file=paste0(set,"_",kmerS,"_",CountType2,"_",names(df_sel)[l],"_barplot_attr.txt"), sep="\t", row.names = F)
  
  # CREAMOS df_ring_heatmap
  # Primero hago una lista con todos los palindromos significativos
  pals_sel <- c()
  for (i in 1:nrow(df_l)) {
    pals_sel <- append(pals_sel,df_l[i,3])
  }
  pals_sel <- unique(pals_sel)
  # Para cada especie agrego el conteo para cada uno de los palindromos significativos (pals_sel)
  df_ring_heatmap_counts <- c()
  for (i in 1:length(df_l[,1])){
    vec_counts <- c(df_l[i,1], df_l[i,3], df_l[i,CountType])
    df_ring_heatmap_counts <- rbind(df_ring_heatmap_counts, vec_counts)
    df_ring_heatmap_counts <- unique(df_ring_heatmap_counts)
  }
  # Para cada especie agrego un CERO al conteo para cada uno de los palindromos significativos (pals_sel)
  df_ring_heatmap_zeros <- c()
  for (i in 1:length(tabla[,1])){
    for (j in 1:length(pals_sel)){
      vec_zeros <- c(tabla[i,1], pals_sel[j], 0)
      df_ring_heatmap_zeros <- rbind(df_ring_heatmap_zeros, vec_zeros)
    }
    df_ring_heatmap_zeros <- unique(df_ring_heatmap_zeros)
  }
  
  # Combino la tabla de CEROS y la de conteo de cada palindromo significativo
  for (i in 1:length(df_ring_heatmap_zeros[,1])){
    for (j in 1:length(df_ring_heatmap_counts[,1])){
      if ((isTRUE(df_ring_heatmap_zeros[i,1]==df_ring_heatmap_counts[j,1])==TRUE)&(isTRUE(df_ring_heatmap_zeros[i,2]==df_ring_heatmap_counts[j,2])==TRUE)){
        df_ring_heatmap_zeros[i,3] = df_ring_heatmap_counts[j,3]
      }
    }
  }
  colnames(df_ring_heatmap_zeros) = c("Spp", "Palindromes", CountType2)
  write.table(df_ring_heatmap_zeros, file=paste0(set,"_",kmerS,"_",CountType2,"_",names(df_sel)[l],"_ring_heatmap.txt"), sep="\t", row.names = F)
  
  ### *** ###
}
}


###
## TREE ANOTATION
###

#setwd('/home/lalibelulalo/R_2/FIGURES_pico/Octanucs_OE/')
#l ="sel256"
if (CountType == 14){
  for (l in dataframes){
    dat2_l<-read.table(file=paste0(set,"_",kmerS,"_",CountType2,"_df_",l,"_ring_heatmap.txt"), header=TRUE, sep="\t")
    dat3_l<-read.table(file=paste0(set,"_",kmerS,"_",CountType2,"_df_",l,"_barplot_attr.txt"), header=TRUE, sep="\t")
  
    ###
    ## palindromos
    pals_l <- c()
    for (i in 1:nrow(dat2_l)) {
      pals_l <- append(pals_l,dat2_l[i,2])
    }
    pals_l
    pals_l <- unique(pals_l)
    pals_l
    ## colores para los palindromos
    cols_l <- rainbow(length(pals_l))
    cols_l
    ##
    if (length(pals_l)==1){
      #
      p <- ggtree(tree, size=0.2, open.angle=5)+ # layout="fan" #, open.angle=5**
        geom_tiplab(size=sppLab, align=TRUE, offset=0.01) +
        ggplot2::expand_limits(x = 1)
      
      p2 <- p
      
      p1 <-p2
      
      p1 <- p1 + new_scale_fill() +
        geom_fruit(data=dat3_l, geom=geom_bar,
                   mapping=aes(y=Spp, x=HigherFrecObs, fill=Palindromes),
                   pwidth=BarPlotPwidths[contador], 
                   orientation="y", 
                   stat="identity",
                   axis.params=list(
                     axis       = "x",
                     text.size  = 3.8*2,
                     title = CountType2, title.size = 3*2, title.height = title.Alt,
                     hjust      = 1,
                     vjust      = 0.5,
                     nbreak     = 3,
                   ),
                   grid.params=list(),offset = BarOffset[contador]) +
        scale_fill_manual(values=cols_l,
                          guide=guide_legend(keywidth = GuideW,#0.3, 
                                             keyheight = GuideH,#0.3,
                                             order=4))+
        geom_treescale(fontsize=2, linesize=0.3, x=0.001, y=0.001) + # X es que tan largas son las ramas
        theme(legend.position=c(0.93, 0.5),
              legend.background=element_rect(fill=NA),
              legend.title=element_text(size=16.5), #size=6.5
              legend.text=element_text(size=14.5), #size=4.5
              legend.spacing.y = unit(0.02, "cm"),
        )
      p1
      p4=p1
      p5=p4 + theme_tree()
      ggsave(paste0(set,"_",kmerS,"_",CountType2,"_",l,"_filogenia.pdf"),p5, width=6, height=6, units="in", scale=3)
      ggsave(paste0(set,"_",kmerS,"_",CountType2,"_",l,"_filogenia_HIG.png"),p1, width=6, height=6, units="in", scale=3)
      contador = contador+1
    }else{
      p <- ggtree(tree, size=0.4, open.angle=5) + # layout="fan" #, open.angle=5**
      geom_tiplab(size=sppLab, align=TRUE, offset=0.01) +
        ggplot2::expand_limits(x = PLOT_OFFSET)
      p2 <- p
      p1 <-p2
      p1 <- p1 + new_scale_fill() +
        
      geom_fruit(data=dat2_l,
                geom=geom_tile,
                mapping=aes(y=Spp,
                            x=Palindromes, 
                            alpha=FrecObs,
                            fill=Palindromes),
                color = "grey50",
                offset = HeatOffset[contador],## Offset QUE TAN LEJOS ESTA EL HEATMAP 0.4 0.02
                size = HeatSize,
                axis.params = list(axis=unique(dat2_l[i,2]), ## etiquetas pal_sel
                text.angle = angles[contador],
                text.size = TextSize[contador],
                hjust=0.5), pwidth = pwidths[contador])+ # pwidth : ancho de celda
      scale_alpha_continuous(range=c(0, 1),                                                   # hjust mueve las etiquetas del heatmap
      guide=guide_legend(keywidth = 0.3, 
      keyheight = 0.3, order=5)) +
      
      geom_fruit(data=dat3_l,
                 geom=geom_bar,
                 mapping=aes(y=Spp,
                             x=HigherFrecObs,
                             fill=Palindromes),
                 pwidth=BarPlotPwidths[contador],
                 orientation="y",
                 stat="identity",
                 axis.params=list(axis = "x",
                                  text.size  = 1.8,
                                  title = CountType2,
                                  title.size = 3*2,
                                  title.height = title.Alt,
                                  hjust = 1,
                                  vjust = 0.5,
                                  nbreak = 3),
                 grid.params=list(),offset = BarOffset[contador]) +
        
      scale_fill_manual(values=cols_l,
                        guide=guide_legend(keywidth = 0.3,
                                           keyheight = 0.3,
                                           order=4))+
      geom_treescale(fontsize=2,
                     linesize=0.3,
                     x=0.001,
                     y=0.001) + # X es que tan largas son las ramas
        
      theme(legend.position=c(0.93, 0.5),
            legend.background=element_rect(fill=NA),
            legend.title=element_text(size=6.5*3.8),
            legend.text=element_text(size=4.5*3.8),
            legend.spacing.y = unit(0.02, "cm")
            )
      
          p1
          p4=p1
          p6 = p4 + theme_tree() 
          p5=p4 + theme_tree() #+
            #geom_hilight(node=c(416), fill="red", alpha=0.2, extend = 4.8)+
            #geom_hilight(node=c(522), fill="orange", alpha=0.2, extend = 4.8)+
            #geom_hilight(node=c(216,488,491,528), fill="yellow", alpha=0.2, extend = 4.8)+
            #geom_hilight(node=c(275,278), fill="green", alpha=0.2, extend = 4.8)+
            #geom_hilight(node=c(281), fill="blue", alpha=0.2, extend = 4.8)
            #geom_hilight(node=c(219), fill="blue", alpha=0.2, extend = 4.8) ## A18-40
          ggsave(paste0(set,"_",kmerS,"_",CountType2,"_",l,"_filogenia_HIG_cases.pdf"),p5, width=6, height=6, units="in", scale=3)
          ggsave(paste0(set,"_",kmerS,"_",CountType2,"_",l,"_filogenia_HIG_cases.png"),p5, width=6, height=6, units="in", scale=3)
          ggsave(paste0(set,"_",kmerS,"_",CountType2,"_",l,"_filogenia_HIG.pdf"),p1, width=6, height=6, units="in", scale=3)
          ggsave(paste0(set,"_",kmerS,"_",CountType2,"_",l,"_filogenia_HIG.png"),p1, width=6, height=6, units="in", scale=3)
          contador = contador+1
          }
  }
}
#-
if (CountType == 15){
  for (l in dataframes){
    dat2_l<-read.table(file=paste0(set,"_",kmerS,"_",CountType2,"_df_",l,"_ring_heatmap.txt"), header=TRUE, sep="\t")
    dat3_l<-read.table(file=paste0(set,"_",kmerS,"_",CountType2,"_df_",l,"_barplot_attr.txt"), header=TRUE, sep="\t")
    
    ###
    ## palindromos
    pals_l <- c()
    for (i in 1:nrow(dat2_l)) {
      pals_l <- append(pals_l,dat2_l[i,2])
    }
    pals_l
    pals_l <- unique(pals_l)
    pals_l
    ## colores para los palindromos
    cols_l <- rainbow(length(pals_l))
    cols_l
    ##
    if (length(pals_l)==1){
      p <- ggtree(tree, size=0.4, open.angle=5)+ # layout="fan" #, open.angle=5**
        geom_tiplab(size=sppLab, align=TRUE, offset=0.01)# +
      p2 <- p
      p1 <-p2
      p1 <- p1 + new_scale_fill() +
        geom_fruit(data=dat3_l,
                   geom=geom_bar,
                   mapping=aes(y=Spp,
                               x=HigherOE,
                               fill=Palindromes),
                   pwidth=BarPlotPwidths[contador], 
                   orientation="y", 
                   stat="identity",
                   axis.params=list(axis = "x",
                                    text.size  = 1.8,
                                    title = CountType2,
                                    title.size = 3,
                                    title.height = title.Alt,
                                    hjust = 1,
                                    vjust = 0.5,
                                    nbreak = 3),
                   grid.params=list(),
                   offset = BarOffset[contador]) +
        scale_fill_manual(values=cols_l,
                          guide=guide_legend(keywidth = 0.3, 
                                             keyheight = 0.3,
                                             order=4))+
        geom_treescale(fontsize=2,
                       linesize=0.3,
                       x=0.001,
                       y=0.001) + # X es que tan largas son las ramas
        theme(legend.position=c(0.93, 0.5),
              legend.background=element_rect(fill=NA),
              legend.title=element_text(size=6.5),
              legend.text=element_text(size=4.5),
              legend.spacing.y = unit(0.02, "cm"),
        )
      p1
      p4=p1
      p5=p4 + theme_tree()
      ggsave(paste0(set,"_",kmerS,"_",CountType2,"_",l,"_filogenia.pdf"),p5, width=6, height=6, units="in", scale=3)
      contador = contador+1
    }else{
      p <- ggtree(tree, size=0.4, open.angle=5,)+ # layout="fan" #, open.angle=5**
      geom_tiplab(size=sppLab, align=TRUE, offset=0.01)+
        ggplot2::expand_limits(x = 1)
      
      p2 <- p
      p1 <-p2
      p1 <- p1 + new_scale_fill() +
        geom_fruit(data=dat2_l,
                   geom=geom_tile,
                   mapping=aes(y=Spp,
                               x=Palindromes, 
                               alpha=OE,
                               fill=Palindromes),
                   color = "grey50",
                   offset = HeatOffset[contador],## Offset QUE TAN LEJOS ESTA EL HEATMAP 0.4 0.02
                   size = HeatSize,
                   axis.params = list(axis=unique(dat2_l[i,2]), ## etiquetas pal_sel
                                      text.angle = angles[contador],
                                      text.size = TextSize[contador],
                                      hjust=0.5),
                   pwidth = pwidths[contador])+ # pwidth : ancho de celda
        scale_alpha_continuous(range=c(0, 1),                                                   # hjust mueve las etiquetas del heatmap
                               guide=guide_legend(keywidth = 0.3, 
                                                  keyheight = 0.3,
                                                  order=5)) +
        geom_fruit(data=dat3_l,
                   geom=geom_bar,
                   mapping=aes(y=Spp,
                               x=HigherOE,
                               fill=Palindromes),
                   pwidth=BarPlotPwidths[contador], 
                   orientation="y", 
                   stat="identity",
                   axis.params=list(axis = "x",
                                    text.size  = 1.8,
                                    title = CountType2,
                                    title.size = 3*2,
                                    title.height = title.Alt,
                                    hjust = 1,
                                    vjust = 0.5,
                                    nbreak = 3),
                   grid.params=list(),
                   offset = BarOffset[contador]) +
        scale_fill_manual(values=cols_l,
                          guide=guide_legend(keywidth = 0.3, 
                                             keyheight = 0.3,
                                             order=4))+
        geom_treescale(fontsize=2,
                       linesize=0.3,
                       x=0.001,
                       y=0.001) + # X es que tan largas son las ramas
        theme(legend.position=c(0.93, 0.5),
              legend.background=element_rect(fill=NA),
              legend.title=element_text(size=6.5*3.8),
              legend.text=element_text(size=4.5*3.8),
              legend.spacing.y = unit(0.02, "cm"),
        )
      p1
      p4=p1
      p6 = p4 + theme_tree()
      p5=p4 + theme_tree()#+
        #geom_hilight(node=c(416), fill="red", alpha=0.2, extend = 4.8)+
        #geom_hilight(node=c(522), fill="orange", alpha=0.2, extend = 4.8)+
        #geom_hilight(node=c(216,488,491,528), fill="yellow", alpha=0.2, extend = 4.8)+
        #geom_hilight(node=c(275,278), fill="green", alpha=0.2, extend = 4.8)+
        #geom_hilight(node=c(281), fill="blue", alpha=0.2, extend = 4.8)
      #geom_hilight(node=c(219), fill="blue", alpha=0.2, extend = 4.8)
      ggsave(paste0(set,"_",kmerS,"_",CountType2,"_",l,"_filogenia_HIG_cases.pdf"),p5, width=6, height=6, units="in", scale=3)
      contador = contador+1
      ggsave(paste0(set,"_",kmerS,"_",CountType2,"_",l,"_filogenia_HIG_cases.png"),p5, width=6, height=6, units="in", scale=3)
      ggsave(paste0(set,"_",kmerS,"_",CountType2,"_",l,"_filogenia_HIG.png"),p6, width=6, height=6, units="in", scale=3)
      #ggsave(benp, file=paste0(names(all_data)[i],"_significative-palindromes.png"),width=7, height=6, units="in", scale=1.5)
      ggsave(paste0(set,"_",kmerS,"_",CountType2,"_",l,"_filogenia_HIG.pdf"),p1, width=6, height=6, units="in", scale=3)
      ggsave(paste0(set,"_",kmerS,"_",CountType2,"_",l,"_filogenia_HIG.png"),p1, width=6, height=6, units="in", scale=3)
    }
  }
}
###



#--------------------------------------------------------------------------------------------------
unique(sel32[,3])
write.table(unique(sel32[,3]),file="palindromes_from_genome_pal_counts_32.txt",sep="\t",row.names = FALSE,col.names = FALSE,quote = FALSE)


setwd('/home/lalibelulalo/R_2/FIGURES_pico/')
tabla = read.table(file="Markov_count_pico_2022_gbff_2022-10-28_10hrs29mins_Octanuc_M3_.txt", header=TRUE, sep="\t")
tabla <- tabla%>%
  filter(spp=='Parasynechococcus_marenigrum_WH_8102'|spp=='Synechococcus_sp_WH_8103'|spp=='Synechococcus_sp_A18-46_1'|spp=='Synechococcus_sp_BOUM118'|spp=='Synechococcus_sp_RS9915'|spp=='Synechococcus_sp_A18-40'|spp=='Synechococcus_sp_A15-24'|spp=='Synechococcus_sp_A15-28')
tabla <- tabla%>%
  filter(palindrome=='CGTTAACG')
write.table(tabla,file='Count_table_calothrix_clade_filtered.txt',sep="\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
