library(dplyr)
x <- tidytree::as_tibble(tree)

p5 <- ggtree::ggtree(tree,branch.length="none",layout='circular') +
  #ggtree::geom_tiplab(color='firebrick', offset = 01.10, hjust=3.0)+
  ggtree::geom_label(ggplot2::aes(label = node),show.legend = TRUE,label.size = 0.05)

ggsave('tree_169_node_labels.pdf',p5, width=6, height=6, units="in", scale=6)

setwd('/home/lalibelulalo/PIPELINES_2023/REPEATS/Calothrix/')
tabla <- read.table("336-3.GCGATCGC.Orthologues_Palindrome_sites.txt", sep = "\t", header = TRUE )

TABLE <-as.data.frame(c())
A=1394:1397
B=2297:2312
C=2948:2962
D=2965:2969
E=3701:3710
F1=4871:4886
INTERVALS = c(A,B,C,D,E,F1)
for (i in INTERVALS){
  vector<-c(tabla[tabla$SiteNumber == i,])
  TABLE <- rbind(TABLE,vector)
}




collapse_rows_dt <- data.frame(C1 = c(rep("a", 10), rep("b", 5)),
                               C2 = c(rep("c", 7), rep("d", 3), rep("c", 2), rep("d", 3)),
                               C3 = 1:15,
                               C4 = sample(c(0,1), 15, replace = TRUE))

kbl(collapse_rows_dt) %>%
  kable_paper(full_width = F) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 1:2, valign = "middle")




tabla_t <- read.table('/home/lalibelulalo/PIPELINES_2023/CODON_MUTATIONS/BLAST_clado_calothrix.csv',sep = "\t",header = TRUE)

library(kableExtra)
kable(tabla_t ) %>% 
  collapse_rows(columns = 1, valign = "middle")



kableExtra::kbl(tabla_t, align = "c") %>%
  kableExtra::kable_paper(full_width = F) %>%
  kableExtra::column_spec(1, bold = T) %>%
  kableExtra::collapse_rows(columns = 1:2, valign = "top")



source("/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/ASR_Orth_Functions/NodeAndEdges.R")


PAL_PATH <- '/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/Thermosynechococcus_clade/PALINDROMES/CAGGCCTG/'
DIRS <- list.files(path=PAL_PATH, pattern=NULL, all.files=FALSE,
           full.names=FALSE)

for (i in 1:length(DIRS)){
  path <- paste0(PAL_PATH,DIRS[i],'/')
  setwd(path)
  getwd()
  Nodes.Edges <- Create_Transition_Table(SitesTable = "Orthologues_Palindrome_sites.txt",
                                         EvolutionModel = "F81",
                                         Method = "bayes",
                                         Phylogeny = "../../../SpeciesTree_rooted.txt",
                                         OrthoPath = "Only_ORTHOLOGUES/")
}




setwd("/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/Thermosynechococcus_clade/")
Counts <- read.table("Markov_count_GCGATCGC_gbff_homologuesThermosynechococcusspaNK55GCF000505665_f0_0taxa_algOMCL_e0__2023-7-16_23hrs7mins_Octanuc_M3_.txt", sep = "\t", header = TRUE)

Spp <- Counts$spp
Spp <- unique(Spp)

i=0
for (n in 1:length(Spp)){ 
  for (m in 1:length(Counts$spp)) {
    if (Spp[n] == Counts[m,2]){
      i=i+Counts[m,5]
    }
  }
  print(paste0(Spp[n]," ",i))
  i=0
}

FileRF1 = "codon_mutations_RF1.rooted.txt"
tableRF1 = read.table(FileRF1,header = TRUE,sep = "\t",row.names = NULL)
