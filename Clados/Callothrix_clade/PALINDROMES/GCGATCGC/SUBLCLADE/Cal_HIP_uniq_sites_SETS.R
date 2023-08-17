source("/home/lalibelulalo/TESIS/ASR_Orth_Functions/NodeAndEdges.R")
setwd("/home/lalibelulalo/TESIS/Clados/Geminocystis_clade/")


clades <- c(#"Geminocystis_clade",
            #"Callothrix_clade",
            "Pseudoanabaena_clade",
            #"Clado_A18-40",
            #"SynechococcusCyanobium_clade",
            #"Cyanobacterium_clade",
            "Thermosynechococcus_clade")
for (j in clades){
  #print(j)
  spps <- system(paste0("ls /home/lalibelulalo/TESIS/Clados/",j,"/PALINDROMES/GCGATCGC/"), intern = TRUE)
  for (i in 1:length(spps)){
    setwd(paste0("/home/lalibelulalo/TESIS/Clados/",j,"/PALINDROMES/GCGATCGC/",spps[i],"/"))
    Create_Transition_Table_No_Fit(SitesTable = "Orthologues_Palindrome_sites.txt",
                                   EvolutionModel = "F81",
                                   Method = "bayes",
                                   Phylogeny = "../../../SpeciesTree_rooted.txt",
                                   OrthoPath = "Only_ORTHOLOGUES/")
  }
}

j = "Thermosynechococcus_clade"
spps <- system(paste0("ls /home/lalibelulalo/TESIS/Clados/",j,"/PALINDROMES/CAGGCCTG/"), intern = TRUE)
for (i in 1:length(spps)){
  setwd(paste0("/home/lalibelulalo/TESIS/Clados/",j,"/PALINDROMES/CAGGCCTG/",spps[i],"/"))
  Create_Transition_Table_No_Fit(SitesTable = "Orthologues_Palindrome_sites.txt",
                                 EvolutionModel = "F81",
                                 Method = "bayes",
                                 Phylogeny = "../../../SpeciesTree_rooted.txt",
                                 OrthoPath = "Only_ORTHOLOGUES/")
}

#------------------------------------

SitesTableA = "/home/lalibelulalo/TESIS/Clados/Callothrix_clade/PALINDROMES/GCGATCGC/336-3/Orthologues_Palindrome_sites.txt"
SitesTableB = "/home/lalibelulalo/TESIS/Clados/Callothrix_clade/PALINDROMES/GCGATCGC/NIES-3974/Orthologues_Palindrome_sites.txt"
SitesTableC = "/home/lalibelulalo/TESIS/Clados/Callothrix_clade/PALINDROMES/GCGATCGC/PCC_6303/Orthologues_Palindrome_sites.txt"


SitesA <- read.table(SitesTableA, sep = "\t", header = TRUE)
SitesA$MergedCoordsOrth <- paste(SitesA$START,"-",SitesA$END,"-",SitesA$FILE)

SitesB <- read.table(SitesTableB, sep = "\t", header = TRUE)
SitesB$MergedCoordsOrth <- paste(SitesB$START,"-",SitesB$END,"-",SitesB$FILE)

SitesC <- read.table(SitesTableC, sep = "\t", header = TRUE)
SitesC$MergedCoordsOrth <- paste(SitesC$START,"-",SitesC$END,"-",SitesC$FILE)


Full_AB <- full_join(SitesA, SitesB, by = "MergedCoordsOrth")
Full_AB <- Full_AB %>% distinct()
Full_AB <- Full_AB[ , c(1,2,3,4,5)]

Full_BC <- full_join(SitesB, SitesC, by = "MergedCoordsOrth")
Full_BC <- Full_BC %>% distinct()
Full_BC <- Full_BC[ , c(1,2,3,4,5)]

Full_AC <- full_join(SitesA, SitesC, by = "MergedCoordsOrth")
Full_AC <- Full_AC %>% distinct()
Full_AC <- Full_AC[ , c(1,2,3,4,5)]

Cal_336 <- SitesA[!(SitesA$MergedCoordsOrth %in% Full_BC$MergedCoordsOrth),]
Cal_3974 <- SitesB[!(SitesB$MergedCoordsOrth %in% Full_AC$MergedCoordsOrth),]
Cal_6303 <- SitesC[!(SitesC$MergedCoordsOrth %in% Full_AB$MergedCoordsOrth),]

Cal_HIP_uniq_sites <- Cal_336 %>%
  bind_rows(Cal_3974) %>%
  bind_rows(Cal_6303)

Cal_HIP_uniq_sites <- Cal_HIP_uniq_sites[ , c(1,2,3,4)]
write.table(Cal_HIP_uniq_sites,file="Cal_HIP_uniq_sites.txt",sep="\t",row.names = FALSE,col.names = TRUE)


source("/home/lalibelulalo/TESIS/ASR_Orth_Functions/NodeAndEdges.R")
setwd("/home/lalibelulalo/TESIS/Clados/Callothrix_clade/PALINDROMES/GCGATCGC/SUBLCLADE/")
Nodes.Edges <- Create_Transition_Table_No_Fit(SitesTable = "Cal_HIP_uniq_sites.txt",
                                              EvolutionModel = "F81",
                                              Method = "bayes",
                                              Phylogeny = "../../../SpeciesTree_rooted.txt",
                                              OrthoPath = "Only_ORTHOLOGUES/")

#------------------------------------------------------------------------------------------------------------------------------------------------
setwd("/home/lalibelulalo/TESIS/Clados/Callothrix_clade/PALINDROMES/GCGATCGC/SUBLCLADE/")
SecondTableA = "/home/lalibelulalo/TESIS/Clados/Callothrix_clade/PALINDROMES/GCGATCGC/336-3/Orthologues_Palindrome_sites.AllFrames.SECOND.txt"
SecondTableB = "/home/lalibelulalo/TESIS/Clados/Callothrix_clade/PALINDROMES/GCGATCGC/NIES-3974/Orthologues_Palindrome_sites.AllFrames.SECOND.txt"
SecondTableC = "/home/lalibelulalo/TESIS/Clados/Callothrix_clade/PALINDROMES/GCGATCGC/PCC_6303/Orthologues_Palindrome_sites.AllFrames.SECOND.txt"

SecondA <- read.table(SecondTableA, sep = "\t", header = TRUE)
SecondA$MergedCoordsOrth <- paste(SecondA$FILE,"-",SecondA$Spp,"-",SecondA$START)

SecondB <- read.table(SecondTableB, sep = "\t", header = TRUE)
SecondB$MergedCoordsOrth <- paste(SecondB$FILE,"-",SecondB$Spp,"-",SecondB$START)

SecondC <- read.table(SecondTableC, sep = "\t", header = TRUE)
SecondC$MergedCoordsOrth <- paste(SecondC$FILE,"-",SecondC$Spp,"-",SecondC$START)


Full_AB <- full_join(SecondA, SecondB, by = "MergedCoordsOrth")
Full_AB <- Full_AB %>% distinct()

Full_BC <- full_join(SecondB, SecondC, by = "MergedCoordsOrth")
Full_BC <- Full_BC %>% distinct()

Full_AC <- full_join(SecondA, SecondC, by = "MergedCoordsOrth")
Full_AC <- Full_AC %>% distinct()


Cal_336 <- SecondA[!(SecondA$MergedCoordsOrth %in% Full_BC$MergedCoordsOrth),]
Cal_336 <- Cal_336[ , c(1:8)]

Cal_3974 <- SecondB[!(SecondB$MergedCoordsOrth %in% Full_AC$MergedCoordsOrth),]
Cal_3974 <- Cal_3974[ , c(1:8)]

Cal_6303 <- SecondC[!(SecondC$MergedCoordsOrth %in% Full_AB$MergedCoordsOrth),]
Cal_6303 <- Cal_6303[ , c(1:8)]

Cal_HIP_uniq_sites <- Cal_336 %>%
  bind_rows(Cal_3974) %>%
  bind_rows(Cal_6303)

#Cal_HIP_uniq_sites <- Cal_HIP_uniq_sites[ , c(1,2,3,4)]
write.table(Cal_HIP_uniq_sites,file="Cal_HIP_uniq_sites.AllFrames.SECOND.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote = FALSE)


RFS <- length(unique(system('awk \'{if(NR!=1) {print $5}}\' Cal_HIP_uniq_sites.AllFrames.SECOND.txt |uniq',intern = TRUE)))

for(RF in 1:RFS){ ## Para los tres marcos de lectura (RF)
  CreateCodonMutationsTable(Tree = '../../../SpeciesTree_rooted.txt',
                            Second = 'Cal_HIP_uniq_sites.AllFrames.SECOND.txt',
                            TaxonNumber = 6,
                            RF = RF)
  system(paste0('python3 ../../../../CodonMutations.py RF',RF,'.RECONSTRUCTION.rooted.parents.txt ',RF,' codon_mutations_RF',RF,'.rooted.txt GCGATCGC'))
}

system('cat *.rooted.txt >codon_mutations_AllF.rooted.txt')



