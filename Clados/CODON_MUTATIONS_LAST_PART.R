library(tidyverse)
setwd('/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/Cyanobacterium_clade/PALINDROMES/GCGATCGC/PCC_8801/')

#Tree <- ggtree::read.tree('../../../SpeciesTree_rooted.txt')
#File = "Orthologues_Palindrome_sites.AllFrames.SECOND.txt"
#table = read.table(File,header = TRUE,sep = "\t",row.names = NULL)
#TaxonNumber = 6
#RF=1

CreateCodonMutationsTable <- function(Tree,Second,TaxonNumber,RF) {
  Tree <- ggtree::read.tree(Tree)
  File = Second
  table = read.table(File,header = TRUE,sep = "\t",row.names = NULL)
  TaxonNumber = TaxonNumber
  RF=RF
  
  Spps <- unique(table[,2])
  ## Primero extraigo el nombre de los nodos parentales del de la filogenia:
  x <- tidytree::as_tibble(Tree)
  #UNROOTED = 7:10
  ROOTED = (TaxonNumber+1):length(Spps)
  for (i in ROOTED){
    x$label[i] <- x$node[i]
  }
  FromTo = as.data.frame(x[ , c(1,4)])
  colnames(FromTo) <- c("from","to")
  
  ## Primero filtro el Marco de lectura 1
  tableRF1<-table%>%
    filter(ReadingFrame==RF)
  
  ## Agrego 3 columnas que contendrán: el nodo parental, secuencia parental y AA parental.
  tableRF1$parent = NA
  tableRF1$parentPAL = NA
  tableRF1$parentAA = NA
  
  ## Quito espacio de cada codon
  for (i in 1:length(tableRF1[,1])){
    #print(paste0('Linea ',i,' de ',length(tableRF1[,1]),'.'))
    tableRF1[i,7] <- gsub(' ','',tableRF1[i,7])
  }
  
  k=0
  for (i in 1:length(tableRF1[,1])){
    k=k+1 ## esta k debe ser IGUAL O MENOR que 11 ya que tengo 11 nodos
    ActualSpp <- tableRF1[i,2]
    ParentalSpp <- FromTo[FromTo$to == ActualSpp,][1,1]
    SubPosInic <- i-k+1 ## es en estos subconjuntos de CADA n nodos (11 en este caso ya que hay 11 especies) donde voya buscar los ANCESTRALES para cada ACTUAL
    ## por lo tanto si este numero es mas grande podria estar buscando el nodo de un conjunto en OTRO conjunto de 11 INCORRECTO
    
    for (j in SubPosInic:(SubPosInic + (length(Spps)-1))){ ## Voy a ir buscando en intervalos de n nodos
      if (tableRF1[j,2] == ParentalSpp){
        ParentalPAL = tableRF1[j,7]
        ParentalAA = tableRF1[j,8]
      }
    }
    tableRF1[i,9] <- ParentalSpp
    tableRF1[i,10] <- ParentalPAL
    tableRF1[i,11] <- ParentalAA
    print(paste0('RENGLON ', i,'.'))
    print (paste0('La especie actual es: ', ActualSpp,'.', ' Su Palindromo es:',tableRF1[i,7],'.'))
    print (paste0('La especie parental es: ', ParentalSpp,'.', ' Su Palindromo es:',ParentalPAL,'.'))
    print ('---------------------')
    
    if(k==length(Spps)){ ## si k es igual al numero de nodos
      k=0
    }
  }
  
  #Quito el nodo 7 ya que es el nodo mas antiguo y no tiene datos “ancestrales”.
  tableRF1<-tableRF1%>%
    filter(Spp!='7')
  write.table(tableRF1,file=paste0('RF',RF,'.RECONSTRUCTION.rooted.parents.txt'),sep="\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
}


CreateCodonMutationsTable(Tree = '../../../SpeciesTree_rooted.txt',
                          Second = 'Orthologues_Palindrome_sites.AllFrames.SECOND.txt',
                          TaxonNumber = 6,
                          RF = 3)

system('python3 ../../../../CodonMutations.py RF1.RECONSTRUCTION.rooted.parents.txt 1 codon_mutations_RF1.rooted.txt GCGATCGC')
system('python3 ../../../../CodonMutations.py RF2.RECONSTRUCTION.rooted.parents.txt 2 codon_mutations_RF2.rooted.txt GCGATCGC')
system('python3 ../../../../CodonMutations.py RF3.RECONSTRUCTION.rooted.parents.txt 3 codon_mutations_RF3.rooted.txt GCGATCGC')

#---------------------------------------------------------------------------------------------------------------------------------------
setwd('/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/Geminocystis_clade/PALINDROMES/GCGATCGC/NIES-4102/')
CreateCodonMutationsTable(Tree = '../../../SpeciesTree_rooted.txt',
                          Second = 'Orthologues_Palindrome_sites.AllFrames.SECOND.txt',
                          TaxonNumber = 11,
                          RF = 3)
system('python3 ../../../../CodonMutations.py RF1.RECONSTRUCTION.rooted.parents.txt 1 codon_mutations_RF1.rooted.txt GCGATCGC')
system('python3 ../../../../CodonMutations.py RF2.RECONSTRUCTION.rooted.parents.txt 2 codon_mutations_RF2.rooted.txt GCGATCGC')
system('python3 ../../../../CodonMutations.py RF3.RECONSTRUCTION.rooted.parents.txt 3 codon_mutations_RF3.rooted.txt GCGATCGC')

#---------------------------------------------------------------------------------------------------------------------------------------
setwd('/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/Pseudoanabaena_clade/PALINDROMES/GCGATCGC/PCC_7336/')
CreateCodonMutationsTable(Tree = '../../../SpeciesTree_rooted.txt',
                          Second = 'Orthologues_Palindrome_sites.AllFrames.SECOND.txt',
                          TaxonNumber = 6,
                          RF = 3)
system('python3 ../../../../CodonMutations.py RF1.RECONSTRUCTION.rooted.parents.txt 1 codon_mutations_RF1.rooted.txt GCGATCGC')
system('python3 ../../../../CodonMutations.py RF2.RECONSTRUCTION.rooted.parents.txt 2 codon_mutations_RF2.rooted.txt GCGATCGC')
system('python3 ../../../../CodonMutations.py RF3.RECONSTRUCTION.rooted.parents.txt 3 codon_mutations_RF3.rooted.txt GCGATCGC')

#---------------------------------------------------------------------------------------------------------------------------------------
setwd('/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/Synechococcus_cyanobium_clade/PALINDROMES/GCGATCGC/CBW1004/')
CreateCodonMutationsTable(Tree = '../../../SpeciesTree_rooted.txt',
                          Second = 'Orthologues_Palindrome_sites.AllFrames.SECOND.txt',
                          TaxonNumber = 8,
                          RF = 3)
system('python3 ../../../../CodonMutations.py RF1.RECONSTRUCTION.rooted.parents.txt 1 codon_mutations_RF1.rooted.txt GCGATCGC')
system('python3 ../../../../CodonMutations.py RF2.RECONSTRUCTION.rooted.parents.txt 2 codon_mutations_RF2.rooted.txt GCGATCGC')
system('python3 ../../../../CodonMutations.py RF3.RECONSTRUCTION.rooted.parents.txt 3 codon_mutations_RF3.rooted.txt GCGATCGC')

#---------------------------------------------------------------------------------------------------------------------------------------
setwd('/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/Thermosynechococcus_clade/PALINDROMES/GCGATCGC/SCTE542/')
CreateCodonMutationsTable(Tree = '../../../SpeciesTree_rooted.txt',
                          Second = 'Orthologues_Palindrome_sites.AllFrames.SECOND.txt',
                          TaxonNumber = 12,
                          RF = 3)
system('python3 ../../../../CodonMutations.py RF1.RECONSTRUCTION.rooted.parents.txt 1 codon_mutations_RF1.rooted.txt GCGATCGC')
system('python3 ../../../../CodonMutations.py RF2.RECONSTRUCTION.rooted.parents.txt 2 codon_mutations_RF2.rooted.txt GCGATCGC')
system('python3 ../../../../CodonMutations.py RF3.RECONSTRUCTION.rooted.parents.txt 3 codon_mutations_RF3.rooted.txt GCGATCGC')
