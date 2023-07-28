library(ggplot2)
library(ggtree)
library(ape)
library(tidyverse)
library(tidytree)
library(phangorn)
library(dplyr)

source("/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/ASR_Orth_Functions/NodeAndEdges.R")
setwd("/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/Cyanobacterium_clade/PALINDROMES/GCGATCGC/")

Sites <- read.table("Orthologues_Palindrome_sites.AllFrames.FIRST.NIES-4071.ACGCGCGT.txt", sep = "\t", header = TRUE)
EvolModel = "F81"
Method = "bayes"
Tree = read.tree("../../SpeciesTree_rooted.txt")
PATH = 'Only_ORTHOLOGUES/'
Sites<-Sites%>% 
  filter(Spp=='NIES-4071')




Create_Reconstruction_Files <- function(SitesTable,EvolutionModel,Method,Phylogeny,OrthoPath) {
  v = 0
  #ORTH = 1
  for (ORTH in 1:length(SitesTable[,1])){
    ## CARGO EL ALINEAMIENTO
    Alignment =paste0(PATH,SitesTable[ORTH,1])
    cyanobacterias <- read.phyDat(Alignment, format = "interleaved")
    cyanobacterias <- subset(cyanobacterias, site.pattern = FALSE) 
    parsimony(Phylogeny, cyanobacterias)
    
    ## CARGO EL MODELO
    fit <- pml(Phylogeny, cyanobacterias,model=EvolutionModel)
    #fit <- optim.pml(fit, model=EvolModel, control = pml.control(trace=0),)
    
    ## CARGO La NUEVA FILOGENIA SIN RAIZ
    #Tree2 <- fit[["tree"]]
    
    ## HAGO LA RECONSTRUCCIÓN
    if(Method == "bayes"){
      anc.bayes <- ancestral.pml(fit, type="bayes", return = "prob")
      Reconstruction <- anc.bayes
    }
    if(Method == "ml"){
      anc.ml <- ancestral.pml(fit, "ml")
      Reconstruction <- anc.ml
    }
    
    ## LEO EL SITIO DE INICIO Y FINAL DE LA SECUENCIA
    ## ES DECIR EL INTERVALO COMPLETO DEL ORTOLOGO
    PalInterval0 = 1:Sites[ORTH,6]
    
    ## Obtengo las posiciones para la posicion PalInterval
    PalindromeNucleotidePositions0 <- attributes(cyanobacterias)$index[PalInterval0]
    
    ## Nombre de los nodos
    Nodos0 <- attributes(Reconstruction)$names
    
    ################## EN ESTA PARTE ARMAREMOS LOS PALINDROMOS DE LAS PUNTAS ###################
    #a.1 <- as_tibble(Tree2)
    a.1 <- as_tibble(Phylogeny)
    ## Extraigo el Boostrap por si se ocupa
    NA.NODE = length(cyanobacterias)+1
    for (i in NA.NODE:length(Nodos0)){ ## ESTO DEPENDE DE LA ESTRUCTURA DEL ARBOL
      a.1$label[i] <- a.1$node[i] 
    }
    ## Agrego los datos de los palindromos al arbol
    b.1 <- tibble::tibble(label = Nodos0)
    c.1 <- full_join(a.1, b.1, by = 'label')
    NodeLabels <- c.1$label
    
    #TipLabels <- Tree2$tip.label
    TipLabels <- Phylogeny$tip.label
    NodeLabels2 <- NodeLabels[-c(1:length(TipLabels))]
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ## Creo una matriz de 8 x length(Nodos) (8 nucleotidos, numero de nodos) para guardar los resultados
    ScoresMtx0 <- matrix(0,
                         nrow = length(Nodos0),
                         ncol = length(PalInterval0))
    colnames(ScoresMtx0) <- PalInterval0
    rownames(ScoresMtx0) <- Nodos0
    ## Relleno la matriz con los scores
    j0 = 1
    i0 = 1
    k0 = 0
    for (Position in PalindromeNucleotidePositions0){
      for (Nodo0 in TipLabels){
        PositionScores0 <- Reconstruction[[Nodo0]][Position,] ## Extraigo los scores para el nodo "NODO" en la posicion "POSITION" de la secuencia
        if (length(unique(PositionScores0)) == 1){
          k0=k0+1
          WinnerNucleotide0 = "-"
          ScoresMtx0[i0,j0] <- WinnerNucleotide0
          i0=i0+1 ## contador para cada Nodo
        }else{
          PositionWinnerScore0 <- max(PositionScores0) ## Extraigo el valor mas grande. Este es el valor de la probabilidad de que el nucleotido ancestral sea ese
          WinnerNucleotide0 <- which(PositionScores0 == PositionWinnerScore0, arr.ind = T) ## pregunto a que numero (Nucleótido) corresponde ese valor
          ScoresMtx0[i0,j0] <- WinnerNucleotide0 ## Guardo dicho valor en la matriz de resultados
          i0=i0+1 ## contador para cada Nodo
        }
      }
      j0=j0+1 ## contador para cada posición
      i0=1 ## Reiniciamoa el contador de Nodo para la siguiente posicion
    }
    ################## EN ESTA PARTE ARMAREMOS LOS PALINDROMOS DE LOS NODOS ###################
    ## Relleno la matriz con los scores
    j0 = 1
    i0 = length(TipLabels)+1
    k0 = 0
    #(length(TipLabels)+1):length(NodeLabels)
    for (Position0 in PalindromeNucleotidePositions0){
      for (Nodo0 in NodeLabels2){
        PositionScores0 <- Reconstruction[[Nodo0]][Position0,] ## Extraigo los scores para el nodo "NODO" en la posicion "POSITION" de la secuencia
        if (length(unique(PositionScores0)) == 1){
          k0=k0+1
          WinnerNucleotide0 = "X"
          ScoresMtx[i0,j0] <- WinnerNucleotide
          i0=i0+1 ## contador para cada Nodo
        }else{
          PositionWinnerScore0 <- max(PositionScores0) ## Extraigo el valor mas grande. Este es el valor de la probabilidad de que el nucleotido ancestral sea ese
          WinnerNucleotide0 <- which(PositionScores0 == PositionWinnerScore0, arr.ind = T) ## pregunto a que numero (Nucleótido) corresponde ese valor
          ScoresMtx0[i0,j0] <- WinnerNucleotide0 ## Guardo dicho valor en la matriz de resultados
          i0=i0+1 ## contador para cada Nodo
        }
      }
      j0=j0+1 ## contador para cada posición
      i0=length(TipLabels)+1 ## Reiniciamoa el contador de Nodo para la siguiente posicion
    }
    ################## EN ESTA PARTE ARMAREMOS LOS PALINDROMOS DE LOS NODOS ###################
    
    ## EN ESTA PARTE ARMAREMOS EL PALINDROMO PARA CADA NODO
    
    ## Creo un matriz de 1x14 para guardar la etiqueta (palindromo) de cada nodo
    PalsMtx0 <- matrix(0,
                       nrow = length(Nodos0),
                       ncol = 1)
    colnames(PalsMtx0) <- c("Sequence")
    rownames(PalsMtx0) <- Nodos0
    
    ## Extraigo cada palindromo, cambio el codigo de numero por el codigo de letras y lo guardo en la matriz 1x14
    i=0
    j=0
    for (Nodo in Nodos0){
      j=j+1
      Palindrome <- as.character(ScoresMtx0[j,]) # Este es el palindromo
      pal =""
      for (nuc in Palindrome){
        if(nuc== 1){nuc = 'A'}
        if(nuc== 2){nuc = 'C'}
        if(nuc== 3){nuc = 'G'}
        if(nuc== 4){nuc = 'T'}
        pal = paste(pal,nuc,sep = "")
      }
      PalsMtx0[j,1] <- pal
    }
    ## CREO EL ARCHIVO FASTA CON LA RECONSTRUCCION CONMPLETA DE SECUENCIAS
    FASTA <-c()
    for (i in 1:length(PalsMtx0[,1])){
      header <- paste0('>',rownames(PalsMtx0)[i])
      seq <- PalsMtx0[i,1]
      FASTA <- rbind(FASTA,header)
      FASTA <- rbind(FASTA,seq)
    }
    
    ## ESCRIBO EL ARCHIVO FASTA
    write.table(FASTA,file=paste0(OrthoPath,SitesTable[ORTH,1],'.RECONSTRUCTION.fasta'),sep="\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    print(paste0("Sitio ",ORTH, " de ",length(SitesTable$FILE),"."))
  }
}