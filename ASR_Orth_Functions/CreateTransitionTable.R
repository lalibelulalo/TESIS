Create_Transition_Table_No_Fit <- function(SitesTable,EvolutionModel,Method,Phylogeny,OrthoPath) {
  ## Primero cargo Orthologues_Palindrome_sites.txt del paso anterior.
  Sites <- read.table(SitesTable, sep = "\t", header = TRUE)
  
  ## Creo la matriz que va a contener unicamente las transiciones
  DF <- matrix(0,
               nrow = 0,
               ncol = 2)
  colnames(DF) <- c("from","to")
  ##Creo una matriz adicional que tendra las transiciones y una columna extra con la direccion de la transición
  LINKS <- matrix(0,
                  nrow = 0,
                  ncol = 3)
  colnames(LINKS) <- c("from","to","direction")
  
  #Creo una matriz adicional que tendra las deleciones:
  DELECIONES <- matrix(0,
                       nrow = 0,
                       ncol = 4)
  colnames(DELECIONES) <- c("ortólogo","sitio","deleciones","species")
  
  ## Cargo el modelo de evolucion,
  ## el metodo de reconstrucción (bayesiano),
  ## la filogenia (previamente calculada con orthofinder).
  EvolModel = EvolutionModel
  #Method = "bayes"
  Tree = ggtree::read.tree(Phylogeny)
  
  for (ORTH in 1:length(Sites[,1])){ ## HAY 2593 ortólogos
    Alignment = paste0(OrthoPath,Sites[ORTH,1])
    SI = Sites[ORTH,3]
    EI = Sites[ORTH,4]
    
    ## Cargo la alineación multiple
    cyanobacterias <- phangorn::read.phyDat(Alignment, format = "interleaved")
    
    ## Cambio el attributo site.patter a FALSE para poder acceder a todos los indices
    cyanobacterias <- subset(cyanobacterias, site.pattern = FALSE) 
    #TreeR = pratchet(cyanobacterias, trace = 0)|>acctran(cyanobacterias)
    phangorn::parsimony(Tree, cyanobacterias)
    
    ## Estimo los parametros para el modelo
    fit <- phangorn::pml(Tree, cyanobacterias, model = EvolModel)
    #fit <- phangorn::optim.pml(fit, model=EvolModel, control = phangorn::pml.control(trace=0))
    
    #Tree2 <- fit[["tree"]]
    
    ## Calculo la probabilidad del árbol filogenético dado el alineamiento y el modelo
    if(Method == "bayes"){
      anc.bayes <- phangorn::ancestral.pml(fit, type="bayes", return = "prob")
      Reconstruction <- anc.bayes
    }
    if(Method == "ml"){
      anc.ml <- phangorn::ancestral.pml(fit, "ml")
      Reconstruction <- anc.ml
    }
    
    #plotAnc(Tree2, anc.bayes, 1)
    PalInterval = SI:EI ## ESTE ES EL INTERVALO DEL SITIO PALINDROMICO
    
    ## Obtengo las posiciones para la posicion PalInterval
    PalindromeNucleotidePositions <- attributes(cyanobacterias)$index[PalInterval]
    
    ## Nombre de los nodos
    Nodos <- attributes(Reconstruction)$names
    
    ########################################################################################################
    ## Transformo el arbol en data frame para poder agregar los palindromos de cada nodo
    a.1 <- tidytree::as_tibble(Tree)
    ## Extraigo el Boostrap por si se ocupa
    NA.NODE = length(cyanobacterias)+1
    for (i in NA.NODE:length(Nodos)){ ## ESTO DEPENDE DE LA ESTRUCTURA DEL ARBOL
      a.1$label[i] <- a.1$node[i] # x$node[i-1]
    }
    ## Agrego los datos de los palindromos al arbol
    b.1 <- tibble::tibble(label = Nodos)
    c.1 <- tidytree::full_join(a.1, b.1, by = 'label')
    NodeLabels <- c.1$label
    
    TipLabels <- Tree$tip.label
    NodeLabels2 <- NodeLabels[-c(1:length(TipLabels))]
    
    ## Creo una matriz de 8 x length(Nodos) (8 nucleotidos, numero de nodos) para guardar los resultados
    ScoresMtx <- matrix(0,
                        nrow = length(Nodos),
                        ncol = 8)
    colnames(ScoresMtx) <- c("1","2","3","4","5","6","7","8")
    rownames(ScoresMtx) <- Nodos
    
    
    ## Relleno la matriz con los scores
    j = 1
    i = 1
    k = 0
    for (Position in PalindromeNucleotidePositions){
      for (Nodo in TipLabels){
        PositionScores <- Reconstruction[[Nodo]][Position,] ## Extraigo los scores para el nodo "NODO" en la posicion "POSITION" de la secuencia
        if (length(unique(PositionScores)) == 1){
          k=k+1
          WinnerNucleotide = "-"
          ScoresMtx[i,j] <- WinnerNucleotide
          i=i+1 ## contador para cada Nodo
        }else{
          PositionWinnerScore <- max(PositionScores) ## Extraigo el valor mas grande. Este es el valor de la probabilidad de que el nucleotido ancestral sea ese
          WinnerNucleotide <- which(PositionScores == PositionWinnerScore, arr.ind = T) ## pregunto a que numero (Nucleótido) corresponde ese valor
          ScoresMtx[i,j] <- WinnerNucleotide ## Guardo dicho valor en la matriz de resultados
          i=i+1 ## contador para cada Nodo
        }
      }
      j=j+1 ## contador para cada posición
      i=1 ## Reiniciamoa el contador de Nodo para la siguiente posicion
    }
    
    ## Relleno la matriz con los scores
    j = 1
    i = length(TipLabels)+1
    k = 0
    #(length(TipLabels)+1):length(NodeLabels)
    for (Position in PalindromeNucleotidePositions){
      for (Nodo in NodeLabels2){
        PositionScores <- Reconstruction[[Nodo]][Position,] ## Extraigo los scores para el nodo "NODO" en la posicion "POSITION" de la secuencia
        if (length(unique(PositionScores)) == 1){
          k=k+1
          WinnerNucleotide = "X"
          ScoresMtx[i,j] <- WinnerNucleotide
          i=i+1 ## contador para cada Nodo
        }else{
          PositionWinnerScore <- max(PositionScores) ## Extraigo el valor mas grande. Este es el valor de la probabilidad de que el nucleotido ancestral sea ese
          WinnerNucleotide <- which(PositionScores == PositionWinnerScore, arr.ind = T) ## pregunto a que numero (Nucleótido) corresponde ese valor
          ScoresMtx[i,j] <- WinnerNucleotide ## Guardo dicho valor en la matriz de resultados
          i=i+1 ## contador para cada Nodo
        }
      }
      j=j+1 ## contador para cada posición
      i=length(TipLabels)+1 ## Reiniciamoa el contador de Nodo para la siguiente posicion
    }
    ########################################################################################################
    
    ## EN ESTA PARTE ARMAREMOS EL PALINDROMO PARA CADA NODO
    
    ## Creo un matriz de 1x14 para guardar la etiqueta (palindromo) de cada nodo
    PalsMtx <- matrix(0,
                      nrow = length(Nodos),
                      ncol = 1)
    colnames(PalsMtx) <- c("Palindrome")
    rownames(PalsMtx) <- Nodos
    
    ## Extraigo cada palindromo, cambio el codigo de numero por el codigo de letras y lo guardo en la matriz 1x14
    i=0
    j=0
    for (Nodo in Nodos){
      j=j+1
      Palindrome <- as.character(ScoresMtx[j,]) # Este es el palindromo
      pal =""
      for (nuc in Palindrome){
        if(nuc== 1){nuc = 'A'}
        if(nuc== 2){nuc = 'C'}
        if(nuc== 3){nuc = 'G'}
        if(nuc== 4){nuc = 'T'}
        pal = paste(pal,nuc,sep = "")
      }
      PalsMtx[j,1] <- pal
    }
    
    ## ESTA PATRTE ES PARA GUARDAR LOS PALINDROMOS DE CADA NODO
    l <- as.list(PalsMtx[,1]) 
    Nodes = names(l)
    NodePals = as.character(PalsMtx[,1])
    
    ## Transformo el arbol en data frame para poder agregar los palindromos de cada nodo
    x <- tidytree::as_tibble(Tree)
    
    ## Extraigo el Boostrap por si se ocupa
    NA.NODE = length(cyanobacterias)+1
    BootStrap <- x$label[NA.NODE:length(Nodes)]
    for (i in NA.NODE:length(Nodes)){ ## ESTO DEPENDE DE LA ESTRUCTURA DEL ARBOL
      x$label[i] <- x$node[i] # x$node[i-1]
    }
    ## Agrego los datos de los palindromos al arbol
    d <- tibble::tibble(label = Nodes,
                        Palindromo = NodePals)
    y <- tidytree::full_join(x, d, by = 'label')
    
    #------------------------
    
    ## Esta es la configuracion del arbol
    from<-y$label#c("336-3","7","7","8","8","9","9","10","10")
    to<-y$parent#c("7","NIES-3974","8","PCC_6303","9","PCC_7716","10","NIES-4071","NIES-4105")
    
    LinksMtx <- matrix(0,
                       nrow = length(from),
                       ncol = 2)
    colnames(LinksMtx) <- c("from","to")
    
    for (i in 1:length(from)){
      LinksMtx[i,1] = from[i]
      LinksMtx[i,2] = to[i]
    }
    
    ## AQUI VOY A CREAR LA MATRIZ QUE ME SERVIRA COMO DICCIONARIO PARA RELLENAR LA MATRIZ DE TRANSICIONES
    Mat <- y[,4:5]
    #Mat <- Mat %>% ## AQJUI QUITO EL RENGLON 7 QUE NO CONTIENE NADA. ESTE CORRESPONDE A LA RAIZ DEL ARBOL
    #  filter(!is.na(Palindromo))
    Mat
    
    Mat2 <- as.data.frame(Mat) ## CONVIERTO A DF
    Mat2
    
    #~~~~~~~~~~~~~~~~~~~~~~~~
    d = 0
    SPP = c()
    for (i in 1:length(Mat2[,1])) {
      if(Mat2[i,2] == '--------'){
        d = d+1
        spp = Mat2[i,1]
        SPP = c(SPP,spp)
      }
    }
    Spp <- paste(SPP, collapse = ",")
    deleciones <-c()
    if(d != 0){
      print(paste0("Hay ",d," delecion(es) en el archivo ",Sites[ORTH,1],". Intervarlo ",Sites[ORTH,3],":",Sites[ORTH,4],"."))
      SITE = paste0(Sites[ORTH,3],":",Sites[ORTH,4])
      deleciones <- c(Sites[ORTH,1],SITE,d,Spp)
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~
    
    ## AQUI CREO LA MATRIZ DE TRANSICIONES.
    LinksMtx2 <- matrix(0,
                        nrow = length(from),
                        ncol = 2)
    colnames(LinksMtx2) <- c("from","to")
    
    ## AQUI RELLENO LA MATRIZ DE TRANSICIONES USANDO COMO DICCIONARIO A LA MATRIZ ANTERIOR (Mat2)
    for (i in 1:length(from)) {
      for (j in 1:length(from)){
        if (LinksMtx[i,1]==Mat2[j,1]){
          LinksMtx2[i,1] = Mat2[j,2]
        }
        if (LinksMtx[i,2]==Mat2[j,1]){
          LinksMtx2[i,2] = Mat2[j,2]
        }
      }
    }
    
    directions<-paste(from,to,sep="--") ## ESTE VECTOR CONTIENE LAS DIRECCIONES
    Links0 = as.data.frame(LinksMtx2) ## HAGO UNA COPIA DE LA MATRIZ DE TRANSICIONES
    Links0$direction<-directions ## AGREGO COMO COLUMNA NUEVA A "directions"
    Links0
    LINKS<-rbind(LINKS,Links0) ## UNO LA MATRIZ ACTUAL DE TRANSICIONES(CON DIRECCIONES) CON LA ANTERIOR PARA CREAR UNA SOLA
    
    df <- as.data.frame(LinksMtx2) ## CONVIERTO A DF LA MATRIZ DE TRANSICIONES
    DF <- rbind(DF,df) ## UNO LA MATRIZ ACTUAL DE TRANSICIONES CON LA ANTERIOR PARA CREAR UNA SOLA
    print(paste0("Sitio ",ORTH, " de ",length(Sites$FILE),"."))
    DELECIONES <- rbind(DELECIONES,deleciones)
  }
  length(DELECIONES[,1])
  DELECIONES = as.data.frame(DELECIONES)
  DELECIONES = DELECIONES[order(DELECIONES$deleciones,decreasing = TRUE),]
  write.table(DELECIONES,file="DELECIONES.txt",sep="\t",row.names = FALSE,col.names = TRUE)
  write.table(DF,file="TRANSICIONES.txt",sep="\t",row.names = FALSE,col.names = TRUE)
  write.table(LINKS,file="TRANSICIONES_DIRECCION.txt",sep="\t",row.names = FALSE,col.names = TRUE)
  #TABLE1 = as.data.frame(DF)
  #TABLE2 = as.data.frame(LINKS)
  #return(list(nodes,edges,TABLE2))
}