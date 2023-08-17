library(dplyr)

Create_Transition_Table <- function(SitesTable,EvolutionModel,Method,Phylogeny,OrthoPath) {
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
    fit <- phangorn::pml(Tree, cyanobacterias)
    fit <- phangorn::optim.pml(fit, model=EvolModel, control = phangorn::pml.control(trace=0))
    
    Tree2 <- fit[["tree"]]
    
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
    a.1 <- tidytree::as_tibble(Tree2)
    ## Extraigo el Boostrap por si se ocupa
    NA.NODE = length(cyanobacterias)+1
    for (i in NA.NODE:length(Nodos)){ ## ESTO DEPENDE DE LA ESTRUCTURA DEL ARBOL
      a.1$label[i] <- a.1$node[i] # x$node[i-1]
    }
    ## Agrego los datos de los palindromos al arbol
    b.1 <- tibble::tibble(label = Nodos)
    c.1 <- tidytree::full_join(a.1, b.1, by = 'label')
    NodeLabels <- c.1$label
    
    TipLabels <- Tree2$tip.label
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
    x <- tidytree::as_tibble(Tree2)
    
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
#---------------------------------------------------------------------------------------------------------------------------------------
Create_Network <- function(Transitions){
  DF = Transitions
  RowCts <- data.table::setDT(DF)[,list(Count=as.numeric(.N)),names(DF)] ## CUENTO EL NUMERO DE OCURRENCIAS DE CADA TRANSICIÓN #7211
  colnames(RowCts) <- c("from","to","weight") 
  RowCts <- transform(RowCts, weight = as.integer(weight)) ## CONVIERTO A ENTEROS LA COLUMNA "weight"
  RowCts <- tidytree::as_tibble(RowCts) ## CONVIERTO A TIBBLE
  RowCts #7088
  
  ## Filtro conteos bajos (<=7) y loops.
  ## Al final me quedo con 63 transiciones.
  RowCts<-RowCts%>% ##QUITO LOS LOOPS. ES DECIR QUITO AQUELLAS TRANSICIONES QUE VAYAN ASI MISMA
    filter(from!=to)
  length(RowCts$from) #5503
  RowCts<-RowCts%>% ## QUITO AQUELLAS TRANSICIONES CON UN PESO MENOR A 10
    filter(weight>=7)
  length(RowCts$from) #98
  
  ## Extraigo todos los nodos posibles y los agrego a una matriz que va a contener los nodos del grafo.
  ## Adiconalmente le agrego una columna score que me indicara que tan parecido es el aplindromo a HIP1
  ## (despues la quito porque no creo que sea un buen índice)
  
  GraphNodes <- c(RowCts$from, RowCts$to) ## OBTENGO LOS PALINDROMOS DE CADA NODO
  length(GraphNodes) #196
  GNMtx <- matrix(0, ## CREO LA MATRIZ QUE CONTENDRÁ LAS NODOS PARA EL GRAFO
                  nrow = length(GraphNodes),
                  ncol = 2)
  colnames(GNMtx) <- c("pal","score") ## ETIQUETO LAS COLUMNAS DE LA MATRIZ
  
  for (i in 1:length(GraphNodes)){ ## RELLENO LA MATRIZ
    GNMtx[i,1] = GraphNodes[i] ## AGREGO LOS NODOS
    GNMtx[i,2] = RecordLinkage::levenshteinSim(GraphNodes[i], 'GCGATCGC')*1000 ## AGREGO UN SCORE QUE SERÁ EL PORCENTAJE DE IDENTIDAD CON HIP1 MULTIPLICADO POR 1000
  }
  
  GNMtx <- as.data.frame(GNMtx) ## CONVIERTO LA MATRIZ DE NODOS A DF
  GNMtx <- GNMtx[!duplicated(GNMtx), ] ## ELIMINO LOS NODOS DUPLICADOS
  
  ## Creo un vector que enumerará los palindromos. Dicha enumeracion correspondera a su ID
  palindromes<-GNMtx[,1] #51  ## EXTRAIGO LOS NODOS DE LA MATRIZ DEL GRAFO. ESTO LO HAGO PARA POSTERIAMENTE ENUMERARLAS Y ASOCIAR CADA NUMERO COMO ID
  
  ids<-c() ## CREO UN VECTOR VACIO QUE CONTENDRÁ LOS ID'S
  for (k in 1:length(palindromes)){ ## ENUMERO LOS PALINDROMOS. ESTOS NUMEROS SERAN LOS IDS
    ids <-append(ids, k)
  }
  
  ## Creo una matriz nueva que contendrá los nodos del grafo pero con su ID para cada palindromo:
  Nodes2 <-cbind(GNMtx,ids) ## CREO UNA MATRIZ NUEVA CON LA MATRIZ DE NODOS Y LOS IDS. ESTA LA VOY A USAR COMO UN DICCIONARIO ID-NODO
  length(Nodes2$pal)
  Nodes2<-Nodes2[,c(3,1,2)] ## REORDENO LAS COLUMNAS DE LA MATRIZ
  colnames(Nodes2)<-c("id","label","score") ## ETIQUETO LAS COLUMNAS
  
  ## Creo la matriz de vertices del grafo y tambien sustituyo el nombre de los nodos por sus respectivos ID'S.
  ## Adicionalmente cambio los scores que tengan 0 por 1's ya que hay errores. Sin embargo, la columna score no la uso mas adelante.
  GEMtx <- RowCts ## ESTA VA A SER LA MATRIZ DE VERTICES
  Edges2 <- as.matrix(RowCts)#RowCts ESTA ES OTRA MATRIZ DE VERTICES DE LA CUAL VOY A SUSTITUIR LOS NODOS POR SUS IDS
  
  for (n in 1:length(Nodes2[,2])){ ##AQUI CAMBIARE LOS PALINDROMOS DE LOS NODOS POR IDS
    Edges2[Edges2==Nodes2[n,2]] <- as.numeric(Nodes2[n,1])
  }
  
  for (n in 1:length(Nodes2[,2])){ ## HANDLE 0's
    Nodes2[Nodes2==0] <- as.numeric(1)
  }
  
  ## Para empezar cargo los **nodos** y los **vertices**. En esta parte quito la columna score de los nodos.
  nodes = tidytree::as_tibble(Nodes2)
  nodes = select(nodes, -score)
  
  edges = tidytree::as_tibble(Edges2)
  edges = transform(edges, from = as.integer(from))
  edges = transform(edges, to = as.integer(to))
  edges = transform(edges, weight = as.integer(weight))
  edges = tidytree::as_tibble(edges)
  return(list(nodes,edges))
}
#---------------------------------------------------------------------------------------------------------------------------------------
Create_Directed_Network <- function(DirectedTransitions,Direction,Weight){
  #DIRECTION = Direction
  RowCtsL <- data.table::setDT(DirectedTransitions)[,list(Count=as.numeric(.N)),names(DirectedTransitions)] ## CUENTO EL NUMERO DE OCURRENCIAS DE CADA TRANSICIÓN #11520 ***
  colnames(RowCtsL) <- c("from","to","direction","weight") ## CAMBIO LA ETIQUETA "counts" A "weight" ***
  RowCtsL = transform(RowCtsL, weight = as.integer(weight)) ## CONVIERTO A ENTEROS LA COLUMNA "weight" ***
  RowCtsL <-tidytree::as_tibble(RowCtsL) ## CONVIERTO A TIBBLE**
  length(RowCtsL$from)
  RowCtsL<-RowCtsL%>% ##QUITO LOS LOOPS. ES DECIR QUITO AQUELLAS TRANSICIONES QUE VAYAN ASI MISMA ***
    filter(from!=to)
  RowCtsL<-RowCtsL%>% ## QUITO AQUELLAS TRANSICIONES CON UN PESO MENOR A 2 ***
    filter(weight>=Weight)
  RowCtsL<-RowCtsL%>% ## filtro las direcciones 9--8
    filter(direction==Direction)
  
  GraphNodesL <- c(RowCtsL$from, RowCtsL$to) ## OBTENGO LOS PALINDROMOS DE CADA NODO ***
  
  GNMtxL <- matrix(0, ## CREO LA MATRIZ QUE CONTENDRÁ LAS NODOS PARA EL GRAFO ***
                   nrow = length(GraphNodesL),
                   ncol = 2)
  colnames(GNMtxL) <- c("pal","score") ## ETIQUETO LAS COLUMNAS DE LA MATRIZ ***
  
  for (i in 1:length(GraphNodesL)){ ## RELLENO LA MATRIZ ***
    GNMtxL[i,1] = GraphNodesL[i] ## AGREGO LOS NODOS ***
    GNMtxL[i,2] = RecordLinkage::levenshteinSim(GraphNodesL[i], 'GCGATCGC')*1000 ## AGREGO UN SCORE QUE SERÁ EL PORCENTAJE DE IDENTIDAD CON HIP1 MULTIPLICADO POR 1000 ***
  }
  
  GNMtxL <- as.data.frame(GNMtxL) ## CONVIERTO LA MATRIZ DE NODOS A DF ***
  GNMtxL <- GNMtxL[!duplicated(GNMtxL), ] ## ELIMINO LOS NODOS DUPLICADOS ***
  palindromesL<-GNMtxL[,1] #25  ## EXTRAIGO LOS NODOS DE LA MATRIZ DEL GRAFO. ESTO LO HAGO PARA POSTERIAMENTE ENUMERARLAS Y ASOCIAR CADA NUMERO COMO ID ***
  
  idsL<-c() ## CREO UN VECTOR VACIO QUE CONTENDRÁ LOS ID'S ***
  for (k in 1:length(palindromesL)){ ## ENUMERO LOS PALINDROMOS. ESTOS NUMEROS SERAN LOS IDS ***
    idsL <-append(idsL, k)
  }
  
  Nodes2L <-cbind(GNMtxL,idsL) ## CREO UNA MATRIZ NUEVA CON LA MATRIZ DE NODOS Y LOS IDS. ESTA LA VOY A USAR COMO UN DICCIONARIO ID-NODO ***
  Nodes2L<-Nodes2L[,c(3,1,2)] ## REORDENO LAS COLUMNAS DE LA MATRIZ ***
  colnames(Nodes2L)<-c("id","label","score") ## ETIQUETO LAS COLUMNAS ***
  
  GEMtxL <- RowCtsL ## ESTA VA A SER LA MATRIZ DE VERTICES ***
  Edges2L <- as.matrix(RowCtsL)#RowCts ESTA ES OTRA MATRIZ DE VERTICES DE LA CUAL VOY A SUSTITUIR LOS NODOS POR SUS IDS ***
  
  for (n in 1:length(Nodes2L[,2])){ ##AQUI CAMBIARE LOS PALINDROMOS DE LOS NODOS POR IDS ***
    Edges2L[Edges2L==Nodes2L[n,2]] <- as.numeric(Nodes2L[n,1])
  }
  
  for (n in 1:length(Nodes2L[,2])){ ## HANDLE 0's
    Nodes2L[Nodes2L==0] <- as.numeric(1)
  }
  
  nodesL = tidytree::as_tibble(Nodes2L)
  nodesL = select(nodesL, -score)
  
  edgesL = tidytree::as_tibble(Edges2L)
  edgesL = transform(edgesL, from = as.integer(from))
  edgesL = transform(edgesL, to = as.integer(to))
  edgesL = transform(edgesL, weight = as.integer(weight))
  edgesL = tidytree::as_tibble(edgesL)
  
  return(list(nodesL,edgesL))
}
#---------------------------------------------------------------------------------------------------------------------------------------
Create_Reconstruction_Files <- function(SitesTable,EvolutionModel,Method,Phylogeny,OrthoPath,TAXON) {
  v = 0
  #ORTH = 1
  Sites<- read.table(SitesTable, sep = "\t", header = TRUE)
  Sites<-Sites%>% 
    filter(Spp==TAXON)
  Tree = ggtree::read.tree(Phylogeny)
  for (ORTH in 1:length(Sites[,1])){
    ## CARGO EL ALINEAMIENTO
    Alignment =paste0(OrthoPath,Sites[ORTH,1])
    cyanobacterias <- phangorn::read.phyDat(Alignment, format = "interleaved")
    cyanobacterias <- subset(cyanobacterias, site.pattern = FALSE) 
   # phangorn::parsimony(Phylogeny, cyanobacterias)
    
    ## CARGO EL MODELO
    fit <- phangorn::pml(Tree, cyanobacterias,model=EvolutionModel)
    #fit <- optim.pml(fit, model=EvolModel, control = pml.control(trace=0),)
    
    ## CARGO La NUEVA FILOGENIA SIN RAIZ
    #Tree2 <- fit[["tree"]]
    
    ## HAGO LA RECONSTRUCCIÓN
    if(Method == "bayes"){
      anc.bayes <- phangorn::ancestral.pml(fit, type="bayes", return = "prob")
      Reconstruction <- anc.bayes
    }
    if(Method == "ml"){
      anc.ml <- phangorn::ancestral.pml(fit, "ml")
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
    a.1 <- tidytree::as_tibble(Tree)
    ## Extraigo el Boostrap por si se ocupa
    NA.NODE = length(cyanobacterias)+1
    for (i in NA.NODE:length(Nodos0)){ ## ESTO DEPENDE DE LA ESTRUCTURA DEL ARBOL
      a.1$label[i] <- a.1$node[i] 
    }
    ## Agrego los datos de los palindromos al arbol
    b.1 <- tibble::tibble(label = Nodos0)
    c.1 <- tidytree::full_join(a.1, b.1, by = 'label')
    NodeLabels <- c.1$label
    
    #TipLabels <- Tree2$tip.label
    TipLabels <- Tree$tip.label
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
    write.table(FASTA,file=paste0(OrthoPath,Sites[ORTH,1],'.RECONSTRUCTION.fasta'),sep="\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    print(paste0("Sitio ",ORTH, " de ",length(Sites$FILE),"."))
  }
}
#---------------------------------------------------------------------------------------------------------------------------------------
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