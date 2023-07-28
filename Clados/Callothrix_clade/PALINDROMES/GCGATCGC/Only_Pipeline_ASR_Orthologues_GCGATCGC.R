setwd("/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/Callothrix_clade/PALINDROMES/GCGATCGC/")

## Filtrado del conteo
## Para empezar escojo la especie con mayor conteo de palindromos segun la filogenia anotada (**PCC_7716**).
## Sin embargo tambien me aseguro de que dicha especie sea la que tenga el número de sitios:
Counts <- read.table("../../Markov_count_GCGATCGC_2023-5-7_18hrs43mins_Octanuc_.txt", sep = "\t", header = TRUE)

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


# GRAFO

## primero cargo las librerias necesarias

#library(ggplot2)
#library(ggtree)
#library(ape)
#library(tidyverse)
#library(tidytree)
#library(phangorn)
library(dplyr)


## Preparacion de datos para el grafo

## Primero cargo Orthologues_Palindrome_sites.txt del paso anterior.
Sites <- read.table("Orthologues_Palindrome_sites.txt", sep = "\t", header = TRUE)

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
EvolModel = "F81"

Method = "bayes"
Tree = ggtree::read.tree("../../SpeciesTree_rooted.txt")

for (ORTH in 1:length(Sites[,1])){ ## HAY 2593 ortólogos
  Alignment = paste0("Only_ORTHOLOGUES/",Sites[ORTH,1])
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


## Filtrado de datos para el GRAFO
## Primero cuento el numero de veces que se repite cada linea de la matriz de transiciones.
## El conteo para cada transicion sera el peso de cada vertice.
## Luego cambio el nombre de la columna *
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

## Creo un grafo dirigido agregando un peso de acuerdo a la columna **weight**. 
routes_tidy <- tidygraph::tbl_graph(nodes = nodes,
                                    edges = edges,
                                    directed = TRUE)

routes_tidy %>% 
  tidygraph::activate(edges) %>% 
  tidygraph::arrange(desc(weight))

## Visualizo el grafo:
ggraph::ggraph(routes_tidy, layout = "graphopt") + 
  ggraph::geom_node_point() +
  ggraph::geom_edge_link(ggplot2::aes(width = weight), alpha = 0.8) + 
  ggraph::scale_edge_width(range = c(0.2, 2)) +
  ggraph::geom_node_text(ggplot2::aes(label = label), repel = TRUE) +
  ggplot2::labs(edge_width = "Times") +
  ggraph::theme_graph()

## GRAFOS INTERACTIVOS
## Este es un grafo sin peso en los vertices.
visNetwork::visNetwork(nodes, edges)

## Este es el mismo grafo pero con peso en sus vertices
edges <- mutate(edges, width = weight/5 + 1)
visNetwork::visNetwork(nodes, edges) %>% 
  visNetwork::visIgraphLayout(layout = "layout_with_fr") %>% 
  visNetwork::visEdges(arrows = "middle")

## Este es otro grafo de fuerzas con peso en sus vertices y te muestra las conceciones de cada nodo de manera mas visual
nodes_d3 <- mutate(nodes, id = id - 1)
edges_d3 <- mutate(edges, from = from - 1, to = to - 1)

networkD3::forceNetwork(Links = edges_d3, Nodes = nodes_d3,
             Source = "from", Target = "to", 
             NodeID = "label", Group = "id", Value = "weight", 
             opacity = 1, fontSize = 16, zoom = TRUE)  

## Este ultimo grafo muestra el grafo anterior pero de una forma mas analizable.
networkD3::sankeyNetwork(Links = edges_d3, Nodes = nodes_d3, Source = "from", Target = "to", 
              NodeID = "label", Value = "weight", fontSize = 16, unit = "TIME(s)")




#---------------------------------------------------------------
# GRAFOS POR DIRECCION

## Las siguientes figuras corresponden a cada una de las direcciones.
WGHT = 30

## 336-3--7
ggtree::ggtree(Tree2) +
  ggtree::geom_tiplab(color='firebrick', offset = .10, hjust=3.0)+
  ggtree::geom_label(ggplot2::aes(label = node),show.legend = TRUE)

DIRECTION = "9--8"

RowCtsL <- data.table::setDT(LINKS)[,list(Count=as.numeric(.N)),names(LINKS)] ## CUENTO EL NUMERO DE OCURRENCIAS DE CADA TRANSICIÓN #11520 ***
colnames(RowCtsL) <- c("from","to","direction","weight") ## CAMBIO LA ETIQUETA "counts" A "weight" ***
RowCtsL = transform(RowCtsL, weight = as.integer(weight)) ## CONVIERTO A ENTEROS LA COLUMNA "weight" ***
RowCtsL <-tidytree::as_tibble(RowCtsL) ## CONVIERTO A TIBBLE**
length(RowCtsL$from)
RowCtsL<-RowCtsL%>% ##QUITO LOS LOOPS. ES DECIR QUITO AQUELLAS TRANSICIONES QUE VAYAN ASI MISMA ***
  filter(from!=to)
RowCtsL<-RowCtsL%>% ## QUITO AQUELLAS TRANSICIONES CON UN PESO MENOR A 2 ***
  filter(weight>=WGHT)
RowCtsL<-RowCtsL%>% ## filtro las direcciones 9--8
  filter(direction==DIRECTION)

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

### ***  
routes_tidyL <- tidygraph::tbl_graph(nodes = nodesL,
                          edges = edgesL,
                          directed = TRUE)

routes_tidyL %>%
  tidygraph::activate(edges) %>% 
  tidygraph::arrange(desc(weight))

edgesL <- mutate(edgesL, width = weight/5 + 1)
visNetwork::visNetwork(nodesL, edgesL) %>% 
  visNetwork::visIgraphLayout(layout = "layout_with_fr") %>% 
  visNetwork::visEdges(arrows = "middle")

nodes_d3L <- mutate(nodesL, id = id - 1)
edges_d3L <- mutate(edgesL, from = from - 1, to = to - 1)

networkD3::sankeyNetwork(Links = edges_d3L, Nodes = nodes_d3L, Source = "from", Target = "to", 
              NodeID = "label", Value = "weight", fontSize = 16, unit = "TIME(s)")

## extraigo las secuencias del nodo 9-->8
sequences<-as.vector(nodes_d3L[,2])
write.table(sequences,file="sequences.txt",sep="\t",row.names = FALSE,col.names = FALSE,quote = FALSE)

## Cambio los saltos de linea por comas
## Este archivo lo voy a usar para buscar dichas secuencias en los ortólogos