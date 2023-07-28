setwd("/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/Callothrix_clade/PALINDROMES/GCGATCGC/")

source("/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/ASR_Orth_Functions/NodeAndEdges.R")

Nodes.Edges <- Create_Transition_Table (SitesTable = "Orthologues_Palindrome_sites.txt",
                                EvolutionModel = "F81",
                                Method = "bayes",
                                Phylogeny = "../../SpeciesTree_rooted.txt",
                                OrthoPath = "Only_ORTHOLOGUES/")


Nodes.Edges <- Create_Transition_Table (SitesTable = "/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/Callothrix_clade/PALINDROMES/GCGATCGC/Orthologues_Palindrome_sites.txt",
                                EvolutionModel = "F81",
                                Method = "bayes",
                                Phylogeny = "/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/Callothrix_clade/SpeciesTree_rooted.txt",
                                OrthoPath = "/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/Callothrix_clade/PALINDROMES/GCGATCGC/Only_ORTHOLOGUES/")

Transition_table <- read.table("/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/Callothrix_clade/PALINDROMES/GCGATCGC/TRANSICIONES.txt", sep = "\t", header = TRUE )
Network <- Create_Network(Transitions = Transition_table)
nodes = Network[[1]]
edges = Network[[2]]

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

#----------------------------------------------------------------------
Directed_Transition_table <- read.table("/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/Callothrix_clade/PALINDROMES/GCGATCGC/TRANSICIONES_DIRECCION.txt", sep = "\t", header = TRUE )
Directed_Network <- Create_Directed_Network(DirectedTransitions = Directed_Transition_table,
                                             Direction = "9--8",
                                             Weight = 30)
nodesL = Directed_Network[[1]]
edgesL = Directed_Network[[2]]

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

table.Transitions <- as.data.frame(edges)
table.nodes <- as.data.frame(nodes)
for (i in 1:length(table.nodes[,1])){
  table.Transitions$from[table.Transitions$from == table.nodes[i,1]] <- table.nodes[i,2]
  table.Transitions$to[table.Transitions$to == table.nodes[i,1]] <- table.nodes[i,2]
}

table.Transitions


Tree = ggtree::read.tree("/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/Callothrix_clade/SpeciesTree_rooted.txt")
ggtree::ggtree(Tree,branch.length='none') +
  ggtree::geom_tiplab(color='firebrick', offset = .14)+
  ggtree::geom_label(ggplot2::aes(label = node),show.legend = FALSE)+
  ggtree::geom_hilight(node=10, fill="steelblue", alpha=.6) +
  ggtree::geom_hilight(node=c(3,2,1), fill="darkgreen", alpha=.6)
  
  
ggtree::ggtree(Tree,branch.length='none')+
  ggtree::geom_label(ggplot2::aes(label = node),show.legend = TRUE)+
  ggtree::geom_cladelab(node=4, label="PCC_7716", align=TRUE, 
                  geom='label', fill='lightblue')+
  ggtree::geom_cladelab(node=5, label="NIES-4071", align=TRUE, 
                        geom='label', fill='lightblue')+
  ggtree::geom_cladelab(node=6, label="NIES-4105", align=TRUE, 
                        geom='label', fill='lightblue')+
  ggtree::geom_cladelab(node=1, label="336-3", align=TRUE, 
                        geom='label', fill='red')+
ggtree::geom_cladelab(node=2, label="PCC_6303", align=TRUE, 
                      geom='label', fill='red')+
  ggtree::geom_cladelab(node=3, label="NIES-3974", align=TRUE, 
                        geom='label', fill='red')


tree2 <- ggtree::groupClade(Tree, c(7,10))
ggtree::ggtree(tree2, ggplot2::aes(color=group),branch.length='none') + 
  ggplot2::theme(legend.position='none') +
  ggtree::scale_color_manual(values=c("red","steelblue"))+
  ggtree::geom_label(ggplot2::aes(label = node),show.legend = TRUE)+
  ggtree::geom_cladelab(node=4, label="PCC_7716", align=TRUE, 
                        geom='label', fill='steelblue')+
  ggtree::geom_cladelab(node=5, label="NIES-4071", align=TRUE, 
                        geom='label', fill='steelblue')+
  ggtree::geom_cladelab(node=6, label="NIES-4105", align=TRUE, 
                        geom='label', fill='steelblue')+
  ggtree::geom_cladelab(node=1, label="336-3", align=TRUE, 
                        geom='label', fill='red')+
  ggtree::geom_cladelab(node=2, label="PCC_6303", align=TRUE, 
                        geom='label', fill='red')+
  ggtree::geom_cladelab(node=3, label="NIES-3974", align=TRUE, 
                        geom='label', fill='red')
  



