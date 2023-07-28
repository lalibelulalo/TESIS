source("/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/ASR_Orth_Functions/NodeAndEdges.R")


PAL_PATH <- '/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/Thermosynechococcus_clade/PALINDROMES/GCGATCGC/'
DIRS <- list.files(path=PAL_PATH, pattern=NULL, all.files=FALSE,
                   full.names=FALSE)
#DIRS[2]
#i=2
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



setwd('/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/Thermosynechococcus_clade/PALINDROMES/GCGATCGC/SCTE542/')
Create_Reconstruction_Files(SitesTable = "Orthologues_Palindrome_sites.AllFrames.FIRST.txt",
                            EvolutionModel = "F81",
                            Method ="bayes",
                            Phylogeny = "../../../SpeciesTree_rooted.txt",
                            OrthoPath = "Only_ORTHOLOGUES/",
                            TAXON = 'PCC_6312')


library(kableExtra)
library(formattable)
#Sys.setenv("OPENSSL_CONF"="/dev/null")

tabla_t <- read.table('/home/lalibelulalo/PIPELINES_2023/CODON_MUTATIONS/BLAST_clado_calothrix.csv',sep = "\t",header = TRUE)
tabla_t <- tabla_t[,c(1,2,3,5,6,8,11,13,14)]
tabla_t <- tabla_t%>%
  filter(Clade!='Calothrix')

tabla_t$RepeatType = cell_spec(tabla_t$RepeatType, color = ifelse(tabla_t$RepeatType == 'TRP', "red", "blue"))

tabla_t[4] <- lapply(tabla_t[4], function(x) {
  cell_spec(x, bold = T, 
            color = spec_color(x, end = 0.9),
            font_size = spec_font_size(x))
})
#tabla_t$RepeatLength <- color_bar("lightgreen")(tabla_t$RepeatLength)

tabla_t$GC <- color_tile("white", "green")(tabla_t$GC)
tabla_t$RepeatFreq <- color_tile("white", "orange")(tabla_t$RepeatFreq)

kbl(tabla_t, align = "c",caption = "Conjuntos de repeticiones con palíndromos hallados en los 6 clados. ",longtable = TRUE,format = "html", escape = F) %>%
  kable_paper(full_width = T) %>%
  column_spec(1, bold = T) %>%
  #row_spec(c(1,3,4,8,10,20,38,51,52), bold = T, color = "white", background = "yellow")%>%
  collapse_rows(columns = 1:2, valign = "top")%>%
  kable_styling(c("striped", "condensed"), 
                latex_options = "striped", 
                full_width = F,
                fixed_thead = T,
                font_size = 12)%>%
  scroll_box(width = "100%", box_css = "border: 0px;")%>%
  save_kable("NOcallothrix_repeats_table.png")

library(dplyr)
tabla2 = read.table(file="/home/lalibelulalo/PIPELINES_2023/CLADO_A18-40/Orthologues_Palindrome_sites.txt", header=TRUE, sep="\t")
tabla2 = tabla2[ , c(1,4,5,6)]
tabla2<-tabla2%>%
  filter(SiteNumber>=7)

knitr::kable(
  head(tabla2[, 1:4], 15), booktabs = TRUE,
  caption = "Ubicación de los sitios CGTTACG. La tabla muestra 15 sitios partiendo del sitio 7 donde comienzan las repeticiones. La primera columna se muestra el numero de sitio. La segunda columna muestra el intervalo en el que se encuentra el palíndromo. La tercera columna muestra a cuantos nucléotidos se encuentra el ultimo sitio. La cuarta columna indica la diferencia entre la distancia del ultimo palindromo y la distancia del siguiente.",
  digits = 2,align = 'r')%>%
  kable_paper(full_width = T) %>%
  column_spec(1, bold = T) %>%
  kable_styling(c("striped", "condensed"), 
                latex_options = "striped", 
                full_width = F,
                fixed_thead = T,
                font_size = 25)%>%
  save_kable("A18-40_repeats_positions.png")


library(dplyr)

tabla = read.table(file="/home/lalibelulalo/R_2/FIGURES_pico/Markov_count_pico_2022_gbff_2022-10-28_10hrs29mins_Octanuc_M3_.txt", header=TRUE, sep="\t")
tabla <- tabla%>%
  filter(spp=='Parasynechococcus_marenigrum_WH_8102'|spp=='Synechococcus_sp_WH_8103'|spp=='Synechococcus_sp_A18-46_1'|spp=='Synechococcus_sp_BOUM118'|spp=='Synechococcus_sp_RS9915'|spp=='Synechococcus_sp_A18-40'|spp=='Synechococcus_sp_A15-24'|spp=='Synechococcus_sp_A15-28')
tabla <- tabla%>%
  filter(palindrome=='ATGCGCAT')#CGTTAACG,ATGCGCAT
tabla = tabla[ , c(1,3,4,5)]

tabla$obs <- color_tile("white", "orange")(tabla$obs)
tabla$spp = cell_spec(tabla$spp, color = ifelse(tabla$spp == 'Synechococcus_sp_A18-40', "red", "black"))

kbl(tabla, align = "c",caption = "Conteo del palíndromo ATGCGCAT en el clado A18-40.",longtable = TRUE,format = "html", escape = F) %>%
  kable_paper(full_width = T) %>%
  column_spec(1, bold = T) %>%
  #row_spec(c(1,3,4,8,10,20,38,51,52), bold = T, color = "white", background = "yellow")%>%
  #collapse_rows(columns = 1:2, valign = "top")%>%
  kable_styling(c("striped", "condensed"), 
                latex_options = "striped", 
                full_width = F,
                fixed_thead = T,
                font_size = 25)%>%
  #scroll_box(width = "100%", box_css = "border: 0px;")%>%
  save_kable("A18-40_Counts.png")



#---------------------------------
source("/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/ASR_Orth_Functions/NodeAndEdges.R")

Transition_table <- read.table("/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/Callothrix_clade/PALINDROMES/GCGATCGC/TRANSICIONES.txt", sep = "\t", header = TRUE )
Network <- Create_Network(Transitions = Transition_table)
nodes = Network[[1]]
edges = Network[[2]]

routes_tidy <- tidygraph::tbl_graph(nodes = nodes,
                         edges = edges,
                         directed = TRUE)

GRAPH <- ggraph::ggraph(routes_tidy, layout = "graphopt") + 
  ggraph::geom_node_point() +
  ggraph::geom_edge_link(ggplot2::aes(width = weight), alpha = 0.8) + 
  ggraph::scale_edge_width(range = c(0.2, 2)) +
  ggraph::geom_node_text(ggplot2::aes(label = label), repel = TRUE) +
  ggplot2::labs(edge_width = "Times") +
  ggraph::theme_graph()

ggplot2::ggsave(GRAPH, file='Calothrix_net.png',width=7, height=6, units="in", scale=1.5)


Directed_Transition_table <- read.table("/home/lalibelulalo/PIPELINES_2023/ASR_ORTHOLOGUES/Callothrix_clade/PALINDROMES/GCGATCGC/TRANSICIONES_DIRECCION.txt", sep = "\t", header = TRUE )
Directed_Network <- Create_Directed_Network(DirectedTransitions = Directed_Transition_table,
                                            Direction = "9--8",
                                            Weight = 30)
nodesL = Directed_Network[[1]]
edgesL = Directed_Network[[2]]

routes_tidy <- tidygraph::tbl_graph(nodes = nodesL,
                                    edges = edgesL,
                                    directed = TRUE)

GRAPH <- ggraph::ggraph(routes_tidy, layout = "graphopt") + 
  ggraph::geom_node_point() +
  ggraph::geom_edge_link(ggplot2::aes(width = weight), alpha = 0.8) + 
  ggraph::scale_edge_width(range = c(0.2, 2)) +
  ggraph::geom_node_text(ggplot2::aes(label = label), repel = TRUE) +
  ggplot2::labs(edge_width = "Times") +
  ggraph::theme_graph()

ggplot2::ggsave(GRAPH, file='Calothrix_net_directed.png',width=7, height=6, units="in", scale=1.5)

