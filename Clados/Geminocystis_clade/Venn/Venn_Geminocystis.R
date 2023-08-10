library(dplyr)
Clado = "Geminocystis_clade"
Spps <- c(
#  "NIES-3708",
          "NIES-3757",
          "PCC_10605",
          "PCC_6605",
          "PCC_7437",
          "Z-M001",
#          "NIES-3709",
          "NIES-4102",
#          "PCC_6308",
          "PCC_7418"
#          "PCC_8305"
)
Spps2 <- c(
  #  "NIES-3708",
  "Stanieria sp NIES-3757",
  "Cyanobacterium aponinum PCC 10605",
  "Chamaesiphon minutus PCC 6605",
  "Stanieria cyanosphaera PCC 7437",
  "Euhalothece natronophila Z-M001",
  #          "NIES-3709",
  "Chondrocystis sp NIES-4102",
  #          "PCC 6308",
  "Halothece sp PCC 7418"
  #          "PCC 8305"
)

z <- list()
for(i in 1:length(Spps)){
  File <- read.table(paste0("/home/lalibelulalo/TESIS/Clados/",Clado,"/PALINDROMES/GCGATCGC/",Spps[i],"/Orthologues_Palindrome_sites.AllFrames.FIRST.txt"),
                     sep = "\t",
                     header = TRUE)
  File <- File%>%
    filter(Spp==Spps[i])
  
  File$MergedCoords <- paste(File$START,"-",File$END)
  File$MergedCoordsOrth <- paste(File$MergedCoords,"-",File$FILE)
  spp=Spps[i]
  Conjunto.spp = File[,10]
  z[[i]] <- Conjunto.spp
}

p <- ggVennDiagram::ggVennDiagram(z,label_alpha = 0.35, color = "black", lwd = 0.8, lty = 1, edge_lty = "solid", edge_size = 1,
                             category.names = Spps2,set_size = 4.5,label_size = 2.3) + 
  ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
  ggplot2::scale_fill_distiller(palette = "YlGnBu")+
  ggplot2::scale_color_brewer(palette = "Set1")+
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.2))+
  ggplot2::labs( fill = "Conteo" )

ggplot2::ggsave(p,file="/home/lalibelulalo/TESIS/Clados/Geminocystis_clade/Venn/Venn_Geminocystis-7.png",width=7, height=6, units="in", scale=1.5)
