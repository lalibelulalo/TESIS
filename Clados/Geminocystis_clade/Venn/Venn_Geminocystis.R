library(dplyr)
Clado = "Geminocystis_clade"
interval = 1:4
Spps <- c("NIES-4102",
          "NIES-3757",
          "PCC_7437",
          "PCC_7418",
          "PCC_10605",
          "Z-M001",
          "PCC_6605"
)
Spps2 <- c("Chondrocystis sp NIES-4102",
           "Stanieria sp NIES-3757",
           "Stanieria cyanosphaera PCC 7437",
           "Halothece sp PCC 7418",
           "Cyanobacterium aponinum PCC 10605",
           "Euhalothece natronophila Z-M001",
           "Chamaesiphon minutus PCC 6605"
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

p2 <- ggVennDiagram::ggVennDiagram(z[interval],label_alpha = 0.35, color = "black", lwd = 0.8, lty = 1, edge_lty = "solid", edge_size = 1,
                                   category.names = Spps2[interval],set_size = 4.5,label_size = 2.3) + 
  ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
  ggplot2::scale_fill_distiller(palette = "YlGnBu")+
  ggplot2::scale_color_brewer(palette = "Set1")+
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.2))+
  ggplot2::labs( fill = "Conteo" )

ggplot2::ggsave(p2,file=paste0("/home/lalibelulalo/TESIS/Clados/",Clado,"/Venn/Venn_",Clado,"-",length(Spps2),"RELEVANT.png"),width=7, height=6, units="in", scale=1.5)
