library(dplyr)
Clado = "Pseudoanabaena_clade"
interval = 1:4
Spps <- c("PCC_7336",
          "ABRG5-3",
          "PCC_7502",
          "PCC_7367",
          "JA-3-3Ab",
          "JA-2-3Ba"
)

Spps2 <- c("Synechococcus sp PCC_7336",
           "Pseudanabaena sp ABRG5-3",
           "Synechococcus sp PCC_7502",
           "Pseudanabaena sp PCC_7367",
           "Synechococcus sp JA-3-3Ab",
           "Synechococcus sp JA-2-3Ba"
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

ggplot2::ggsave(p,file=paste0("/home/lalibelulalo/TESIS/Clados/",Clado,"/Venn/Venn_",Clado,"-",length(Spps2),".png"),width=7, height=6, units="in", scale=1.5)

p2 <- ggVennDiagram::ggVennDiagram(z[interval],label_alpha = 0.35, color = "black", lwd = 0.8, lty = 1, edge_lty = "solid", edge_size = 1,
                                   category.names = Spps2[interval],set_size = 4.5,label_size = 2.3) + 
  ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
  ggplot2::scale_fill_distiller(palette = "YlGnBu")+
  ggplot2::scale_color_brewer(palette = "Set1")+
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.2))+
  ggplot2::labs( fill = "Conteo" )

ggplot2::ggsave(p2,file=paste0("/home/lalibelulalo/TESIS/Clados/",Clado,"/Venn/Venn_",Clado,"-",length(Spps2),"RELEVANT.png"),width=7, height=6, units="in", scale=1.5)
