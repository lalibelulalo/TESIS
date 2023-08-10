library(dplyr)
Clado = "Pseudoanabaena_clade"
Spps <- c("ABRG5-3",
          "JA-2-3Ba",
          "JA-3-3Ab",
          "PCC_7336",
          "PCC_7367",
          "PCC_7502"
)

Spps2 <- c("Pseudanabaena sp ABRG5-3",
           "Synechococcus sp JA-2-3Ba",
           "Synechococcus sp JA-3-3Ab",
           "Synechococcus sp PCC_7336",
           "Pseudanabaena sp PCC_7367",
           "Synechococcus sp PCC_7502"
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
