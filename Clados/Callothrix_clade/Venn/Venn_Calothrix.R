library(dplyr)
Clado = "Callothrix_clade"
interval = 1:3
Spps <- c("336-3",
          "NIES-3974",
          "PCC_6303",
          "NIES-4105",
          "NIES-4071",
          "PCC_7716")
Spps2 <- c("Calothrix 336-3",
           "Calothrix NIES-3974",
           "Calothrix PCC_6303",
           "Calothrix NIES-4105",
           "Calothrix NIES-4071",
           "Calothrix PCC_7716")

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

