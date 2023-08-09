File1 <- read.table("Clados/Callothrix_clade/PALINDROMES/GCGATCGC/336-3/Orthologues_Palindrome_sites.AllFrames.FIRST.txt", sep = "\t", header = TRUE)
File1 <- File1%>%
  filter(Spp=='336-3')

File2 <- read.table("Clados/Callothrix_clade/PALINDROMES/GCGATCGC/NIES-3974/Orthologues_Palindrome_sites.AllFrames.FIRST.txt", sep = "\t", header = TRUE)
File2 <- File2%>%
  filter(Spp=='NIES-3974')

File3 <- read.table("Clados/Callothrix_clade/PALINDROMES/GCGATCGC/PCC_6303/Orthologues_Palindrome_sites.AllFrames.FIRST.txt", sep = "\t", header = TRUE)
File3 <- File3%>%
  filter(Spp=='PCC_6303')

File1$MergedCoords <- paste(File1$START,"-",File1$END)
File1$MergedCoordsOrth <- paste(File1$MergedCoords,"-",File1$FILE)

File2$MergedCoords <- paste(File2$START,"-",File2$END)
File2$MergedCoordsOrth <- paste(File2$MergedCoords,"-",File2$FILE)

File3$MergedCoords <- paste(File3$START,"-",File3$END)
File3$MergedCoordsOrth <- paste(File3$MergedCoords,"-",File3$FILE)

length(paste(File1$START,"-",File1$END)) ## 2591 coordenadas de sitios
length(unique(paste(File1$START,"-",File1$END))) ## 1209  coordenadas de sitios UNICAS

length(File1[,10]) ## 2591
length(unique(File1[,10])) ## 2591

y <- list(
  Calothrix_336.3 = File1[,10], 
  Calothrix_NIES.3974 = File2[,10], 
  Calothrix_PCC_6303 = File3[,10]
)

ggVennDiagram::ggVennDiagram(y,label_alpha = 0.5, color = "black", lwd = 0.8, lty = 1, edge_lty = "dashed", edge_size = 1,
                             category.names = c("Calothrix 336-3","Calothrix NIES-3974","Calothrix PCC_6303"),set_size = 3.5) + 
  ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
  ggplot2::scale_fill_distiller(palette = "YlGnBu")+
  ggplot2::scale_color_brewer(palette = "Greys")+
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.2))+
  ggplot2::labs( fill = "Conteo" )

#--------------------------------------
File1 <- read.table("Clados/Cyanobacterium_clade/PALINDROMES/GCGATCGC/ATCC_51142/Orthologues_Palindrome_sites.AllFrames.FIRST.txt", sep = "\t", header = TRUE)
File1 <- File1%>%
  filter(Spp=='ATCC_51142')

File2 <- read.table("Clados/Cyanobacterium_clade/PALINDROMES/GCGATCGC/PCC_8801/Orthologues_Palindrome_sites.AllFrames.FIRST.txt", sep = "\t", header = TRUE)
File2 <- File2%>%
  filter(Spp=='PCC_8801')

File1$MergedCoords <- paste(File1$START,"-",File1$END)
File1$MergedCoordsOrth <- paste(File1$MergedCoords,"-",File1$FILE)

File2$MergedCoords <- paste(File2$START,"-",File2$END)
File2$MergedCoordsOrth <- paste(File2$MergedCoords,"-",File2$FILE)

y <- list(
  A = File1[,10], 
  B = File2[,10]
)

ggVennDiagram::ggVennDiagram(y,label_alpha = 0.5, color = "black", lwd = 0.8, lty = 1, edge_lty = "dashed", edge_size = 1,
                             category.names = c("ATCC_51142","PCC_8801"),set_size = 4.5) + 
  ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
  ggplot2::scale_fill_distiller(palette = "YlGnBu")+
  ggplot2::scale_color_brewer(palette = "Greys")+
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.2))+
  ggplot2::labs( fill = "Conteo" )

