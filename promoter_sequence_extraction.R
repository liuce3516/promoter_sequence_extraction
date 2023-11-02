#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("Biostrings")

library(Biostrings)
library(tidyverse)

load("/Users/liuce/Desktop/R/Broccoli.genome.seq.RData")#Bo.seq
load("/Users/liuce/Desktop/R/HDEM.gene.info.RData")#HDEM.info

find.Bol.pro.seq <- function(gene.id,pro.length=2000,whether_gene.seq=TRUE){
  gene.id.i <- gene.id
  chr.i <- str_match(gene.id.i,"(?<=Bol).+(?=t)")[1,1]
  strand.i <- HDEM.info$mRNA %>% filter(gene.id==gene.id.i) %>% .$strand
  
  if(strand.i=="+"){
    mRNA.i <- HDEM.info$mRNA %>% filter(gene.id==gene.id.i)
    #promoter position
    Pro.start.i <- mRNA.i$start[1]-pro.length
    Pro.end.i <- mRNA.i$start[1]-1
    Pro.seq.i <- Bo.seq[[chr.i]][Pro.start.i:Pro.end.i]
    gene.seq.i <- Bo.seq[[chr.i]][mRNA.i$start:mRNA.i$end]#DNA
    
    if(whether_gene.seq==TRUE){
      filename <- str_c("Pro.",gene.id.i,".fa",sep = "")
      sink(filename)
      cat(paste0(">Pro.",gene.id.i,"\n"))
      cat(paste0(Pro.seq.i,"\n"))
      cat(paste0(">Gene.",gene.id.i,"\n"))
      cat(paste0(gene.seq.i,"\n"))
      sink()
    }else if(whether_gene.seq==FALSE){
      filename <- str_c("Pro.",gene.id.i,".fa",sep = "")
      sink(filename)
      cat(paste0(">Pro.",gene.id.i,"\n"))
      cat(paste0(Pro.seq.i,"\n"))
      sink()
    }
    
  }else if(strand.i=="-"){
    mRNA.i <- HDEM.info$mRNA %>% filter(gene.id==gene.id.i)
    #promoter position
    Pro.start.i <- mRNA.i$end[1]+1
    Pro.end.i <- mRNA.i$end[1]+pro.length
    Pro.seq.i <- reverseComplement(Bo.seq[[chr.i]][Pro.start.i:Pro.end.i])
    gene.seq.i <- reverseComplement(Bo.seq[[chr.i]][mRNA.i$start:mRNA.i$end])#DNA
    if(whether_gene.seq==TRUE){
      filename <- str_c("Pro.",gene.id.i,".fa",sep = "")
      sink(filename)
      cat(paste0(">Pro.",gene.id.i,"\n"))
      cat(paste0(Pro.seq.i,"\n"))
      cat(paste0(">Gene.",gene.id.i,"\n"))
      cat(paste0(gene.seq.i,"\n"))
      sink()
    }else if(whether_gene.seq==FALSE){
      filename <- str_c("Pro.",gene.id.i,".fa",sep = "")
      sink(filename)
      cat(paste0(">Pro.",gene.id.i,"\n"))
      cat(paste0(Pro.seq.i,"\n"))
      sink()
    }
  }
  return(as.character(Pro.seq.i))
  print(filename,"file has been generated,\n which is located in folder",getwd())
}

#example
#Pro.seq <- find.Bol.pro.seq("BolC1t01574H")
#Pro.seq

#Pro.seq <- find.Bol.pro.seq("BolC4t27245H")
