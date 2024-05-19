setwd("E:/MSc Bioinformatics - MLU/4. Semester/Angewandte Bioinformatik/Übungen/übung04/anbio_ueb04_3")
library(rhdf5)
library(tximport)
library(readr)
library(DESeq2)

# Aufgabe 4.3a) -----------------------------------------------------------

# Alle Abundances einlesen und mit DEseq2 einen Datensatz anlegen
samples <- read.table("sample_table.tsv",header=T)
dir <- "E:/MSc Bioinformatics - MLU/4. Semester/Angewandte Bioinformatik/Übungen/übung04/anbio_ueb04_3/kallisto"
files <- file.path(dir,samples$experiment,"abundance.tsv")
names(files) <- samples$experiment
txi <- tximport(files,type="kallisto",txOut = T)
dds <- DESeqDataSetFromTximport(txi,
                                colData = samples,
                                design = ~condition)

# Aufgabe 4.3b) -----------------------------------------------------------

# Normierung der Zählwerte (norm/r-log/vst) und Basis-Stufe zu Mock
dds$condition <- relevel(dds$condition,"Mock")
n <- list()

# counts
a <- DESeq(dds)
norm <- counts(a,normalized=T); lognorm <- log(norm+1); n[[1]] <- lognorm

# R-log (RLD)
RLD <- rlog(dds); n[[2]] <- assay(RLD) 

# Variance stabilizing transformation (VST)
VST <- vst(dds); n[[3]] <- assay(VST)

# Aufgabe 4.3c) -----------------------------------------------------------

# Vergleich der Verteilung der normierten Zählwerte

# Normiert und nicht-logarithmiert
tbl <- as_tibble(norm, rownames="Gene") %>% 
  pivot_longer(names_to = "Experiment", values_to="Count",-Gene)
ggplot(tbl,aes(x=Count,color=Experiment)) +
  geom_density()

# Normiert und logarithmiert
ggplot(tbl,aes(x=log(Count+1),color=Experiment)) +
  geom_density()

## Vergleich verschiedener Normalisierungsfunktionen (counts/RLD/VST)
for(i in seq_along(n)){
  a <- n[[i]]; fun <- c("lognorm","RLD","VST")
  png(paste0("Hist_density_",fun[i],".png"),1920,1080)
  par(mfrow=c(3,3))
  for(j in 1:ncol(a)){
    hist(a[,j],xlab="Zählwerte",ylab="Dichte", probability = T, col=i+2,
         main=paste("Verteilung",fun[i],colnames(a)[j]), cex.main=2, cex.lab=1.7)
    lines(density(a[,j]),lwd=2,col="red")
  }
  par(mfrow=c(1,1))
  dev.off()
  
  png(paste0("Boxplot_",fun[i],".png"),1920,1080)
  par(mar=c(4.5,7,3,2))
  boxplot(a,horizontal=T,las=1, col=i+1, cex.main=2, cex.lab=1.7,
          main=paste("Boxplot",fun[i]))
  par(mar=c(5,4,2,2)+0.1)
  dev.off()
}


# Aufgabe 4.3d) -----------------------------------------------------------
for(i in seq_along(n)) {
  a <- n[[i]]; fun <- c("lognorm","RLD","VST")
  for(j in 1:ncol(a)){
    png(paste0("Pairs_",fun[i],".png"),1920,1080)
    pairs(a,upper.panel = NULL, cex.main=2, cex.lab=1.7,
          main=paste("Verteilung",fun[i],colnames(a)[j]))
    dev.off()
  }
}
