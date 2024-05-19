setwd("E:/MSc Bioinformatics - MLU/4. Semester/Angewandte Bioinformatik/Übungen/übung04/anbio_ueb04_4")
library(rhdf5)
library(tximport)
library(readr)
library(DESeq2)
library(lattice)
library(UpSetR)
library(gridExtra)
library(ggVennDiagram)

# Aufgabe 4.4a) -----------------------------------------------------------

# Bestimmung der p-Werte/LFC mit results() 
a <- DESeq(dds); res <- results(a); res.df <- as.data.frame(res)

# BLS vs. Mock
BLS_Mock <- results(a,name=resultsNames(a)[2]); res.BLS <- as.data.frame(BLS_Mock)

# RS vs. Mock
RS_Mock <- results(a,name=resultsNames(a)[3]); res.RS <- as.data.frame(RS_Mock)


# Aufgabe 4.4b) -----------------------------------------------------------

# Bestimmung der IDs der DE-Gene nach p < 0.01 & L2FC > 4
ID_BLS <- rownames(res.BLS[!is.na(res.BLS$padj) & 
                             res.BLS$padj < 0.01 & 
                             res.BLS$log2FoldChange > 4,])

ID_RS <- rownames(res.RS[!is.na(res.RS$padj) & 
                           res.RS$padj < 0.01 & 
                           res.RS$log2FoldChange > 4,])

# norm: counts(dds,normalized=T) aus Aufgabe 4.3b). Falls "norm" nicht gefunden,
# bitte Skript anbio_ueb04_3 ausführen

# DEGs und Heatmap von BLS_Mock und RS_Mock
deg_bls <- norm[ID_BLS,]; heatmap(deg_bls)
deg_rs <- norm[ID_RS,]; heatmap(deg_rs)

# Aufgabe 4.4c) -----------------------------------------------------------

# Venn-Diagramme und Upset-plots der DE-Transkripte für BLS und RS

# BLS
degs.pval <- rownames(res.BLS[!is.na(res.BLS$padj) &res.BLS$padj < 0.01,])
degs.up <- rownames(res.BLS[!is.na(res.BLS$padj) & res.BLS$log2FoldChange > 4,])
degs.down <- rownames(res.BLS[!is.na(res.BLS$padj) & res.BLS$log2FoldChange < 4,])
degs.list <- list("pval < 0.01"=degs.pval, "up (lfc > 4)"=degs.up, "down (lfc < 4)"=degs.down)

png("upset_venn_BLS.png",1920,1080)
par(mfrow=c(1,2))
upset <- upset(fromList(degs.list))
venn <- ggVennDiagram(degs.list)
grid.arrange(venn,upset,ncol=2)
par(mfrow=c(1,1))
dev.off()
