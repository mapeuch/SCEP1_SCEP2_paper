####### Parameters ------
path<-"Z:/EQUIPE/meiorec/PROJECTS/TransciptomeScreen/DOUBLE HYBRIDE SCEPs/Results"

####### Libraries ------

library(ComplexHeatmap)
library(circlize)
library(rjson)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(EnvStats)
library(tinytex)
library(tibble)
library(gmodels)
library(reshape2)
library(RColorBrewer)
library(Matrix)
library(grid)
library(sjPlot)

####### Code ------

raw_df<-read.csv(paste(path,"Y2H_Raw_Results_paper1.csv", sep="/"), sep=";")
raw_df["LWH"][is.na(raw_df["LWH"])] <- 0
raw_df["LWHA"][is.na(raw_df["LWHA"])] <- 0
raw_df$test<- with(raw_df, paste(pGAD,pGBK, sep="_"))
raw_df

LWH_df<-aggregate(raw_df$LWH, by=list(test=raw_df$test), FUN=sum)
colnames(LWH_df)<-c("test","LWH")
LWHA_df<-aggregate(raw_df$LWHA, by=list(test=raw_df$test), FUN=sum)
colnames(LWHA_df)<-c("test","LWHA")

n_df<-raw_df %>%
  group_by(test) %>%
  summarise(Rep=n())

InterRatio_df<-data.frame(n_df$test, n_df$Rep, LWH_df$LWH, LWHA_df$LWHA)
colnames(InterRatio_df)<-c("test", "rep", "LWH", "LWHA")
InterRatio_df$InterRatio<-(InterRatio_df$LWHA+InterRatio_df$LWH*0.5)/InterRatio_df$rep

InterRatio_df[c('pGAD', 'pGBK')] <- str_split_fixed(InterRatio_df$test, '_', 2)

InterRatio_df = subset(InterRatio_df, select = c(pGAD, pGBK, InterRatio, rep) )
InterRatio_df

IR_wide_df=subset(InterRatio_df, select=-c(rep))
IR_wide_df<-spread(IR_wide_df, key = pGBK, value = InterRatio)

order <- c("pGAD", "SCEP1", "SCEP1 Nter",
           "SCEP1 Cter",	"SCEP1 Cter2",
           "SCEP1 Dctd",	"SCEP2",
           "SCEP2 Nter",	"SCEP2 Cter",
           "ZYP1 Nter","ZYP1 Cter",
           "ZIP4 Cter")
str_sort(colnames(IR_wide_df))
str_sort(order)
IR_wide_df=subset(IR_wide_df, select=order)
IR_wide_df<-IR_wide_df %>% arrange(match(pGAD, order))

IR_wide_df

IR_wide_df<-IR_wide_df %>% remove_rownames %>% column_to_rownames(var="pGAD")

IR_mat <- as.matrix(IR_wide_df)
matrix<-IR_mat

anno_df<-n_df
anno_df[c('pGAD', 'pGBK')] <- str_split_fixed(n_df$test, '_', 2)
anno_df = subset(anno_df, select = c(pGAD, pGBK, Rep) )

anno_wide_df<-spread(anno_df, key = pGBK, value = Rep)
anno_wide_df=subset(anno_wide_df, select=order)
anno_wide_df<-anno_wide_df %>% arrange(match(pGAD, order))
anno_wide_df<-anno_wide_df %>% remove_rownames %>% column_to_rownames(var="pGAD")
anno_wide_df[is.na(anno_wide_df)] <- 0

annotations<-as.matrix(anno_wide_df)

##### select parts of matrix ----

matrix_unchanged<-matrix
annotations_unchanged<-annotations


##### draw matrix -----

col_fun = colorRamp2(c(0.2,0.5,1), c("beige","orange", "red"))

colnames(matrix)<-paste("BD",colnames(matrix), sep="-")
rownames(matrix)<-paste("AD",rownames(matrix), sep="-")

map<-Heatmap(matrix,
             width = ncol(matrix)*unit(12, "mm"),
             height = nrow(matrix)*unit(12, "mm"),
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             show_heatmap_legend=TRUE,
             heatmap_legend_param = list(title = "InterRatio", color_bar = "continuous",
                                         labels_gp = gpar(fontsize = 14),
                                         title_gp=gpar(fontsize=14, fontface='bold')),
             column_names_gp=grid::gpar(fontsize=16),
             row_names_gp=grid::gpar(fontsize=16),
             col=col_fun,
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.text(sprintf("%.0f", annotations[i, j]), x, y, gp = gpar(fontsize = 14))
             }
)

draw(map, row_title = "                  Fusion to Activating Domain", row_title_gp = gpar(fontsize = 8),
     column_title = "Fusion to DNA-Binding Domain                    ", column_title_gp = gpar(fontsize = 8))


setwd("Z:/EQUIPE/meiorec/PROJECTS/TransciptomeScreen/DOUBLE HYBRIDE SCEPs/Results/Heatmaps")

png("Heatmap.png", units="cm", width=30, height=30, res=400)
draw(map, row_title = "                  Fusion to Activating Domain", row_title_gp = gpar(fontsize = 20),
     column_title = "Fusion to DNA-Binding Domain                  ", column_title_gp = gpar(fontsize = 20))
dev.off()

png("Heatmap.png", units="cm", width=30, height=30, res=400)
draw(map)
dev.off()
