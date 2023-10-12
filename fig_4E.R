##### load libraries ----

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(rstatix)
library(dunn.test)
library(ggtext)
library(xfun)
library(ggthemes)
library(tidyverse)
library(svglite)
library(multcompView)

##### get and prepare data ----

folder<-"Z:/EQUIPE/meiorec/PROJECTS/TransciptomeScreen/Immuno/STED/Cologne_Avril_Mesure_SCEP1"

setwd(folder)

df<-read.csv2("fig_4E_data.csv")
df$antibody<-as.factor(df$antibody)
df$cell<-as.factor(df$cell)

df$conc<-paste(df$antibody, df$cell)
df$conc<-as.character(df$conc)

df

##### draw graph ----
###### pub theme -----

theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(0.9)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size=10),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))

}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

theme_set(theme_pubclean())
theme_set(theme_Publication())

###### sort colors out ----

cp<-c(rep("grey",2),
      rep("magenta",2),
      rep("#7fc97f",3))

###### prepare stat values ----

k_test <- kruskal.test(distances ~ antibody, data=df)
k_test # says there are differences between the groups
d_test <- dunn.test(df$distances, df$antibody, method="bonferroni")  # gives details on these differences

# although we should not use a wilcox test here, I used it to annotate the final graph
# as it gives the same significance levels as the dunn test

stat.test <- df %>% wilcox_test(distances ~ antibody)
stat.test <- stat.test %>% add_xy_position(x = "antibody")
stat.test
print(stat.test, n=50)

###### draw plot ------

df_raw<-df
df$names<-df$antibody

# need to do this because the order of the antibodies in the graph
stat.test[1,13]<-3
stat.test[2,12]<-2
stat.test[2,13]<-3
stat.test[3,12]<-1
stat.test[3,13]<-2

v_plot<-ggplot(df, aes(x=factor(names, level=c('SCEP1', 'ZYP1 Cter', 'REC8')), y = distances))+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.01, dotsize=1.5, aes(fill=conc), color="grey42")+# ,color = "magenta", fill = "magenta") +
  xlab("")+
  scale_fill_manual(values=cp)+
  scale_colour_manual(values=cp)+
  theme(axis.text.x = element_text(angle = 40, vjust=1, hjust = 1))+
  theme(axis.text.x = ggtext::element_markdown())+
  ylab("Distance between antibody lines (μm)")+
  ylim(0,0.35)+
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult=1),
               geom = "pointrange", color = "black", size=0.5)+
  stat_pvalue_manual( stat.test, label = "p.adj.signif", tip.length = 0.01, size=6,
                      y.position = c(0.3,0.27,0.27), bracket.shorten = 0.05)

v_plot


###### save plot ----

png("fig_4E_paper1.png", units="cm", width=15, height=12.5, res=2000)
v_plot
dev.off()


###### plot with letters per cell -----

df$names<-df$conc

anova <- aov(distances ~ names, data = df)
tukey <- TukeyHSD(anova)
cld <- multcompLetters4(anova, tukey)
cld <- as.data.frame.list(cld$names)
cld$names<-row.names(cld)
df <- merge(df, cld, by.x = "names")

label_pos <- df %>%
  group_by(conc) %>%
  summarize(mean_dist = mean(distances))

# changed this manually because otherwise c was always the first letter and thats not so pretty
df["Letters"][df["conc"] == "SCEP1 cell 1"] <- "a"
df["Letters"][df["conc"] == "SCEP1 cell 2"] <- "a"
df["Letters"][df["conc"] == "ZYP1 Cter cell 1"] <- "b"
df["Letters"][df["conc"] == "ZYP1 Cter cell 2"] <- "bc"
df["Letters"][df["conc"] == "ZYP1 Cter cell 3"] <- "bc"
df["Letters"][df["conc"] == "REC8 cell 1"] <- "c"
df["Letters"][df["conc"] == "REC8 cell 2"] <- "bc"



v_plot<-ggplot(df, aes(x=factor(conc, level=c('SCEP1 cell 1', 'SCEP1 cell 2', 'ZYP1 Cter cell 1', 'ZYP1 Cter cell 2', 'ZYP1 Cter cell 3',
                                              'REC8 cell 1', 'REC8 cell 2')), y = distances))+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.01, dotsize=1.3, aes(fill=conc), color="grey42")+# ,color = "magenta", fill = "magenta") +
  xlab("")+
  scale_fill_manual(values=cp)+
  scale_colour_manual(values=cp)+
  theme(axis.text.x = element_text(angle = 30, vjust=1, hjust = 1))+
  theme(axis.text.x = ggtext::element_markdown())+
  ylab("Distance between REC8 axes (μm)")+
  geom_text(aes(label = Letters), y=0.3, size=6, family="Roboto", color="grey42") +
  ylim(0,0.35)+stat_summary(fun.data = "mean_sdl", fun.args = list(mult=1),
             geom = "pointrange", color = "black", size=0.5)

v_plot

###### plot with letters per antibody -----
df<-df_raw
df$names<-df$antibody

anova <- aov(distances ~ names, data = df)
tukey <- TukeyHSD(anova)
cld <- multcompLetters4(anova, tukey)
cld <- as.data.frame.list(cld$names)
cld$names<-row.names(cld)
df <- merge(df, cld, by.x = "names")

label_pos <- df %>%
  group_by(names) %>%
  summarize(mean_dist = mean(distances))

tukey_results <- aov(distances ~ names, data = df) %>% TukeyHSD()
pvals <- data.frame(
  comparison = names(tukey_results[[1]][, "p adj"]),
  pvalue = tukey_results[[1]][, "p adj"]
)

# changed this manually because otherwise c was always the first letter and thats not so pretty
df["Letters"][df["names"] == "SCEP1"] <- "a"
df["Letters"][df["names"] == "ZYP1 Cter"] <- "b"
df["Letters"][df["names"] == "REC8"] <- "c"
df$names<-as.factor(df$names)

v_plot<-ggplot(df, aes(x=factor(names, level=c('SCEP1', 'ZYP1 Cter', 'REC8')), y = distances))+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.01, dotsize=1.3, aes(fill=conc), color="grey42")+# ,color = "magenta", fill = "magenta") +
  xlab("")+
  scale_fill_manual(values=cp)+
  scale_colour_manual(values=cp)+
  theme(axis.text.x = element_text(angle = 30, vjust=1, hjust = 1))+
  theme(axis.text.x = ggtext::element_markdown())+
  ylab("Distance between respective antibody axes (μm)")+
  geom_text(aes(label = Letters), y=0.3, size=6, family="Roboto", color="grey42") +
  ylim(0,0.36)+stat_summary(fun.data = "mean_sdl", fun.args = list(mult=1),
                            geom = "pointrange", color = "black", size=0.5)+
  geom_text(data = pvals, aes(x = c(1.5,2,2.5), y = c(0.33, 0.36,0.33),
                              label = paste0("p = ", signif(pvalue, 2))))+
  geom_segment(aes(x = 1.05 , y = 0.32, xend = 1.95, yend = 0.32)) +
  geom_segment(aes(x = 1.05 , y = 0.35, xend = 2.95, yend = 0.35)) +
  geom_segment(aes(x = 2.05 , y = 0.32, xend = 2.95, yend = 0.32))

v_plot



###### save plot ----

png("fig_4E_paper1.png", units="cm", width=15, height=12.5, res=2000)
v_plot
dev.off()

ggsave("fig_4E_paper1.svg", plot = v_plot, device = "svg")

