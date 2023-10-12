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
library(multcompView)

##### get and prepare data ----

folder<-"Z:/EQUIPE/meiorec/PROJECTS/TransciptomeScreen/Immuno/STED/Cologne_Janvier/Results/comparison"

setwd(folder)

files<-list.files(folder, pattern = "\\.csv$")
files_raw<-files

full_df<-data.frame(file_ID=character(),
                      distance=double())
for (i in files){
  raw_data<-read.csv(i)
  full_df<-rbind(full_df,raw_data)
}

raw_data<-read.csv(files[1])
raw_data

full_df[c('Mutant', 'Image', 'Line')] <- str_split_fixed(full_df$file_ID, '_', 3)
full_df = subset(full_df, select = c(file_ID, Mutant, Image, Line, distances) )

full_df$conc<-paste(full_df$Mutant, full_df$Image)
full_df$conc<-as.character(full_df$conc)

full_df["conc"][full_df["conc"] == "scep1 1"] <- "*scep1-1* cell 1"
full_df["conc"][full_df["conc"] == "scep1 2"] <- "*scep1-1* cell 2"
full_df["conc"][full_df["conc"] == "scep1 3"] <- "*scep1-1* cell 3"
full_df["conc"][full_df["conc"] == "scep1 4"] <- "*scep1-1* cell 4"
full_df["conc"][full_df["conc"] == "scep2 1"] <- "*scep2-1* cell 1"
full_df["conc"][full_df["conc"] == "scep2 2"] <- "*scep2-1* cell 2"
full_df["conc"][full_df["conc"] == "scep2 3"] <- "*scep2-1* cell 3"
full_df["conc"][full_df["conc"] == "scep2 4"] <- "*scep2-1* cell 4"
full_df["conc"][full_df["conc"] == "wt 1"] <- "wt cell 1"
full_df["conc"][full_df["conc"] == "wt 2"] <- "wt cell 2"

full_df$conc<-as.factor(full_df$conc)

full_df
#write.csv2(full_df, paste(folder,"/distances_full_df.csv", sep=""))
tail(full_df, n=50)

##### draw plot -----

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
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size=15),
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


###### prepare stat values ----

k_test <- kruskal.test(distances ~ conc, data=full_df)
k_test # says there are differences between the groups
d_test <- dunn.test(full_df$distances, full_df$conc, method="bonferroni")  # gives details on these differences

# although we should not use a wilcox test here, I used it to annotate the final graph
# as it gives the same significance levels as the dunn test

stat.test <- full_df %>% wilcox_test(distances ~ conc)
stat.test <- stat.test %>% add_xy_position(x = "conc")
stat.test
print(stat.test, n=50)

###### sort colors out ----

cp<-c(rep("magenta",4),rep("deepskyblue2",4), rep("grey",2))

###### draw plot ------

v_plot<-ggplot(full_df, aes(x = conc, y = distances))+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.01, dotsize=2, aes(fill=conc), color="grey42")+# ,color = "magenta", fill = "magenta") +
  xlab("")+
  scale_fill_manual(values=cp)+
  scale_colour_manual(values=cp)+
  theme(axis.text.x = element_text(angle = 30, vjust=1, hjust = 1))+
  theme(axis.text.x = ggtext::element_markdown())+
  ylab("Distance between REC8 axes (μm)")+
  ylim(0,0.8)+
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult=1),
               geom = "pointrange", color = "black", size=0.5)+
  stat_pvalue_manual( stat.test, label = "p.adj.signif", tip.length = 0, size=6,
                      y.position = c(rep(-2,3), 0.79, rep(-2,3), 0.726, rep(-2,25),
                                     0.656, rep(-2,10), 0.59), bracket.shorten = 0.05)
v_plot

###### plot with letters ----

full_df$names<- str_replace(full_df$conc, "-", "_")

anova <- aov(distances ~ names, data = full_df)
tukey <- TukeyHSD(anova)
cld <- multcompLetters4(anova, tukey)
cld <- as.data.frame.list(cld$names)
cld$names<-row.names(cld)
full_df <- merge(full_df, cld, by.x = "names")

label_pos <- full_df %>%
  group_by(conc) %>%
  summarize(mean_dist = mean(distances))

v_plot<-ggplot(full_df, aes(x = conc, y = distances))+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.01, dotsize=2, aes(fill=conc), color="grey42")+# ,color = "magenta", fill = "magenta") +
  xlab("")+
  scale_fill_manual(values=cp)+
  scale_colour_manual(values=cp)+
  theme(axis.text.x = element_text(angle = 30, vjust=1, hjust = 1))+
  theme(axis.text.x = ggtext::element_markdown())+
  ylab("Distance between REC8 axes (μm)")+
  geom_text(aes(label = Letters), y=0.65, size=6, family="Roboto", color="grey42") +
  ylim(0,0.8)

v_plot

###### save plot ----

ggsave("fig_2C_paper1.svg", plot = v_plot, device = "svg")

png("per_image_distances.png", units="cm", width=15, height=12.5, res=2000)
v_plot
dev.off()

##### run additional tests -----

lvls<-levels(full_df$conc)

###### check normality of each image -----

for (i in lvls){
  sub<-subset(full_df,conc==i)
  print(i)
  test<-shapiro.test(sub$distances)
  print(test)
  if (test$p.value > 0.05){
    message(paste(i,"is normal"," "))
  }
  else {
    warning(paste(i, "is not normal", " "))
  }
}

#### Get mean and sd ----
full_df
means<-aggregate(full_df$distances, list(full_df$Mutant), FUN=mean)
sds<-aggregate(full_df$distances, list(full_df$Mutant), FUN=sd)
sum<-aggregate(full_df$distances, list(full_df$Mutant), FUN=sum)
effectifs<-sum$x/means$x


