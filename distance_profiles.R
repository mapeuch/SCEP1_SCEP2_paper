
folder<-"Z:/EQUIPE/meiorec/PROJECTS/TransciptomeScreen/Immuno/STED/Cologne_Janvier/scep1/SD23010-scep1-o1i1-ASY1-HEI10-REC8_data"

image_name="scep1_1"

library(ggplot2)
library(ggridges)

setwd(folder)

files<-list.files(folder, pattern = "\\.csv$")
files_raw<-files

read.csv(files[1])

full_df_rec_8<-data.frame(Distance=double(),
                    Gray_Value=double())

full_df_scep1<-data.frame(Distance=double(),
                          Gray_Value=double())

# add more conditions if other proteins than rec8 and scep1
for (i in files){
  firstCharacter = substr(i ,1,4)
  raw_data<-read.csv(i)
  if (firstCharacter == "rec8"){
    full_df_rec_8<-rbind(full_df_rec_8,raw_data)
  }
    else {
      full_df_scep1<-rbind(full_df_scep1,raw_data)
    }
}

full_df_rec_8$protein<-"rec8"
full_df_scep1$protein<-"scep1"

max_rec8<-max(full_df_rec_8$Gray_Value)
max_scep1<-max(full_df_scep1$Gray_Value)

full_df_rec_8$normalized_gray_value<-full_df_rec_8$Gray_Value/max_rec8
full_df_scep1$normalized_gray_value<-full_df_scep1$Gray_Value/max_scep1

full_df<-rbind(full_df_rec_8,full_df_scep1)

plot(full_df$Distance_.microns.,full_df$normalized_gray_value)

e <- ggplot(full_df, aes(x = Distance_.microns., y = Gray_Value))+
  xlab("gray value")+
  xlab("Distance (μm)")

e+geom_point(aes(col=protein))

p <- ggplot(full_df, aes(x=Distance_.microns.)) +
  geom_density(aes(col=protein))
p

x1<-full_df_rec_8$Distance_.microns.
y1<-full_df_rec_8$normalized_gray_value
lo1 <- loess(y1~x1, span=0.4)
xl1 <- seq(min(x1),max(x1), (max(x1) - min(x1))/1000)

x2<-full_df_scep1$Distance_.microns.
y2<-full_df_scep1$normalized_gray_value
lo2 <- loess(y2~x2, span=0.4)
xl2 <- seq(min(x2),max(x2), (max(x2) - min(x2))/1000)

color_pallete_function <- colorRampPalette(
  colors = c("green3", "magenta"),
  space = "Lab" # Option used when colors do not represent a quantitative scale
)
full_df$protein<-as.factor(full_df$protein)
num_colors <- nlevels(full_df$protein)
protein_color_colors <- color_pallete_function(num_colors)


x<-full_df$Distance_.microns.
y<-full_df$normalized_gray_value
plot(x,y, title(image_name),
     col = protein_color_colors[full_df$protein],
     ylab="Normalized Gray Value",
     xlab="Distance (μm)",
     pch=20
)

legend(
  x ="topright",
  legend = levels(full_df$protein), # for readability of legend
  col = protein_color_colors,
  pch = 19, # same as pch=20, just smaller
  cex = .7 # scale the legend to look attractively sized
)

lines(xl1, predict(lo1,xl1), col='green3', lwd=2)
lines(xl2, predict(lo2,xl2), col='magenta', lwd=2)

##### Save plot ----
jpeg(paste(image_name, "plot.jpg", sep="_"), quality = 500)
plot(x,y, title(image_name),
     col = protein_color_colors[full_df$protein],
     ylab="Normalized Gray Value",
     xlab="Distance (μm)",
     pch=20
)
legend(
  x ="topright",
  legend = levels(full_df$protein), # for readability of legend
  col = protein_color_colors,
  pch = 19, # same as pch=20, just smaller
  cex = 1 # scale the legend to look attractively sized
)
lines(xl1, predict(lo1,xl1), col='magenta', lwd=2)
lines(xl2, predict(lo2,xl2), col='green3', lwd=2)
dev.off()

##### SCEP1 normality --------
scep1_normality_test<-shapiro.test(full_df_scep1$Gray_Value)


if (scep1_normality_test$p.value>0.01){
  warning("SCEP1 does not show normal distribution")
  print(scep1_normality_test)
}

if (scep1_normality_test$p.value<0.01){
  print(scep1_normality_test)
  message("SCEP1 shows normal distribution")
  }

##### REC8 normality ----
rec8_normality_test<-shapiro.test(full_df_rec_8$Gray_Value)

if (rec8_normality_test$p.value>0.01){
  warning("REC8 does not show normal distribution")
  print(rec8_normality_test)
}

if (rec8_normality_test$p.value<0.01){
  print(rec8_normality_test)
  message("REC8 shows normal distribution")
}
