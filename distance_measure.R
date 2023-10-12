# Entry files for this script are csv files produced by imageJ when plotting the profile along lines. Run this for each image analyzed.
##### Parameters ----

folder<-"Z:/EQUIPE/meiorec/PROJECTS/TransciptomeScreen/Immuno/STED/WT_zyp1_from_steph/SD128-1_O3i4_decon"
image_ID<-"SD128-1_O3i4_decon"

##### get file names -----

setwd(folder)
getwd()
files<-list.files(folder, pattern = "\\.csv$")
files_raw<-files

a=0
for (i in files){
  a=a+1
  print(a)
  files[a]=substring(files[a],1, nchar(files[a])-4)
}


##### run on each file and make a df ------

distances<-c(1:length(files))
maxi1<-c(1:length(files))
maxi2<-c(1:length(files))

#depending on how many lines were drawn, use either of these two following lines
#file_numbers<-c(1,10,2,3,4,5,6,7,8,9)
file_numbers<-c(1,10,11,12,13,14,15,16,17,18,19,2,20,3,4,5,6,7,8,9)

a=0
c=0
#i="line1.csv"
for (i in files_raw){
  a=c
  a=a+1
  c=a
  a=file_numbers[a]
  print(a)
  img<-read.csv(i)
  x<-img$Distance_.microns.
  y<-img$Gray_Value
  lo <- loess(y~x, span=0.5)
  jpeg(paste(substring(i,1, nchar(i)-4),"jpg", sep="."))
  plot(x,y, title(i))
  xl <- seq(min(x),max(x), (max(x) - min(x))/1000)
  lines(xl, predict(lo,xl), col='red', lwd=2)
  dev.off()

  model<-loess(y~x, span=0.5)

  xpredicted<-seq(from=x[1], to=max(x), by=0.0001)
  df <- data.frame(xpredicted)
  ypredicted<-predict(model,xpredicted)
  interpolated<-data.frame(xpredicted,ypredicted)

  index<-which(diff(sign(diff(interpolated$ypredicted)))==-2)+1

  if (length(index)>2) {
    stop("MORE THAN TWO LOCAL MAXIMA: CHOOSE ANOTHER PROFILE")
  }
  max_x1<-interpolated$xpredicted[index[1]]
  maxi1[a]<-max_x1
  max_x2<-interpolated$xpredicted[index[2]]
  maxi2[a]<-max_x2
  distance<-max_x2-max_x1
  print(distance)
  distances[a]<-distance
}

maxis<-data.frame(maxi1,maxi2)
maxis

##### create df and save as csv ----

file_ID<-paste(image_ID,files,sep = "_")
output_df<-data.frame(file_ID,distances)
print(output_df)

write.csv(output_df,paste(paste(folder, paste("distances", paste(image_ID,"csv",sep="."),sep="_"), sep="/"), sep="/"), row.names = FALSE)

warning("REMOVE THE RESULT FOLDER BEFORE RERUNNING IN SAME FOLDER")

