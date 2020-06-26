##CONTEXT
# We are going to analyse the results of a BeadChip Array experiment (Infinium HumanMethylation450 BeadChip)

#STEP 1
rm(list=ls()) # Removing any previous jobs made in R
setwd("~/Desktop/Bioinformatics/DRD/Report") # Set the working directory in which we have our dataset
suppressMessages(library(minfi)) # We load the minfi library. The supress message avoids showing the output of loading this package.

SampleSheet <- read.table("~/Desktop/Bioinformatics/DRD/Report/Input_data/Samplesheet_report_2020.csv",sep=",",header=T) # We read the samplesheet of the experiment.
SampleSheet #Checking that the object SampleSheet is correct. We have 4 wild-type samples and 4 DS samples.
baseDir <- ("~/Desktop/Bioinformatics/DRD/Report/Input_data")
targets <- read.metharray.sheet(baseDir) # Read the sample sheet of the Illumina experiment.
targets #Checking that the object targets is correct.
RGset <- read.metharray.exp(targets = targets) # From the read sheet we read each idat file
save(RGset,file="RGset.RData") # We save the RGset object in the working directory with the name RGset.RData. This is only necessary if we are not going to do the whole analysis in only one session
RGset
#STEP 2
Red <- data.frame(getRed(RGset)) # We extract and store the red fluorescence in a new object, named Red
dim(Red) # To check the dimension of the Red object.
Green <- data.frame(getGreen(RGset)) # We extract and store the green fluorescence in a new object, named green
dim(Green) # To check the dimension of the Green object.

#STEP 3

load('~/Desktop/Bioinformatics/DRD/Lesson_2/Illumina450Manifest_clean.RData') # We upload the information of the Manifest of this Array (that contains information about the annotation of the probes) in order to use later this information.
# Actually the previous step would require uploading and cleaning the manifest, but this was previously done in a practical lesson of DRD, specifically in the Lesson 2, Chunk 2. Nevertheless I will put the commands required for obtaining the clean manifest herebelow.
# We need to download the Manifest from Internet before reading it with R, from the following website: http://support.illumina.com/array/array_kits/infinium_humanmethylation450_beadchip_kit/downloads.html

## COMMENTED CODE
# Illumina450Manifest <- read.table("HumanMethylation450_15017482_v1-2.csv",sep=",",header=T, nrows=20)
# Illumina450Manifest_clean <- Illumina450Manifest[Illumina450Manifest$CHR!="",]
# Illumina450Manifest_clean <- droplevels(Illumina450Manifest_clean)
# save(Illumina450Manifest_clean,file="Illumina450Manifest_clean.RData")
## END OF COMMENTED CODE

ls()

Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID=="59609499",] #We check, for the adress that has been given to us, the infomation related to the probe, like the ID.
Red[rownames(Red)=="59609499",] #Red fluorescences of the Adress that corresponds to us in the samples of the SampleSheet
Green[rownames(Green)=="59609499",] #Green fluorescences of the Adress that corresponds to us in the samples of the SampleSheet

#The filled table of this step can be seen in the Results file of the following GitHub repository: https://github.com/Sanjo23Chaval/DRD

#STEP 4

MSet.raw <- preprocessRaw(RGset) #Upload the methylation set, that stores unmethylathed and methylated values for each Illumina ID. This is the raw data, without normalization or background correction.
MSet.raw
save(MSet.raw,file="MSet_raw.RData") #We save the raw Mset in the file MSet_raw.RData

#STEP 5

qc <- getQC(MSet.raw)
qc
png(file="QC.png") # We save the QC plot in the QC.png file
plotQC(qc) # As we can see from the plot, saved in the Results file in the Github repository previously mentioned, all the samples have a good quality.
dev.off()

df_TypeControl <- data.frame(getProbeInfo(RGset, type = "Control")) # We select the control probes to check their names and numbers related to, and we store the information in the df_TypeControl object
png(file="Intensities.png") # We save the Intensities plot in the Intensities.png file
controlStripPlot(RGset, controls="NEGATIVE") # Here we plot the intensity values of each type of control probe; nevertheless in the version 4 of R there is a bug that plots wrongly this plot, so the plot is not correctly done, as the Red table indicates the Green fluorescence and viceversa.
dev.off()
detP <- detectionP(RGset) # We calculate the detection P-value for the SampleSheet samples and probes. Depending on its value the probes are reliable or not (the lower the p-value the more reliable they are).
save(detP,file="detP.RData") # We save the p-values in the detP.RData file
load("detP.RData") # or detP <- readRDS("detP.rds")
str(detP)
failed <- detP>0.01 # We set that a p-value higher than 0.01 indicates a not reliable probe, and therefore probes with these values should be filtered out.
dim(failed)
table(failed) # Here we can see how many probes have a p-value higher than 0.01 (TRUE) and lower than 0.01 (FALSE). The TRUE probes are reported as failed positions in the Results file. Take into account that this reports the total number of probes summing all the samples. In the next line we can plot the failed and non failed positions for each sample.
summary(failed) # Here we observe the distribution of probes according to the p-value threshold for each sample.

#STEP 6

beta <- getBeta(MSet.raw) # We calculate the beta values of the samples
beta
beta_wt <- beta[, c(3:5,8)]
beta_ds <- beta[, c(1,2,6,7)]
# With the previous steps the have splitted the beta values into two groups, the WT and the DS.
mean_of_beta_wt <- apply(beta_wt,1,mean,na.rm=T)
mean_of_beta_ds <- apply(beta_ds,1,mean,na.rm=T)
d_mean_of_beta_wt <- density(mean_of_beta_wt)
d_mean_of_beta_ds <- density(mean_of_beta_ds)
# With the previous steps we have calculated the mean of beta of the groups and the density distribution.
M <- getM(MSet.raw) # We calculate the M values of the samples
M_wt <- M[, c(3:5,8)]
M_ds <- M[, c(1,2,6,7)]
# With the two previous lines we have divided our M results into WT and DS respectively.
mean_of_M_wt <- apply(M_wt,1,mean,na.rm=T)
mean_of_M_ds <- apply(M_ds,1,mean,na.rm=T)
d_mean_of_M_ds <- density(mean_of_M_ds)
d_mean_of_M_wt <- density(mean_of_M_wt)
# With the previous steps we have calculated the mean of M of the groups and the density distribution.

png(file="Beta.png")
plot(d_mean_of_beta_ds,main="Density of Beta Values",col="red")
lines(d_mean_of_beta_wt,col="blue")
legend("topright",
  legend = c("WT", "DS"),
  col = c("blue", "red"),
  lwd = c(3,3),
  bty = "n",
  pt.cex = 2,
  cex = 1.2,
  text.col = "black",
  horiz = F ,
  inset = c(0.1, 0.1))
dev.off()
# With the previous block we have represented the density of beta values in a plot that can be seen in the Results file
png(file="M.png")
plot(d_mean_of_M_ds,main="Density of M Values",col="red")
lines(d_mean_of_M_wt,col="blue")
legend("topright",
    legend = c("WT", "DS"),
    col = c("blue", "red"),
    lwd = c(3,3),
    bty = "n",
    pt.cex = 2,
    cex = 1.2,
    text.col = "black",
    horiz = F ,
    inset = c(0.1, 0.1))
dev.off()
# With the previous block we have represented the density of M values in a plot that can be seen in the Results file

dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",] # We identify the Type I probes and store them in the dfI object
dfI <- droplevels(dfI)

dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",] # We identify the Type I probes and store them in the dfII object
dfII <- droplevels(dfII)

# STEP 7

# First of all we will plot the mean beta distribution of all the samples, as we did before but dividing the dataset in two groups.
  mean_of_beta <- apply(beta,1,mean,na.rm=T)
  d_mean_of_beta <- density(mean_of_beta)
  png(file="Beta_all.png")
  plot(d_mean_of_beta,main="Density of Beta Values",col="red")
  dev.off()

beta_I <- beta[rownames(beta) %in% dfI$IlmnID,]
beta_II <- beta[rownames(beta) %in% dfII$IlmnID, ]
mean_of_beta_I <- apply(beta_I,1,mean)
mean_of_beta_II <- apply(beta_II,1,mean)
d_mean_of_beta_I <- density(mean_of_beta_I,na.rm=T)
d_mean_of_beta_II <- density(mean_of_beta_II,na.rm=T)

sd_of_beta_I <- apply(beta_I,1,sd,na.rm=T) # Calculate the standard deviation of the beta values of probes Type I
sd_of_beta_II <- apply(beta_II,1,sd,na.rm=T) # Calculate the standard deviation of the beta values of probes Type II
d_sd_of_beta_I <- density(sd_of_beta_I,) # Calculate the density of standard deviation of the beta values of probes Type I
d_sd_of_beta_II <- density(sd_of_beta_II) # Calculate the density of standard deviation of the beta values of probes Type II


preprocessNoob_results <- preprocessNoob(RGset) #We have produced the normalization using the corresponding function assigned to me.
save(preprocessNoob_results,file="preprocessNoob_results.RData")
load('~/Desktop/Bioinformatics/DRD/Report/preprocessNoob_results.RData')

beta_preprocessNoob <- getBeta(preprocessNoob_results) # We have calculated the Beta values of the normalized dataset

beta_preprocessNoob_I <- beta_preprocessNoob[rownames(beta_preprocessNoob) %in% dfI$IlmnID,]
beta_preprocessNoob_II <- beta_preprocessNoob[rownames(beta_preprocessNoob) %in% dfII$IlmnID,] # With this line and the previous one we have divided the Type I and Type II probes
mean_of_beta_preprocessNoob_I <- apply(beta_preprocessNoob_I,1,mean)
mean_of_beta_preprocessNoob_II <- apply(beta_preprocessNoob_II,1,mean) # With this line and the previous one we have calculated the mean of beta values for each of the two types of probes
d_mean_of_beta_preprocessNoob_I <- density(mean_of_beta_preprocessNoob_I,na.rm=T)
d_mean_of_beta_preprocessNoob_II <- density(mean_of_beta_preprocessNoob_II,na.rm=T) # We have calculated the density of the means of beta values for each type of probe
sd_of_beta_preprocessNoob_II <- apply(beta_preprocessNoob_II,1,sd)
sd_of_beta_preprocessNoob_I <- apply(beta_preprocessNoob_I,1,sd) # We have calculated the standard deviation of the beta values of the two types of probes.
d_sd_of_beta_preprocessNoob_I <- density(sd_of_beta_preprocessNoob_I,na.rm=T)
d_sd_of_beta_preprocessNoob_II <- density(sd_of_beta_preprocessNoob_II,na.rm=T) # Calculate the density of standard deviation of the beta values the two types of probes.


png(file="Step7.png")
par(mfrow=c(2,3))
plot(d_mean_of_beta_I,col="blue",main="Density of raw beta values",xlim=c(0,1),ylim=c(0,6))
lines(d_mean_of_beta_II,col="red")
legend("topright",
    legend = c("TypeI", "TypeII"),
    col = c("blue", "red"),
    lwd = c(3,3),
    bty = "n",
    pt.cex = 2,
    cex = 1.2,
    text.col = "black",
    horiz = F ,
    inset = c(0.1, 0.1))
plot(d_sd_of_beta_I,col="blue",main="Density of raw sd",xlim=c(0,0.6),ylim=c(0,90))
lines(d_sd_of_beta_II,col="red")
legend("topright",
    legend = c("TypeI", "TypeII"),
    col = c("blue", "red"),
    lwd = c(3,3),
    bty = "n",
    pt.cex = 2,
    cex = 1.2,
    text.col = "black",
    horiz = F ,
    inset = c(0.1, 0.1))
boxplot(beta,col=c('red', 'red', 'green', 'green', 'green', 'red', 'red', 'green'),ylim=c(0,1))
plot(d_mean_of_beta_preprocessNoob_I,col="blue",main="Normalized beta density (Noob)",xlim=c(0,1),ylim=c(0,6))
lines(d_mean_of_beta_preprocessNoob_II,col="red")
legend("topright",
    legend = c("TypeI", "TypeII"),
    col = c("blue", "red"),
    lwd = c(3,3),
    bty = "n",
    pt.cex = 2,
    cex = 1.2,
    text.col = "black",
    horiz = F ,
    inset = c(0.1, 0.1))
plot(d_sd_of_beta_preprocessNoob_I,col="blue",main="Normalized sd density (Noob)",xlim=c(0,0.6),ylim=c(0,130))
lines(d_sd_of_beta_preprocessNoob_II,col="red")
legend("topright",
    legend = c("TypeI", "TypeII"),
    col = c("blue", "red"),
    lwd = c(3,3),
    bty = "n",
    pt.cex = 2,
    cex = 1.2,
    text.col = "black",
    horiz = F ,
    inset = c(0.1, 0.1))
boxplot(beta_preprocessNoob,col=c('red', 'red', 'green', 'green', 'green', 'red', 'red', 'green'),ylim=c(0,1))
dev.off()
# As we can see Type II probes are more represented by methylated probes while the Type I probes show the opposite behaviour. The standard deviation of the Type I probes is very high compared to the sd of Type I probes.
# From the boxplot we can see that the means of the different samples are not equal, nevertheless there are no big differences between the two groups at least if we only observe the box-plots
# There are also big diferrences between the raw and normalized plots. Mainly we can observe that in the normalized beta values the two types of probes are more proximal. Also in the normalized box plot the samples have similar distributions while in the raw boxplot disttibutions are not similar.
# The plots can be seen in the Results file

# STEP 8
png(file="Step8.png")
pca_results <- prcomp(t(beta_preprocessNoob), scale = T) #We apply a PCA analysis on the matrix of beta values obtained after the normalization.
plot(pca_results$x[,1], pca_results$x[,2],cex=2,pch=2,col=c("red","red","green","green","green","red", "red","green"),xlab="PC1",ylab="PC2",xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_results$x[,1], pca_results$x[,2],labels=rownames(pca_results$x),cex=0.5,pos=1)
legend("bottomright",legend=c("WT","DS"),col=c("green","red"),pch=2)

dev.off()
# According to the plot, seems that WT samples have low PC1 and PC2 values (except one outlier observed in the right side) while the DS samples are divided into two sub-clusters, one with High PC1 values and the other with high PC2 values.

# STEP 9

My_ttest_function <- function(x) {
  t_test <- t.test(x~   SampleSheet$Group)
  return(t_test$p.value)
}

pValues_ttest <- apply(beta_preprocessNoob,1, My_ttest_function) # We obtain the p-value of the t-test between the two groups for each probe.

final_ttest <- data.frame(beta_preprocessNoob, pValues_ttest) # We create a new data.frame with the beta values along with the p-values of the t-t_test
final_ttest <- final_ttest[order(final_ttest$pValues_ttest),]
final_ttest

save(final_ttest,file="final_ttest.RData")
load('~/Desktop/Bioinformatics/DRD/Report/final_ttest.RData')
#Finally we can check the probes that are significantly differentially expressed. The significance threshold chosen is 0.05 as it is indicated in the Report pipeline.

final_ttest_pvalueth <- final_ttest[final_ttest$pValues_ttest<=0.05,]
dim(final_ttest_pvalueth)

#There are 35396 probes differentially expressed (significantly). This means that the methylation in these probes between the two conditions is diffent and significant.

# STEP 10

raw_pValues <- final_ttest[,9] # We store all the p-values in a vector.

corrected_pValues_BH <- p.adjust(raw_pValues,"BH") #We apply a Benjamini & Hochberg correction
corrected_pValues_Bonf <- p.adjust(raw_pValues,"bonferroni")  #We apply a Bonferroni correction
final_ttest_corrected <- data.frame(final_ttest, corrected_pValues_BH, corrected_pValues_Bonf)

dim(final_ttest[final_ttest$corrected_pValues_BH<=0.05,])
dim(final_ttest[final_ttest$corrected_pValues_Bonf<=0.05,])

#With the corrections we do not find any significant probe. This is normal as we have many probes so the significance threshold is very stringent.

#STEP 11

library(gplots) # We need the library gplots to produce heat maps.

input_heatmap=as.matrix(final_ttest_corrected[1:100,1:8]) #We create the heatmap using the top 100 probes significantly expressed

colorbar <- c("red","red","green","green","green","red","red","green") # We colour the WT samples with green and the DS samples with Red

png(file="heatmap.png")
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'average'),dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F)
dev.off()

#From the heatmap, visible in the results file, we can see that the clustarization is perfect, as we can see the two groups totally divided.

# STEP 12

beta <- final_ttest[,1:8]
beta_groupA <- beta[,SampleSheet$Group=="WT"]
mean_beta_groupA <- apply(beta_groupA,1,mean)
beta_groupB <- beta[,SampleSheet$Group=="DS"]
mean_beta_groupB <- apply(beta_groupB,1,mean)
# With the previous lines we have calculated two matrices, one for the beta values of the Group WT and the other one the beta values of the Group DS.

difference <- mean_beta_groupB-mean_beta_groupA
head(difference)
toVolcPlot <- data.frame(difference, -log10(final_ttest_corrected$pValues_ttest))
head(toVolcPlot)
# We create a dataframe that will store the -log10 of the p-values of the probes as well as the difference between the means of the two groups. These values correspond to the Y and X axis of a Volcano plot respectively.

png(file="volcano.png")
toHighlight <- toVolcPlot[abs(toVolcPlot[,1])>0.1 & toVolcPlot[,2]>(-log10(0.01)),] # WIth this line we are "storing" the probes (points) that are significantly over or underexpressd (with a p-value lower than 0.01 and with a mean difference higher than 0.1 (absolute value))
plot(toVolcPlot[,1], toVolcPlot[,2], main ="Volcano plot", xlab="Mean difference", ylab="-log10(p-value)",pch=16,cex=0.5)
points(toHighlight[,1], toHighlight[,2],pch=16,cex=0.7,col="red")
dev.off()

# The volcano plot can be observed in the Results file. We can see highlighted in red all the probes that are differentially expressed that are significant.

library(gap) #We need to load this library to produce the Manhattan Plot.

final_ttest_corrected_man <- data.frame(rownames(final_ttest_corrected),final_ttest_corrected) #In order to merge on the basis of CpG islands.
colnames(final_ttest_corrected_man)[1] <- "IlmnID"

final_ttest_corrected_man_annotated <- merge(final_ttest_corrected_man, Illumina450Manifest_clean,by="IlmnID")

input_Manhattan <- data.frame(final_ttest_corrected_man_annotated$CHR, final_ttest_corrected_man_annotated$MAPINFO, final_ttest_corrected_man_annotated$pValues_ttest)

input_Manhattan$final_ttest_corrected_man_annotated.CHR <- factor(input_Manhattan$final_ttest_corrected_man_annotated.CHR,levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))
levels(input_Manhattan$final_ttest_corrected_man_annotated.CHR)

palette <- rainbow(24)
png(file="manhattan.png")
mhtplot(input_Manhattan,control=mht.control(colors=palette))
dev.off()

# OPTIONAL PART

chr21 <- Illumina450Manifest_clean[Illumina450Manifest_clean$CHR=="21",] # We choose from the Manifest only the probes associated with the chromosome 21.
beta_21 <- beta[rownames(beta) %in% chr21$IlmnID,]
beta_wt_21 <- beta_21[, c(3:5,8)]
beta_ds_21 <- beta_21[, c(1,2,6,7)]
# With the previous steps the have splitted the beta values into two groups, the WT and the DS.
mean_of_beta_wt_21 <- apply(beta_wt_21,1,mean,na.rm=T)
mean_of_beta_ds_21 <- apply(beta_ds_21,1,mean,na.rm=T)
d_mean_of_beta_wt_21 <- density(mean_of_beta_wt_21)
d_mean_of_beta_ds_21 <- density(mean_of_beta_ds_21)

beta_normal_21


png(file="Beta_chr21.png")
plot(d_mean_of_beta_ds_21,main="Density of Beta Values (Chr21)",col="red")
lines(d_mean_of_beta_wt_21,col="blue")
legend("topright",
  legend = c("WT", "DS"),
  col = c("blue", "red"),
  lwd = c(3,3),
  bty = "n",
  pt.cex = 2,
  cex = 1.2,
  text.col = "black",
  horiz = F ,
  inset = c(0.1, 0.1))
dev.off()

#We can see in the Results file that there are no significant differences between the two groups in terms of total methylation. The beta densities are similar between the two groups. Nevertheless we can observe a small increase in the methylation in the WT groups, as the peak near to the value of beta = 1 is higher than the DS group

beta_21_n <- beta_preprocessNoob[rownames(beta_preprocessNoob) %in% chr21$IlmnID,]
beta_wt_21_n <- beta_21_n[, c(3:5,8)]
beta_ds_21_n <- beta_21_n[, c(1,2,6,7)]
mean_of_beta_wt_21_n <- apply(beta_wt_21_n,1,mean,na.rm=T)
mean_of_beta_ds_21_n <- apply(beta_ds_21_n,1,mean,na.rm=T)
d_mean_of_beta_wt_21_n <- density(mean_of_beta_wt_21_n)
d_mean_of_beta_ds_21_n <- density(mean_of_beta_ds_21_n)

png(file="Beta_chr21_n.png")
plot(d_mean_of_beta_ds_21_n,main="Density of Beta Values (Chr21)",col="red")
lines(d_mean_of_beta_wt_21_n,col="blue")
legend("topright",
  legend = c("WT", "DS"),
  col = c("blue", "red"),
  lwd = c(3,3),
  bty = "n",
  pt.cex = 2,
  cex = 1.2,
  text.col = "black",
  horiz = F ,
  inset = c(0.1, 0.1))
dev.off()
# If we observe the normalized distributions of beta values we can observe again that small difference in the peak close to 1, so we can affirm that there is more methylation in the Chr21 probes in the WT condition.

pValues_ttest_21 <- apply(beta_21_n,1, My_ttest_function)
final_ttest_21 <- data.frame(beta_21_n, pValues_ttest_21)

final_ttest_pvalueth_21 <- final_ttest_21[final_ttest_21$pValues_ttest_21<=0.05,]
dim(final_ttest_pvalueth_21)

# There are 587 differentially expressed without doing any correction
raw_pValues_21 <- final_ttest_21[,9]
corrected_pValues_BH_21 <- p.adjust(raw_pValues_21,"BH") #We apply a Benjamini & Hochberg correction
corrected_pValues_Bonf_21 <- p.adjust(raw_pValues_21,"bonferroni")  #We apply a Bonferroni correction

dim(final_ttest_21[final_ttest_21$corrected_pValues_BH_21<=0.05,])
dim(final_ttest_21[final_ttest_21$corrected_pValues_Bonf_21<=0.05,])

# With both corrections there are no differentially expressed probes.
