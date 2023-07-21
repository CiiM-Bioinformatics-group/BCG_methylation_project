# Step 1: DNA methylation data QC 
  # Step 1.1: load in data 
  # Step 1.2: removes samples with a bad call rate
  # Step 1.3: Bisulfite Conversion Control
  # Step 1.4: check sample mixup
  # Step 1.5: check sex 
  # Step 1.6: probes QC 
  # Step 1.7: preprocessing
  # Step 1.8: estimate the cell type proportions

# Step 2: Epigenome wide association analysis of trained immunity
  # Step 2.1: DAN methylation trimming
  # Step 2.2: Linear mixed model to explore the BCG effect on DNA methylation
  # Step 2.3: Robust linear model to perform the EWAS of trained immunity
  # Step 2.4: associations between urate change and trained immunity 

# Step 3: Variance explanation and prediction model 

# Step 4: Mediation analysis



########################################## Step 1 #######################################################

######################### Step 1: DNA methylation data QC ################### 
# reads in the targets and the idat files
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
library(missMethyl)
library(Gviz)
library(DMRcate)
library(stringr)
library(ggplot2)
library(wateRmelon) # We use this function to calculate the number of beads per cpg
library(future)
library(gtools)
library(ggplot2)

############# Step 1.1: Load in data #################
# Reading in the sample sheet
targets <- read.metharray.sheet(ConfigList$SAMPLE_SHEET, pattern="csv$")
targets$Basename <- paste0(ConfigList$IDAT_PATH, targets$Slide, "/", targets$Slide, "_", targets$Array)
# Some sample names use - instead of _
# We rename these samples
targets$Sample_Name <- gsub("-", "_", targets$Sample_Name)
# We re-order the targets	
targets	<- targets[order(targets$Sample_Name), ]
# We load the meta data and add it to the targets object
metadata <- read.csv(ConfigList$SAMPLE_INFO)
targets <- merge(targets,metadata, all.x = TRUE, by = "Sample_Name")
# We have some samples that are technical replicates. These replicates get a d after their name.
# 300BCG007_v1 --> 300BCG007d_v1
duplicates <- targets[duplicated(targets$Sample_Name),1]
targets[duplicated(targets$Sample_Name),1] <- paste0(substr(duplicates, start = 1, stop = 9), "d", substr(duplicates, start = 10, stop = nchar(duplicates)))
# read in the raw data from the IDAT files
rgSet <- read.metharray.exp(targets = targets, extended = TRUE)
##### We save the intermediate outputs
# We write the targets file as a csv
print("We save the partially filtered targets files as a csv")
saveRDS(targets, paste0(ConfigList$OUTPUT_DIR, "targets_partial_filtered_01.rds"))
# We save the partially filtered red and green channel set object
print("We save the partially filtered extended Channel set")
saveRDS(rgSet, paste0(ConfigList$OUTPUT_DIR, "rgSet_partial_filtered_01.rds"))



######################### Step 1.2: remove samples with bad call rate #####################
# For each sample, we calculate the percentage of probes that have a detection p-value samller than the detection threshold.
# We then check if this percentage is larger than 0.96, which is our condition to include a sample. 
# calculate the detection p-values for each position
detP <- detectionP(rgSet)
# We check which probes have a detection p-value smaller then the detection threshold
detP.thresholded <- detP <= as.numeric(ConfigList$DETECTION_THRESHOLD)
keep <- colMeans(detP.thresholded) > 0.96
# We remove all the poor quality samples from the channel set
rgSet <- rgSet[,keep]
# We remove the poor quality samples from the detection p-value table and the thresholded table
detP <- detP[,keep]
detP.thresholded <- detP.thresholded[,keep]
# We remove the poor quality samples from the targets object
targets <- targets[keep,]
# We note how many samples were removed because they had a bad call rate
#filtered.samples.ls <- read.csv(paste0(ConfigList$OUTPUT_DIR, "List_removed_SamplesReasons.csv"), stringsAsFactors = FALSE)
if(length(names(keep)[!keep]) > 0){
  samples.remove <- data.frame(Sample = names(keep)[!keep],
                               Reason = "Bad Call Rate")
  filtered.samples.ls <- samples.remove
} else{
  filtered.samples.ls <- data.frame(Sample = "None",
                                    Reason = "Bad Call Rate")
}
# We save the list of samples we removed and for what reason
write.csv(filtered.samples.ls, paste0(ConfigList$OUTPUT_DIR, "List_removed_SamplesReasons.csv"), row.names = FALSE, quote = FALSE)
# We write how many samples were removed
system(paste0("echo how many samples were removed because they have a bad call rate >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo Samples: ", sum(!keep), " >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo Probes: ", dim(detP)[1], " >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
# We save the intermediate results
# We save the partially filtered targets file
print("We save the partially filtered targets file after call rate removal")
saveRDS(targets, paste0(ConfigList$OUTPUT_DIR, "targets_partial_filtered_02.rds"))
# We save the partially filtered red and green channel set object
print("We save the partially filtered extended Channel set after call rate removal")
saveRDS(rgSet, paste0(ConfigList$OUTPUT_DIR, "rgSet_partial_filtered_02.rds"))
print("We save the partially filtered detection p-value matrix call rate removal")
saveRDS(detP, paste0(ConfigList$OUTPUT_DIR, "detP_partial_filtered_02.rds"))
print("We save the partially filtered thresholded detection p-value matrix call rate removal")
saveRDS(detP.thresholded, paste0(ConfigList$OUTPUT_DIR, "detP_partial_filtered_thresholded_02.rds"))



#################### Step 1.3: Bisulfite Conversion Control ###########################
# We calculate the mean intensity of the Bisulfite Conversion Control probes for all positions and samples.
# We do this for the Type I and Type II control probes separately.
# We then compare determine for each sample, if their mean intensity is more than 3 SD away from the mean intensity.
# If this is true for the Type I or Type II probes, we will remove the sample.
# We read in the output from the previous file, but only if the variable does not already exist
# This way, we save time when we run the script in order, but we can run them indepently.
if(!exists("targets")){
  print("Reading the targets")
  targets <- readRDS(paste0(ConfigList$OUTPUT_DIR, "targets_partial_filtered_02.rds"))
}

if(!exists("rgSet")){
  print("Reading the rgSet")
  rgSet <- readRDS(paste0(ConfigList$OUTPUT_DIR, "rgSet_partial_filtered_02.rds"))
}

if(!exists("detP")){
  print("Reading the detP matrix")
  detP <- readRDS(paste0(ConfigList$OUTPUT_DIR, "detP_partial_filtered_02.rds"))
}

if(!exists("detP.thresholded")){
  print("Reading the detP thresholded matrix")
  detP.thresholded <- readRDS(paste0(ConfigList$OUTPUT_DIR, "detP_partial_filtered_thresholded_02.rds"))
}

# Annotation
#if(!exists("anno850k")){
#  print("Reading the EPIC annotation")
#  anno850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
#}

print("We get the intensities of the control probes")
# We use the probes recommended by Zhou et al. to calculate the GCT (Green CpC to TpC) score
# Zhou et al. do not give a threshold, but most of their samples have a GCT value below 1.5 
GCT <- control.intensities(rgSet)
# Which samples have a GCT larger than 1.5?
samples.exclude <- names(GCT[GCT > 1.5])
# We check if we even have samples with a GCT below 
if(length(samples.exclude) > 0){
  # We load the samples we already removed
  filtered.samples.ls <- read.csv(paste0(ConfigList$OUTPUT_DIR, "List_removed_SamplesReasons.csv"), stringsAsFactors = FALSE)

  # We add the samples that are removed to the list of removed samples
  samples.exclude <- data.frame(Sample = samples.exclude,
                                Reason = "GCT larger 1.5")
  filtered.samples.ls <- rbind(filtered.samples.ls, samples.exclude)
  write.table(filtered.samples.ls, paste0(ConfigList$OUTPUT_DIR, "List_removed_SamplesReasons.csv"), row.names = FALSE, sep = ",", quote = FALSE)

  print("Remove the samples from the various objects")
  keep <- !targets$Sample_Name %in% samples.exclude$Sample
  rgSet <- rgSet[,keep]
  targets <- targets[keep,]
  detP <- detP[,keep]
  detP.thresholded <- detP.thresholded[,keep]

  print("We write how many samples were removed")
  system(paste0("echo How many samples were removed because they had a low Bisulfite Control Intensity >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
  system(paste0("echo Samples: ", sum(!keep), " >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
  system(paste0("echo Probes: ", dim(detP)[1], " >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
  system(paste0("echo >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
} else {
  print("We write how many samples were removed")
  system(paste0("echo How many samples were removed because they had a low Bisulfite Control Intensity >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
  system(paste0("echo No samples wre removed because they had a GCT larger than 1.5. >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
  system(paste0("echo Samples: 0 >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
  system(paste0("echo Probes: ", dim(detP)[1], " >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
  system(paste0("echo >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
}
# We save our partially filtered data
saveRDS(targets, paste0(ConfigList$OUTPUT_DIR, "targets_partial_filtered_03.rds"))
write.csv(targets, paste0(ConfigList$OUTPUT_DIR, "targets_premixup.csv"), row.names = FALSE)
saveRDS(rgSet, paste0(ConfigList$OUTPUT_DIR, "rgSet_partial_filtered_03.rds"))



######################### Step 1.4 check samples mixup ####################################
# The script takes the output from 01_LodingData.R and calculates the correlations of the SNPs
# After running the script you have to check the result by hand and determine which samples are mixed up
# We get the beta values for the SNPs of the array.
# Using this we can determine if a sample has been mixed up.
snps <- getSnpBeta(rgSet)
snps <- as.data.frame(snps)

# We write the SNPs to easily recreate the plots
saveRDS(snps, paste0(ConfigList$OUTPUT_DIR_CORR, "snpsbeta.rds"))

colnames(snps) <- targets$Sample_Name
sample.prefices <- gsub("_.*", "", colnames(snps))
sample.prefices <- substr(sample.prefices, start = 1, stop = 9)
sample.prefices <- unique(sample.prefices)

# We iterate over all samples
print("We generate the correlation plots per sample combination")
correlations <- c()
for(i in 1:length(sample.prefices)){
  # We get the position of the samples in the snps matrix
  positions <- grep(sample.prefices[i], colnames(snps))
  permutations <- t(combn(positions,2))
  if(length(positions) > 1){
    for(j in 1:dim(permutations)[1]){
      result <- cor(snps[,permutations[j,1]], snps[,permutations[j,2]])
      result <- data.frame(x = colnames(snps)[permutations[j,1]],
                           y = colnames(snps)[permutations[j,2]],
                           Cor = result)
      correlations <- rbind(correlations,result)

      if(result$Cor <= 0.9){
        pdf(paste0(ConfigList$OUTPUT_DIR_CORR, "Correlation_", colnames(snps)[permutations[j,1]], "_", colnames(snps)[permutations[j,2]], ".pdf"), width = 6, height = 6)
        print(ggplot(snps, aes(x = snps[,permutations[j,1]], y = snps[,permutations[j,2]])) +
                geom_point() +
                ggtitle(label = paste0(colnames(snps)[permutations[j,1]], " vs ", colnames(snps)[permutations[j,2]])) +
                xlab(label = colnames(snps)[permutations[j,1]]) +
                ylab(label = colnames(snps)[permutations[j,2]]))
        dev.off()
      }
    }
  }
}

# We save the correlation values
print("We save the correlation values")
write.csv(correlations, paste0(ConfigList$OUTPUT_DIR, "SNP_correlations.csv"), row.names = FALSE, quote = FALSE)
#write.csv(subset(correlations, Cor < as.numeric(ConfigList$SNP_COR_THRESHOLD)), paste0(ConfigList$OUTPUT_DIR, "SNP_correlations_belowthreshold.csv"), row.names = FALSE, quote = FALSE)

correlations.relevant <- subset(correlations, Cor < as.numeric(ConfigList$SNP_COR_THRESHOLD))
correlations.relevant <- unique(correlations.relevant$x)
correlations.relevant <- substr(correlations.relevant, start = 1, stop = 9)
correlations.relevant <- correlations[substr(correlations$x, start = 1, stop = 9) %in% correlations.relevant,]
write.csv(correlations.relevant, paste0(ConfigList$OUTPUT_DIR, "SNP_correlations_belowthreshold.csv"), row.names = FALSE, quote = FALSE)

# We genrate a scatter plot indicating which samples have failed
print("We generate a scatter plot indicating which samples have been mixed up")
correlations$Comparison <- paste0(correlations$x, "_vs_", correlations$y)
correlations_t1 <- correlations[grep("_v1", correlations$x),]
correlations_t1t2 <- correlations_t1[grep("_v2", correlations_t1$y),]
correlations_t1t3 <- correlations_t1[grep("_v3", correlations_t1$y),]

correlations_t2 <- correlations[grep("_v2", correlations$x),]
correlations_t2t3 <- correlations[grep("_v3", correlations_t2$y),]

reverse.scatterplot <- function(data, title = "Insert Title", ylim = c(min(correlations$Cor),1)){
  ggplot(data, aes(x = Comparison, y = Cor)) +
    geom_point() +
    geom_hline(yintercept = 0.9) +
    ylim(ylim) +
    ggtitle(label = title) +
    xlab(label = "Compared samples") +
    ylab(label = "Pearson Correlation") +
    theme(axis.text.x = element_text(size = 10,
                                     angle = 90,
                                     vjust = 0.5))
}

parts <- ceiling(dim(correlations_t1t2)[1]/4)
starts <- 0:3*parts+1
stops <- c(1:3*parts,dim(correlations_t1t2)[1])
p1 <- reverse.scatterplot(data = correlations_t1t2[starts[1]:stops[1],], title = "Correlation of SNP beta values T1 vs T2, part 1")
p2 <- reverse.scatterplot(data = correlations_t1t2[starts[2]:stops[2],], title = "Correlation of SNP beta values T1 vs T2, part 2")
p3 <- reverse.scatterplot(data = correlations_t1t2[starts[3]:stops[3],], title = "Correlation of SNP beta values T1 vs T2, part 3")
p4 <- reverse.scatterplot(data = correlations_t1t2[starts[4]:stops[4],], title = "Correlation of SNP beta values T1 vs T2, part 4")

pdf(paste0(ConfigList$OUTPUT_DIR, "SNP_correlations_t1_t2.pdf"), width = 30, height = 10)
grid.arrange(p1,p2,p3,p4, layout_matrix = matrix(1:4, ncol = 2))
dev.off()

parts <- ceiling(dim(correlations_t1t3)[1]/4)
starts <- 0:3*parts+1
stops <- c(1:3*parts,dim(correlations_t1t3)[1])
p1 <- reverse.scatterplot(data = correlations_t1t3[starts[1]:stops[1],], title = "Correlation of SNP beta values T1 vs T3, part 1")
p2 <- reverse.scatterplot(data = correlations_t1t3[starts[2]:stops[2],], title = "Correlation of SNP beta values T1 vs T3, part 2")
p3 <- reverse.scatterplot(data = correlations_t1t3[starts[3]:stops[3],], title = "Correlation of SNP beta values T1 vs T3, part 3")
p4 <- reverse.scatterplot(data = correlations_t1t3[starts[4]:stops[4],], title = "Correlation of SNP beta values T1 vs T3, part 4")

pdf(paste0(ConfigList$OUTPUT_DIR, "SNP_correlations_t1_t3.pdf"), width = 30, height = 10)
grid.arrange(p1,p2,p3,p4, layout_matrix = matrix(1:4, ncol = 2))
dev.off()
parts <- ceiling(dim(correlations_t2t3)[1]/4)
starts <- 0:3*parts+1
stops <- c(1:3*parts,dim(correlations_t2t3)[1])
p1 <- reverse.scatterplot(data = correlations_t2t3[starts[1]:stops[1],], title = "Correlation of SNP beta values T2 vs T3, part 1")
p2 <- reverse.scatterplot(data = correlations_t2t3[starts[2]:stops[2],], title = "Correlation of SNP beta values T2 vs T3, part 2")
p3 <- reverse.scatterplot(data = correlations_t2t3[starts[3]:stops[3],], title = "Correlation of SNP beta values T2 vs T3, part 3")
p4 <- reverse.scatterplot(data = correlations_t2t3[starts[4]:stops[4],], title = "Correlation of SNP beta values T2 vs T3, part 4")

pdf(paste0(ConfigList$OUTPUT_DIR, "SNP_correlations_t2_t3.pdf"), width = 30, height = 10)
grid.arrange(p1,p2,p3,p4, layout_matrix = matrix(1:4, ncol = 2))
dev.off()




########################## Step 1.5: check sex #########################################
####################################################################################################
# We check, if the predicted sex for each sample is equal to the sex we have
# We create a GenomicMethylSet object, which we need to perform the sex check on
# The preprocessFunnorm function also estimates the sex using the getSex function,
# if we don't supply the sex for each sample.
# According to the code male is 1 and female is 2.
mSetSex <- preprocessQuantile(rgSet)

SexEstimatePre <- mSetSex$predictedSex

# We note how many and which samples have a different sex according to the prediction
samples.sex.ls <- data.frame(Sample       = mSetSex@colData$Sample_Name,
                             Sex          = mSetSex@colData$sex,
                             PredictedSex = SexEstimatePre)
samples.sex.ls$Sex <- ifelse(samples.sex.ls$Sex == "female", "F", "M")
samples.sex.ls$Congruent <- samples.sex.ls$Sex == samples.sex.ls$PredictedSex

write.csv(samples.sex.ls, paste0(ConfigList$OUTPUT_DIR, "SexCheck.csv"))

filtered.samples.ls <- read.csv(paste0(ConfigList$OUTPUT_DIR, "List_removed_SamplesReasons.csv"))
# We note which samples were removed because their predicted sex and reported sex is not congruent
samples.remove <- subset(samples.sex.ls, !samples.sex.ls$Congruent)
samples.remove <- data.frame(Sample = samples.remove$Sample,
                             Reason = "wrong predicted sex")
filtered.samples.ls <- rbind(filtered.samples.ls, samples.remove)
write.csv(filtered.samples.ls, paste0(ConfigList$OUTPUT_DIR, "List_removed_SamplesReasons.csv"), row.names = FALSE, quote = FALSE)

# We plot the sex of the samples
pdf(paste0(ConfigList$OUTPUT_DIR, "Plot_Sex.pdf"), height = 5, width = 5)
plotSex(mSetSex)
dev.off()

# We remove the samples with the wrong predicted sex from our data
keep <- !targets$Sample_Name %in% samples.remove$Sample
targets <- subset(targets, !Sample_Name %in% samples.remove$Sample)
rgSet <- rgSet[,keep]

# We calculate the detection p-values for each position
detP <- detectionP(rgSet)

# We check which probes have a detection p-value smaller then the detection threshold
detP.thresholded <- detP <= as.numeric(ConfigList$DETECTION_THRESHOLD)
# We write which samples are removed, because their predicted sex and their actual sex are not congruent
system(paste0("echo How many samples are removed, because their reported and predicted sex are not congruent >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo Samples: ", sum(!samples.sex.ls$Congruent), " >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))

##### We save the intermediate results again
# We save the partially filtered targets file
print("We save the partially filtered targets file after wrong predicted sex removal")
saveRDS(targets, paste0(ConfigList$OUTPUT_DIR, "targets_partial_filtered_06.rds"))

# We save the partially filtered red and green channel set object
print("We save the partially filtered extended Channel set after wrong predicted sex removal")
saveRDS(rgSet, paste0(ConfigList$OUTPUT_DIR, "rgSet_partial_filtered_06.rds"))

print("We save the partially filtered detection p-value matrix call rate removal")
saveRDS(detP, paste0(ConfigList$OUTPUT_DIR, "detP_partial_filtered_06.rds"))

print("We save the partially filtered thresholded detection p-value matrix call rate removal")
saveRDS(detP.thresholded, paste0(ConfigList$OUTPUT_DIR, "detP_partial_filtered_thresholded_06.rds"))





################################# Step 1.6: probes QC ######################################
# We now perform the call rate test for the probes our self
# We remove all probes that have a detection p-value larger than 0.01 in 10% of the samples
keep <- !rowMeans(detP.thresholded) < 0.9

detP <- detP[keep,]
detP.thresholded <- detP.thresholded[keep,]
rgSet <- subsetByLoci(rgSet, includeLoci = rownames(detP))

# We note which probes have a bad call rate
probes.remove <- data.frame(Probe_Id = names(keep)[!keep],
                            Reason = "Call Rate")

filtered.probes.ls <- probes.remove

# We note how many probes have a bad call rate
system(paste0("echo How many probes have been removed because they have a bad call rate >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo Samples: ", dim(detP)[2], " >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo Probes: ", sum(!keep), " >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))

# We now remove probes that have a bead number of < 3 in 10% of the samples
beads_per_cpg <- beadcount(rgSet)

# Calculating how many CpGs have less then 3 beads in 10% of the samples
keep <- !rowSums(is.na(beads_per_cpg)) > ceiling(0.1 * ncol(beads_per_cpg))

detP <- detP[keep,]
detP.thresholded <- detP.thresholded[keep,]
rgSet <- subsetByLoci(rgSet, includeLoci = rownames(detP))

# We note which probes have a bead number less than 3 in 10 percent of the samples
probes.remove <- data.frame(Probe_Id = names(keep)[!keep],
                            Reason = "Bead Number")

filtered.probes.ls <- rbind(filtered.probes.ls, probes.remove)

# How many probes have been removed because of their bead number
system(paste0("echo How many probes have been removed because they have a bead number less than 3 in 10 percent of the samples >> ",ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo Samples: ", dim(rgSet)[2], " >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo Probes: ", sum(!keep), " >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))


# We remove probes that are cross-reactive or are significantly influenced by SNPs
# We remove the cross-reactive probes from Pidsley et al. 2016 and Zhou et al 2017
# From Zhou we only remove the probes for the Europeans (EPIC.hg19.manifest.pop.tsv, Column MASK_general_EUR, MASK_snp5_EUR)
cross <- read.csv(ConfigList$PROBES_CROSS, stringsAsFactors = FALSE)
snp <- read.csv(ConfigList$PROBES_SNP, stringsAsFactors = FALSE)

# We remove the cross reactive probes
keep <- !rownames(detP) %in% cross$probes
cross <- cross[cross$probes %in% rownames(detP),]

detP <- detP[keep,]
detP.thresholded <- detP.thresholded[keep,]
rgSet <- subsetByLoci(rgSet, includeLoci = rownames(detP))

# We note which probes are cross reactive
probes.remove <- data.frame(Probe_Id = cross,
                            Reason = "Cross-Reactive")

filtered.probes.ls <- rbind(filtered.probes.ls, probes.remove)

# How many probes have been removed because they are cross-reactive
system(paste0("echo How many probes have been removed because they are cross-reactive >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo Samples: ", dim(rgSet)[2], " >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo Probes: ", sum(!keep), " >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))

# We remove the probes affected by an SNP
keep <- !rownames(detP) %in% snp$probes
snp <- snp[snp$probes %in% rownames(detP),]

detP <- detP[keep,]
detP.thresholded <- detP.thresholded[keep,]
rgSet <- subsetByLoci(rgSet, includeLoci = rownames(detP))

# We note which probes have a bad call rate
probes.remove <- data.frame(Probe_Id = snp,
                            Reason = "SNP")

filtered.probes.ls <- rbind(filtered.probes.ls, probes.remove)

# How many probes have been removed because they are affected by an SNP
system(paste0("echo How many probes have been removed because they are affected by an SNP >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo Samples: ", dim(rgSet)[2], " >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo Probes: ", sum(!keep), " >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))


# We write which probes have been removed
write.csv(filtered.probes.ls, paste0(ConfigList$OUTPUT_DIR, "List_removed_ProbesReasons.csv"), row.names = FALSE)

# We write how many samples and probes we are left with after QC
system(paste0("echo How many samples and probes we are left with after QC >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo Samples: ", dim(detP)[2], " >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo Probes: ", dim(detP)[1], " >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))

# We save the intermediate results
# We save the partially filtered targets file
print("We save the filtered targets file after probe QC")
saveRDS(targets, paste0(ConfigList$OUTPUT_DIR, "targets_filtered_bulk.rds"))

# We save the partially filtered red and green channel set object
print("We save the filtered extended Channel set after probe QC")
saveRDS(rgSet, paste0(ConfigList$OUTPUT_DIR, "rgSet_filtered_bulk.rds"))

print("We save the filtered detection p-value matrix probe QC")
saveRDS(detP, paste0(ConfigList$OUTPUT_DIR, "detP_filtered_bulk.rds"))

print("We save the filtered thresholded detection p-value matrix probe QC")
saveRDS(detP.thresholded, paste0(ConfigList$OUTPUT_DIR, "detP_filtered_thresholded_bulk.rds"))




############################## Step 1.7: preprocessing ######################################
# We preprocess the data 
# The preprocessFunnorm function produces a genomic ratio set
sex <- as.character(rgSet@colData$sex)
sex <- ifelse(sex == "female", 2, 1)
GRSet<- preprocessFunnorm(rgSet, sex = sex)

# We subset to the sex chromosomes
keep <- (featureNames(GRSet) %in% anno850k$Name[anno850k$chr %in% c("chrX","chrY")])
GRSetSex <- GRSet[keep,]
rgSetSex <- subsetByLoci(rgSet, includeLoci = rownames(GRSetSex))

sex <- ifelse(targets$sex == "female", "F", "M")
GRSetSex <- preprocessQuantile(rgSetSex, sex = sex)

print("We save the filtered rgSet object, that contains only the sex chromosomes")
#saveRDS(rgSetSex, paste0(ConfigList$OUTPUT_DIR, "rgSet_filtered_bulk_XY.rds"))

print("We save the filtered GRSet object, that contains only the sex chromosomes")
#saveRDS(GRSetSex, paste0(ConfigList$OUTPUT_DIR, "GRSet_filtered_bulk_XY_quantile.rds"))

# We remove the sex chromosomes from the GRSet object
keep <- !keep
GRSet <- GRSet[keep,]
rgSet <- subsetByLoci(rgSet, includeLoci = rownames(GRSet))

# We write which probes have been removed because they are on the sex chromosomes
print("We write which probes have been removed because they are on the sex chromosomes")
filtered.probes.ls <- read.csv(paste0(ConfigList$OUTPUT_DIR, "List_removed_ProbesReasons.csv"))
x.probes <- GRSetSex[(featureNames(GRSetSex) %in% anno850k$Name[anno850k$chr %in% c("chrX")]),]
x.probes <- featureNames(x.probes)
x.probes <- data.frame(Probe_Id = x.probes, Reason = "chrX")
filtered.probes.ls <- rbind(filtered.probes.ls,x.probes)
y.probes <- GRSetSex[(featureNames(GRSetSex) %in% anno850k$Name[anno850k$chr %in% c("chrY")]),]
y.probes <- featureNames(y.probes)
y.probes <- data.frame(Probe_Id = y.probes, Reason = "chrY")
filtered.probes.ls <- rbind(filtered.probes.ls,y.probes)
write.csv(filtered.probes.ls, paste0(ConfigList$OUTPUT_DIR, "List_removed_ProbesReasons.csv"), row.names = FALSE)

# How many probes were removed, because they are located on the sex chromosome
print("How many probes were removed, because they are located on the sex chromosomes")
system(paste0("echo How many probes were removed because they are located on the sex chromosomes >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo Samples: ", dim(GRSet)[2], " >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo Probes: ", sum(!keep), " >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo This leaves us with this many probes >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo Probes: ", sum(keep), " >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))
system(paste0("echo >> ", ConfigList$OUTPUT_DIR, "Probes_Samples"))

# We perform Noob normalisation
GRSet.noob <- preprocessNoob(rgSet)
saveRDS(GRSet.noob, paste0(ConfigList$OUTPUT_DIR, "GRSet_filtered_bulk_nosex_noob.rds"))




################################ Step 1.8: estimate the cell type proportions ################
# We estimate the cell types using the estimateCellCounts2 function
# This function uses a newer reference set and produces more accurate results.
# http://52.71.54.154/packages/devel/data/experiment/vignettes/FlowSorted.Blood.EPIC/inst/doc/FlowSorted.Blood.EPIC.html

cells2 <- estimateCellCounts2(rgset, sex = targets$sex, referencePlatform = "IlluminaHumanMethylationEPIC")
cells2 <- cells2$counts
cells2[cells2 < 0] <- 0



########################################## Step 2 #######################################################


############################## Step 2.1: DAN methylation trimming #############################
M.val <- readRDS("mVals_filtered.rds")
removeOutliers<-function(probes){
  require(matrixStats)
  if(nrow(probes) < ncol(probes)) warning("expecting probes are rows (long dataset)")
  rowIQR <- rowIQRs(probes, na.rm = T)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
  maskL <- probes < row2575[,1] - 3 * rowIQR
  maskU <- probes > row2575[,2] + 3 * rowIQR
  initial_NAs<-rowSums(is.na(probes))
  probes[maskL] <- NA
  removed_lower <- rowSums(is.na(probes))-initial_NAs
  probes[maskU] <- NA
  removed_upper <- rowSums(is.na(probes))-removed_lower-initial_NAs
  N_for_probe<-rowSums(!is.na(probes))
  Log<-data.frame(initial_NAs,removed_lower,removed_upper,N_for_probe)
  return(list(probes, Log))
}
system.time(OutlierResults<-removeOutliers(M.val))
M.val<-OutlierResults[[1]]
Log<-OutlierResults[[2]]
save(M.val, file="mVals_filtered_trimmed.Rdata")



############################## Step 2.2: Linear mixed model to explore the BCG effect on DNA methylation ###############

pheno <- readRDS("pheno_filtered.rds")
load("mVals_filtered_trimmed.Rdata")
sum(colnames(M.val) == rownames(pheno))
remove <- rownames(pheno %>% filter(Sample_Plate== "plate11"))
pheno <- pheno[-which(rownames(pheno)%in% remove),]
mVals <- M.val[,-which(colnames(M.val)%in% remove)]
sum(colnames(mVals) == rownames(pheno))
mVals.t <- t(mVals)
pheno$Sample_Plate <- droplevels(pheno$Sample_Plate)
pheno$time <- as.factor(pheno$time)
pheno$age <- as.numeric(pheno$age)
str(pheno)
lmer.test <- function(column, pheno, mVals, comparison){
  cpg <- colnames(mVals)[column]
  percentage <- column/dim(mVals)[2]
  print(paste0("Comparison: ", comparison, ", CpG: ", cpg, ", ", percentage))
  data <- cbind(mVals[,column], pheno)
  colnames(data)[1] <- "methylation"
  ML <- lmer(methylation ~ time + sex + age + (1|id) + Sample_Plate + CD8T + CD4T + NK + Bcell + Mono + Neu, data = data)
  result <- Anova(ML)
  if(is.null(ML@optinfo$conv$lme4$messages)){
    message <- "Nothing"
  } else {
    message <- ML@optinfo$conv$lme4$messages
    message <- paste0(message, collapse = " & ")
  }
  result <- data.table(CpG = cpg, result[1,], Message = message)
  colnames(result)[4] <- "P"
  result
}

columns <- dim(mVals.t)[2]
num.cores = as.numeric(10)
print("Start Modelling")
system.time(Comparison <- mclapply(1:columns, lmer.test, pheno = pheno, mVals = mVals.t, comparison = "LMEM", mc.cores = num.cores))

print("We rbind the results")
Comparison <- rbindlist(Comparison)
Comparison$P.adj <- p.adjust(Comparison$P, "fdr")
Comparison <- Comparison[order(Comparison$P.adj, decreasing = FALSE),]

save(Comparison, file = "output/lme/all_lme.rdata")




############################## Step 2.3: Robust linear model to perform the EWAS of trained immunity ###############
# we take IFN-gamma as an example #
TI <- read_excel("input/ti_value.xlsx") ### 278 samples 
load("mVals_filtered_trimmed.Rdata")
pd <- readRDS("pheno_filtered.rds")
pdv1 <- pd %>% filter(time == 0)
mv1.t <- t(M.val)
pdv1$id2 <- paste0("X", rownames(pdv1))

remove <- pdv1 %>% filter(Sample_Plate== "plate11") %>% pull(id2)
TI <- TI[-which(paste0("X", TI$ID2) %in% remove),]   ## 276 samples left   
mv1.t <- mv1.t[-which(paste0("X", rownames(mv1.t)) %in% remove),]           
pdv1 <- pdv1[-which(pdv1$id2 %in% remove),]

female <- pdv1 %>% filter(sex %in% c("female", "male")) %>% pull(id2) # 123
TI.female <- TI[which(paste0("X", TI$ID2) %in% female),] ## 100 
mv1.t <- mv1.t[which(rownames(mv1.t) %in% TI.female$ID2),] # 100
pdv1 <- pdv1[which(pdv1$id2 %in% paste0("X",TI.female$ID2)),]

shapiro.test(na.omit(TI.female$FC_IL1b_T3_13)) 
shapiro.test(na.omit(TI.female$FC_IL6_T3_13)) 
shapiro.test(na.omit(TI.female$FC_TNF_T3_13))
shapiro.test(na.omit(TI.female$FC_IFNg_W3_13)) 

ti2 <- data.frame(id = TI.female$ID2, value = TI.female$FC_IFNg_W3_13)
ti2 <- na.omit(ti2) 
mv1.t <- mv1.t[which(rownames(mv1.t) %in% ti2$id),] 
pdv1 <- pdv1[which(pdv1$id2 %in% paste0("X",ti2$id)),]

pdv1$Sample_Plate <- droplevels(pdv1$Sample_Plate)
str(pdv1)

sum(pdv1$id2 == paste0("X",ti2$id))
sum(pdv1$id2 == paste0("X", rownames(mv1.t)))

ti2$value_nor <- qnorm((rank(ti2$value, na.last="keep") - 0.5) / sum(!is.na(ti2$value)))
shapiro.test(ti2$value_nor)

RLMtest = function(meth_matrix, methcol, TI, age, gender, plate, CD8, CD4, NK, B, Mono, Neu) {mod = try(rlm(TI ~ meth_matrix[, methcol]+age+gender+plate+CD8+CD4+NK+B+Mono+Neu, maxit=200))
	cf = try(coeftest(mod, vcov=vcovHC(mod, type="HC0")))

if (class(cf)=="try-error") {
  bad <- as.numeric(rep(NA, 3))
  names(bad)<- c("Estimate", "Std. Error", "Pr(>|z|)")
  bad
}
else{
  cf[2, c("Estimate", "Std. Error", "Pr(>|z|)")]
}
}

res <- lapply(setNames(seq_len(ncol(mv1.t)), dimnames(mv1.t)[[2]]), RLMtest, meth_matrix=mv1.t, TI=ti2$value_nor, age=pdv1$age, gender = pdv1$sex, plate=pdv1$Sample_Plate, CD8=pdv1$CD8T, CD4=pdv1$CD4T, NK=pdv1$NK, B=pdv1$Bcell, Mono=pdv1$Mono, Neu=pdv1$Neu)

setattr(res, 'class', 'data.frame')
setattr(res, "row.names", c(NA_integer_,4))
setattr(res, "names", make.names(names(res), unique=TRUE))
probelistnamesB <- names(res)
result <- t(data.table(res))
result<-data.table(result)
result[, probeID := probelistnamesB]
setnames(result, c("BETA","SE", "P_VAL", "probeID")) # rename columns
setcolorder(result, c("probeID","BETA","SE", "P_VAL"))
result$padj <- p.adjust(result$P_VAL, method ="BH")
write.xlsx(result, file = "output/IFN_v1_methy_rlm_mval_trimmed.xlsx")





############################## Step 2.4: associations between urate change and trained immunity ###############

# we calculated the DNA methylation change between every two time points ########
load("mv1.t.rata") ## time 1 DNA methylation M value matrix
load("pdv1.rdata")

remove <- pdv1 %>% filter(Sample_Plate== "plate11") %>% pull(id2)
mv1.t <- mv1.t[-which(rownames(mv1.t) %in% remove),]  
pdv1 <- pdv1[-which(pdv1$id2 %in% remove),]

phe.use <- pdv1 %>% select(plate, CD8T, CD4T, NK, Bcell, Mono, Neu)
sum(substr(rownames(phe.use),1,9) == substr(rownames(mv1.t), 2, 10))

LMtest = function(meth_matrix, methcol, plate, CD8, CD4, NK, B, Mono, Neu) {
	mod = lm(meth_matrix[, methcol] ~ plate+CD8+CD4+NK+B+Mono+Neu)
	mod$residuals	
	}

res <- lapply(setNames(seq_len(ncol(mv1.t)), dimnames(mv1.t)[[2]]), LMtest, meth_matrix=mv1.t, plate=pdv1$Sample_Plate, CD8=pdv1$CD8T, CD4=pdv1$CD4T, NK=pdv1$NK, B=pdv1$Bcell, Mono=pdv1$Mono, Neu=pdv1$Neu)

setattr(res, 'class', 'data.frame')
setattr(res, "row.names", c(NA_integer_,283))
setattr(res, "names", make.names(names(res), unique=TRUE))
probelistnamesB <- names(res)
result <- t(data.table(res))
result<-data.table(result)
colnames(result) <- rownames(mv1.t)
result[, probeID := probelistnamesB]

t1_raw_regress <- result
save(t1_raw_regress, file = "t1_raw_regress.rdata")

##### Same process was done for time 2 and time 3. then we will get the t2_raw_regress.rdata and t3_raw_regress.rdata
load("t1_raw_regress.rdata")
load("t2_raw_regress.rdata")
load("t3_raw_regress.rdata")

sum(substr(colnames(t1_raw_regress),2,10)== substr(colnames(t2_raw_regress),2,10))
sum(t1_raw_regress$probeID == t2_raw_regress$probeID)
t1_raw_regress <- data.frame(t1_raw_regress)
t2_raw_regress <- data.frame(t2_raw_regress)
rownames(t1_raw_regress) <- t1_raw_regress$probeID
rownames(t2_raw_regress) <- t2_raw_regress$probeID
t1_raw_regress <- t1_raw_regress[,-282]
t2_raw_regress <- t2_raw_regress[,-282]
mv2_v1_regress <- t2_raw_regress - t1_raw_regress
save(mv2_v1_regress, file = "mv2_v1_regress.rdata")

#### same process for the time 3 vs time 1 and time 3 vs time 2, we got mv3_v1_regress.rdata and mv3_v2_regress.rdata

mv3_v1_regress <- as.matrix(mv3_v1_regress)

removeOutliers<-function(probes){
  require(matrixStats)
  if(nrow(probes) < ncol(probes)) warning("expecting probes are rows (long dataset)")
  rowIQR <- rowIQRs(probes, na.rm = T)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
  maskL <- probes < row2575[,1] - 3 * rowIQR 
  maskU <- probes > row2575[,2] + 3 * rowIQR 
  initial_NAs<-rowSums(is.na(probes))
  probes[maskL] <- NA
  removed_lower <- rowSums(is.na(probes))-initial_NAs
  probes[maskU] <- NA
  removed_upper <- rowSums(is.na(probes))-removed_lower-initial_NAs
  N_for_probe<-rowSums(!is.na(probes))
  Log<-data.frame(initial_NAs,removed_lower,removed_upper,N_for_probe)
  return(list(probes, Log))
}
system.time(OutlierResults<-removeOutliers(mv3_v1_regress))
M.val<-OutlierResults[[1]]

mv1.t <- t(M.val)

remove <- pdv1 %>% filter(Sample_Plate== "plate11") %>% pull(id2)
remove2 <- pdv1 %>% filter(Sample_Plate== "plate11") %>% pull(id)

pdv1 <- pdv1[-which(pdv1$id2 %in% remove),]

female <- pdv1 %>% filter(sex %in% c("female", "male")) %>% pull(id)
uv2_minus_uv1_female <- ti2[which(ti2$id %in% female),]
uv2_minus_uv1_female$minus <- uv2_minus_uv1_female$value
shapiro.test(uv2_minus_uv1_female$minus)
uv2_minus_uv1_female$minus <- qnorm((rank(uv2_minus_uv1_female$minus, na.last="keep") - 0.5) / sum(!is.na(uv2_minus_uv1_female$minus)))
shapiro.test(uv2_minus_uv1_female$minus)

idx <- intersect(uv2_minus_uv1_female$id, intersect(rownames(mv1.t), pdv1$id))
pdv1 <- pdv1[which(pdv1$id %in% idx),]
mv1.t <- mv1.t[which(rownames(mv1.t) %in% idx),]
uv2_minus_uv1_female <- uv2_minus_uv1_female[which(uv2_minus_uv1_female$id %in% idx),]

pdv1$Sample_Plate <- droplevels(pdv1$Sample_Plate)
str(pdv1)

sum(pdv1$id == uv2_minus_uv1_female$id)
sum(pdv1$id == rownames(mv1.t))
pdv1$age <- as.numeric(pdv1$age)

mv1.t <- as.data.frame(mv1.t)

is.recursive(mv1.t)
is.recursive(uv2_minus_uv1_female)


RLMtest = function(meth_matrix, methcol, urate, age, sex) {mod = try(rlm(urate ~ meth_matrix[, methcol]+age+sex, maxit=200))
	cf = try(coeftest(mod, vcov=vcovHC(mod, type="HC0")))

if (class(cf)=="try-error") {
  bad <- as.numeric(rep(NA, 3))
  names(bad)<- c("Estimate", "Std. Error", "Pr(>|z|)")
  bad
}
else{
  cf[2, c("Estimate", "Std. Error", "Pr(>|z|)")]
}
}

res <- mclapply(setNames(seq_len(ncol(mv1.t)), dimnames(mv1.t)[[2]]), RLMtest, meth_matrix=mv1.t, urate=uv2_minus_uv1_female$minus, age=pdv1$age, sex = pdv1$gender)

setattr(res, 'class', 'data.frame')
setattr(res, "row.names", c(NA_integer_,4))
setattr(res, "names", make.names(names(res), unique=TRUE))
probelistnamesB <- names(res)
result <- t(data.table(res))
result<-data.table(result)
result[, probeID := probelistnamesB]
setnames(result, c("BETA","SE", "P_VAL", "probeID")) # rename columns
setcolorder(result, c("probeID","BETA","SE", "P_VAL"))
result$padj <- p.adjust(result$P_VAL, method ="BH")

lambda = median(qchisq(as.numeric(as.character(result$P_VAL)),df=1,lower.tail = F),na.rm=T)/qchisq(0.5,1)
lambda 
result$lambda <- lambda 
save(result, file = "IFNg_all_v3_v1_methy_regress_trim.rdata")





############################################### Step 4 #######################################################

############################# Step 4: Mediation analysis #####################################
# we take IFN-gamma as an examples, other three markers were done with same process ###
set.seed(123)
geno <- read_delim("BCG_dosage.txt",delim = "\t", col_names = T)
colnames(geno) <- paste0("300BCG", str_pad(colnames(geno), 3, side="left", "0"))
### Phenotype 
pd <- read.csv("sample_info.csv")
colnames(pd)[1] <- "Sample_Name"
pd$id <- substr(pd$Sample_Name, 1,9)
pdv1 <- pd %>% filter(time == "v1")
pdv1$gender <- ifelse(pdv1$sex == 'female', 0, 1)
# choose the age and gender from the phenotype file 
rownames(pdv1) <- pdv1$id
cova <- pdv1 %>% dplyr::select(gender, age)
### DNAm change 
load("mv3_v1_regress.rdata")
mv3_v1_regress <- as.matrix(mv3_v1_regress)
removeOutliers<-function(probes){
  require(matrixStats)
  if(nrow(probes) < ncol(probes)) warning("expecting probes are rows (long dataset)")
  rowIQR <- rowIQRs(probes, na.rm = T)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
  maskL <- probes < row2575[,1] - 3 * rowIQR 
  maskU <- probes > row2575[,2] + 3 * rowIQR 
  initial_NAs<-rowSums(is.na(probes))
  probes[maskL] <- NA
  removed_lower <- rowSums(is.na(probes))-initial_NAs
  probes[maskU] <- NA
  removed_upper <- rowSums(is.na(probes))-removed_lower-initial_NAs
  N_for_probe<-rowSums(!is.na(probes))
  Log<-data.frame(initial_NAs,removed_lower,removed_upper,N_for_probe)
  return(list(probes, Log))
}
system.time(OutlierResults<-removeOutliers(mv3_v1_regress))
M.val<-OutlierResults[[1]]

### we only do for the sig. sites identified from DNAm change TI change EWAS studies 
load("IFNg_all_v3_v1_methy_regress_trim.rdata")
IFN.sig = result %>% filter(padj < 0.05) %>% pull(probeID)
sig =  IFN.sig## 47 sites 
DNAm = M.val[which(rownames(M.val) %in% sig),] ## 14, 280 

overlap <- intersect(intersect(colnames(geno),colnames(DNAm)),rownames(cova))
length(overlap) ## 244

DNAm2 <- data.frame(DNAm[,which(colnames(DNAm) %in% overlap)])
head(DNAm2)
colnames(DNAm2) = substr(colnames(DNAm2), 2, 10) ## 244
## match the genotype id with uv and cova
geno2 <- as.matrix(geno[,which(colnames(geno) %in% overlap)])
rownames(geno2) <- geno$'300BCG0ID' ## 4296841     244
cova2 <- cova[which(rownames(cova) %in% overlap),]
head(cova2)

sum(rownames(cova2) == colnames(DNAm2))
sum(rownames(cova2) == colnames(geno2))

extract_mediation_summary <- function (x) { 

  clp <- 100 * x$conf.level
  isLinear.y <- ((class(x$model.y)[1] %in% c("lm", "rq")) || 
                   (inherits(x$model.y, "glm") && x$model.y$family$family == 
                      "gaussian" && x$model.y$family$link == "identity") || 
                   (inherits(x$model.y, "survreg") && x$model.y$dist == 
                      "gaussian"))

  printone <- !x$INT && isLinear.y

  if (printone) {

    smat <- c(x$d1, x$d1.ci, x$d1.p)
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))

    rownames(smat) <- c("ACME", "ADE", "Total Effect", "Prop. Mediated")

  } else {
    smat <- c(x$d0, x$d0.ci, x$d0.p)
    smat <- rbind(smat, c(x$d1, x$d1.ci, x$d1.p))
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$z1, x$z1.ci, x$z1.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    smat <- rbind(smat, c(x$n1, x$n1.ci, x$n1.p))
    smat <- rbind(smat, c(x$d.avg, x$d.avg.ci, x$d.avg.p))
    smat <- rbind(smat, c(x$z.avg, x$z.avg.ci, x$z.avg.p))
    smat <- rbind(smat, c(x$n.avg, x$n.avg.ci, x$n.avg.p))

    rownames(smat) <- c("ACME (control)", "ACME (treated)", 
                        "ADE (control)", "ADE (treated)", "Total Effect", 
                        "Prop. Mediated (control)", "Prop. Mediated (treated)", 
                        "ACME (average)", "ADE (average)", "Prop. Mediated (average)")

  }

  colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep = ""), 
                      paste(clp, "% CI Upper", sep = ""), "p-value")
  smat

} ## this is the function we used for extract the summary results

statistics = read.table("TI_rank_GWAS_p0.05.tsv", head = T)
cis = read.table("DNAmchange_cis_GWAS_p0.05.tsv", head = T)

statistics2 = statistics %>% filter(gene == "FC_IFNg_W3_13")
cis2 = cis %>% filter(gene %in% IFN.sig) # 4102
snp = intersect(statistics2$snps, cis2$snps) 
length(snp)
### 4 snps with sig. result for both DNAm change and TI

TI = statistics2 %>% filter(snps %in% snp)
DNAmchange = cis2 %>% filter(snps %in% snp)

### 244 samples used for the mediation analysis 
# genotype: geno2; DNAm change: DNAm2; TI: trained_immunity; phenotype: cova2

trained_immunity <- data.frame(read_excel("ti_value.xlsx")) ### 278 samples 
rownames(trained_immunity) <- substr(trained_immunity$ID2, 1, 9)
trained_immunity <- trained_immunity %>% dplyr::select(FC_IL1b_T3_13, FC_IL6_T3_13, FC_TNF_T3_13, FC_IFNg_W3_13)
dim(trained_immunity)
shapiro.test(trained_immunity$FC_IFNg_W3_13)
trained_immunity$FC_IFNg_W3_13 = qnorm((rank(trained_immunity$FC_IFNg_W3_13, na.last="keep") - 0.5) / sum(!is.na(trained_immunity$FC_IFNg_W3_13)))

sum(rownames(cova2) == colnames(geno2))
sum(rownames(trained_immunity) == colnames(geno2))
sum(colnames(DNAm2) == colnames(geno2))

#### DNAm change as the mediator, TI as the outcome
result = NULL
a.res = NULL
for(i in 1:269){
	  treat = DNAmchange[i,1]
    mediator = DNAmchange[i,2]
    outcome = "FC_IFNg_W3_13"
    data <- na.omit(data.frame(treat = geno2[treat,], mediator = as.numeric(DNAm2[mediator,]), outcome = trained_immunity[,outcome], age = cova2[,"age"], gender = cova2[,"gender"]))
    samplesize = dim(data)[1]
    med.fit <- lm(mediator ~ treat + age + gender, data = data)
    out.fit <- lm(outcome ~ treat + mediator + age + gender, data = data)
    med.out <- mediate(med.fit, out.fit, treat = "treat", mediator = "mediator", sims = 1000)
    a = extract_mediation_summary(med.out)
    ACME=a[1,4]
    ADE=a[2,4]
    total=a[3,4]
    prop=a[4,4]
    a = data.frame(a)
    a$treat = treat
    a$mediator = mediator
    a$outcome = outcome
    a.res = data.frame(rbind(a.res, a))
    res = data.frame(treat = treat, mediator = mediator, outcome = outcome, ACME = ACME, ADE = ADE, total = total, prop = prop, samplesize = samplesize)
    result = rbind(result, res)
}

saveRDS(result, file = "mediation_IFN_direction1.rds")
saveRDS(a.res, file = "all_output_mediation_IFN_direction1.rds")



#### TI as the mediator, DNA methylation change as the outcome
result = NULL
a.res = NULL
for(i in 1:269){
	treat = DNAmchange[i,1]
    outcome = DNAmchange[i,2]
    mediator = "FC_IFNg_W3_13"
    data <- na.omit(data.frame(treat = geno2[treat,], outcome = as.numeric(DNAm2[outcome,]),mediator = trained_immunity[,mediator], age = cova2[,"age"], gender = cova2[,"gender"]))
    samplesize = dim(data)[1]
    med.fit <- lm(mediator ~ treat + age + gender, data = data)
    out.fit <- lm(outcome ~ treat + mediator + age + gender, data = data)
    med.out <- mediate(med.fit, out.fit, treat = "treat", mediator = "mediator", sims = 1000)
    a = extract_mediation_summary(med.out)
    ACME=a[1,4]
    ADE=a[2,4]
    total=a[3,4]
    prop=a[4,4]
    a = data.frame(a)
    a$treat = treat
    a$mediator = mediator
    a$outcome = outcome
    a.res = data.frame(rbind(a.res, a))
    res = data.frame(treat = treat, mediator = mediator, outcome = outcome, ACME = ACME, ADE = ADE, total = total, prop = prop, samplesize = samplesize)
    result = rbind(result, res)
}


saveRDS(result, file = "mediation_IFN_direction2.rds")
saveRDS(a.res, file = "all_output_mediation_IFN_direction2.rds")








