################################################################################
## This program calculates the Differential expressed genes through limma
## Only work for cases and controls the first row of your input matrix should 
## depict your controls (with ones "1") and cases with zeros "0"
## Author: Raúl Alejandro Mejía Pedroza github: https://github.com/raulmejia
################################################################################
### Data given by the user
################################################################################
#args <- commandArgs(trailingOnly = TRUE)
#exp_mat_path <-args[1]
# exp_mat_path <-c("../Data/toyMETABRIC_Controls_LumA.txt") ; Path_of_Code<-c("./") ; Path_of_Results<-c("../Results/DEG/METABRIC/") ; mypvalue = 0.05 ; mylfc = 0.6 ; Label_for_results <-"METABRIC_some_LumAs_vs_Controls_"  
# exp_mat_path <- c("../Data/joined_indicator_METABRIC.txt")  
#Path_of_Code <-args[2]
#Path_of_Results <-args[3]
#mylfc <-args[4]
#mypvalue <- args[5]
#Label_for_results <- args[6]
#Filter_value <- args[ ] # Filter low value genes # METABRIC = 3

###############################################################################
### Installing and/or loading required packages
###############################################################################
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("illuminaHumanv3.db")) {
  BiocManager::install("illuminaHumanv3.db", ask =FALSE)
  library("illuminaHumanv3.db")
}
if (!require("hgu133plus2.db")) {
  BiocManager::install("hgu133plus2.db", ask =FALSE)
  library("hgu133plus2.db")
}
if (!require("annotate")) {
  BiocManager::install("annotate", ask =FALSE)
  library("annotate")
}
if (!require("limma")) {
  BiocManager::install("limma", ask =FALSE)
  library("limma")
}
if (!require("gProfileR")) {
  BiocManager::install("gProfileR", ask =FALSE)
  library("gProfileR")
}
if (!require("org.Hs.eg.db")) {
  BiocManager::install("org.Hs.eg.db", ask =FALSE)
  library("org.Hs.eg.db")}

#Exp_Mat  
#MyLabels  
#Path_of_Code
#Path_of_Results <- Path_of_Results_DEG 
#Label_for_results 
#mylfc <- mylfctreshold
#mypvalue <- myp.adj

DGE_limma <- function(Exp_Mat , MyLabels ,Path_of_Code, Path_of_Results, Label_for_results,
                      mylfc, mypvalue  ){
##########################################
###### Reading the data 
##########################################
# Load the data frame (expression matrix)
#Exp_Mat <- read.table(exp_mat_path ,header=TRUE , sep="\t")

############################################
### Preprocessing some data given y the user
############################################
dir.create(Path_of_Results, recursive = TRUE)
mypvalue <- as.numeric(mypvalue) # coercing to numeric 
mylfc <- as.numeric( mylfc) # coercing to numeric 

###########################################
####### Hadling the labels
###########################################
dfLabels <- data.frame(Exp_Mat[1,])
dfLabels <- gsub(1, "Control", dfLabels)
dfLabels <- gsub(0, "Case", dfLabels)
samples <- relevel(as.factor(dfLabels), "Control")

#############################################
# Erasing the row of normals
#############################################
normals_pos <- which(rownames(Exp_Mat) %in% "NORMALS")
Exp_Mat <- Exp_Mat[-normals_pos ,]

######################
 # set up the experimental design
 design <- model.matrix(~0 + samples)
 colnames(design) <- c("Control", "Case") # Como level el segundo 
 
 # fit the linear model to the your expression set "eset"
 eset<-Exp_Mat
 fit <- lmFit(eset, design)
 # set up a contrast matrix to compare tissues v cell line
 contrast.matrix <- makeContrasts( Case_Controls = Case - Control,  levels=design)
 # Now the contrast matrix is combined with the per-probeset linear model fit.
 huvec_fits <- contrasts.fit(fit, contrast.matrix)
 huvec_ebFit <- eBayes(huvec_fits)
 
#####################################
##DGEs with  Affy ids no anotation###
#####################################
 results_only_affys<-list(length(colnames(contrast.matrix)))
 for(i in 1:length(colnames(contrast.matrix))){
   DGEs<-topTable(huvec_ebFit, coef=colnames(contrast.matrix)[i], number=10000, p.value = mypvalue, lfc = mylfc)
   DGEs<-DGEs[order(rownames(DGEs)),]
   # Saving the matrix
   write.table(DGEs[,c(1,5)], paste(Path_of_Results ,"DGE_",Label_for_results, mypvalue, mylfc, colnames(contrast.matrix)[i], "_ID_with_ILMN.txt", sep=""), sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)
   # Saving with no ILMN_ identifiers
   write.table(DGEs[,c(1,5)][!grepl("ILMN_",rownames(DGEs[,c(1,5)])),], paste(Path_of_Results ,"DGE_",Label_for_results, mypvalue, mylfc, colnames(contrast.matrix)[i], "_ID_no_ILMN_.txt", sep=""), sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)
   results_only_affys[[i]]<-DGEs
 }
 saveRDS(results_only_affys,file= paste(Path_of_Results ,"DGE_",Label_for_results, mypvalue, mylfc, colnames(contrast.matrix)[i], "_ID_with_ILMN.RDS", sep=""))
 names(results_only_affys) <- Label_for_results
 return(results_only_affys[[1]])
}