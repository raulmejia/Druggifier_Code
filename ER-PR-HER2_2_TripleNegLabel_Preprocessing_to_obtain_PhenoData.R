# This program receives a data frame and add a "triple-neg" Label according
# to Her2,ER,PR information inside that data frame
#################################################
### Loading and installing required libraries
#################################################
args <- commandArgs(trailingOnly = TRUE)
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
###########################################
### Data given by the user
###########################################
Path_to_PhenoData_METABRIC <-c("../Data/METABRIC/PhenoDataDiscovSet.tsv") ; Path_of_Results_METABRIC <-c("../Data/METABRIC/") ; Title_of_Results_METABRIC <- c("METABRIC-PhenoD-triple-neg-Label") ; PAM50_labels_METABRIC<- c("../Data/METABRIC/Labels_Ctrl_and_NL_Recal_separated_METABRIC.txt")
Path_to_PhenoData_TCGA <-c("../Data/TCGA/brca_tcga_pub2015_clinical_data.tsv") ; Path_of_Results_TCGA <-c("../Data/TCGA/") ; Title_of_Results_TCGA <- c("TCGA-PhenoD-triple-neg-Label") ; PAM50_labels_TCGA_Path<- c("../Data/TCGA/Labels_Controls_and_Normal_separated_TCGA.txt")
# Path_to_your_PhenoData <- args[1]
# Path_of_Results <- args[2]
# Title_of_Results <- args[3]
############################################
### Reading the data METABRIC
############################################
PhenoDataMETABRIC <- read.table(Path_to_PhenoData_METABRIC , header=TRUE, sep= "\t", quote = "")
head(PhenoDataMETABRIC)

########################################
###### Labeling triple Negatives
########################################
PhenoDataMETABRIC$TripleNeg <- rep(NA,dim(PhenoDataMETABRIC)[1],)
ERstatus <- PhenoDataMETABRIC$ER.Status == "Negative"
HER2status <- PhenoDataMETABRIC$HER2.Status == "Negative"
TripleNeg <- (ERstatus & HER2status) 
PhenoDataMETABRIC$TripleNeg <- TripleNeg 

########################################
###### Adding PAM 50
########################################
PAMLabMET<-read.table( PAM50_labels_METABRIC )

rownames(PAMLabMET)[PAMLabMET == "Control"]
dim(PhenoDataMETABRIC)
PhenocontrolsMET <- matrix(rep(NA, length(rownames(PAMLabMET)[PAMLabMET == "Control"]) * dim(PhenoDataMETABRIC)[2]),
                      nrow= length(rownames(PAMLabMET)[PAMLabMET == "Control"]), 
                      dimnames = list(rownames(PAMLabMET)[PAMLabMET == "Control"], colnames(PhenoDataMETABRIC)))

PhenoDataMETABRIC_Ctrls <- rbind(PhenocontrolsMET, PhenoDataMETABRIC)
PhenoDataMETABRIC_Ctrls_ordered <- PhenoDataMETABRIC_Ctrls[order(match( rownames(PhenoDataMETABRIC_Ctrls)  ,
                                                                        rownames(PAMLabMET))),] # ordering according Labels PAM50 MET

PhenoDataMETABRIC_PAM50added <- cbind(PhenoDataMETABRIC_Ctrls_ordered,PAMLabMET)

write.table( PhenoDataMETABRIC_PAM50added , file = paste0( Path_of_Results_METABRIC,"METABRIC_PhenoData_TripleNeg_PAM50.tsv"), sep="\t", quote=FALSE , row.names= TRUE, col.names= TRUE )


## -------------TCGA -----------------##
############################################
### Reading the data TCGA
############################################
PhenoDataTCGA <- read.table(Path_to_PhenoData_TCGA , header=TRUE, sep= "\t", quote = "")
head(PhenoDataTCGA)
row.names(PhenoDataTCGA) <- gsub("-",".",PhenoDataTCGA$Sample.ID)

PAMLabTCGA<-read.table( PAM50_labels_TCGA_Path )
###########################################
### Same number of samples (remember that you erase the males)
##########################################
positions <- which( rownames(PhenoDataTCGA) %in% rownames(PAMLabTCGA))
PhenoDataTCGA <- PhenoDataTCGA[positions,]

########################################
###### Labeling triple Negatives
########################################
PhenoDataTCGA$TripleNeg <- rep(NA,dim(PhenoDataTCGA)[1],)
HER2status_TCGA <- PhenoDataTCGA$IHC.HER2 == "Negative"
ERstatus_TCGA <- PhenoDataTCGA$ER.Status.By.IHC == "Negative"
TripleNeg_TCGA <- (HER2status_TCGA & ERstatus_TCGA)
PhenoDataTCGA$TripleNeg <- TripleNeg_TCGA

########################################
###### Adding PAM 50
########################################


rownames(PAMLabTCGA)[PAMLabTCGA == "Control"]
dim(PhenoDataTCGA)
PhenocontrolsTCGA <- matrix(rep(NA, length(rownames(PAMLabTCGA)[PAMLabTCGA == "Control"]) * dim(PhenoDataTCGA)[2]),
                           nrow= length(rownames(PAMLabTCGA)[PAMLabTCGA == "Control"]), 
                           dimnames = list(rownames(PAMLabTCGA)[PAMLabTCGA == "Control"], colnames(PhenoDataTCGA)))

PhenoDataTCGA_Ctrls <- rbind(PhenocontrolsTCGA, PhenoDataTCGA)
PhenoDataTCGA_Ctrls_ordered <- PhenoDataTCGA_Ctrls[order(match( rownames(PhenoDataTCGA_Ctrls)  ,
                                                                        rownames(PAMLabTCGA))),] # ordering according Labels PAM50 TCGA
rownames(PhenoDataTCGA_Ctrls_ordered) == rownames(PAMLabTCGA)

PhenoDataTCGA_PAM50added <- cbind(PhenoDataTCGA_Ctrls_ordered,PAMLabTCGA)

write.table(PhenoDataTCGA_PAM50added , file = paste0( Path_of_Results_TCGA,"TCGA_PhenoData_TripleNeg_PAM50.tsv"), sep="\t", quote=FALSE , row.names= TRUE, col.names= TRUE )
