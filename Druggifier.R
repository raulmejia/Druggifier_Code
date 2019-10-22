##############################################################################
### Druggifier 
## Author: Raúl Alejandro Mejía Pedroza github: https://github.com/raulmejia
##############################################################################
###########################################
## Data given by the user
###########################################
args <- commandArgs(trailingOnly = TRUE)
Path_to_your_Matrix<-args[1] # The path to your matrix
Path_to_myPhenoData <- args[2] # The path to your Feature data
Path_of_Code<-args[3] # The path to your code
Path_of_Results <-args[4] # # where do you want to save your results?
Labels <-args[5] # Label for your results, usually where the expression data came from
log2transformation <- args[6] # Do you want log2transformation for your expression matrix ?
Name_of_the_column_that_contain_your_control_or_reference <- args[7]
Name_of_your_control_or_reference <-args[8]
filter_by_this_characteristic_of_the_pheatmap <- args[9]
character_value_to_use_to_filter <- args[10]
Labels_to_your_values <- args[11]
path_to_KEGG <- args[12]
Stabilizing_pathifier <- args[13]
Filter_value_pathifier <- args[14] 
Matrix_Top_N_pathifier <- args[15]
DiffExp_method <- args[16] # "limma" or "DESeq2"
# Path_to_your_Matrix<-c("../Data/METABRIC/joined_indicator_METABRIC.txt") ; Path_to_myPhenoData <- c("../Results/Clusterig/METABRIC/Heatmaps/Lehmann_Only_ABCS/PhenoData_Clusters_Heatmap_METABRIC_Lehmann_subtypes_ABC_transporters.tsv") ; Path_of_Code<-c("./") ; Path_of_Results<-c("../Results/Druggifier/METABRIC/") ; Labels <- "METABRIC_Cluster_ABCs_on_Lehmann_subtypes" ; log2transformation <- "log2-transformation-no" ;Name_of_the_column_that_contain_your_control_or_reference <-"Labels"; Name_of_your_control_or_reference <-"Control"; filter_by_this_characteristic_of_the_pheatmap <- "cutree_rowcenter_Heatmap_METABRIC_Lehmann_subtypes_ABC_transporters" ; character_value_to_use_to_filter<-c("6","5","1"); Labels_to_your_values <-c("CL1","CL2","CL3") ; path_to_KEGG <- "../KEGG/" ; Stabilizing_pathifier <- 5; Filter_value_pathifier <- 3 ; Matrix_Top_N_pathifier <- 50 ; DiffExp_method <- "limma" ; mylfctreshold <- 0.6; myp.adj <- 0.05 ; mycoresavaiable <- 7
# Path_to_your_Matrix<-c("../Data/TCGA/expMatrix_TCGA_cBioPortal_no_males_withindicator.txt") ; Path_to_myPhenoData <- c("../Results/Clusterig/TCGA/Heatmaps/Lehmann_Only_ABCS/PhenoData_Clusters_Heatmap_TCGA_Lehmann_subtypes_ABC_transporters.tsv") ; Path_of_Code<-c("./") ; Path_of_Results<-c("../Results/Druggifier/TCGA/") ; Labels <- "TCGA_Cluster_ABCs_on_Lehmann_subtypes" ; log2transformation <- "log2-transformation-yes" ;Name_of_the_column_that_contain_your_control_or_reference <-"Labels"; Name_of_your_control_or_reference <-"Control"; filter_by_this_characteristic_of_the_pheatmap <- "cutree_rowcenter_Heatmap_TCGA_Lehmann_subtypes_ABC_transporters" ; character_value_to_use_to_filter<-c("6","4","2","3","1"); Labels_to_your_values <-c("CL1","CL2","CL2-weak","CL3","CL3-weak") ; path_to_KEGG <- "../KEGG/" ; Stabilizing_pathifier <- 5; Filter_value_pathifier <- 3.75 ; Matrix_Top_N_pathifier <- 50 ; DiffExp_method <- "DESeq2" ; mylfctreshold <- 0.6; myp.adj <- 0.05 ; mycoresavaiable <- 7 

# Path_to_your_Matrix<-c("../Data/TCGA/expMatrix_TCGA_cBioPortal_no_males_withindicator.txt") ; Path_to_myPhenoData <- c("../Results/Lehmann-STROMA4/TCGA/PhenoData_with_Lehmann_Subt_and_properties_TCGA-TripleNeg.tsv") ; Path_of_Code<-c("./") ; Path_of_Results<-c("../Results/Clusterig/TCGA/Heatmaps/Lehmann_Only_ABCS/") ; Labels <- "Heatmap_TCGA_Lehmann_subtypes_ABC_transporters" ; log2transformation <- "log2-transformation-yes"; Path_to_your_gene_list <- c("../Data/List_of_Genes/all_myABC.tsv") ; filter_by_this_characteristic_of_the_pheatmap <-;
# Path_to_your_Matrix <- choose.files() # choose.dir()
###############################################################################
### Installing and/or loading required packages
###############################################################################
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("Hiiragi2013")) {
  BiocManager::install("Hiiragi2013", dependencies = TRUE)
  library(Hiiragi2013)
}
if (!require("pheatmap")) {
  BiocManager::install("pheatmap", dependencies = TRUE)
  library(pheatmap)
}
if (!require("dplyr")) {
  BiocManager::install("dplyr", dependencies = TRUE)
  library(dplyr)
}
if (!require("som")) {
  install.packages("som", dependencies = TRUE)
  library(som)
}

###########################################################
##### Converting the input variables to the rigth type
###########################################################
Stabilizing_pathifier <- as.numeric(Stabilizing_pathifier)
Filter_value_pathifier <- as.numeric(Filter_value_pathifier)
Matrix_Top_N_pathifier <- as.numeric(Matrix_Top_N_pathifier)

#############################################################
### Reading the data
#############################################################
Mymatrix_raw <-read.table(Path_to_your_Matrix,header=TRUE,row.names = 1)
allPhenoData <-read.table(Path_to_myPhenoData ,header=TRUE, sep = "\t", quote = "")

#############################################################
### Log 2 transformation or not 
#############################################################
Mymatrix <- Mymatrix_raw

if(log2transformation == "log2-transformation-yes"){
  if( length(grep("NORMAL",rownames(Mymatrix_raw))) > 0   ){
    Mymatrix <- Mymatrix_raw[-grep("NORMAL",rownames(Mymatrix_raw)) , ] # Erasing the indicator 1 (Normal) for the log2 kahorik 
  }
  Mymatrix <- log2(Mymatrix+1)  
  Mymatrix <- rbind(Mymatrix_raw[grep("NORMAL",rownames(Mymatrix_raw)) , ] , Mymatrix)
  print("your data were log-2 transformed by this program")
  Mymatrix[1:6,1:6]
}else{
  if(log2transformation == "log2-transformation-no"  ){
    print("your data will NOT log-2 transformed by this program")
    Mymatrix[1:6,1:6]
  }else{
    print("Please for the argument [log2transformation] type 'log2-transformation-yes' or 'log2-transformation-no'")
  }
}
# Take this out of here and make a new module
###############################################
## Creating paths to organize your results according to your characteristics and Labels
###############################################
Path_of_Results <- normalizePath(Path_of_Results)
dir.create(Path_of_Results , recursive = TRUE )
list_of_Path_of_Results <- list()
for(w in  Labels_to_your_values){
  list_of_Path_of_Results[[w]] <- paste0(Path_of_Results,"/",w)
}

#############################################################
### Submatrix of controls
#############################################################
PhenoData_NoNAs_inReferencescolumn <- allPhenoData[ !is.na(allPhenoData[,Name_of_the_column_that_contain_your_control_or_reference]),] # PhenoData without NA in column that contain the controls
indicator_reference <- PhenoData_NoNAs_inReferencescolumn[,Name_of_the_column_that_contain_your_control_or_reference ] == Name_of_your_control_or_reference # indicator to extract your controls
Mymatrix_reforCtrl <- Mymatrix[ ,rownames(PhenoData_NoNAs_inReferencescolumn)[indicator_reference] ] # Submatrix of controls
# Pending to use eSet for this task

#############################################################
### Submatrix of your election
#############################################################
PhenoDataNoNas_inthat_characteristic <- allPhenoData[ !is.na(allPhenoData[,filter_by_this_characteristic_of_the_pheatmap]),] # Phenodata with no NAs in your characteristic to filter
list_of_submatrices_withcontrols <- list() 
for(k in 1:length(character_value_to_use_to_filter)){ # Creating submatrices Controls pasted with the submatrix of samples that fulfill each our your values to filter
  indicator_mycurrent_value <- as.character(PhenoDataNoNas_inthat_characteristic[,filter_by_this_characteristic_of_the_pheatmap]) == character_value_to_use_to_filter[k]  
  Rownames_with_your_feature <- rownames(PhenoDataNoNas_inthat_characteristic)[indicator_mycurrent_value]
  list_of_submatrices_withcontrols[[k]] <- cbind(Mymatrix_reforCtrl,Mymatrix[,Rownames_with_your_feature]) # pasting controls submatrix with the submatrix of samples of that value
}
names(list_of_submatrices_withcontrols) <- paste0( filter_by_this_characteristic_of_the_pheatmap,"-",  Labels_to_your_values )

#############################################################
### Pathifier
#############################################################
###### Download KEGG
source(paste0(Path_of_Code,"UptodateKEGG_R_function.R"))
KEGGdb_path <- UptodateKEGG_R_function(paste0(Path_of_Results,"/KEGG/"))

###### Runing_Pathifier
source(paste0(Path_of_Code,"Pathifier_As_R_function.R"))
list_Pathifiers<- list()
ptm <- proc.time()
for( k in 1:length( names( list_of_submatrices_withcontrols))){
  dir.create(paste0(list_of_Path_of_Results[[k]],"/","Pathifier"), recursive = TRUE)
  list_Pathifiers[[k]] <- Pathifier_As_R_Function(as.matrix(list_of_submatrices_withcontrols[[k]]) , KEGGdb_path , Path_of_Code,
                                                  paste0(list_of_Path_of_Results[[k]],"/","Pathifier"), 
                                                  paste0(Labels,"---",Labels_to_your_values[k]) ,
                                                  Stabilizing_pathifier, Filter_value_pathifier, Matrix_Top_N_pathifier )
  print(k)
}
Pathifier_list_time <-proc.time() - ptm
names( list_Pathifiers ) <-  paste0( "Pathifier",filter_by_this_characteristic_of_the_pheatmap,"-", character_value_to_use_to_filter)

################################################################
######### DGEs through DESeq2
################################################################
# "limma/DESeq2"

Path_of_Results <- normalizePath(Path_of_Results)
Path_of_Results_DEG <- paste0(Path_of_Results,"/","DEG")

Exp_Mat <- list_of_submatrices_withcontrols[[1]]

pos_NORMAL_indicator <-  grep("NORMAL",rownames(list_of_submatrices_withcontrols[[1]]))
Thelabels <- gsub(1,"Control",list_of_submatrices_withcontrols[[1]][pos_NORMAL_indicator,]) ; Thelabels <-  gsub("0",Labels_to_your_values[1],Thelabels)
MyLabels <- data.frame( Thelabels ,
                        row.names= colnames(list_of_submatrices_withcontrols[[1]]),
                        stringsAsFactors = FALSE); colnames(MyLabels)<- "Labels"
MyLabels[,1] <- as.character(MyLabels[,1])

Labels_DGE <-paste0("DGE_",Labels)
Label_for_results <-  Labels_DGE

mylfctreshold ; myp.adj ; mycoresavaiable

if(DiffExp_method == "limma"){
  
}else{
  if(DiffExp_method == "DESeq2"  ){
    DESeq2_list <- DESeq2_with_log2_transformed_data_only(
      Exp_Mat , MyLabels , Path_of_Code,
      Path_of_Results_DEG , Label_for_results, mylfctreshold, 
      myp.adj , mycoresavaiable)
  }else{
    print("Please for the argument [DiffExp_method] type limma or DESeq2")
  }
}

################################################################
######### DGEs through DESeq2
################################################################

