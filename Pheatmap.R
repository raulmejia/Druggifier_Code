################################################################################
# 
## Author: Raúl Alejandro Mejía Pedroza github: https://github.com/raulmejia
################################################################################
###########################################
## Data given by the user
###########################################
args <- commandArgs(trailingOnly = TRUE)
Path_to_your_Matrix<-args[1] # The path to your matrix
Path_to_myPhenoData <- args[2] # The path to your Feature data
Path_of_Code<-args[3] # The path to your code
Path_of_Results <-args[4] # # where do you want to save your results?
Labels <-args[5] # Label for your results
log2transformation <- args[6] # Do you want log2transformation for your expression matrix ?
# Path_to_your_Matrix<-c("../Data/METABRIC/joined_indicator_METABRIC.txt") ; Path_to_myPhenoData <- c("../Results/Lehmann-STROMA4/METABRIC/PhenoData_with_Lehmann_Subt_and_properties_METABRIC-TripleNeg.tsv") ; Path_of_Code<-c("./") ; Path_of_Results<-c("../Results/Clusterig/METABRIC/Heatmaps/Lehmann_Only_ABCS/") ; Labels <- "Heatmap_METABRIC_Lehmann_subtypes_ABC_transporters" ; log2transformation <- "log2-transformation-no"; Path_to_your_gene_list <- c("../Data/List_of_Genes/all_myABC.tsv") ; filter_by_this_characteristic_of_the_pheatmap <-"TripleNeg"; character_value_to_use_to_filter<-"TRUE"
# Path_to_your_Matrix<-c("../Data/TCGA/expMatrix_TCGA_cBioPortal_no_males_withindicator.txt") ; Path_to_myPhenoData <- c("../Results/Lehmann-STROMA4/TCGA/PhenoData_with_Lehmann_Subt_and_properties_TCGA-TripleNeg.tsv") ; Path_of_Code<-c("./") ; Path_of_Results<-c("../Results/Clusterig/TCGA/Heatmaps/Lehmann_Only_ABCS/") ; Labels <- "Heatmap_TCGA_Lehmann_subtypes_ABC_transporters" ; log2transformation <- "log2-transformation-yes"; Path_to_your_gene_list <- c("../Data/List_of_Genes/all_myABC.tsv") ; filter_by_this_characteristic_of_the_pheatmap <-"TripleNeg"; character_value_to_use_to_filter<-"TRUE"
# Path_to_your_Matrix <- choose.files() # choose.dir()
Path_to_your_gene_list <- args[7]
filter_by_this_characteristic_of_the_pheatmap <- args[8]
character_value_to_use_to_filter <- args[9]
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


#####################################################
### Reading your list of genes
######################################################
gene_list_raw <- read.table(Path_to_your_gene_list, header = FALSE, sep="\t", stringsAsFactors = FALSE) 
gene_list <- unique(gene_list_raw)
rownames(gene_list) <-  gene_list[,1]; colnames(gene_list) = "genes"

#####################################################
### Reading your Pheatmap
######################################################
allPhenoData <-read.table(Path_to_myPhenoData ,header=TRUE, sep = "\t", quote = "")
indicator_noNas_PhenoD <- !is.na(allPhenoData[,filter_by_this_characteristic_of_the_pheatmap])
myPhenoData <- allPhenoData[indicator_noNas_PhenoD,]
indicator_myvalue_and_characteristic <- myPhenoData[,filter_by_this_characteristic_of_the_pheatmap] == character_value_to_use_to_filter
myPhenoData <- myPhenoData[indicator_myvalue_and_characteristic,]

#############################################################
### Reading the data and creating the folder for the results
#############################################################
dir.create(Path_of_Results , recursive = TRUE )
Fullmatrix <-read.table(Path_to_your_Matrix,header=TRUE,row.names = 1)
indicator <- which(rownames(Fullmatrix) %in% rownames(gene_list))
M.matrix <- Fullmatrix[indicator,] # Filtering by genes
M.matrix <- M.matrix[,rownames(myPhenoData)] # Filtering by your condition

# transforming your matrix to log2
if( log2transformation == "log2-transformation-no"){
  print("Your expression matrix has NOT being log2 Transformed by this program")
  M.matrix[1:5,1:5]
}else{
  if( log2transformation == "log2-transformation-yes" ){
    print("Your expression matrix has being log2 Transformed by this program")
    M.matrix <- log(M.matrix+1)
    M.matrix[1:5,1:5]
  }else{
    print("ERROR: Please type one option 'log2 transformation no' or 'log2 transformation no' ")
    exit("no")
  }  
}


##################################
### Building ExpressionSet Object
##################################
## PhenoData 
#myPhenoData <-read.table(Path_to_myPhenoData ,header=TRUE,row.names = 1)
#myPhenoData <- as.data.frame(myPhenoData )
My_pData_Anotdf <- AnnotatedDataFrame( data= myPhenoData)
validObject( My_pData_Anotdf )
myExpressionSet <- ExpressionSet(as.matrix( M.matrix ), phenoData= My_pData_Anotdf  ) ## ExpresisionSet

########################
### Pheatmap
########################
rowCenter = function(x) { x - rowMeans(x) } # Defining a function that only make the 
PDSmatrix_zscores <- som::normalize( Biobase::exprs( myExpressionSet ) , byrow = TRUE)
colnames(PDSmatrix_zscores) <- colnames( Biobase::exprs( myExpressionSet ))

mymin <- min(exprs( myExpressionSet ))
mymax <- max(exprs( myExpressionSet ))
mydim <- dim(rowCenter(Biobase::exprs( myExpressionSet )))[2]

pdf( file = paste0(Path_of_Results,Labels,"_Z-scores_Euclidean_wardD2.pdf"))
pheatmap(  PDSmatrix_zscores,
           show_rownames = TRUE, show_colnames = FALSE,
           method = c("euclidean"),
           clustering_method = "ward.D2",
           fontsize = 5,
           main = gsub("_"," ",Labels),
           #breaks = seq(mymin, mymax, length = mydim+1 ),
           annotation_col =
             pData(myExpressionSet)[,c("Labels","D.stroma.property","B.stroma.property","T.stroma.property","E.stroma.property","MSL.property","M.property","LAR.property","IM.property","BL1.property","BL2.property")],
           annotation_colors = list(
             #D.stroma.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #B.stroma.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #T.stroma.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #E.stroma.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #MSL.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #M.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #LAR.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #IM.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #BL1.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #BL2.property = c(`1` = "black", `0` = "gray", `-1` = "white")
           ),
           cutree_rows = 6,
           cutree_cols = 6
)
dev.off()


pdf( file = paste0(Path_of_Results,Labels,"_Row_Center_Euclidean_wardD2.pdf") )
#png(file = paste0(Path_of_Results,Labels,"_Row_Center_Euclidean_wardD2.png"),width = 480, height = 480, pointsize = 2)
pheatmap(  rowCenter(Biobase::exprs( myExpressionSet )) ,
           main = gsub("_"," ",Labels),
           show_rownames = TRUE, show_colnames = FALSE,
           method = c("euclidean"),
           clustering_method = "ward.D2",
           fontsize = 5,
           #breaks = seq(mymin, mymax, length = mydim+1 ),
           annotation_col =
             pData(myExpressionSet)[,c("Labels","D.stroma.property","B.stroma.property","T.stroma.property","E.stroma.property","MSL.property","M.property","LAR.property","IM.property","BL1.property","BL2.property")],
           annotation_colors = list(
             #D.stroma.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #B.stroma.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #T.stroma.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #E.stroma.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #MSL.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #M.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #LAR.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #IM.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #BL1.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #BL2.property = c(`1` = "black", `0` = "gray", `-1` = "white")
           ),
           cutree_rows = 6,
           cutree_cols = 6
)
dev.off()

########################################################
## Gettting the clusters and adding up to the PhenoData
########################################################
dist_zscores <- dist(t(PDSmatrix_zscores) , "euclidean")
myhclust_zscores <- hclust( dist_zscores , method="ward.D2")
cutree_zscores <- cutree(myhclust_zscores,k=6)
cutree_zscoresDf <- data.frame(cutree_zscores,Sample_ID=names(cutree_zscores))
colnames(cutree_zscoresDf)[1] <- paste0("cutree_zscores_",Labels)
myPhenoData_cutree_zscores <- merge( x=myPhenoData, y=cutree_zscoresDf , by="Sample_ID", all=TRUE)
rownames(myPhenoData_cutree_zscores) <- rownames(myPhenoData)

dist_rowcenter <-  dist(t(rowCenter(Biobase::exprs( myExpressionSet ))), "euclidean")
myhclust_rowcenter <- hclust(dist_rowcenter, method="ward.D2")
cutree_rowcenter <- cutree(myhclust_rowcenter , k=6)
cutree_rowcenterDf <- data.frame( cutree_rowcenter, Sample_ID=names(cutree_rowcenter))
colnames(cutree_rowcenterDf)[1] <- paste0("cutree_rowcenter_",Labels)
myPhenoData_cutree_zscores_rowcenter <- merge( x=myPhenoData_cutree_zscores, y=cutree_rowcenterDf , by="Sample_ID", all=TRUE)
rownames(myPhenoData_cutree_zscores_rowcenter) <-rownames(myPhenoData_cutree_zscores) 

########################################################
## Pasting the new variables to the big PhenoData
########################################################
complement <- matrix( rep(NA, length(setdiff(rownames(allPhenoData),rownames(myPhenoData_cutree_zscores_rowcenter))) * 2),
        ncol=2, dimnames = list(setdiff(rownames(allPhenoData),rownames(myPhenoData_cutree_zscores_rowcenter)),
                                c(paste0("cutree_zscores_",Labels),paste0("cutree_rowcenter_",Labels)) ) )

PhenoNoclusters <- allPhenoData[setdiff(rownames(allPhenoData),rownames(myPhenoData_cutree_zscores_rowcenter)),]
PhenoNoclusters <- cbind(PhenoNoclusters, complement)

allPhenowithClusters <- rbind(PhenoNoclusters, myPhenoData_cutree_zscores_rowcenter)

allPhenowithClusters_ordered <- allPhenowithClusters[order(match( rownames(allPhenowithClusters)  ,
                                                      rownames(allPhenoData))),] # Ordering according the initial PhenoData (which conserve the same order as the exprresion matrix)

write.table( allPhenowithClusters_ordered , file = paste0(Path_of_Results,"PhenoData_Clusters_",Labels,".tsv" ), sep="\t", quote=FALSE , row.names= TRUE, col.names= TRUE )
saveRDS( allPhenowithClusters_ordered, file = paste0(Path_of_Results,"PhenoData_Clusters_",Labels,".RDS" ))

