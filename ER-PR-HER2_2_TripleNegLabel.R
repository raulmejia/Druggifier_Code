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
# Path_to_your_ExpressionMatrix <-c("") ; Path_of_Results<-c("../Results/") ; Title_of_Results <- c("METABRIC-triple-neg-Label") 
Path_to_your_ExpressionMatrix <- args[1]
Path_of_Results <- args[2]
Title_of_Results <- args[3]
############################################
### Reading the data
############################################
MyEMatrix <- read.table(Path_to_your_ExpressionMatrix, header=FALSE, sep= "\t", quote = "", row.names = FALSE)
?read.table
myPhenoData <-read.table(Path_to_your_PhenoData)

############################################
### Building an expressionSet
############################################
My_pData_Anotdf <- AnnotatedDataFrame( data= myPhenoData) # Creating the PhenoData as an annotated data frame  
validObject(My_pData_Anotdf)

myExpressionSet <- ExpressionSet(as.matrix(MyEMatrix), phenoData=My_pData_Anotdf )

############################################
### Quality control 
############################################
if (length(colnames(MyEMatrix)) / ( sum(colnames(MyEMatrix) == rownames(myPhenoData)) ) == 1 ){
  print("Ok: Your expression matrix and Labels have the same names in the same order")
} else{
  print("ERROR: Your expression matrix and Labels DOESN'T have the same names in the same order")
} # Checking consistency between the given data


###########################################
########## Saving the results
############################################
Path_of_Results_Normalized<- normalizePath(Path_of_Results)
