# This program Split an expression matrix (EM) in submatrices that include
# the controls submatrix with the submatrix of certain label (in that order)
# according to a given data frame that contains the column names of the EM in the same order that the labels. 

#################################################
### Loading and installing required libraries
#################################################
args <- commandArgs(trailingOnly = TRUE)
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("GOexpress")) {
  BiocManager::install("GOexpress", ask =FALSE)
  library("GOexpress")
}

###########################################
### Data given by the user
###########################################
# Path_to_your_ExpressionMatrix <-c("../Data/expMatrix_TCGA_cBioPortal_no_males_withindicator.txt") ; Path_to_your_PhenoData <-c("../Data/Labels_Controls_and_Normal_separated_TCGA.txt") ; Path_of_Results<-c("../Results/Splited_SubMatrices_with_controls/") ; Title_of_Results <- c("TCGA"); colname_of_PhenoData_to_split<-c("Labels") ; values_to_use_for_split<- c("Control") 
Path_to_your_ExpressionMatrix <- args[1]
Path_to_your_PhenoData <- args[2]
Path_of_Results <- args[3]
Title_of_Results <- args[4]
colname_of_PhenoData_to_split <- args[5]
values_to_use_for_split <- args[6] # Houston we have a problem
############################################
### Reading the data
############################################
MyEMatrix <- read.table(Path_to_your_ExpressionMatrix,header=TRUE, sep= "\t", quote = "")
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

################################################################
### Splitting the original expressionSet in subexpressionSets
################################################################
pData(myExpressionSet)

sub.myExpressionSet <- subEset(
  eSet = myExpressionSet,
  subset = list(
    Labels = c("Control")
  )
)

sub.myExpressionSet <- subEset(
  eSet = myExpressionSet,
  subset = list(
    get(colname_of_PhenoData_to_split)
    )
)



myfun <- function(var, env = globalenv()) {
  assign(eval(substitute(var)),values_to_use_for_split, envir = env)
  print(get(var))
}
  
myfun(colname_of_PhenoData_to_split)
Labels

sub.AlvMac <- subEset(
  eSet=AlvMac,
  subset=list(
    Treatment=c("CN","MB"),
    Time=c("2H","6H")
  )
)
hola <- ls()
catch_my_var1 <- get(colname_of_PhenoData_to_split)

subset <- list(get(colname_of_PhenoData_to_split))
colname_of_PhenoData_to_split
Labels
########################################################
### Splitting the original data frame in subdataframes
########################################################
### Preparing the data for the splitting function
EMLabeled <- cbind(as.character(PhenoData$Labels), as.data.frame(t(MyEMatrix) , StringAsFactors = FALSE) ) # We coerced to a dataframe because the cutting tool only work in that way 
colnames(EMLabeled)[1] <- "Labels" 
list_of_submatrices <- split( EMLabeled , EMLabeled$Labels )

### Deleting the first column "Labels"
for( k in 1:length(list_of_submatrices)){
  DeleteThisPositions <- which(colnames(list_of_submatrices[[k]]) %in% c("Labels"))
  list_of_submatrices[[k]] <- list_of_submatrices[[k]][, -DeleteThisPositions]
}

### Transposing matrices to their original orientation
list_of_submatrices <- lapply( list_of_submatrices, t)

##########################################################################
#########  Pasting the submatrix of "Controls" to the other submatrices
##########################################################################
Position_of_controls <- which(names(list_of_submatrices) %in% "Control")
if (  (names(list_of_submatrices)[Position_of_controls] == "Control") == TRUE ){
  print("Ok: You have a Label of -Control- that will be used for further processing")
} else{
  print("ERROR: You DON'T have a Label of -Control- (Exact match)")
} # Checking consistency between the given data

### Pasting the submatrices with the submatrix of Control
list_of_submatrices_pasted_with_Control <- list()
for( k in setdiff(names(list_of_submatrices),"Control")  ){
  list_of_submatrices_pasted_with_Control[[k]] <- cbind(list_of_submatrices[["Control"]],list_of_submatrices[[k]])
}

###########################################
######### Save the results
###########################################
dir.create( Path_of_Results, recursive = TRUE)
saveRDS(list_of_submatrices_pasted_with_Control , file = paste0(Path_of_Results,"List_of_Control_with_subexpression_matrix_of_",Title_of_Results,".RDS" ))
for( k in names(list_of_submatrices_pasted_with_Control)  ){
  write.table( list_of_submatrices_pasted_with_Control[[k]], col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE,
              file = paste0(Path_of_Results,"Control_with_subexpression_matrix_of_",k,"_from_",Title_of_Results,"_.tsv" ))
}

data(AlvMac)
str(pData(AlvMac))
assayData(AlvMac)
experimentData(AlvMac)
str(exprs(AlvMac))

# Subset it to only samples of "CN" and "MB" treatments, and also only "2H",
# "6H", and "24H" time-points
sub.AlvMac <- subEset(
  eSet=AlvMac,
  subset=list(
    Treatment=c("CN","MB"),
    Time=c("2H","6H")
  )
)

str(pData(sub.AlvMac))
str(exprs(sub.AlvMac))
