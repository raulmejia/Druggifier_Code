# clinical Data (Pheno Data) of metabric downloaded from https://www.cbioportal.org/study/clinicalData?id=brca_metabric 
# 
Path_to_your_ExpressionMatrix <-c("/home/raulmejia/Documentos/Druggifier/Data/METABRIC/joined_indicator_METABRIC.txt")
Path_to_your_PhenoData <-c("/home/raulmejia/Documentos/Druggifier/Data/METABRIC/brca_metabric_clinical_data.tsv")
Path_of_Results <-c("/home/raulmejia/Documentos/Druggifier/Data/METABRIC/")
ExpressionMatrix <- read.table(Path_to_your_ExpressionMatrix, header=TRUE, sep= "\t", quote = "")

PhenoData <- read.table(Path_to_your_PhenoData, header=TRUE, sep= "\t", quote = "")

rownames(PhenoData)<- gsub("-",".",PhenoData[,1])
positions <- which( rownames(PhenoData) %in% colnames(ExpressionMatrix))
PhenoDiscovSet <- PhenoData[positions,]
#head(PhenoData)
write.table(PhenoDiscovSet , file = paste0(Path_of_Results,"PhenoDataDiscovSet.tsv"), sep="\t", quote=FALSE , row.names= TRUE, col.names= TRUE )
