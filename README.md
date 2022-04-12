This repo describes the approach used in https://www.frontiersin.org/articles/10.3389/fphar.2018.00905/full  

Modules:

## 1) Subtyping
Obtaining the subtypes for your gene expression datasets

## 2) Selecting genes
#### 2.A) DEG  
Differential expressed genes in the datasets
#### 2.B) Pathway enrichment
Pathway enrichment analysis (in this case through Pathifier) see the paper
#### 2.C) Intersecting DEG n Pathway enrichment
Integrating the above mentioned methods to select targets 

## 3) Drug-target_Annotating
Considering the dataframe of integrated information you can identify which of your genes are druggable

## 4) Ranking
Once obtained the previous dataframe you can apply several criteria to choose your candidate targets

## 5) Multiporpuse libraries
General purpuse libraries required for the programs 
