This repo describes the approach used in https://www.frontiersin.org/articles/10.3389/fphar.2018.00905/full  



Modules:

## 1) Subtyping (B)
Obtaining the subtypes for your gene expression datasets

## 2) Selecting genes
#### 2.1) DEG (E)  
Differential expressed genes in the datasets
#### 2.2) Pathway enrichment (C)
Pathway enrichment analysis (in this case through Pathifier) see the paper
#### 2.3) Intersecting DEG n Pathway enrichment (F)
Integrating the above mentioned methods to select targets 

## 3) Drug-target_Annotating (D)
Considering the dataframe of integrated information you can identify which of your genes are druggable

## 4) Ranking
Once obtained the previous dataframe you can apply several criteria to choose your candidate targets

## 5) Multiporpuse libraries
General purpuse libraries required for the programs 

![alt text](https://www.frontiersin.org/files/Articles/390753/fphar-09-00905-HTML/image_m/fphar-09-00905-g001.jpg)
