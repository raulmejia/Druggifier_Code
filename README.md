This repo describes the approach used in https://www.frontiersin.org/articles/10.3389/fphar.2018.00905/full  


---
Modules:
---
### 0) Multiporpuse libraries
General purpose libraries required for the programs

### 1) Subtyping 
Obtaining the subtypes for your gene expression datasets (B)

### 2) Selecting genes
#### 2.1) DEG   
Differential expressed genes in the datasets (E)
#### 2.2) Pathway enrichment 
Pathway enrichment analysis (in this case through Pathifier) see the paper (C)
#### 2.3) Intersecting DEG n Pathway enrichment 
Integrating the above-mentioned methods to select targets (F)

### 3) Drug-target_Annotating 
Considering the data frame of integrated information you can identify which of your genes are druggable (D)

### 4) Ranking
Once obtained the previous dataframe you are free to apply several criteria to choose your candidate targets

 

![alt text](https://www.frontiersin.org/files/Articles/390753/fphar-09-00905-HTML/image_m/fphar-09-00905-g001.jpg)
