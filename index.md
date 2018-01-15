# Muhammad Arif's BB2490 Individual Diary

## Table of contents
1. [2017-11-28](#1)
2. [2017-11-29](#2)
3. [2017-11-30](#3)
4. [2017-12-02](#4)
5. [2017-12-03](#5)
6. [2017-12-05 (DESeq2)](#6)
7. [2017-12-05 (PIANO)](#7)
8. [2017-12-09](#8)
9. [2017-12-09](#9)
10. [2017-12-28](#10)
11. [2018-1-11](#11)
12. [2018-1-12](#12)

## Diary
### 2018-1-12<a name="12"></a>:
Our last team meeting before the preso: Finalizing our paper! Ready for the presentation

### 2018-1-11<a name="11"></a>:
We met again to finalize our results and to create the poster. We have decided what to include in the posters (you will see it during the presentation :D). My next task for today: Create a network figure for the metabolic pathways with number of shared significant genes as the edge. This figure will help us to emphasize our finding that fatty acid and amino acid metabolic pathways are disturbed.

Here's the figure:

![net](https://muharif.github.io/BB2490/Network_done_v3.png "net")



### 2017-12-28<a name="10"></a>:
I was away for several days and these past few days I was trying to search for more information to validate our result. Collected several papers, the top 3 that I think it's good to mention:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5572395/
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5541415/
https://www.ncbi.nlm.nih.gov/pubmed/28992545

These papers support most of our findings.

### 2017-12-11<a name="9"></a>:
3rd Group Meeting. We decided that now we have to carry on to the next step, validation of the result. I'm in charge of reviewing the functional analysis result, especially the progession from Healthy to Steathosis.

### 2017-12-09<a name="8"></a>:
To get more sensitivity in the analysis, we decided to do one-way ANOVA test to find the significantly changing genes. We understood that it would not give much information about in which stage is the changes occur. After we did this, we found out that the number of significantly different genes in terms of expression was quite similar to the DESeq2 result.

### 2017-12-05 (PIANO) <a name="7"></a>:
Functional or gene enrichment analysis was performed using PIANO. The input for piano was:
1. Gene Set Collection (retrieved from MSigDB): Mapping from pathway/biological process to gene
2. Fold Changes and P-Values from DESeq2

As the result of this, we found several pathways and processes that are significantly affected by the disease.

```R
library(DESeq2)
conds=as.factor(substr(colnames(count_sub), 1, 1))
subject=colnames(count_sub)
coldata <- data.frame(row.names=subject,conds)
dds <- DESeqDataSetFromMatrix(
  countData=count_sub,
  colData=coldata,
  design=~conds)
dds <- DESeq(dds)
res=results(dds,contrast=c('conds',cond2,cond1))
cuk=mcols(res, use.names=TRUE)[2,2]
res=data.frame(res)
res=cbind(res,hpa[rownames(res),])
write.table(res,paste0('../Result/DifferentialExpression/1_',strsplit(cuk,': ')[[1]][2],'.txt'),sep='\t')
```

### 2017-12-05<a name="6"></a>:
Differential Gene Expression analysis was performed using DESeq2 result and we found very less number of significant gene in the first stage 1 vs 2 and stage 2 vs 3 (~5). But in the last stage we found quite a number of significantly (~4000 genes). The result from DESeq2 will be used for the functional analysis using PIANO.

```R
library(DESeq2)
conds=as.factor(substr(colnames(count_sub), 1, 1))
subject=colnames(count_sub)
coldata <- data.frame(row.names=subject,conds)
dds <- DESeqDataSetFromMatrix(
  countData=count_sub,
  colData=coldata,
  design=~conds)
dds <- DESeq(dds)
res=results(dds,contrast=c('conds',cond2,cond1))
cuk=mcols(res, use.names=TRUE)[2,2]
res=data.frame(res)
res=cbind(res,hpa[rownames(res),])
write.table(res,paste0('../Result/DifferentialExpression/1_',strsplit(cuk,': ')[[1]][2],'.txt'),sep='\t')
```

### 2017-12-03<a name="5"></a>:
Filtering the protein coding genes only. The list of protein coding genes were taken from http://proteinatlas.org/

As the result, we have 87 samples with ~18000 genes. This result will be used for all future analysis.

### 2017-12-02<a name="4"></a>:
Starting the data pre-processing. All codes were implemented in Python with the help of Pandas module for data manipulation. 

As the result: a clean count file, with 87 samples and ~40000 genes (will be filtered more for only protein coding genes)
```python
import pandas as pd

metadata=pd.read_table('../Data/Metadata_Processed.txt')
count=pd.read_table('../Data/Count_Processed.txt')

count_gene=pd.concat([count,gene['Gene']],1).groupby('Gene').sum()
metadata2=metadata
metadata2.index=metadata2['Sample_ID']
combi=pd.concat([count_gene.T,metadata2['NASH_IDENTIFIER']],1).groupby('NASH_IDENTIFIER').sum().T
combi=combi.T[~(combi.columns.str.contains('X') | combi.columns.str.contains('U'))].T
combi.dropna().to_csv('../Data/Done_deal.txt',sep='\t')
```

### 2017-11-30<a name="3"></a>:
Group meeting to kick off the work: We decided that everyone will take part in the quality control and data normalization to be able to understand the data better

### 2017-11-29<a name="2"></a>:
Preparing libraries for the data analysis, for discussion tomorrow with the group:
* Gene To Transcript data
* Gene Length data (for calculating normalized count)

### 2017-11-28<a name="1"></a>:
First entry of the diary!
We had our first group meeting at KTH Library at 14.00 today. Several things was discussed:
* Overview about the project
* Overview of the data
* Discussion about the expected timeline
* Preparation for first seminar presentation




