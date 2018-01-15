# Muhammad Arif's BB2490 Individual Diary

## Table of contents
1. [2017-11-28](#1)
2. [2017-11-29](#2)
3. [2017-11-30](#3)
4. [2017-12-02](#4)
5. [2017-12-03](#5)
6. [2017-12-05 (DESeq2)](#6)
7. [2017-12-05 (PIANO)](#7)

## Diary
### 2017-12-05 (PIANO) <a name="7"></a>:
Functional or gene enrichment analysis was performed using PIANO. The input for piano was:
1. Gene Set Collection (retrieved from MSigDB): Mapping from pathway/biological process to gene
2. Fold Changes and P-Values from DESeq2

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




