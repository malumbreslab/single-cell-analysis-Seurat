# Single-cell RNA sequencing analysis

## Description

Title: Single-cell analysis from CTC in culture from CDK4/6 inhibitor-treated metastatic breast cancer patients
Analysis: Single-cell RNA sequencing
Date: 04/04/23
Author: Agustín Sánchez Belmonte (asanchezb@cnio.es)
Institution: Spanish National Research Cancer Centre (CNIO)

## Table of contents

- [Workflow](#workflow)
- [Contents of the repository](#contents-of-the-repository)
- [Pipeline](#pipeline)
- [Recomendations](#recomendations)

## Workflow

![This is an image](/images/workflow.png)

## Contents of the repository

- The R script `analysis_seurat.R` that can be used for doing Single-cell RNA sequencing analysis.

## Pipeline

### Step 0: Set up

Set working directory

```
setwd('<path_working_directory>')
```

Required packages

```
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(SCINA)
library(patchwork)
library(ggplot2)
library(readxl)
library(ggrepel)
``` 

Theme

```
theme_set(theme_bw(base_size = 14))
```

Results path

```
PLOTS <- '<PATH_RESULTS'
```
  
### Step 1: Read Input files

Read STARsolo output raw data

```
counts <- Read10X(data.dir = '/Users/asanchezb/Desktop/gcbonel_scGEX_220921/CTC/CTC.out/Gene/raw/')  # Seurat function to read in 10x count data
dim(counts)
```

Create Seurat object

```
CTC_seurat <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200, project = "CTC")
CTC_seurat
```

`min.cells` and `min.features` is to filter genes and cells



### Step 2: Filtering

#### A) Mitochondrial

```
grep("^MT",rownames(CTC_seurat@assays$RNA@counts),value = TRUE)
CTC_seurat[["percent.mt"]] <- PercentageFeatureSet(CTC_seurat, pattern = "^MT")

#### B) Ribosomal

```
grep("^RP[LS]",rownames(CTC_seurat@assays$RNA@counts),value = TRUE)
CTC_seurat[["percent.ribo"]] <- PercentageFeatureSet(CTC_seurat, pattern = "RP[LS]")

#### C) Counts

```
feat.max <- round(mean(CTC_seurat$nFeature_RNA) + 2 * sd(CTC_seurat$nFeature_RNA), digits = -2)
feat.min <- round(mean(CTC_seurat$nFeature_RNA) - 1 * sd(CTC_seurat$nFeature_RNA), digits = -2)
CTC_seurat <- subset(x = CTC_seurat, 
                          subset = nFeature_RNA > feat.min & nFeature_RNA < feat.max & percent.mt < 20 & percent.ribo < 40)

#### D) Largest gene

```
CTC_seurat.nomalat <- CTC_seurat[rownames(CTC_seurat) != "MALAT1",]
CTC_seurat.nomalat$largest_count <- apply(CTC_seurat.nomalat@assays$RNA@counts,2,max)
CTC_seurat.nomalat$largest_index <- apply(CTC_seurat.nomalat@assays$RNA@counts,2,which.max)

CTC_seurat.nomalat$largest_gene <- rownames(CTC_seurat.nomalat)[CTC_seurat.nomalat$largest_index]
CTC_seurat.nomalat$percent.Largest.Gene <- 100 * CTC_seurat.nomalat$largest_count / CTC_seurat.nomalat$nCount_RNA

CTC_seurat$largest_gene <- CTC_seurat.nomalat$largest_gene
CTC_seurat$percent.Largest.Gene <- CTC_seurat.nomalat$percent.Largest.Gene 

rm(CTC_seurat.nomalat)
```


#### Visualization

```
png(paste0(PLOTS,'Distribution.png'), width = 1600, height = 900)
VlnPlot(CTC_seurat, features=c("nFeature_RNA","nCount_RNA","percent.mt", "percent.ribo","percent.Largest.Gene"),ncol = 5)+ scale_y_log10()
dev.off()

png(paste0(PLOTS,'Distribution_scatter.png'), width = 1600, height = 900)
plot1 <- FeatureScatter(CTC_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CTC_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

qc.metrics <- as_tibble(CTC_seurat[[]],rownames="Cell.Barcode") 
head(qc.metrics)

qc.metrics %>%
  arrange(percent.mt) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("Example of plotting QC metrics") +
  geom_hline(yintercept = 750) +
  geom_hline(yintercept = 2000) + scale_x_log10() + scale_y_log10()

qc.metrics %>%
  arrange(percent.ribo) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.ribo)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("Example of plotting QC metrics") +
  geom_hline(yintercept = 750) +
  geom_hline(yintercept = 2000) + scale_x_log10() + scale_y_log10()

qc.metrics %>%
  group_by(largest_gene) %>%
  count() %>%
  arrange(desc(n)) -> largest_gene_list

largest_gene_list

largest_gene_list %>%
  filter(n>140) %>%
  pull(largest_gene) -> largest_genes_to_plot

qc.metrics %>%
  filter(largest_gene %in% largest_genes_to_plot) %>%
  mutate(largest_gene=factor(largest_gene, levels=largest_genes_to_plot)) %>%
  arrange(largest_gene) %>%
  ggplot(aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), colour=largest_gene)) +
  geom_point(size=1) +
  scale_colour_manual(values=c("grey",RColorBrewer::brewer.pal(9,"Set1")))

qc.metrics %>%
  ggplot(aes(percent.Largest.Gene)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Largest Gene") +
  geom_vline(xintercept = 10)

qc.metrics %>%
  ggplot(aes(percent.mt)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage Mitochondrion") +
  geom_vline(xintercept = 10)

qc.metrics %>%
  ggplot(aes(percent.ribo)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage Ribosomal") +
  geom_vline(xintercept = 10)


### Step 3: Create table

`chopchop.py` will need a [table](http://genome.ucsc.edu/cgi-bin/hgTables?command=start) to look up genomic coordinates if you want to supply names of the genes rather than coordinates. To get example genePred table:

- Select organism and assembly
- Select group: Genes and Gene Predictions
- Select track: RefSeq Genes or Ensemble Genes
- Select table: refGene or ensGene
- Select region: genome
- Select output format: all fields from selected table
- Fill name with extension ".gene_table' e.g. danRer10.gene_table
- Get output

```
mkdir genePred_folder
```
Save `file.gene_table` inside of `genePred_folder`.

### Step 4: Download genome

Download *.2bit compressed [genome](http://hgdownload.soe.ucsc.edu/downloads.html):

- Select organism in complete annotation sets section
- Select Full data set
- download *.2bit file

```
mkdir 2bit_folder
wget -P 2bit_folder http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.2bit
```

### Step 5: Download genome

Create fasta version of genome by running twoBitToFa on *.2bit file

``` 
./twoBitToFa 2bit_folder/hg38.2bit hg38.fasta
```

### Step 6: Create bowtie version

Make [bowtie](http://bowtie-bio.sourceforge.net/manual.shtml#the-bowtie-build-indexer) compressed version of genome using your new *.fasta file

```
mkdir ebwt_folder
./bowtie/bowtie-build hg38.fasta ebwt_folder/hg38
```

### Step 6: Create bowtie version

Change `config.json` file, replace paths with your own for .2bit genome files, bowtie (.ewbt) genome files and *.gene_table files

Observe `config.json` in order to see an example.

### Step 7: Permissions

Make sure all these files and programs have proper access rights. You can use the `chmod` command in order to change permissions. Maybe some packages may require compilation for your operating system.

### Step 8a: Run KO pipeline (single gen) 

You must run this in your terminal shell and in gen must type the name of the interest gen (be carefull, you must write gene name correctly, some genes have several names, but it is only in one way).

```
./chopchop.py -G hg38 -o results -Target <gen> --scoringMethod DOENCH_2016 -consensusUnion -t CODING > results/<gen>.txt
```
- -G is the genome to search
- -o output folder
- -Target Target genes or regions
- -t Target the whole gene CODING/WHOLE/UTR5/UTR3/SPLICE
- -consensusUnion this option specifies union of isoforms

When the gene is very small, the design the guides will fail and -t WHOLE is recommended.

### Step 8b: Run Activation pipeline (single gen)

You must run this in your terminal shell and in gen must type the name of the interest gen (be carefull, you must write gene name correctly, some genes have several names, but it is only in one way).

```
./chopchop.py -G hg38 -o results -Target <gen> --scoringMethod DOENCH_2016 -consensusUnion -t PROMOTER -TDP 0 -TUP 300 > results/<gen>.txt
```
- -t Promoter
- -TDP how many bp to target downstream of TSS
- -TUP how many bp to target upstream of TSS

### Step 9: Run pipeline (several genes)

You must run this in your terminal shell and type interest genes separated by spaces.

```
bash chop_pipeline.sh <gen1> <gen2> <gen3> <gen4>
```
## Recomendations

`chopchop.py` has a lot of funtionalities and arguments that you can change, it would be well for you observe this in the chochop [link](https://bitbucket.org/valenlab/chopchop/src/master/) or doing this:

```
./chopchop.py --help
```
