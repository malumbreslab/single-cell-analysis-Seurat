# Single-cell RNA sequencing analysis

## Description

- Analysis: Single-cell RNA sequencing
- Date: 04/04/23
- Author: Agustín Sánchez Belmonte (asanchezb@cnio.es)
- Institution: Spanish National Research Cancer Centre (CNIO)

## Table of contents

- [Step 0: Set up](#step-0-set-up)
- [Step 1: Read Input files](#step-1-read-input-files)
- [Step 2: Filtering](#step-2-filtering)
- [Step 3: Normalization](#step-3-normalization)
- [Step 4: Add aneuploidy information](#step-4-add-aneuploidy-information)
- [Step 5: Higly variable genes](#step-5-higly-variable-genes)
- [Step 6: Scale data](#step-6-scale-data)
- [Step 7: Dimensionality reduction](#step-7-dimensionality-reduction)
- [Step 8: Cluster markers](#step-8-cluster-markers)
- [Step 9: Scores](#step-9-scores)
- [Step 10: Annotation](#step-10-annotation)
- [Step 11: Subsets](#step-11-subsets)
- [Step 12: Save files](#step-12-save-files)


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
PLOTS <- '<path_results>'
```
  
### Step 1: Read Input files

Read STARsolo output raw data

```
counts <- Read10X(data.dir = 'path_to_STARsolo_output>/CTC.out/Gene/raw/')  # Seurat function to read in 10x count data
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
```

#### B) Ribosomal

```
grep("^RP[LS]",rownames(CTC_seurat@assays$RNA@counts),value = TRUE)
CTC_seurat[["percent.ribo"]] <- PercentageFeatureSet(CTC_seurat, pattern = "RP[LS]")
```

#### C) Counts

```
feat.max <- round(mean(CTC_seurat$nFeature_RNA) + 2 * sd(CTC_seurat$nFeature_RNA), digits = -2)
feat.min <- round(mean(CTC_seurat$nFeature_RNA) - 1 * sd(CTC_seurat$nFeature_RNA), digits = -2)
CTC_seurat <- subset(x = CTC_seurat, 
                          subset = nFeature_RNA > feat.min & nFeature_RNA < feat.max & percent.mt < 20 & percent.ribo < 40)
```

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
```

### Step 3: Normalization

Choose one of them

#### A) LogNormalize method

```
CTC_seurat <- NormalizeData(CTC_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
apply(CTC_seurat@assays$RNA@data,1,mean) -> gene.expression
sort(gene.expression, decreasing = TRUE) -> gene.expression
head(gene.expression, n=50)

ggplot(mapping = aes(CTC_seurat@assays$RNA@data["GAPDH",])) + 
  geom_histogram(binwidth = 0.05, fill="yellow", colour="black") + 
  ggtitle("GAPDH expression")
```

#### B) CLR method

```
CTC_seurat <- NormalizeData(CTC_seurat, normalization.method = "CLR", margin = 2) 
ggplot(mapping = aes(CTC_seurat@assays$RNA@data["GAPDH",])) + 
  geom_histogram(binwidth = 0.05, fill="yellow", colour="black") + 
  ggtitle("GAPDH expression")
```
#### C) SCTransform method (recommended)

```
CTC_seurat <- SCTransform(CTC_seurat, vst.flavor = "v2")
ggplot(mapping = aes(CTC_seurat@assays$RNA@data["GAPDH",])) + 
  geom_histogram(binwidth = 0.05, fill="yellow", colour="black") + 
  ggtitle("GAPDH expression")
```

### Step 4: Add aneuploidy information

Read copykat prediction

```
pred.test <- read.table('path_copykat_predictions', sep = '\t', header = TRUE)
```

Remove undefined cells

```
pred.test <- pred.test[-which(pred.test$copykat.pred=="not.defined"),]
```

Merge with metadata and add to Seurat object

```
new_data <- merge(CTC_seurat@meta.data, pred.test, by.x = 0, by.y = 'cell.names', all.x = TRUE)
CTC_seurat@meta.data$copykat.pred <- new_data$copykat.pred
```

### Step 5: Higly variable genes

Identification of highly variable features (feature selection)

```
CTC_seurat <- FindVariableFeatures(CTC_seurat, selection.method = "vst", nfeatures = 2000)
```

Identify the 10 most highly variable genes

```
top10 <- head(VariableFeatures(CTC_seurat), 10)
```

Plot variable features with and without labels

```
png(paste0(PLOTS,'Variable_features.png'), width = 1600, height = 900)
plot1 <- VariableFeaturePlot(CTC_seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2 
dev.off()
```

### Step 6: Scale data

Please skip this step if you have used `SCTransform` normalization method

Scaling data

```
all.genes <- rownames(CTC_seurat)
CTC_seurat <- ScaleData(CTC_seurat, features = all.genes) # vars.to.regress = c("S.Score", "G2M.Score")
```


### Step 7: Dimensionality reduction

#### A) Principal Component Analysis (PCA)

Perform linear dimensional reduction

```
CTC_seurat <- RunPCA(CTC_seurat, features = VariableFeatures(object = CTC_seurat))
```
Assign Cell-Cycle Scores (for human)

```
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
CTC_seurat <- CellCycleScoring(CTC_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
CTC_seurat[[]]
as_tibble(CTC_seurat[[]]) %>%
  ggplot(aes(Phase)) + geom_bar()

as_tibble(CTC_seurat[[]]) %>%
  ggplot(aes(x=S.Score, y=G2M.Score, color=Phase)) + 
  geom_point() +
  coord_cartesian(xlim=c(-0.15,0.15), ylim=c(-0.15,0.15))

RidgePlot(CTC_seurat, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
```

Performing PCA analysis
```
CTC_seurat <- RunPCA(CTC_seurat, features = c(s.genes, g2m.genes))
DimPlot(CTC_seurat)
```

Examine and visualize PCA results a few different ways

```
print(CTC_seurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(CTC_seurat, dims = 1:7, reduction = "pca")

png(paste0(PLOTS,'PCA_orign.png'), width = 1400, height = 1200)
DimPlot(CTC_seurat, reduction = "pca", group.by = 'orig.ident')
dev.off()

png(paste0(PLOTS,'heatmap_PCAs.png'), width = 1400, height = 1200)
DimHeatmap(CTC_seurat, dims = 1:6, balanced = TRUE)
dev.off()
```

Selection

```
CTC_seurat <- JackStraw(CTC_seurat, num.replicate = 100)
CTC_seurat <- ScoreJackStraw(CTC_seurat, dims = 1:20)

png(paste0(PLOTS,'JackStrawPlot.png'), width = 1400, height = 1200)
JackStrawPlot(CTC_seurat, dims = 1:10)
dev.off()

png(paste0(PLOTS,'Elbow_plot.png'), width = 1400, height = 1200)
ElbowPlot(CTC_seurat)
dev.off()
```

#### B) UMAP

```
CTC_seurat <- FindNeighbors(CTC_seurat, dims = 1:10)
CTC_seurat <- FindClusters(CTC_seurat, resolution = 0.5)
```

Look at cluster IDs of the first 5 cells

```
head(Idents(CTC_seurat), 5)
```
```
CTC_seurat <- RunUMAP(CTC_seurat, dims = 1:10)
DimPlot(CTC_seurat, reduction = "umap")
DimPlot(CTC_seurat, reduction = "umap", group.by = 'copykat.pred')

png(paste0(PLOTS,'UMAP_clusters.png'), width = 1200, height = 1200)
DimPlot(CTC_seurat, reduction = "umap")
dev.off()
```

Quality visualization

```
FeaturePlot(CTC_seurat, features = 'nCount_RNA')
FeaturePlot(CTC_seurat, features = 'percent.mt')
FeaturePlot(CTC_seurat, features = 'percent.ribo')
```

### Step 8: Cluster markers

```
CTC.markers <- FindAllMarkers(CTC_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(CTC.markers, 'markers_clusters.txt', sep = "\t")
top5 <- CTC.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

CTC.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

png(paste0(PLOTS,'Heatmap.png'), width = 1400, height = 1200)
DoHeatmap(CTC_seurat, features = top10$gene) + NoLegend()
dev.off()
```

### Step 9: Scores

Scoring Epithelial markers

```
CTC_features <- list(c("EPCAM", "KRT8", "KRT4", "KRT5","KRT18","KRT19","KRT9","KRT14","KRT15","ESR1","PGR","AR","ERBB2","ERBB1","MKI67"))

CTC_seurat <- AddModuleScore(
  object = CTC_seurat,
  features = CTC_features,
  ctrl = 5,
  name = 'CTC_Features'
)
FeaturePlot(CTC_seurat, features = "CTC_Features1")
```
### Step 10: Annotation

#### A) Manual annotation

1) CTCs 

```
CTC <- c("EPCAM", "KRT8", "KRT4", "KRT5","KRT18","KRT19","KRT9","KRT14","KRT15","ESR1","PGR","AR","ERBB2","ERBB1","MKI67")

png(paste0(PLOTS,'CTC_markers.png'), width = 1400, height = 1000)
FeaturePlot(CTC_seurat, features = CTC)
dev.off()
```

2) Immune system

```
FeaturePlot(CTC_seurat, features = 'PTPRC')

Macrophages <- c("ITGAL","ITGAM","ITGAX","CD14","FUT4","FCGR3A","CD33","FCGR1A","CD68","CD80","CCR5","TLR2","MSR1","MRC1","TLR4","HLA-DRA","CD163")
Tcells <- c("CD4","CD3D","CD3E","CD3G","FOXP3","PTPRC","IL2RA","CD6","CCR6","CD8A","CCR7","IL7R","SELL","CD27","CTLA4","GZMB")
NK <- c("NCAM1","CD2","FCGR3A","KLRD1","NCR1","KLRB1","KLRK1","CD69","NKG7","IL2RB","NKG2D","GZMB","GZMA","GZMM","FCGR3A","GNLY","COX6A2","ZMAT4","KIR2DL4")
Bcells <- c("CD19","CD24","MS4A1","PAX5","JCHAIN","CD37","CD74","CD79A","CD79B","HLA-DPA1","HLA-DRA")
Platelets <- c("CD9","ITGB1","PECAM1","FCGR2A","CD36","ITGA2B","GP1BA","SPN","CD46","CD47","CD48","ITGA2","ITGA6","PPBP","PF4","GNG11")
DC <- c("CD4","CD1A","CD1B","CD1C","ITGAM","ITGAX","CD40","CCR7","IL3RA","CST3","NRP1","FCER1A")
Fibroblast <- c("MME", "ITGB1", "CD47", "CD81","LRP1","IL1R1")
Erythrocytes <- c("CD24","GYPA","RUVBL1","HBB","HBD")
```


A) Tcells

```
png(paste0(PLOTS,'Tcells_markers.png'), width = 1400, height = 1200)
FeaturePlot(CTC_seurat, features = Tcells)
dev.off()
```

B) NK

```
png(paste0(PLOTS,'NK_markers.png'), width = 1400, height = 1200)
FeaturePlot(CTC_seurat, features = NK)
dev.off()
```



C) Macrophages/Monocytes

```
png(paste0(PLOTS,'Macrophages_markers.png'), width = 1400, height = 1200)
FeaturePlot(CTC_seurat, features = Macrophages)
dev.off()
```

D) Bcells

```
png(paste0(PLOTS,'Bcells_markers.png'), width = 1400, height = 1200)
FeaturePlot(CTC_seurat, features = Bcells)
dev.off()
```

E) Platelets
```
png(paste0(PLOTS,'Platelets_markers.png'), width = 1400, height = 1200)
FeaturePlot(CTC_seurat, features = Platelets)
dev.off()
```

F) DC

```
png(paste0(PLOTS,'DC_markers.png'), width = 1400, height = 1200)
FeaturePlot(CTC_seurat, features = DC)
dev.off()
```

G) Fibroblast

```
png(paste0(PLOTS,'Fibroblast_markers.png'), width = 1400, height = 1200)
FeaturePlot(CTC_seurat, features = Fibroblast)
dev.off()
```

Combination

```
Comb <- c("CD3D","IL7R","NKG7","GNLY","CD14","CD68","CD79A","MS4A1")
png(paste0(PLOTS,'Comb_markers.png'), width = 1400, height = 1200)
FeaturePlot(CTC_seurat, features = Comb)
dev.off()
```

#### B) Automatic annotation

1. SCINA

```
as.data.frame(CTC_seurat@assays$RNA[,]) -> scina.data
load(system.file('extdata','example_signatures.RData', package = "SCINA"))

names(signatures)[1] <- 'Monocytes'
names(signatures)[2] <- 'Bcells'
names(signatures)[3] <- 'NK'

signatures$CTC <- CTC
signatures$Fibroblast <- Fibroblast
signatures$Platelets <- Platelets
signatures$Erythrocytes <- Erythrocytes
signatures$Macrophages <- Macrophages
signatures$Tcells <- Tcells
signatures$Bcells <- Bcells
signatures$DC <- DC
signatures$NK <- NK


SCINA(
  scina.data,
  signatures, 
  max_iter = 100, 
  convergence_n = 10, 
  convergence_rate = 0.999, 
  sensitivity_cutoff = 0.9, 
  rm_overlap=TRUE, 
  allow_unknown=FALSE

) -> scina.results

CTC_seurat$scina_labels <- scina.results$cell_labels

png(paste0(PLOTS,'UMAP_scina_annotations_less_unknown.png'), width = 1200, height = 1200)
DimPlot(CTC_seurat,reduction = "umap", pt.size = 0.8, label = TRUE, group.by = "scina_labels", label.size = 5)
dev.off()

table(CTC_seurat$scina_labels)


Idents(object = CTC_seurat) <- "scina_labels"
CTC.markers.scina <- FindAllMarkers(CTC_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
```

2. AZIMUV

Azimuth predictions

```
predictions <- read.delim('<path_azimuv_predictions>', row.names = 1)
CTC_seurat <- AddMetaData(
  object = CTC_seurat,
  metadata = predictions)

DimPlot(CTC_seurat,reduction = "umap", pt.size = 0.8, label = TRUE, group.by = "predicted.celltype.l2", label.size = 5)
DimPlot(CTC_seurat,reduction = "umap", pt.size = 0.8, label = TRUE, group.by = "predicted.celltype.l1", label.size = 5)
```

### Step 11: Subset

Subset CTCs

```
CTC_seurat_CTC <- subset(x = CTC_seurat, 
                          subset = scina_labels == 'CTC')
```


### Step 12: Save files

Save H5 object

```
SaveH5Seurat(CTC_seurat, overwrite = TRUE)
CTC_seurat <- LoadH5Seurat("CTC/CTC.h5Seurat")
```

Save an object to a file

```
saveRDS(CTC_seurat, file = "/Users/asanchezb/Desktop/gcbonel_scGEX_220921/CTC/CTC_seurat.rds")
CTC_seurat <- readRDS(file = "/Users/asanchezb/Desktop/gcbonel_scGEX_220921/CTC/CTC_seurat.rds")
```




