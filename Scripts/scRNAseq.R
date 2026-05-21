#This is a script to analyse single cell RNA sequencing data

#load in libraries 
library(tidyverse)
library(Seurat)
library(patchwork)
library(hdf5r)
library(dplyr)

#Data analysis 
#Separate script

# 1. Load in the count matrix and create a Seurat object 

#Load in count matrix 

pbmc.n <- Read10X_h5("raw_data/Data020/sc5p_v2_hs_PBMC_10k_raw_feature_bc_matrix.h5")
#Check count matrix 
str(pbmc.n)
#Only choose Gene Expression 
pbmc_counts <- pbmc.n$`Gene Expression`

#Load data in Seurat object 
#It is possible to add the functions min.cells and min.features for filtering 
pbmc.seurat.obj <- CreateSeuratObject(pbmc_counts, project="PBMC", min.cells = 3, min.features = 200)
pbmc.seurat.obj
#36601 features across 737280 samples 

#Check object 
str(pbmc.seurat.obj)


#General analysis 
# 2. QC on data 
View(pbmc.seurat.obj@meta.data)
#Check for poor quality cells or doublets 

#Check mitochondrial transcript (MT) percentage. dying/low quality cells have hgher MT concentrations 
pbmc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(pbmc.seurat.obj, pattern = "^MT[-\\.]")

#Check if MTT percentage are added 
View(pbmc.seurat.obj@meta.data)

#Look at distribution data 
VlnPlot(pbmc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#Other output 
VlnPlot(pbmc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, 
        pt.size=0)


#View correlation of the data 
library(patchwork)
plot1 <- FeatureScatter(pbmc.seurat.obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Filter out low quality cells based on number of genes and MT percentage 
pbmc.seurat.obj <- subset(pbmc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                   percent.mt < 5)
#Check cell amount 
pbmc.seurat.obj
#20233 features across 7574 samples 

#Can filter out for doublets -> see if data gets better 
#Package: DoubletFinder (2 cells being sequenced together as 1 cell)


#3. normalize data
#Divide gene measurement in each cell by the total expression multiplied by a scaling factor, log trasnsform 
pbmc.seurat.obj <- NormalizeData(pbmc.seurat.obj)
#Find slot command, post filtering the information has been saved 
str(pbmc.seurat.obj)

#4. Identify highly variable features 
#Select a subset of features that exhibit high cell to cell variation 

#Most of the time, the value is between 2000 and 5000 
#Check the other options there are. 
pbmc.seurat.obj <- FindVariableFeatures(pbmc.seurat.obj, selection.method = "vst", nfeatures = 2000)


#Check top 1 highly variable genes
top10_HVG <- head(VariableFeatures(pbmc.seurat.obj), 10)

#Plot features 
plot1 <- VariableFeaturePlot(pbmc.seurat.obj)
LabelPoints(plot = plot1, points = top10_HVG, repel = TRUE)

#Plot variable features 
top_features <- head(VariableFeatures(pbmc.seurat.obj), 20)
plot1 <- VariableFeaturePlot(pbmc.seurat.obj)
plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)
plot1 + plot2

#5. Data scaling 
#Acount for sources of variance 
#Biologival sources: Difference in cell cycle 
#Technical noise: Bacth effect 

all.genes <- rownames(pbmc.seurat.obj)
pbmc.seurat.obj <- ScaleData(pbmc.seurat.obj, features = all.genes)
#You can add regress out more variation. First look at the results and then regress out more variation of needed. 

#See slots in Seurat object 
str(pbmc.seurat.obj)


#6. Perform linear dimension reduction (PCA)
#Identify sources of heterogenity in dataset 

#Run PCA
#default option 
pbmc.seurat.obj <- RunPCA(pbmc.seurat.obj, features = VariableFeatures(object = pbmc.seurat.obj))

#Visualize PCA results 
print(pbmc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)

#Which principle components needs to be added in the downstream analysis?
#Heatmap 
DimHeatmap(pbmc.seurat.obj, dims = 1:5, cells = 500, balanced = TRUE)

#Determine dimensionality of the data
#Determine only significant principle components which capture the signal in Down stream analysis 
#Plot default PC
ElbowPlot(pbmc.seurat.obj)

#Plot all PC
ElbowPlot(pbmc.seurat.obj, ndims = ncol(Embeddings(pbmc.seurat.obj, "pca")))
#Look at where there is high standard deviation 

#7. Cluster cells
#Cells with similair feature expression patterns are clustered together 
pbmc.seurat.obj <- FindNeighbors(pbmc.seurat.obj, dims = 1:15)
View(pbmc.seurat.obj@meta.data)

#Resolution
#Defines the resolution of the clusters. The higher the number, there are more clusters
#See which resolution works best to distinct the cells into different clusters 
pbmc.seurat.obj <- FindClusters(pbmc.seurat.obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1))

#Plot to see which resolution works best 
DimPlot(pbmc.seurat.obj, group.by = "RNA_snn_res.0.3", label = TRUE)
colnames(pbmc.seurat.obj@meta.data)

#Set identity of cells to cluster 

#See default resolution 
Idents(pbmc.seurat.obj)

#Change resolution 
Idents(pbmc.seurat.obj) <- "RNA_snn_res.0.3"

#8. non-linear dimensionality reduction for display 
#UMAP or t-SNE can be used 
pbmc.seurat.obj <- RunUMAP(pbmc.seurat.obj, dims = 1:15)
pbmc.seurat.obj <- RunTSNE(pbmc.seurat.obj, dims = 1:15)

plot1 <- TSNEPlot(pbmc.seurat.obj, label = TRUE)
plot2 <- UMAPPlot(pbmc.seurat.obj, label = TRUE)
plot1 + plot2
#Individual clusters 
DimPlot(pbmc.seurat.obj, reduction = "umap")


#Annotate each cluster 
# Option 1, identify signature genes and compate those with literature 
cluster_markers <- FindAllMarkers(pbmc.seurat.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
cluster_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)




#Exploratory analysis

#DEG analysis 

#GO/KEGG enrichment 

#GSVA 

#TF identification 

#Cell trajectory 

#Cell-cell interaction 

#Cell Cycle 

#Spatial transcriptome 



