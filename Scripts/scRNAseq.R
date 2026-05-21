#This is a script to analyse single cell RNA sequencing data

#load in libraries 
library(tidyverse)
library(Seurat)
library(patchwork)
library(hdf5r)

#Load in the count matrix and create a Seurat object 

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

# 1. QC on data 
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

#