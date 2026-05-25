#This is a script to analyse single cell RNA sequencing data

#load in libraries 
library(tidyverse)
library(Seurat)
library(patchwork)
library(hdf5r)
library(dplyr)
library(DoubletFinder)
library(ggplot2)
library(celldex)
library(SingleR)
library(viridis)
library(pheatmap)

#Data analysis 
#Separate script

# 1. Load in the count matrix and create a Seurat object 

#Load in count matrix 

pbmc.n <- Read10X_h5(filename = "raw_data/Data020/sc5p_v2_hs_PBMC_10k_raw_feature_bc_matrix.h5")
#Check count matrix 
str(pbmc.n)
#Only choose Gene Expression 
pbmc_counts <- pbmc.n$`Gene Expression`

#Load data in Seurat object 
#It is possible to add the functions min.cells and min.features for filtering 
pbmc.seurat.obj <- CreateSeuratObject(pbmc_counts, project="PBMC", min.cells = 3, min.features = 200)
pbmc.seurat.obj
#20233 features across 10690 samples 

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
pbmc.seurat.filtered <- subset(pbmc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                   percent.mt < 5)
#Check cell amount 
pbmc.seurat.filtered
#20233 features across 7574 samples 

#3. normalize data
#Divide gene measurement in each cell by the total expression multiplied by a scaling factor, log trasnsform 
pbmc.seurat.filtered <- NormalizeData(pbmc.seurat.obj)
#Find slot command, post filtering the information has been saved 
str(pbmc.seurat.filtered)

#4. Identify highly variable features 
#Select a subset of features that exhibit high cell to cell variation 

#Most of the time, the value is between 2000 and 5000 
#Check the other options there are. 
pbmc.seurat.filtered <- FindVariableFeatures(pbmc.seurat.filtered, selection.method = "vst", nfeatures = 2000)


#Check top 1 highly variable genes
top10_HVG <- head(VariableFeatures(pbmc.seurat.filtered), 10)

#Plot variable features 
plot1 <- VariableFeaturePlot(pbmc.seurat.filtered)
LabelPoints(plot = plot1, points = top10_HVG, repel = TRUE)

#5. Data scaling 
#Acount for sources of variance 
#Biologival sources: Difference in cell cycle 
#Technical noise: Bacth effect 

all.genes <- rownames(pbmc.seurat.filtered)
pbmc.seurat.filtered <- ScaleData(pbmc.seurat.filtered, features = all.genes)
#You can add regress out more variation. First look at the results and then regress out more variation of needed. 

#See slots in Seurat object 
str(pbmc.seurat.filtered)


#6. Perform linear dimension reduction (PCA)
#Identify sources of heterogenity in dataset 

#Run PCA
#default option 
pbmc.seurat.filtered <- RunPCA(pbmc.seurat.filtered, features = VariableFeatures(object = pbmc.seurat.filtered))

#Visualize PCA results 
print(pbmc.seurat.filtered[["pca"]], dims = 1:5, nfeatures = 5)

#Which principle components needs to be added in the downstream analysis?
#Heatmap 
DimHeatmap(pbmc.seurat.filtered, dims = 1:5, cells = 500, balanced = TRUE)

#Determine dimensionality of the data
#Determine only significant principle components which capture the signal in Down stream analysis 
#Plot default PC
ElbowPlot(pbmc.seurat.filtered)

#Plot all PC
ElbowPlot(pbmc.seurat.filtered, ndims = ncol(Embeddings(pbmc.seurat.filtered, "pca")))
#Look at where there is high standard deviation 

#7. Cluster cells
#Cells with similair feature expression patterns are clustered together 
pbmc.seurat.filtered <- FindNeighbors(pbmc.seurat.filtered, dims = 1:15)
View(pbmc.seurat.filtered@meta.data)

#Resolution
#Defines the resolution of the clusters. The higher the number, there are more clusters
#See which resolution works best to distinct the cells into different clusters 
pbmc.seurat.filtered <- FindClusters(pbmc.seurat.filtered, resolution = c(0.1, 0.3, 0.5, 0.7, 1))

#Plot to see which resolution works best 
DimPlot(pbmc.seurat.filtered, group.by = "RNA_snn_res.0.3", label = TRUE)
colnames(pbmc.seurat.filtered@meta.data)

#Set identity of cells to cluster 

#See default resolution 
Idents(pbmc.seurat.filtered)

#Change resolution 
Idents(pbmc.seurat.filtered) <- "RNA_snn_res.0.3"

#8. non-linear dimensionality reduction for display 
#UMAP or t-SNE can be used 
pbmc.seurat.filtered <- RunUMAP(pbmc.seurat.filtered, dims = 1:15)
pbmc.seurat.filtered <- RunTSNE(pbmc.seurat.filtered, dims = 1:15)

#Both tests
plot1 <- TSNEPlot(pbmc.seurat.filtered, label = TRUE)
plot2 <- UMAPPlot(pbmc.seurat.filtered, label = TRUE)
plot1 + plot2

#Only UMAP 
DimPlot(pbmc.seurat.filtered, reduction = "umap")

#DoubletFinder
#Filter out doublets in dataset 

#Test different Pk values 

sweep.res.list_pbmc <- paramSweep(pbmc.seurat.filtered, PCs = 1:15, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT= FALSE)
#Find optimal pK value 
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)


#Plot Pk against BCmetric 
ggplot(bcmvn_pbmc, aes(pK, BCmetric, group = 1)) +
  geom_point()+
  geom_line()

#Store optimal pK value to pK variable 
pK <- bcmvn_pbmc %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

#Calculate expexted number of doublets 
#See for each cell what cluster they belong to 
annotations <- pbmc.seurat.filtered@meta.data$seurat_clusters

#Model proportions homotypic doublets 
homotypic.prop <- modelHomotypic(annotations)
#Proportion homotypic doublets 
homotypic.prop 

#Estimate number of doublets expected 
multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
                                  'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
                                  'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))

multiplet_rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < nrow(pbmc.seurat.filtered@meta.data)) %>% 
  dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
  dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells
#Expected number of doublets 
nExp.poi <- round(multiplet_rate * nrow(pbmc.seurat.filtered@meta.data))
#Expected number of doublets after adjustment of homotypic doublets 
nExp.poi_adj <- round(nExp.poi*(1-homotypic.prop))

#Run DoubletFinder 
pbmc.seurat.filtered <- doubletFinder(pbmc.seurat.filtered, 
                                      PCs = 1:15,
                                      pN = 0.25, 
                                      pK = pK, 
                                      nExp = nExp.poi_adj, 
                                      reuse.pANN = FALSE, sct = FALSE)

#9. Annotate each cluster 
#One reference based approach


#Pick a reference that you expect to be in the dataset 
ref <- celldex::HumanPrimaryCellAtlasData()
#Check celltypes that are available in the reference 
View(as.data.frame(colData(ref)))

pbmc <- GetAssayData(pbmc.seurat.filtered, layer = 'counts')

prediction <- SingleR(test = pbmc, 
        ref = ref, 
        labels = ref$label.main)

pbmc.seurat.filtered$singleR.labels <- prediction$labels[match(rownames(pbmc.seurat.filtered@meta.data), rownames(prediction))]
DimPlot(pbmc.seurat.filtered, reduction = "umap", group.by = "singleR.labels")

#Calculate performance 
prediction$scores

plotScoreHeatmap(prediction)

#Plot delta values, to assess poor quality 
#Low delta values, assigning is uncertain 
#Detect cells that are outliers

#Which of the labels have unclear assignment, the cells in that label can be further investigated 
plotDeltaDistribution(prediction)

#Compare to unsupervised clustering 
tab <- table(Assigned=prediction$labels, Clusters=pbmc.seurat.filtered$seurat_clusters)
pheatmap(log10(tab+10), color = colorRampPalette(c("white", "blue"))(10))
View(pbmc.seurat.filtered@meta.data)

#Marker based approach 



#10. Pseudotemporal cell ordering 

#Trajectory analysis 
View(pbmc.seurat.filtered)



#Exploratory analysis

#DEG analysis 

#GO/KEGG enrichment 

#GSVA 

#TF identification 

#Cell trajectory 

#Cell-cell interaction 

#Cell Cycle 

#Spatial transcriptome 



