# DESCRIPTION ####
# The following code will install all dependencies in R allowing for analysis to be performed in the paper Xu and Qadir et al.
# This code was written and developed by Dr. Fahd Qadir PhD (mqadir@tulane.edu)
# Install and run R v4 (x64 bit) and RStudio v1.2.1335 (x64 bit) for running this code
# Seurat is a multimodal single Cell RNA seq analysis algorithm created by
# The Satija Lab. For more information please see: https://satijalab.org/seurat/

# INSTALLATION ####
# Please download and install Rtools 3.5 from http://cran.r-project.org/bin/windows/Rtools/
install.packages("devtools")
install.packages("pkgbuild")

library(devtools)
library(pkgbuild)
find_rtools() # This should come back as true

# Install additional packages
install.packages("Matrix")
install.packages("ggridges")
install.packages("cowplot")
install.packages('ggrepel')
install.packages("R.utils")
install.packages("gridExtra")
install.packages("Seurat")
install.packages("plotly")
install.packages("clustree")
install.packages('multtest')

# Install the scRNAseq analysis package Seurat
install.packages('Seurat') #This installs Seurat v3.0.0 as of 5/20/2019, if Seurat is updated direct installation to v3.0.0

# Install the scRNAseq pseudotemporal analysis package Monocle
# source("http://bioconductor.org/biocLite.R") #Monocle v2.8.0 is used
# biocLite()
# biocLite("monocle")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("monocle")

# Install multtest for Seurat
BiocManager::install("multtest")

# LOAD LIBRARIES ####
# Restart Rstudio or R
# Run the following code once you have Seurat installed

# LOAD LIBRARIES ####
# Restart Rstudio or R
# Run the following code once you have Seurat installed

install.packages('clustree')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

suppressWarnings(
  {
    library(ggplot2)
    library(cowplot)
    library(Matrix)
    library(ggridges)
    library(ggrepel)
    library(dplyr)
    library(Seurat)
    library(monocle3)
    library(plotly)
    library(clustree)
    library(patchwork)
    library(future)
    library(DoubletFinder)
    library(EnhancedVolcano)
    }
  )

# CONFIRM CORRECT INSTALL ####
# Confirm package version of Seurat and Monocle
packageVersion("Seurat")
packageVersion("monocle3")

# Set global environment parameter
options(future.globals.maxSize = 8000 * 1024^2)

# OBJECT SETUP AND NORMALIZATION ####
# STEP 1: Load 10X data ####
HP2107001_ctrl.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Raw Data\F7a GEX_HP21070_01\filtered_feature_bc_matrix)")
HP2107001_DHT.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Raw Data\F7b GEX_HP21070_01_DHT\filtered_feature_bc_matrix)")
HP2107701_ctrl.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Raw Data\F8a_GEX_HP21077_01\filtered_feature_bc_matrix)")
HP2107701_DHT.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Raw Data\F8b_GEX_HP21077_01_DHT\filtered_feature_bc_matrix)")
HP2107901_ctrl.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Raw Data\F9a_GEX_HP21079_01\filtered_feature_bc_matrix)")
HP2107901_DHT.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Raw Data\F9b_GEX_HP21079_01_DHT\filtered_feature_bc_matrix)")

# STEP 1: Load 10X data home PC####
#HP2107001_ctrl.data <- Read10X(data.dir = r"(C:\Users\Qonos\Box\Lab 2301\RNAseq DHT data\Raw Data\F7a GEX_HP21070_01\filtered_feature_bc_matrix)")
#HP2107001_DHT.data <- Read10X(data.dir = r"(C:\Users\Qonos\Box\Lab 2301\RNAseq DHT data\Raw Data\F7b GEX_HP21070_01_DHT\filtered_feature_bc_matrix)")
#HP2107701_ctrl.data <- Read10X(data.dir = r"(C:\Users\Qonos\Box\Lab 2301\RNAseq DHT data\Raw Data\F8a_GEX_HP21077_01\filtered_feature_bc_matrix)")
#HP2107701_DHT.data <- Read10X(data.dir = r"(C:\Users\Qonos\Box\Lab 2301\RNAseq DHT data\Raw Data\F8b_GEX_HP21077_01_DHT\filtered_feature_bc_matrix)")
#HP2107901_ctrl.data <- Read10X(data.dir = r"(C:\Users\Qonos\Box\Lab 2301\RNAseq DHT data\Raw Data\F9a_GEX_HP21079_01\filtered_feature_bc_matrix)")
#HP2107901_DHT.data <- Read10X(data.dir = r"(C:\Users\Qonos\Box\Lab 2301\RNAseq DHT data\Raw Data\F9b_GEX_HP21079_01_DHT\filtered_feature_bc_matrix)")

# STEP 2: Create Seurat objects ####
HP2107001_ctrl <- CreateSeuratObject(counts = HP2107001_ctrl.data, min.features = 500)
HP2107001_DHT <- CreateSeuratObject(counts = HP2107001_DHT.data, min.features = 500)
HP2107701_ctrl <- CreateSeuratObject(counts = HP2107701_ctrl.data, min.features = 500)
HP2107701_DHT <- CreateSeuratObject(counts = HP2107701_DHT.data, min.features = 500)
HP2107901_ctrl <- CreateSeuratObject(counts = HP2107901_ctrl.data, min.features = 500)
HP2107901_DHT <- CreateSeuratObject(counts = HP2107901_DHT.data, min.features = 500)

# Sample specific Metadata addition
HP2107001_ctrl$sample <- "HP2107001_ctrl"
HP2107001_DHT$sample <- "HP2107001_DHT"
HP2107701_ctrl$sample <- "HP2107701_ctrl"
HP2107701_DHT$sample <- "HP2107701_DHT"
HP2107901_ctrl$sample <- "HP2107901_ctrl"
HP2107901_DHT$sample <- "HP2107901_DHT"

# Sex specific Metadata addition
HP2107001_ctrl$sex <- "Male"
HP2107001_DHT$sex <- "Male"
HP2107701_ctrl$sex <- "Male"
HP2107701_DHT$sex <- "Male"
HP2107901_ctrl$sex <- "Male"
HP2107901_DHT$sex <- "Male"

# Treatment specific Metadata addition
HP2107001_ctrl$treatment <- "EtOH"
HP2107001_DHT$treatment <- "DHT[10nM]"
HP2107701_ctrl$treatment <- "EtOH"
HP2107701_DHT$treatment <- "DHT[10nM]"
HP2107901_ctrl$treatment <- "EtOH"
HP2107901_DHT$treatment <- "DHT[10nM]"

# STEP 3: Thresholding ####
# The operator can add columns to object metadata. This is a great place to stash QC stats
HP2107001_ctrl[["percent.mt"]] <- PercentageFeatureSet(object = HP2107001_ctrl, pattern = "^MT-")
HP2107001_DHT[["percent.mt"]] <- PercentageFeatureSet(object = HP2107001_DHT, pattern = "^MT-")
HP2107701_ctrl[["percent.mt"]] <- PercentageFeatureSet(object = HP2107701_ctrl, pattern = "^MT-")
HP2107701_DHT[["percent.mt"]] <- PercentageFeatureSet(object = HP2107701_DHT, pattern = "^MT-")
HP2107901_ctrl[["percent.mt"]] <- PercentageFeatureSet(object = HP2107901_ctrl, pattern = "^MT-")
HP2107901_DHT[["percent.mt"]] <- PercentageFeatureSet(object = HP2107901_DHT, pattern = "^MT-")

# QC information before thresholding
summary(head(HP2107001_ctrl@meta.data))
summary(head(HP2107001_DHT@meta.data))
summary(head(HP2107701_ctrl@meta.data))
summary(head(HP2107701_DHT@meta.data))
summary(head(HP2107901_ctrl@meta.data))
summary(head(HP2107901_DHT@meta.data))

# Visualize QC metrics as a violin plot
VlnPlot(object = HP2107001_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107001_DHT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107701_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107701_DHT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107901_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107901_DHT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# RNA based cell thresholding
HP2107001_ctrl <- subset(x = HP2107001_ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
HP2107001_DHT <- subset(x = HP2107001_DHT, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
HP2107701_ctrl <- subset(x = HP2107701_ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
HP2107701_DHT <- subset(x = HP2107701_DHT, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
HP2107901_ctrl <- subset(x = HP2107901_ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
HP2107901_DHT <- subset(x = HP2107901_DHT, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)

# QC information after thresholding
summary(head(HP2107001_ctrl@meta.data))
summary(head(HP2107001_DHT@meta.data))
summary(head(HP2107701_ctrl@meta.data))
summary(head(HP2107701_DHT@meta.data))
summary(head(HP2107901_ctrl@meta.data))
summary(head(HP2107901_DHT@meta.data))

# Visualize QC metrics post thresholding as a violin plot
VlnPlot(object = HP2107001_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107001_DHT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107701_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107701_DHT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107901_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107901_DHT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = pancreas.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, cols = c("red", "blue"))
# Step 4: Add cell IDs ####
# Add cell IDs
HP2107001_ctrl <- RenameCells(HP2107001_ctrl, add.cell.id = "HP2107001_ctrl")
HP2107001_DHT <- RenameCells(HP2107001_DHT, add.cell.id = "HP2107001_DHT")
HP2107701_ctrl <- RenameCells(HP2107701_ctrl, add.cell.id = "HP2107701_ctrl")
HP2107701_DHT <- RenameCells(HP2107701_DHT, add.cell.id = "HP2107701_DHT")
HP2107901_ctrl <- RenameCells(HP2107901_ctrl, add.cell.id = "HP2107901_ctrl")
HP2107901_DHT <- RenameCells(HP2107901_DHT, add.cell.id = "HP2107901_DHT")

# Step 5: Merge Datasets
# Based on comment to Issue #4753 https://github.com/satijalab/seurat/issues/4753
# We use RPCA to yield conserved mapping and set Tx as control refrence samples
# Merge panc_sex datasets
pancreas.list <- list("HP2107001_ctrl" = HP2107001_ctrl, "HP2107001_DHT" = HP2107001_DHT,
                      "HP2107701_ctrl" = HP2107701_ctrl, "HP2107701_DHT" = HP2107701_DHT,
                      "HP2107901_ctrl" = HP2107901_ctrl, "HP2107901_DHT" = HP2107901_DHT)

# Step 6: Data normalization
#Normalise data
pancreas.list <- lapply(X = pancreas.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

# Step 7: Feature selection
# Select features for downstream integration
pancreas.features <- SelectIntegrationFeatures(object.list = pancreas.list)
pancreas.list <- lapply(X = pancreas.list, FUN = function(x) {
  x <- ScaleData(x, features = pancreas.features, verbose = FALSE)
  x <- RunPCA(x, features = pancreas.features, verbose = FALSE)
})

# Step 8: Anchor identification and data integration
# Identify anchors and integrate dataset
pancreas.list[c(1, 3, 5)] # check that you are correctly picking up control datasets for refrence integration
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, reference = c(1,3,5), # Takes 18min 13sec to run when using cca as reduction
                                           reduction = "rpca", dims = 1:30, verbose = TRUE)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30, verbose = TRUE)

# Step 9: Linear dimensionality assessment
# Look at your default assay
DefaultAssay(object = pancreas.integrated)

# Change default assay to integrated, to view dimensionality
DefaultAssay(object = pancreas.integrated) <- "integrated"

# Scaling this is weird, but as done in https://satijalab.org/seurat/articles/integration_large_datasets.html
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)

# Dimensionality assessment using PCA analysis
pancreas.integrated <- RunPCA(pancreas.integrated, features = VariableFeatures(object = pancreas.integrated))

# Examine data dimensionality
ElbowPlot(pancreas.integrated)

# Step 9a: CLUSTER ANALYSIS ####
# Clustering, we run multiple permutations to allow clustree to analyze optimal clustering resolution.
pancreas.integrated <- FindNeighbors(object = pancreas.integrated, dims = 1:30)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.1)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.2)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.3)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.4)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.5)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.6)

# Cluster-tree analysis, looking appropriate non-anomalous clustering resolution
clustree(pancreas.integrated, prefix = "integrated_snn_res.")

# Based of clustree assessment choose res = 0.5
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 0.5)

# Alternatively build a cluster tree
DefaultAssay(object = pancreas.integrated) <- "integrated"
pancreas.integrated=BuildClusterTree(pancreas.integrated, slot = "scale.data")
PlotClusterTree(pancreas.integrated)

# Step 10: non-linear dimensionality assessment ####
# Run PCA and UMAP calculations
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:30)

# Change default assay to integrated, to view dimensionality
Idents(pancreas.integrated) <- "treatment"
DimPlot(pancreas.integrated, reduction = "umap", 
        cols = c('black', 'red'), 
        label = FALSE,
        order = FALSE)
Idents(pancreas.integrated) <- "sex"
Idents(pancreas.integrated) <- "sample"
Idents(pancreas.integrated) <- "seurat_clusters"
#Idents(pancreas.integrated) <- "integrated_snn_res.0.3"
DimPlot(pancreas.integrated, reduction = "umap", label = FALSE)

#Visualize gene expression
DefaultAssay(object = pancreas.integrated) <- "RNA"
DefaultAssay(object = pancreas.integrated)
FeaturePlot(object = pancreas.integrated,
            features = c("INS", "GCG", "SST", "PPY", "GHRL",
                         "KRT19", "CPA1",
                         "COL1A1", "VWF", "SOX10",
                         "TPSAB1", "SDS", "TRAC"),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            #min.cutoff = 0,
            #max.cutoff = 1,
            slot = 'counts',
            order = TRUE)

FeaturePlot(object = pancreas.integrated,
            features = c("GLP1R"),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            #min.cutoff = 0,
            max.cutoff = 2,
            slot = 'counts',
            order = TRUE)


VlnPlot(
  object = pancreas.integrated,
  features = c("AR"),
  assay = 'RNA',
  slot = 'counts',
  cols = c('red',
           'red4',
           'orange',
           'lightgoldenrod3',
           'sienna',
           'indianred',
           'orangered1',
           'black',
           'darkturquoise',
           'paleturquoise',
           'lightgreen',
           'springgreen4',
           'darkolivegreen',
           'purple4',
           'purple',
           'deeppink',
           'violetred',
           'violet'),
  #y.max = 3,
  pt.size = 1
)

#Rename Idents
pancreas.integrated <- RenameIdents(pancreas.integrated, 
                                    "0" = "Beta INS-hi", 
                                    "1" = "Alpha GCG-hi",
                                    "2" = "Ductal", 
                                    "3" = "Transdifferentiating Beta",
                                    "4" = "Beta INS-low", 
                                    "5" = "Acinar",
                                    "6" = "Activated Stellate", 
                                    "7" = "Ductal",
                                    "8" = "Endothelial", 
                                    "9" = "Delta",
                                    "10" = "Alpha GCG-hi", 
                                    "11" = "Quiescent Stellate",
                                    "12" = "Alpha GCG-low",
                                    "13" = "Ductal",
                                    "14" = "Quiescent Stellate",
                                    "15" = "Gamma",
                                    "16" = "Macrophage",
                                    "17" = "Proliferating Stellate",
                                    "18" = "Schwann",
                                    "19" = "Mast",
                                    "20" = "T-Lymphocyte"
                                    )

#plot <- DimPlot(pancreas.integrated, reduction = "umap")
DefaultAssay(object = pancreas.integrated) <- "RNA"
Idents(pancreas.integrated, WhichCells(object = pancreas.integrated, expression = GHRL > 1, slot = 'counts')) <- 'Epsilon'
#pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Epsilon")

DimPlot(pancreas.integrated, reduction = "umap", label = TRUE)

# Saving this information in the metadata slot
table(Idents(pancreas.integrated))
pancreas.integrated$celltype <- Idents(pancreas.integrated)
head(pancreas.integrated@meta.data)

# Define an order of cluster identities remember after this step-
# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("Beta INS-hi", "Beta INS-low", "Transdifferentiating Beta", "Alpha GCG-hi", "Alpha GCG-low", "Delta", "Gamma", "Epsilon",
               "Ductal", "Acinar", 
               "Quiescent Stellate", "Activated Stellate", "Proliferating Stellate",
               "Macrophage", "T-Lymphocyte", "Mast",
               "Schwann", "Endothelial")
head(pancreas.integrated@meta.data$celltype)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
pancreas.integrated@meta.data$celltype <- factor(x = pancreas.integrated@meta.data$celltype, levels = my_levels)
Idents(pancreas.integrated) <- "celltype"

# Observing cells
DimPlot(pancreas.integrated, split.by = "sample", group.by = "celltype", label = FALSE, ncol = 2,  cols = c("red",
                                                                                                            "red4",
                                                                                                            "orange",
                                                                                                            "lightgoldenrod3",
                                                                                                            "sienna",
                                                                                                            "indianred",
                                                                                                            "orangered1",
                                                                                                            "black",
                                                                                                            "darkturquoise",
                                                                                                            "paleturquoise",
                                                                                                            "lightgreen",
                                                                                                            "springgreen4",
                                                                                                            "darkolivegreen",
                                                                                                            "purple4",
                                                                                                            "purple",
                                                                                                            "deeppink",
                                                                                                            "violetred",
                                                                                                            "violet"
))

Idents(pancreas.integrated) <- "treatment"
DHT <- subset(pancreas.integrated, idents = "DHT[10nM]")
DimPlot(DHT, group.by = "treatment", cols = "red")
EtOH <- subset(pancreas.integrated, idents = "EtOH")
DimPlot(EtOH, group.by = "treatment", cols = "blue")

DimPlot(pancreas.integrated, group.by = "treatment")
UMAPPlot(pancreas.integrated, reduction = "umap",
        pt.size = .75,
        cols = c("red",
                 "red4",
                 "orange",
                 "lightgoldenrod3",
                 "sienna",
                 "indianred",
                 "orangered1",
                 "black",
                 "darkturquoise",
                 "paleturquoise",
                 "lightgreen",
                 "springgreen4",
                 "darkolivegreen",
                 "purple4",
                 "purple",
                 "deeppink",
                 "violetred",
                 "violet"
                 ),
                 label = FALSE)
#beta.hi <- subset(pancreas.integrated, idents = "Beta INS-hi")

#table(beta.hi$treatment)

DotPlot(pancreas.integrated,
        group.by = "treatment",
        #split.by = "treatment",
        features = c("AR"), 
        cols = c("yellow", "red"), 
        col.min = -10, 
        col.max = 10)

# Advanced coding for ggplot2
# Create a new metadata slot containing combined info, segregating clusters and samples
Idents(object = pancreas.integrated) <- "celltype"

# Select only beta cells
beta.hi <- subset(pancreas.integrated, idents = "Beta INS-hi")

pancreas.integrated$celltype.sample <- paste(Idents(pancreas.integrated),pancreas.integrated$treatment, sep = "_")
table(pancreas.integrated@meta.data$celltype.sample)

# New metadata column is not paired, so we need to pair
my_levels2 <- c("Beta INS-hi_EtOH", "Beta INS-hi_DHT[10nM]", "Beta INS-low_EtOH", "Beta INS-low_DHT[10nM]", "Transdifferentiating Beta_EtOH", "Transdifferentiating Beta_DHT[10nM]", "Alpha GCG-hi_EtOH", "Alpha GCG-hi_DHT[10nM]", "Alpha GCG-low_EtOH", "Alpha GCG-low_DHT[10nM]", "Delta_EtOH", "Delta_DHT[10nM]", "Gamma_EtOH", "Gamma_DHT[10nM]", "Epsilon_EtOH", "Epsilon_DHT[10nM]",
                "Ductal_EtOH", "Ductal_DHT[10nM]", "Acinar_EtOH", "Acinar_DHT[10nM]", 
                "Quiescent Stellate_EtOH", "Quiescent Stellate_DHT[10nM]", "Activated Stellate_EtOH", "Activated Stellate_DHT[10nM]", "Proliferating Stellate_EtOH", "Proliferating Stellate_DHT[10nM]",
                "Macrophage_EtOH", "Macrophage_DHT[10nM]", "T-Lymphocyte_EtOH", "T-Lymphocyte_DHT[10nM]", "Mast_EtOH", "Mast_DHT[10nM]", "Schwann_EtOH", "Schwann_DHT[10nM]", "Endothelial_EtOH", "Endothelial_DHT[10nM]")
head(pancreas.integrated@meta.data$celltype)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
pancreas.integrated@meta.data$celltype.sample <- factor(x = pancreas.integrated@meta.data$celltype.sample, levels = my_levels2)
table(pancreas.integrated@meta.data$celltype.sample)

# Re select organized idents
Idents(pancreas.integrated) <- "celltype.sample"
DefaultAssay(object = pancreas.integrated) <- "RNA"

# Selected genes
markers.to.plot <- c("MT-ATP8", "MT-ATP6",
                     "MT-CO1", "MT-CO2", "MT-CO3",
                     "MT-CYB",
                     "MT-ND6", "MT-ND5", "MT-ND4", "MT-ND4L", "MT-ND3", "MT-ND2", "MT-ND1")

# Dotplot
DotPlot(pancreas.integrated,  
        dot.scale = 8,
        col.min = -1, #minimum level
        col.max = 1,  #maximum level
        features = rev(markers.to.plot)) + 
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
  theme_light() +
  #facet_wrap(~??? what metadata should be here??)
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  scale_colour_gradient2(low =c("dodgerblue"), mid = c("white"), high =c("red3")) +
  guides(color = guide_colorbar(title = 'Average Expression'))



# Select only beta cells
Idents(pancreas.integrated) <- "celltype.sample"
pathway <- subset(pancreas.integrated, idents = c("Beta INS-hi_EtOH", "Beta INS-hi_DHT[10nM]", "Beta INS-low_EtOH", "Beta INS-low_DHT[10nM]", 
                                                  "Alpha GCG-hi_EtOH", "Alpha GCG-hi_DHT[10nM]", "Alpha GCG-low_EtOH", "Alpha GCG-low_DHT[10nM]"))

# Selected genes
Idents(pancreas.integrated) <- "celltype.sample"
markers.to.plot <- c("MT-ATP8", "MT-ATP6",
                     "MT-CO1", "MT-CO2", "MT-CO3",
                     "MT-CYB",
                     "MT-ND6", "MT-ND5", "MT-ND4", "MT-ND4L", "MT-ND3", "MT-ND2", "MT-ND1")

markers.to.plot <- c("GCK", "ALDOA", "BPGM",
                     "ENO1", "ENO2", "GAPDH",
                     "GPI",
                     "PFKL", "PFKM", "PGAM1", "PGK1", "PKM", "TPI1")
markers.to.plot <- c("AR", "SLC25A4")
# Dotplot
DotPlot(pathway,
        dot.scale = 8,
        col.min = -1, #minimum level
        col.max = 1,  #maximum level
        features = rev(markers.to.plot)) + 
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
  theme_light() +
  #facet_wrap(~??? what metadata should be here??)
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  scale_colour_gradient2(low =c("dodgerblue"), mid = c("white"), high =c("red3")) +
  guides(color = guide_colorbar(title = 'Average Expression'))

VlnPlot(pathway, features = markers.to.plot, slot = "counts", #"data" 
        split.by = 'treatment', ncol = 1)

# Visualize co-expression of two features simultaneously
FeaturePlot(pancreas.integrated, features = c("GLP1R", "AR"), blend = TRUE, order = TRUE, blend.threshold = 0.01, pt.size = 2)
FeatureScatter(pancreas.integrated, feature1 = "GLP1R", feature2 = "AR", pt.size = 3,
               cols = c("red",
                        "red4",
                        "orange",
                        "lightgoldenrod3",
                        "sienna",
                        "indianred",
                        "orangered1",
                        "black",
                        "darkturquoise",
                        "paleturquoise",
                        "lightgreen",
                        "springgreen4",
                        "darkolivegreen",
                        "purple4",
                        "purple",
                        "deeppink",
                        "violetred",
                        "violet"
               ))

# Calculate %of GLP1R+ Beta-hi cells expressing AR
Idents(pancreas.integrated) <- "celltype"
beta.hi <- subset(pancreas.integrated, idents = "Beta INS-hi")

Idents(beta.hi) <- "sample"
HP2107001_ctrl_beta.hi <- subset(beta.hi, idents = "HP2107001_ctrl")
Idents(HP2107001_ctrl_beta.hi) <- "celltype"
HP2107001_ctrl_beta.hi.GLP1R <- subset(HP2107001_ctrl_beta.hi, subset = GLP1R > 0)
length(WhichCells(object = HP2107001_ctrl_beta.hi.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107001_ctrl_beta.hi@meta.data)*100

HP2107001_DHT_beta.hi <- subset(beta.hi, idents = "HP2107001_DHT")
Idents(HP2107001_DHT_beta.hi) <- "celltype"
HP2107001_DHT_beta.hi.GLP1R <- subset(HP2107001_DHT_beta.hi, subset = GLP1R > 0)
length(WhichCells(object = HP2107001_DHT_beta.hi.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107001_DHT_beta.hi@meta.data)*100

HP2107701_ctrl_beta.hi <- subset(beta.hi, idents = "HP2107701_ctrl")
Idents(HP2107701_ctrl_beta.hi) <- "celltype"
HP2107701_ctrl_beta.hi.GLP1R <- subset(HP2107701_ctrl_beta.hi, subset = GLP1R > 0)
length(WhichCells(object = HP2107701_ctrl_beta.hi.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107701_ctrl_beta.hi@meta.data)*100

HP2107701_DHT_beta.hi <- subset(beta.hi, idents = "HP2107701_DHT")
Idents(HP2107701_DHT_beta.hi) <- "celltype"
HP2107701_DHT_beta.hi.GLP1R <- subset(HP2107701_DHT_beta.hi, subset = GLP1R > 0)
length(WhichCells(object = HP2107701_DHT_beta.hi.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107701_DHT_beta.hi@meta.data)*100

HP2107901_ctrl_beta.hi <- subset(beta.hi, idents = "HP2107901_ctrl")
Idents(HP2107901_ctrl_beta.hi) <- "celltype"
HP2107901_ctrl_beta.hi.GLP1R <- subset(HP2107901_ctrl_beta.hi, subset = GLP1R > 0)
length(WhichCells(object = HP2107901_ctrl_beta.hi.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107901_ctrl_beta.hi@meta.data)*100

HP2107901_DHT_beta.hi <- subset(beta.hi, idents = "HP2107901_DHT")
Idents(HP2107901_DHT_beta.hi) <- "celltype"
HP2107901_DHT_beta.hi.GLP1R <- subset(HP2107901_DHT_beta.hi, subset = GLP1R > 0)
length(WhichCells(object = HP2107901_DHT_beta.hi.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107901_DHT_beta.hi@meta.data)*100

# Calculate %of GLP1R+ Beta-low cells expressing AR
Idents(pancreas.integrated) <- "celltype"
beta.low <- subset(pancreas.integrated, idents = "Beta INS-low")

Idents(beta.low) <- "sample"
HP2107001_ctrl_beta.low <- subset(beta.low, idents = "HP2107001_ctrl")
Idents(HP2107001_ctrl_beta.low) <- "celltype"
HP2107001_ctrl_beta.low.GLP1R <- subset(HP2107001_ctrl_beta.low, subset = GLP1R > 0)
length(WhichCells(object = HP2107001_ctrl_beta.low.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107001_ctrl_beta.low@meta.data)*100

HP2107001_DHT_beta.low <- subset(beta.low, idents = "HP2107001_DHT")
Idents(HP2107001_DHT_beta.low) <- "celltype"
HP2107001_DHT_beta.low.GLP1R <- subset(HP2107001_DHT_beta.low, subset = GLP1R > 0)
length(WhichCells(object = HP2107001_DHT_beta.low.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107001_DHT_beta.low@meta.data)*100

HP2107701_ctrl_beta.low <- subset(beta.low, idents = "HP2107701_ctrl")
Idents(HP2107701_ctrl_beta.low) <- "celltype"
HP2107701_ctrl_beta.low.GLP1R <- subset(HP2107701_ctrl_beta.low, subset = GLP1R > 0)
length(WhichCells(object = HP2107701_ctrl_beta.low.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107701_ctrl_beta.low@meta.data)*100

HP2107701_DHT_beta.low <- subset(beta.low, idents = "HP2107701_DHT")
Idents(HP2107701_DHT_beta.low) <- "celltype"
HP2107701_DHT_beta.low.GLP1R <- subset(HP2107701_DHT_beta.low, subset = GLP1R > 0)
length(WhichCells(object = HP2107701_DHT_beta.low.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107701_DHT_beta.low@meta.data)*100

HP2107901_ctrl_beta.low <- subset(beta.low, idents = "HP2107901_ctrl")
Idents(HP2107901_ctrl_beta.low) <- "celltype"
HP2107901_ctrl_beta.low.GLP1R <- subset(HP2107901_ctrl_beta.low, subset = GLP1R > 0)
length(WhichCells(object = HP2107901_ctrl_beta.low.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107901_ctrl_beta.low@meta.data)*100

HP2107901_DHT_beta.low <- subset(beta.low, idents = "HP2107901_DHT")
Idents(HP2107901_DHT_beta.low) <- "celltype"
HP2107901_DHT_beta.low.GLP1R <- subset(HP2107901_DHT_beta.low, subset = GLP1R > 0)
length(WhichCells(object = HP2107901_DHT_beta.low.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107901_DHT_beta.low@meta.data)*100

# Calculate %of GLP1R+ Trans-beta cells expressing AR
Idents(pancreas.integrated) <- "celltype"
trans.beta <- subset(pancreas.integrated, idents = "Transdifferentiating Beta")

Idents(trans.beta) <- "sample"
HP2107001_ctrl_trans.beta <- subset(trans.beta, idents = "HP2107001_ctrl")
Idents(HP2107001_ctrl_trans.beta) <- "celltype"
HP2107001_ctrl_trans.beta.GLP1R <- subset(HP2107001_ctrl_trans.beta, subset = GLP1R > 0)
length(WhichCells(object = HP2107001_ctrl_trans.beta.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107001_ctrl_trans.beta@meta.data)*100

HP2107001_DHT_trans.beta <- subset(trans.beta, idents = "HP2107001_DHT")
Idents(HP2107001_DHT_trans.beta) <- "celltype"
HP2107001_DHT_trans.beta.GLP1R <- subset(HP2107001_DHT_trans.beta, subset = GLP1R > 0)
length(WhichCells(object = HP2107001_DHT_trans.beta.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107001_DHT_trans.beta@meta.data)*100

HP2107701_ctrl_trans.beta <- subset(trans.beta, idents = "HP2107701_ctrl")
Idents(HP2107701_ctrl_trans.beta) <- "celltype"
HP2107701_ctrl_trans.beta.GLP1R <- subset(HP2107701_ctrl_trans.beta, subset = GLP1R > 0)
length(WhichCells(object = HP2107701_ctrl_trans.beta.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107701_ctrl_trans.beta@meta.data)*100

HP2107701_DHT_trans.beta <- subset(trans.beta, idents = "HP2107701_DHT")
Idents(HP2107701_DHT_trans.beta) <- "celltype"
HP2107701_DHT_trans.beta.GLP1R <- subset(HP2107701_DHT_trans.beta, subset = GLP1R > 0)
length(WhichCells(object = HP2107701_DHT_trans.beta.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107701_DHT_trans.beta@meta.data)*100

HP2107901_ctrl_trans.beta <- subset(trans.beta, idents = "HP2107901_ctrl")
Idents(HP2107901_ctrl_trans.beta) <- "celltype"
HP2107901_ctrl_trans.beta.GLP1R <- subset(HP2107901_ctrl_trans.beta, subset = GLP1R > 0)
length(WhichCells(object = HP2107901_ctrl_trans.beta.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107901_ctrl_trans.beta@meta.data)*100

HP2107901_DHT_trans.beta <- subset(trans.beta, idents = "HP2107901_DHT")
Idents(HP2107901_DHT_trans.beta) <- "celltype"
HP2107901_DHT_trans.beta.GLP1R <- subset(HP2107901_DHT_trans.beta, subset = GLP1R > 0)
length(WhichCells(object = HP2107901_DHT_trans.beta.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107901_DHT_trans.beta@meta.data)*100

FeatureScatter(beta.hi, feature1 = "GLP1R", feature2 = "AR", pt.size = 2,
               cols = c("red",
                        "red4",
                        "orange",
                        "lightgoldenrod3",
                        "sienna",
                        "indianred",
                        "orangered1",
                        "black",
                        "darkturquoise",
                        "paleturquoise",
                        "lightgreen",
                        "springgreen4",
                        "darkolivegreen",
                        "purple4",
                        "purple",
                        "deeppink",
                        "violetred",
                        "violet"
               ))

FeatureScatter(beta.low, feature1 = "GLP1R", feature2 = "AR", pt.size = 2,
               cols = c("red",
                        "red4",
                        "orange",
                        "lightgoldenrod3",
                        "sienna",
                        "indianred",
                        "orangered1",
                        "black",
                        "darkturquoise",
                        "paleturquoise",
                        "lightgreen",
                        "springgreen4",
                        "darkolivegreen",
                        "purple4",
                        "purple",
                        "deeppink",
                        "violetred",
                        "violet"
               ))

tranbeta <- subset(pancreas.integrated, idents = "Transdifferentiating Beta")
FeatureScatter(tranbeta, feature1 = "GLP1R", feature2 = "AR", pt.size = 2,
               cols = c("red",
                        "red4",
                        "orange",
                        "lightgoldenrod3",
                        "sienna",
                        "indianred",
                        "orangered1",
                        "black",
                        "darkturquoise",
                        "paleturquoise",
                        "lightgreen",
                        "springgreen4",
                        "darkolivegreen",
                        "purple4",
                        "purple",
                        "deeppink",
                        "violetred",
                        "violet"
               ))
# Save file this will change but for showing them on 07132021 its fine
#saveRDS(pancreas.integrated, file = r"(C:/Users/mqadir/Box/Lab 2301/RNAseq DHT data/wkdir/pancreas.integrated.rds)")
#pancreas.integrated <- readRDS(r"(C:/Users/mqadir/Box/Lab 2301/RNAseq DHT data/wkdir/pancreas.integrated.rds)")

# Identify conserved cell markers
DefaultAssay(pancreas.integrated) <- "RNA"
markers.beta.hi <- FindConservedMarkers(pancreas.integrated, ident.1 = "Beta INS-hi", grouping.var = "treatment", verbose = TRUE)
head(markers.beta.hi)

markers.beta.low <- FindConservedMarkers(pancreas.integrated, ident.1 = "Beta INS-low", grouping.var = "treatment", verbose = TRUE)
head(markers.beta.low)

markers.alpha.hi <- FindConservedMarkers(pancreas.integrated, ident.1 = "Alpha GCG-hi", grouping.var = "treatment", verbose = TRUE)
head(markers.alpha.hi)

markers.alpha.low <- FindConservedMarkers(pancreas.integrated, ident.1 = "Alpha GCG-low", grouping.var = "treatment", verbose = TRUE)
head(markers.alpha.low)

markers.transbeta <- FindConservedMarkers(pancreas.integrated, ident.1 = "Transdifferentiating Beta", grouping.var = "treatment", verbose = TRUE)
head(markers.transbeta)

markers.delta <- FindConservedMarkers(pancreas.integrated, ident.1 = "Delta", grouping.var = "treatment", verbose = TRUE)
head(markers.delta)

markers.gamma <- FindConservedMarkers(pancreas.integrated, ident.1 = "Gamma", grouping.var = "treatment", verbose = TRUE)
head(markers.gamma)

markers.epsilon <- FindConservedMarkers(pancreas.integrated, ident.1 = "Epsilon", grouping.var = "treatment", verbose = TRUE)
head(markers.epsilon)

markers.ductal <- FindConservedMarkers(pancreas.integrated, ident.1 = "Ductal", grouping.var = "treatment", verbose = TRUE)
head(markers.ductal)

markers.acinar <- FindConservedMarkers(pancreas.integrated, ident.1 = "Acinar", grouping.var = "treatment", verbose = TRUE)
head(markers.acinar)

markers.quiescentstellate <- FindConservedMarkers(pancreas.integrated, ident.1 = "Quiescent Stellate", grouping.var = "treatment", verbose = TRUE)
head(markers.quiescentstellate)

markers.activatedstellate <- FindConservedMarkers(pancreas.integrated, ident.1 = "Activated Stellate", grouping.var = "treatment", verbose = TRUE)
head(markers.activatedstellate)

markers.prolifstellate <- FindConservedMarkers(pancreas.integrated, ident.1 = "Proliferating Stellate", grouping.var = "treatment", verbose = TRUE)
head(markers.prolifstellate)

markers.macrophage <- FindConservedMarkers(pancreas.integrated, ident.1 = "Macrophage", grouping.var = "treatment", verbose = TRUE)
head(markers.macrophage)

markers.tlympho <- FindConservedMarkers(pancreas.integrated, ident.1 = "T-Lymphocyte", grouping.var = "treatment", verbose = TRUE)
head(markers.tlympho)

markers.mast <- FindConservedMarkers(pancreas.integrated, ident.1 = "Mast", grouping.var = "treatment", verbose = TRUE)
head(markers.mast)

markers.schwann <- FindConservedMarkers(pancreas.integrated, ident.1 = "Schwann", grouping.var = "treatment", verbose = TRUE)
head(markers.schwann)

markers.endothelial <- FindConservedMarkers(pancreas.integrated, ident.1 = "Endothelial", grouping.var = "treatment", verbose = TRUE)
head(markers.endothelial)

# Now over-write the SCT assay with new-analyzed data from this subsetted data
pancreas.integrated <- SCTransform(pancreas.integrated, assay = "RNA", new.assay.name = "SCT", verbose = TRUE, return.only.var.genes = TRUE)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
# Here we define a DE gene as a gene which has:
# Fold Change of >1.1x
# Atleast 10% of cells express that gene
Idents(object = pancreas.integrated) <- "celltype"
DefaultAssay(object = pancreas.integrated) <- "RNA"
pancreas.integrated.markers <- FindAllMarkers(object = pancreas.integrated, 
                                                  features = VariableFeatures(pancreas.integrated, assay = 'integrated'), 
                                                  only.pos = TRUE, 
                                                  min.pct = 0.1, 
                                                  logfc.threshold = 0.137504, 
                                                  assay = 'RNA',
                                                  slot = c('data'))

pancreas.integrated.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.csv(pancreas.integrated.markers, r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\pancreas.integrated.markers.RNA.csv)")

# Look at your default assay
DefaultAssay(object = pancreas.integrated)

# Change default assay to SCT, save information in the "SCT" assay
# You can toggle between integrated, SCT and RNA to see different expression profiles/different normalizations
DefaultAssay(object = pancreas.integrated) <- "integrated"

# Create heatmap using doheatmap
top10.nomes <- pancreas.integrated.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(object = pancreas.integrated, 
          features = top10.nomes$gene, 
          disp.min = -1, 
          disp.max = 1,
          label = FALSE) + scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", 
                                                                            "#fbfcbd", 
                                                                            "#ff0000"))(256))
# Identify conserved cell markers
Idents(pancreas.integrated) <- factor(Idents(pancreas.integrated), levels = c("Beta INS-hi", "Beta INS-low", "Transdifferentiating Beta", "Alpha GCG-hi", "Alpha GCG-low", "Delta", "Gamma", "Epsilon",
                                                                              "Ductal", "Acinar", 
                                                                              "Quiescent Stellate", "Activated Stellate", "Proliferating Stellate", 
                                                                              "Macrophage", "T-Lymphocyte", "Mast", "Schwann", "Endothelial"))
markers.to.plot <- c("INS", "IAPP", "NKX6-1", "MAFA", "MAFB", "GCG", "DPP4", "GC", "LEPR", "SST", "FRZB", "PPY", "CALB1", "THSD7A", "GHRL", "PHGR1",
                     "CFTR", "KRT19", "MMP7", "CELA2A", "CELA2B", "CELA3A", "RGS5", "CSRP2", "FABP4", "COL3A1", "FMOD", "PDGFRB", "MKI67", "HIST1H4C", "STMN1", 
                     "CD86", "CSF1R", "SDS", "NKG7", "IL2RB", "CCL5", "RGS13", "TPSB2", "TPSAB1", "SOX10", "CDH19", "NGFR", "CD34", "ENG", "VWF", "UCN3")
markers.to.plot <- c("INS", "IAPP", "PDX1", "MAFA", "MAFB", "GCG", "DPP4", "GC", "LEPR", "SST", "PPY", "THSD7A", "GHRL", "FRZB",
                     "CFTR", "MMP7", "CELA2A", "CELA3A", "RGS5", "FABP4", "COL3A1", "FMOD", "MKI67", "STMN1", 
                     "CSF1R", "SDS", "NKG7", "CCL5", "TPSB2", "TPSAB1", "SOX10", "NGFR", "ENG", "VWF")

# Advanced coding for ggplot2
# Create a new metadata slot containing combined info, segregating clusters and samples
Idents(object = pancreas.integrated) <- "celltype"
pancreas.integrated$celltype.sample <- paste(Idents(pancreas.integrated),pancreas.integrated$treatment, sep = "_")
table(pancreas.integrated@meta.data$celltype.sample)

# New metadata column is not paired, so we need to pair
my_levels2 <- c("Beta INS-hi_EtOH", "Beta INS-hi_DHT[10nM]", "Beta INS-low_EtOH", "Beta INS-low_DHT[10nM]", "Transdifferentiating Beta_EtOH", "Transdifferentiating Beta_DHT[10nM]", "Alpha GCG-hi_EtOH", "Alpha GCG-hi_DHT[10nM]", "Alpha GCG-low_EtOH", "Alpha GCG-low_DHT[10nM]", "Delta_EtOH", "Delta_DHT[10nM]", "Gamma_EtOH", "Gamma_DHT[10nM]", "Epsilon_EtOH", "Epsilon_DHT[10nM]",
                "Ductal_EtOH", "Ductal_DHT[10nM]", "Acinar_EtOH", "Acinar_DHT[10nM]", 
                "Quiescent Stellate_EtOH", "Quiescent Stellate_DHT[10nM]", "Activated Stellate_EtOH", "Activated Stellate_DHT[10nM]", "Proliferating Stellate_EtOH", "Proliferating Stellate_DHT[10nM]",
                "Macrophage_EtOH", "Macrophage_DHT[10nM]", "T-Lymphocyte_EtOH", "T-Lymphocyte_DHT[10nM]", "Mast_EtOH", "Mast_DHT[10nM]", "Schwann_EtOH", "Schwann_DHT[10nM]", "Endothelial_EtOH", "Endothelial_DHT[10nM]")
head(pancreas.integrated@meta.data$celltype)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
pancreas.integrated@meta.data$celltype.sample <- factor(x = pancreas.integrated@meta.data$celltype.sample, levels = my_levels2)
table(pancreas.integrated@meta.data$celltype.sample)

# Re select organized idents
Idents(pancreas.integrated) <- "celltype.sample"
Idents(pancreas.integrated) <- "celltype"
DefaultAssay(object = pancreas.integrated) <- "RNA"
DotPlot(pancreas.integrated,  
        dot.scale = 8,
        scale = TRUE,
        col.min = -2, #minimum level
        col.max = 3,  #maximum level
        features = rev(markers.to.plot)) + 
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
  theme_light() +
  #facet_wrap(~??? what metadata should be here??)
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  scale_colour_gradient2(low =c("dodgerblue"), mid =c("white"), high =c("red")) +
  guides(color = guide_colorbar(title = 'Average Expression'))

# Older dotplot configuration, shows only percentage not expression
Idents(pancreas.integrated) <- "celltype"
DotPlot(pancreas.integrated, features = rev(markers.to.plot), 
        cols = c("blue", "red"), 
        dot.scale = 8, 
        split.by = "treatment") + 
  RotatedAxis() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  theme_light() + 
  #coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =10, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =8, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=10), 
        legend.text=element_text(size=10))

# Diff gene testing across conditions
# pancreas.integrated$treatment.dht <- paste(Idents(pancreas.integrated), pancreas.integrated$treatment, sep = "_")
# pancreas.integrated$celltype.split <- Idents(pancreas.integrated)
# choosing only those genes which are differentially expressed
# Optimise idents and assay
Idents(pancreas.integrated) <- "celltype.sample"
DefaultAssay(object = pancreas.integrated) <- "RNA"

# 1.Beta-cells (INS Hi)
beta.INSHi.DHT.response <- FindMarkers(pancreas.integrated, 
                                 ident.1 = "Beta INS-hi_DHT[10nM]", ident.2 = "Beta INS-hi_EtOH", 
                                 test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                 min.pct = 0.1,
                                 logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                 pseudocount.use = 1,
                                 verbose = TRUE)
head(beta.INSHi.DHT.response, n = 15)
write.csv(beta.INSHi.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\beta.INShi.DHT.response.csv)")

# 2.Beta-cells (INS low)
beta.INSLow.DHT.response <- FindMarkers(pancreas.integrated, 
                                       ident.1 = "Beta INS-low_DHT[10nM]", ident.2 = "Beta INS-low_EtOH", 
                                       test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                       min.pct = 0.1,
                                       logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                       pseudocount.use = 1,
                                       verbose = TRUE)
head(beta.INSLow.DHT.response, n = 15)
write.csv(beta.INSLow.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\beta.INSLow.DHT.response.csv)")

# 3.Alpha-cells (GCG hi)
alpha.GCGHi.DHT.response <- FindMarkers(pancreas.integrated, 
                                 ident.1 = "Alpha GCG-hi_DHT[10nM]", ident.2 = "Alpha GCG-hi_EtOH", 
                                 test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                 min.pct = 0.1,
                                 logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                 pseudocount.use = 1,
                                 verbose = TRUE)
head(alpha.GCGHi.DHT.response, n = 15)
write.csv(alpha.GCGHi.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\alpha.GCGHi.DHT.response.csv)")

# 4.Alpha-cells (GCG low)
alpha.GCGLow.DHT.response <- FindMarkers(pancreas.integrated, 
                                        ident.1 = "Alpha GCG-low_DHT[10nM]", ident.2 = "Alpha GCG-low_EtOH", 
                                        test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                        min.pct = 0.1,
                                        logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                        pseudocount.use = 1,
                                        verbose = TRUE)
head(alpha.GCGLow.DHT.response, n = 15)
write.csv(alpha.GCGLow.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\alpha.GCGLow.DHT.response.csv)")

# 5.Trandifferentiating Endocrine-Cells
tranbeta.DHT.response <- FindMarkers(pancreas.integrated, 
                                  ident.1 = "Transdifferentiating Beta_DHT[10nM]", ident.2 = "Transdifferentiating Beta_EtOH", 
                                  test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                  min.pct = 0.1,
                                  logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                  pseudocount.use = 1,
                                  verbose = TRUE)
head(tranbeta.DHT.response, n = 15)
write.csv(tranbeta.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\tranbeta.DHT.response.csv)")

# 6.Delta-Cells
delta.DHT.response <- FindMarkers(pancreas.integrated, 
                                  ident.1 = "Delta_DHT[10nM]", ident.2 = "Delta_EtOH", 
                                  test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                  min.pct = 0.1,
                                  logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                  pseudocount.use = 1,
                                  verbose = TRUE)
head(delta.DHT.response, n = 15)
write.csv(delta.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\delta.DHT.response.csv)")

# 7.Gamma-Cells
gamma.DHT.response <- FindMarkers(pancreas.integrated, 
                                  ident.1 = "Gamma_DHT[10nM]", ident.2 = "Gamma_EtOH", 
                                  test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                  min.pct = 0.1,
                                  logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                  pseudocount.use = 1,
                                  verbose = TRUE)
head(gamma.DHT.response, n = 15)
write.csv(gamma.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\gamma.DHT.response.csv)")

# 7.Epsilon-Cells
epsilon.DHT.response <- FindMarkers(pancreas.integrated, 
                                  ident.1 = "Epsilon_DHT[10nM]", ident.2 = "Epsilon_EtOH", 
                                  test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                  min.pct = 0.1,
                                  logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                  pseudocount.use = 1,
                                  verbose = TRUE)
head(epsilon.DHT.response, n = 15)
write.csv(epsilon.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\epsilon.DHT.response.csv)")

# 8.Ductal-Cells
ductal.DHT.response <- FindMarkers(pancreas.integrated, 
                                  ident.1 = "Ductal_DHT[10nM]", ident.2 = "Ductal_EtOH", 
                                  test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                  min.pct = 0.1,
                                  logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                  pseudocount.use = 1,
                                  verbose = TRUE)
head(ductal.DHT.response, n = 15)
write.csv(ductal.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\ductal.DHT.response.csv)")

# 9.Acinar-Cells
acinar.DHT.response <- FindMarkers(pancreas.integrated, 
                                   ident.1 = "Acinar_DHT[10nM]", ident.2 = "Acinar_EtOH", 
                                   test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                   min.pct = 0.1,
                                   logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                   pseudocount.use = 1,
                                   verbose = TRUE)
head(acinar.DHT.response, n = 15)
write.csv(acinar.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\acinar.DHT.response.csv)")

# 10.Quiescent Stellate-Cells
qstellate.DHT.response <- FindMarkers(pancreas.integrated, 
                                   ident.1 = "Quiescent Stellate_DHT[10nM]", ident.2 = "Quiescent Stellate_EtOH", 
                                   test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                   min.pct = 0.1,
                                   logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                   pseudocount.use = 1,
                                   verbose = TRUE)
head(qstellate.DHT.response, n = 15)
write.csv(qstellate.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\qstellate.DHT.response.csv)")

# 11.Activated Stellate-Cells
astellate.DHT.response <- FindMarkers(pancreas.integrated, 
                                      ident.1 = "Activated Stellate_DHT[10nM]", ident.2 = "Activated Stellate_EtOH", 
                                      test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                      min.pct = 0.1,
                                      logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                      pseudocount.use = 1,
                                      verbose = TRUE)
head(astellate.DHT.response, n = 15)
write.csv(astellate.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\astellate.DHT.response.csv)")

# 12.Proliferating Stellate-Cells
pstellate.DHT.response <- FindMarkers(pancreas.integrated, 
                                      ident.1 = "Proliferating Stellate_DHT[10nM]", ident.2 = "Proliferating Stellate_EtOH", 
                                      test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                      min.pct = 0.1,
                                      logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                      pseudocount.use = 1,
                                      verbose = TRUE)
head(pstellate.DHT.response, n = 15)
write.csv(pstellate.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\pstellate.DHT.response.csv)")

# 13.Macrophage-Cells
macrophage.DHT.response <- FindMarkers(pancreas.integrated, 
                                      ident.1 = "Macrophage_DHT[10nM]", ident.2 = "Macrophage_EtOH", 
                                      test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                      min.pct = 0.1,
                                      logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                      pseudocount.use = 1,
                                      verbose = TRUE)
head(macrophage.DHT.response, n = 15)
write.csv(macrophage.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\macrophage.DHT.response.csv)")

# 14.T Lymphocyte-Cells
tlympho.DHT.response <- FindMarkers(pancreas.integrated, 
                                      ident.1 = "T-Lymphocyte_DHT[10nM]", ident.2 = "T-Lymphocyte_EtOH", 
                                      test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                      min.pct = 0.1,
                                      logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                      pseudocount.use = 1,
                                      verbose = TRUE)
head(tlympho.DHT.response, n = 15)
write.csv(tlympho.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\tlympho.DHT.response.csv)")

# 15.Mast-Cells
mast.DHT.response <- FindMarkers(pancreas.integrated, 
                                      ident.1 = "Mast_DHT[10nM]", ident.2 = "Mast_EtOH", 
                                      test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                      min.pct = 0.1,
                                      logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                      pseudocount.use = 1,
                                      verbose = TRUE)
head(mast.DHT.response, n = 15)
write.csv(mast.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\mast.DHT.response.csv)")

# 16.Schwann-Cells
schwann.DHT.response <- FindMarkers(pancreas.integrated, 
                                 ident.1 = "Schwann_DHT[10nM]", ident.2 = "Schwann_EtOH", 
                                 test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                 min.pct = 0.1,
                                 logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                 pseudocount.use = 1,
                                 verbose = TRUE)
head(schwann.DHT.response, n = 15)
write.csv(schwann.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\schwann.DHT.response.csv)")

# 17.Endothelial-Cells
endothelial.DHT.response <- FindMarkers(pancreas.integrated, 
                                 ident.1 = "Endothelial_DHT[10nM]", ident.2 = "Endothelial_EtOH", 
                                 test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                 min.pct = 0.1,
                                 logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                 pseudocount.use = 1,
                                 verbose = TRUE)
head(endothelial.DHT.response, n = 15)
write.csv(endothelial.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\endothelial.DHT.response.csv)")

#
#
#
#
#
# Running DE for all genes irrlevant of FC and PCT filtering
# Optimise idents
Idents(pancreas.integrated) <- "celltype.sample"
DefaultAssay(object = pancreas.integrated) <- "RNA"

# 1.Beta-cells (INS Hi)
beta.INSHi.DHT.response <- FindMarkers(pancreas.integrated, 
                                       ident.1 = "Beta INS-hi_DHT[10nM]", ident.2 = "Beta INS-hi_EtOH", 
                                       test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                       min.pct = 0,
                                       logfc.threshold = 0, 
                                       pseudocount.use = 1,
                                       verbose = TRUE)
head(beta.INSHi.DHT.response, n = 15)
write.csv(beta.INSHi.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1beta.INSHi.DHT.response.csv)")

# 2.Beta-cells (INS low)
beta.INSLow.DHT.response <- FindMarkers(pancreas.integrated, 
                                        ident.1 = "Beta INS-low_DHT[10nM]", ident.2 = "Beta INS-low_EtOH", 
                                        test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                        min.pct = 0,
                                        logfc.threshold = 0, 
                                        pseudocount.use = 1,
                                        verbose = TRUE)
head(beta.INSLow.DHT.response, n = 15)
write.csv(beta.INSLow.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1beta.INSLow.DHT.response.csv)")

# 3.Alpha-cells (GCG hi)
alpha.GCGHi.DHT.response <- FindMarkers(pancreas.integrated, 
                                        ident.1 = "Alpha GCG-hi_DHT[10nM]", ident.2 = "Alpha GCG-hi_EtOH", 
                                        test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                        min.pct = 0,
                                        logfc.threshold = 0,
                                        pseudocount.use = 1,
                                        verbose = TRUE)
head(alpha.GCGHi.DHT.response, n = 15)
write.csv(alpha.GCGHi.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1alpha.GCGHi.DHT.response.csv)")

# 4.Alpha-cells (GCG low)
alpha.GCGLow.DHT.response <- FindMarkers(pancreas.integrated, 
                                         ident.1 = "Alpha GCG-low_DHT[10nM]", ident.2 = "Alpha GCG-low_EtOH", 
                                         test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                         min.pct = 0,
                                         logfc.threshold = 0, 
                                         pseudocount.use = 1,
                                         verbose = TRUE)
head(alpha.GCGLow.DHT.response, n = 15)
write.csv(alpha.GCGLow.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1alpha.GCGLow.DHT.response.csv)")

# 5.Trandifferentiating Endocrine-Cells
tranbeta.DHT.response <- FindMarkers(pancreas.integrated, 
                                     ident.1 = "Transdifferentiating Beta_DHT[10nM]", ident.2 = "Transdifferentiating Beta_EtOH", 
                                     test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                     min.pct = 0,
                                     logfc.threshold = 0,
                                     pseudocount.use = 1,
                                     verbose = TRUE)
head(tranbeta.DHT.response, n = 15)
write.csv(tranbeta.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1tranbeta.DHT.response.csv)")

# 6.Delta-Cells
delta.DHT.response <- FindMarkers(pancreas.integrated, 
                                  ident.1 = "Delta_DHT[10nM]", ident.2 = "Delta_EtOH", 
                                  test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                  min.pct = 0,
                                  logfc.threshold = 0, 
                                  pseudocount.use = 1,
                                  verbose = TRUE)
head(delta.DHT.response, n = 15)
write.csv(delta.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1delta.DHT.response.csv)")

# 7.Gamma-Cells
gamma.DHT.response <- FindMarkers(pancreas.integrated, 
                                  ident.1 = "Gamma_DHT[10nM]", ident.2 = "Gamma_EtOH", 
                                  test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                  min.pct = 0,
                                  logfc.threshold = 0, 
                                  pseudocount.use = 1,
                                  verbose = TRUE)
head(gamma.DHT.response, n = 15)
write.csv(gamma.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1gamma.DHT.response.csv)")

# 7.Epsilon-Cells
epsilon.DHT.response <- FindMarkers(pancreas.integrated, 
                                    ident.1 = "Epsilon_DHT[10nM]", ident.2 = "Epsilon_EtOH", 
                                    test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                    min.pct = 0,
                                    logfc.threshold = 0, # based on output log2 so 0.137504 is ~1.1 FC
                                    pseudocount.use = 1,
                                    verbose = TRUE)
head(epsilon.DHT.response, n = 15)
write.csv(epsilon.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\epsilon.DHT.response.csv)")

# 8.Ductal-Cells
ductal.DHT.response <- FindMarkers(pancreas.integrated, 
                                   ident.1 = "Ductal_DHT[10nM]", ident.2 = "Ductal_EtOH", 
                                   test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                   min.pct = 0,
                                   logfc.threshold = 0, 
                                   pseudocount.use = 1,
                                   verbose = TRUE)
head(ductal.DHT.response, n = 15)
write.csv(ductal.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1ductal.DHT.response.csv)")

# 9.Acinar-Cells
acinar.DHT.response <- FindMarkers(pancreas.integrated, 
                                   ident.1 = "Acinar_DHT[10nM]", ident.2 = "Acinar_EtOH", 
                                   test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                   min.pct = 0,
                                   logfc.threshold = 0, 
                                   pseudocount.use = 1,
                                   verbose = TRUE)
head(acinar.DHT.response, n = 15)
write.csv(acinar.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1acinar.DHT.response.csv)")

# 10.Quiescent Stellate-Cells
qstellate.DHT.response <- FindMarkers(pancreas.integrated, 
                                      ident.1 = "Quiescent Stellate_DHT[10nM]", ident.2 = "Quiescent Stellate_EtOH", 
                                      test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                      min.pct = 0,
                                      logfc.threshold = 0, 
                                      pseudocount.use = 1,
                                      verbose = TRUE)
head(qstellate.DHT.response, n = 15)
write.csv(qstellate.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1qstellate.DHT.response.csv)")

# 11.Activated Stellate-Cells
astellate.DHT.response <- FindMarkers(pancreas.integrated, 
                                      ident.1 = "Activated Stellate_DHT[10nM]", ident.2 = "Activated Stellate_EtOH", 
                                      test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                      min.pct = 0,
                                      logfc.threshold = 0, 
                                      pseudocount.use = 1,
                                      verbose = TRUE)
head(astellate.DHT.response, n = 15)
write.csv(astellate.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1astellate.DHT.response.csv)")

# 12.Proliferating Stellate-Cells
pstellate.DHT.response <- FindMarkers(pancreas.integrated, 
                                      ident.1 = "Proliferating Stellate_DHT[10nM]", ident.2 = "Proliferating Stellate_EtOH", 
                                      test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                      min.pct = 0,
                                      logfc.threshold = 0, 
                                      pseudocount.use = 1,
                                      verbose = TRUE)
head(pstellate.DHT.response, n = 15)
write.csv(pstellate.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1pstellate.DHT.response.csv)")

# 13.Macrophage-Cells
macrophage.DHT.response <- FindMarkers(pancreas.integrated, 
                                       ident.1 = "Macrophage_DHT[10nM]", ident.2 = "Macrophage_EtOH", 
                                       test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                       min.pct = 0,
                                       logfc.threshold = 0, # based on output log2 so 0.137504 is ~1.1 FC
                                       pseudocount.use = 1,
                                       verbose = TRUE)
head(macrophage.DHT.response, n = 15)
write.csv(macrophage.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1macrophage.DHT.response.csv)")

# 14.T Lymphocyte-Cells
tlympho.DHT.response <- FindMarkers(pancreas.integrated, 
                                    ident.1 = "T-Lymphocyte_DHT[10nM]", ident.2 = "T-Lymphocyte_EtOH", 
                                    test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                    min.pct = 0,
                                    logfc.threshold = 0, # based on output log2 so 0.137504 is ~1.1 FC
                                    pseudocount.use = 1,
                                    verbose = TRUE)
head(tlympho.DHT.response, n = 15)
write.csv(tlympho.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1tlympho.DHT.response.csv)")

# 15.Mast-Cells
mast.DHT.response <- FindMarkers(pancreas.integrated, 
                                 ident.1 = "Mast_DHT[10nM]", ident.2 = "Mast_EtOH", 
                                 test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                 min.pct = 0,
                                 logfc.threshold = 0, # based on output log2 so 0.137504 is ~1.1 FC
                                 pseudocount.use = 1,
                                 verbose = TRUE)
head(mast.DHT.response, n = 15)
write.csv(mast.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1mast.DHT.response.csv)")

# 16.Schwann-Cells
schwann.DHT.response <- FindMarkers(pancreas.integrated, 
                                    ident.1 = "Schwann_DHT[10nM]", ident.2 = "Schwann_EtOH", 
                                    test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                    min.pct = 0,
                                    logfc.threshold = 0, # based on output log2 so 0.137504 is ~1.1 FC
                                    pseudocount.use = 1,
                                    verbose = TRUE)
head(schwann.DHT.response, n = 15)
write.csv(schwann.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1schwann.DHT.response.csv)")

# 17.Endothelial-Cells
endothelial.DHT.response <- FindMarkers(pancreas.integrated, 
                                        ident.1 = "Endothelial_DHT[10nM]", ident.2 = "Endothelial_EtOH", 
                                        test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                        min.pct = 0,
                                        logfc.threshold = 0, # based on output log2 so 0.137504 is ~1.1 FC
                                        pseudocount.use = 1,
                                        verbose = TRUE)
head(endothelial.DHT.response, n = 15)
write.csv(endothelial.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1endothelial.DHT.response.csv)")


# Plotting DE genes ###
Idents(pancreas.integrated) <- "celltype"
beta.cells <- subset(pancreas.integrated, idents = "Beta INS-hi")
plots <- VlnPlot(beta.cells, features = c("INS", "DDIT3", "MIF", "DEPP1", "PLCG2", "IAPP"), group.by = "treatment", 
                 pt.size = 0, combine = TRUE)
plots <- VlnPlot(beta.cells, features = c("MT-CO3", "MT-ND1", "MT-ND4", "MT-ATP6", "MT-CO1", "MT-CYB"), group.by = "treatment", 
                 pt.size = 0, combine = TRUE)
wrap_plots(plots = plots, nrow = 1, ncol = 1)

# Load data
volcanodat <- read.csv(r"(C:\Users\mqadir\Box\!FAHD\1. AR-DHT Project\DHT_scRNAseq_Islets\1. DGE_analysis\1. All Genes\1alpha.GCGHi.DHT.response.csv)",
                                    header = TRUE, sep = ",", row.names = 1)
#volcanodat <- epsilon.DHT.response

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# set the base colour as 'black'
keyvals <- rep('black', nrow(volcanodat))

# set the base name/label as 'Mid'
names(keyvals) <- rep('Mid', nrow(volcanodat))

# modify keyvals for variables with fold change > 1.1
keyvals[which(volcanodat$avg_log2FC > 0.137504 & volcanodat$p_val < 0.05)] <- 'red'
names(keyvals)[which(volcanodat$avg_log2FC > 0.137504 & volcanodat$p_val < 0.05)] <- 'high'

# modify keyvals for variables with fold change < -1.1
keyvals[which(volcanodat$avg_log2FC < -0.137504 & volcanodat$p_val < 0.05)] <- 'royalblue'
names(keyvals)[which(volcanodat$avg_log2FC < -0.137504 & volcanodat$p_val < 0.05)] <- 'low'

unique(names(keyvals))

unique(keyvals)
keyvals[1:20]
options(ggrepel.max.overlaps = Inf)
EnhancedVolcano(volcanodat,
                lab = rownames(volcanodat),
                x = 'avg_log2FC',
                y = 'p_val',
                #selectLab = FALSE,
                #selectLab = rownames(volcanodat)[which(names(keyvals) %in% c('high', 'low'))],
                selectLab = c('ATP5F1E', 'ATP1B1', 'COX17', 'MT-ND4L',
                              'ATP5MC1', 'ATP5MD', 'COX7A1', 'ATP5ME', 'UQCR10',
                              'COX7A2', 'NDUFC1', 'LDHA', 'COX8A', 'COX7B', 'COX6C', 'NDUFC2', 
                              'NDUFAB1', 'UQCRQ',
                              'MT-CO3', 'MT-ND1','MT-ND4', 'MT-Co1', 'MT-ATP6', 'MT-CYB', 'MT-ND2',
                              'MT-CO2', 'MT-ND5', 'MT-ND3', 'PCSK', 'SOD2', 'MT-ND6', 'ATP9A'), # use this for labelling genes on plot
                #boxedLabels = TRUE,
                xlim = c(-1,0.8),
                ylim = c(0,300),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Custom colour over-ride',
                pCutoff = 0.05,
                FCcutoff = 0.137504,
                #pointSize = c(ifelse(volcanodat$avg_log2FC < -1 | volcanodat$avg_log2FC > 1, 6, 4)),
                pointSize = 3,
                labSize = 5,
                labFace = 'bold',
                #boxedLabels = TRUE,
                labCol = 'black',
                shape = c(20, 20, 20, 20),
                colCustom = keyvals,
                colAlpha = 2,
                legendPosition = 'right',
                legendLabSize = 15,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 1,
                typeConnectors = 'closed', 
                colConnectors = 'black',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                border = 'partial',
                borderWidth = 1.5,
                borderColour = 'black') + theme(axis.text.x = element_text(colour = "black"),
                                                axis.text.y = element_text(colour = "black"),
                                                axis.title.x = element_text(colour = "black"),
                                                axis.title.y = element_text(colour = "black"))


# Load data
volcanodat <- read.csv(r"(C:\Users\mqadir\Box\!FAHD\1. AR-DHT Project\DHT_scRNAseq_Islets\1. DGE_analysis\1. All Genes\1alpha.GCGLow.DHT.response.csv)",
                       header = TRUE, sep = ",", row.names = 1)
#volcanodat <- epsilon.DHT.response

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# set the base colour as 'black'
keyvals <- rep('black', nrow(volcanodat))

# set the base name/label as 'Mid'
names(keyvals) <- rep('Mid', nrow(volcanodat))

# modify keyvals for variables with fold change > 1.1
keyvals[which(volcanodat$avg_log2FC > 0.137504 & volcanodat$p_val < 0.05)] <- 'red'
names(keyvals)[which(volcanodat$avg_log2FC > 0.137504 & volcanodat$p_val < 0.05)] <- 'high'

# modify keyvals for variables with fold change < -1.1
keyvals[which(volcanodat$avg_log2FC < -0.137504 & volcanodat$p_val < 0.05)] <- 'royalblue'
names(keyvals)[which(volcanodat$avg_log2FC < -0.137504 & volcanodat$p_val < 0.05)] <- 'low'

unique(names(keyvals))

unique(keyvals)
keyvals[1:20]
options(ggrepel.max.overlaps = Inf)
EnhancedVolcano(volcanodat,
                lab = rownames(volcanodat),
                x = 'avg_log2FC',
                y = 'p_val',
                #selectLab = FALSE,
                #selectLab = rownames(volcanodat)[which(names(keyvals) %in% c('high', 'low'))],
                selectLab = c('PCSK1N', 'COX4I1', 'ATP5F1B', 'MT-CO3',
                              'MT-ND4L', 'MT-ATP8', 'MT-ND5', 'COX20', 'MT-ND6',
                              'ABCC8', 'ATP1B1', 'ATP5MC1'), # use this for labelling genes on plot
                #boxedLabels = TRUE,
                xlim = c(-1,1.6),
                ylim = c(0,8),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Custom colour over-ride',
                pCutoff = 0.05,
                FCcutoff = 0.137504,
                #pointSize = c(ifelse(volcanodat$avg_log2FC < -1 | volcanodat$avg_log2FC > 1, 6, 4)),
                pointSize = 3,
                labSize = 5,
                labFace = 'bold',
                #boxedLabels = TRUE,
                labCol = 'black',
                shape = c(20, 20, 20, 20),
                colCustom = keyvals,
                colAlpha = 2,
                legendPosition = 'right',
                legendLabSize = 15,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 1,
                typeConnectors = 'closed', 
                colConnectors = 'black',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                border = 'partial',
                borderWidth = 1.5,
                borderColour = 'black') + theme(axis.text.x = element_text(colour = "black"),
                                                axis.text.y = element_text(colour = "black"),
                                                axis.title.x = element_text(colour = "black"),
                                                axis.title.y = element_text(colour = "black"))

# Load data
volcanodat <- read.csv(r"(C:\Users\mqadir\Box\!FAHD\1. AR-DHT Project\DHT_scRNAseq_Islets\1. DGE_analysis\1. All Genes\1beta.INSHi.DHT.response.csv)",
                       header = TRUE, sep = ",", row.names = 1)
#volcanodat <- epsilon.DHT.response

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# set the base colour as 'black'
keyvals <- rep('black', nrow(volcanodat))

# set the base name/label as 'Mid'
names(keyvals) <- rep('Mid', nrow(volcanodat))

# modify keyvals for variables with fold change > 1.1
keyvals[which(volcanodat$avg_log2FC > 0.137504 & volcanodat$p_val < 0.05)] <- 'red'
names(keyvals)[which(volcanodat$avg_log2FC > 0.137504 & volcanodat$p_val < 0.05)] <- 'high'

# modify keyvals for variables with fold change < -1.1
keyvals[which(volcanodat$avg_log2FC < -0.137504 & volcanodat$p_val < 0.05)] <- 'royalblue'
names(keyvals)[which(volcanodat$avg_log2FC < -0.137504 & volcanodat$p_val < 0.05)] <- 'low'

unique(names(keyvals))

unique(keyvals)
keyvals[1:20]
options(ggrepel.max.overlaps = Inf)
EnhancedVolcano(volcanodat,
                lab = rownames(volcanodat),
                x = 'avg_log2FC',
                y = 'p_val',
                #selectLab = FALSE,
                #selectLab = rownames(volcanodat)[which(names(keyvals) %in% c('high', 'low'))],
                selectLab = c('TPI1', 'IAPP', 'PGK1', 'NPY', 'CTNNB1', 'VEGFA', 'PKM',
                              'MT-CO3', 'MT-ND1', 'MT-ND4', 'MT-CO1', 'MT-ATP6', 'MT-CYB', 'MT-ND2', 'MT-CO2', 'MT-ND3', 'MT-ND6', 'ABCC8',
                              'MAFA', 'PDX1', 'PCSK1', 'MAFB', 'ACTB', 'GSN'), # use this for labelling genes on plot
                #boxedLabels = TRUE,
                xlim = c(-1.5,1.5),
                ylim = c(0,300),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Custom colour over-ride',
                pCutoff = 0.05,
                FCcutoff = 0.137504,
                #pointSize = c(ifelse(volcanodat$avg_log2FC < -1 | volcanodat$avg_log2FC > 1, 6, 4)),
                pointSize = 3,
                labSize = 5,
                labFace = 'bold',
                #boxedLabels = TRUE,
                labCol = 'black',
                shape = c(20, 20, 20, 20),
                colCustom = keyvals,
                colAlpha = 2,
                legendPosition = 'right',
                legendLabSize = 15,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 1,
                typeConnectors = 'closed', 
                colConnectors = 'black',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                border = 'partial',
                borderWidth = 1.5,
                borderColour = 'black') + theme(axis.text.x = element_text(colour = "black"),
                                                axis.text.y = element_text(colour = "black"),
                                                axis.title.x = element_text(colour = "black"),
                                                axis.title.y = element_text(colour = "black"))

# Load data
volcanodat <- read.csv(r"(C:\Users\mqadir\Box\!FAHD\1. AR-DHT Project\DHT_scRNAseq_Islets\1. DGE_analysis\1. All Genes\1beta.INSLow.DHT.response.csv)",
                       header = TRUE, sep = ",", row.names = 1)
#volcanodat <- epsilon.DHT.response

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# set the base colour as 'black'
keyvals <- rep('black', nrow(volcanodat))

# set the base name/label as 'Mid'
names(keyvals) <- rep('Mid', nrow(volcanodat))

# modify keyvals for variables with fold change > 1.1
keyvals[which(volcanodat$avg_log2FC > 0.137504 & volcanodat$p_val < 0.05)] <- 'red'
names(keyvals)[which(volcanodat$avg_log2FC > 0.137504 & volcanodat$p_val < 0.05)] <- 'high'

# modify keyvals for variables with fold change < -1.1
keyvals[which(volcanodat$avg_log2FC < -0.137504 & volcanodat$p_val < 0.05)] <- 'royalblue'
names(keyvals)[which(volcanodat$avg_log2FC < -0.137504 & volcanodat$p_val < 0.05)] <- 'low'

unique(names(keyvals))

unique(keyvals)
keyvals[1:20]
options(ggrepel.max.overlaps = Inf)
EnhancedVolcano(volcanodat,
                lab = rownames(volcanodat),
                x = 'avg_log2FC',
                y = 'p_val',
                #selectLab = FALSE,
                #selectLab = rownames(volcanodat)[which(names(keyvals) %in% c('high', 'low'))],
                selectLab = c('ABCC8', 'NKX6-1', 'VEGFA', 'INS', 'PGK1', 'KIF1A', 'NKX2-2', 'IAPP', 'PDK', 'KCNK1',
                              'GCK', 'MYO6', 'GCGR', 'ENO2', 'MT-CO3', 'MT-ND1', 'MT-ATP6', 'MT-ND4', 'SOD2', 'MT-CYB', 'MT-CO1',
                              'COX4I1', 'NDUFB6', 'UQCRB', 'ATP5F1E', 'ATP6V0E2', 'NDUFS8'), # use this for labelling genes on plot
                #boxedLabels = TRUE,
                xlim = c(-1.2,1.6),
                ylim = c(0,45),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Custom colour over-ride',
                pCutoff = 0.05,
                FCcutoff = 0.137504,
                #pointSize = c(ifelse(volcanodat$avg_log2FC < -1 | volcanodat$avg_log2FC > 1, 6, 4)),
                pointSize = 3,
                labSize = 5,
                labFace = 'bold',
                #boxedLabels = TRUE,
                labCol = 'black',
                shape = c(20, 20, 20, 20),
                colCustom = keyvals,
                colAlpha = 2,
                legendPosition = 'right',
                legendLabSize = 15,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 1,
                typeConnectors = 'closed', 
                colConnectors = 'black',
                lengthConnectors = unit(0.02, 'npc'),
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                border = 'partial',
                borderWidth = 1.5,
                borderColour = 'black') + theme(axis.text.x = element_text(colour = "black"),
                                                axis.text.y = element_text(colour = "black"),
                                                axis.title.x = element_text(colour = "black"),
                                                axis.title.y = element_text(colour = "black"))
# Calculating percentages
X1 <- NULL
table(x = FetchData(pancreas.integrated, vars = c('celltype', 'sample')))
x1 <- subset(pancreas.integrated, subset = (celltype == c("beta", "alpha")) & (sex == "Male"))
table(x = FetchData(x1, vars = c('celltype', 'sex')))
x2 <- subset(pancreas.integrated, subset = (celltype != c("beta")) & (sex != "Male")) # wont run because you cant subset a vector with no cells which is what is left once all male cells are removed :)
table(x = FetchData(x2, vars = c('celltype', 'sex')))






theme_set(theme_cowplot())
beta.cells <- subset(pancreas.integrated, idents = "beta")
Idents(beta.cells) <- "treatment"
avg.beta.cells <- log1p(AverageExpression(beta.cells, verbose = FALSE)$RNA)
avg.beta.cells$gene <- rownames(avg.beta.cells)

alpha.cells <- subset(pancreas.integrated, idents = "alpha")
Idents(alpha.cells) <- "treatment"
avg.alpha.cells <- log1p(AverageExpression(alpha.cells, verbose = FALSE)$RNA)
avg.alpha.cells$gene <- rownames(avg.alpha.cells)

genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
p1 <- ggplot(avg.beta.cells, aes("ctrl", "DHT[10nM]")) + geom_point() + ggtitle("Beta Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
plot_grid(p1, p2)



# Extra
# Installation of the latest released version
install.packages('GOplot')
library(GOplot)
packageVersion("GOplot")

gene.data <- read.csv("D:/R-Projects/DHT/Data output/beta.DHT.de.data.csv")
up.go <- read.csv("D:/R-Projects/DHT/Data output/up/DHTGOup.csv")
down.go <- read.csv("D:/R-Projects/DHT/Data output/Down/DHTGOdown.csv")

head(gene.data)
head(up.go)
circ <- circle_dat(up.go, gene.data)
circ
process <- List('cellular response to decreased oxygen levels', "cellular response to hypoxia",
                'response to unfolded protein', 'canonical glycolysis',
                'glucose catabolic process to pyruvate', 'glycolytic process through glucose-6-phosphate',
                'gluconeogenesis', 'cellular response to oxidative stress',
                'protein stabilization ', 'amino acid transport ')

process
chord <- chord_dat(data = circ, genes = gene.data)
chord <- chord_dat(data = circ, process = process)
chord <- chord_dat(data = circ, genes = gene.data, process = process)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)


GOBubble(circ, labels = 1)


