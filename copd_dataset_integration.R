library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
options(future.globals.maxSize = 40000 * 1024^2)
library(RColorBrewer)
library(viridis)
library(hrbrthemes)

#control data
hd65.data <- Read10X_h5("~/Desktop/sa_control/vuhd065_filtered_feature_bc_matrix.h5")
hd66.data <- Read10X_h5("~/Desktop/sa_control/vuhd066_filtered_feature_bc_matrix.h5")
hd67.data <- Read10X_h5("~/Desktop/sa_control/vuhd067_filtered_feature_bc_matrix.h5")
hd68.data <- Read10X_h5("~/Desktop/sa_control/vuhd068_filtered_feature_bc_matrix.h5")
hd69.data <- Read10X_h5("~/Desktop/sa_control/vuhd069_filtered_feature_bc_matrix.h5")
hd70.data <- Read10X_h5("~/Desktop/sa_control/vuhd070_filtered_feature_bc_matrix.h5")
hd71.data <- Read10X_h5("~/Desktop/sa_control/vuhd071_filtered_feature_bc_matrix.h5")
thd001.data <- Read10X_h5("~/Desktop/sa_control/thd001_filtered_feature_bc_matrix.h5")
thd002.data <- Read10X_h5("~/Desktop/sa_control/thd002_filtered_feature_bc_matrix.h5")
thd005a.data <- Read10X_h5("~/Desktop/sa_control/thd005a_filtered_feature_bc_matrix.h5")
thd005b.data <- Read10X_h5("~/Desktop/sa_control/thd005b_filtered_feature_bc_matrix.h5")
thd005c.data <- Read10X_h5("~/Desktop/sa_control/thd005c_filtered_feature_bc_matrix.h5")

hd65 <- CreateSeuratObject(counts = hd65.data, min.features = 750, project = "VU_HD_65")
hd66 <- CreateSeuratObject(counts = hd66.data, min.features = 750, project = "VU_HD_66")
hd67 <- CreateSeuratObject(counts = hd67.data, min.features = 750, project = "VU_HD_67")
hd68 <- CreateSeuratObject(counts = hd68.data, min.features = 750, project = "VU_HD_68")
hd69 <- CreateSeuratObject(counts = hd69.data, min.features = 750, project = "VU_HD_69")
hd70 <- CreateSeuratObject(counts = hd70.data, min.features = 750, project = "VU_HD_70")
hd71 <- CreateSeuratObject(counts = hd71.data, min.features = 750, project = "VU_HD_71")
thd001 <- CreateSeuratObject(counts = thd001.data, min.features = 750, project = "T_HD_001")
thd002 <- CreateSeuratObject(counts = thd002.data, min.features = 750, project = "T_HD_002")
thd005a <- CreateSeuratObject(counts = thd005a.data, min.features = 750, project = "T_HD_005a")
thd005b <- CreateSeuratObject(counts = thd005b.data, min.features = 750, project = "T_HD_005b")
thd005c <- CreateSeuratObject(counts = thd005c.data, min.features = 750, project = "T_HD_005c")

#copd data
tcopd13.data <- Read10X_h5("~/Desktop/new_copd/t_copd_13.h5")
tcopd50.data <- Read10X_h5("~/Desktop/new_copd/t_copd_50.h5")
tcopd63.data <- Read10X_h5("~/Desktop/new_copd/t_copd_63.h5")
vu_copd_29.data <- Read10X_h5("~/Desktop/new_copd/copd29v3.h5")
vu_copd_30.data <- Read10X_h5("~/Desktop/new_copd/vu_copd_30.h5")
vu_copd_34.data <- Read10X_h5("~/Desktop/new_copd/vu_copd_34.h5")
vu_copd_38.data <- Read10X_h5("~/Desktop/new_copd/vu_copd_38.h5")
vu_copd_40.data <- Read10X_h5("~/Desktop/new_copd/vu_copd_40.h5")
vu_copd_41.data <- Read10X_h5("~/Desktop/new_copd/vu_copd_41.h5")
vu_copd_42.data <- Read10X_h5("~/Desktop/new_copd/vu_copd_42.h5")
vu_copd_44.data <- Read10X_h5("~/Desktop/new_copd/vu_copd_44.h5")
vu_copd_45.data <- Read10X_h5("~/Desktop/new_copd/vu_copd_45.h5")

tcopd13 <- CreateSeuratObject(counts = tcopd13.data, min.features = 750, project = "T_COPD_13")
tcopd50 <- CreateSeuratObject(counts = tcopd50.data, min.features = 750, project = "T_COPD_50")
tcopd63 <- CreateSeuratObject(counts = tcopd63.data, min.features = 750, project = "T_COPD_63")
vucopd29 <- CreateSeuratObject(counts = vu_copd_29.data, min.features = 750, project = "VU_COPD_29")
vucopd30 <- CreateSeuratObject(counts = vu_copd_30.data, min.features = 750, project = "VU_COPD_30")
vucopd34 <- CreateSeuratObject(counts = vu_copd_34.data, min.features = 750, project = "VU_COPD_34")
vucopd38 <- CreateSeuratObject(counts = vu_copd_38.data, min.features = 750, project = "VU_COPD_38")
vucopd40 <- CreateSeuratObject(counts = vu_copd_40.data, min.features = 750, project = "VU_COPD_40")
vucopd41 <- CreateSeuratObject(counts = vu_copd_41.data, min.features = 750, project = "VU_COPD_41")
vucopd42 <- CreateSeuratObject(counts = vu_copd_42.data, min.features = 750, project = "VU_COPD_42")
vucopd44 <- CreateSeuratObject(counts = vu_copd_44.data, min.features = 750, project = "VU_COPD_44")
vucopd45 <- CreateSeuratObject(counts = vu_copd_45.data, min.features = 750, project = "VU_COPD_45")

#merge copd and control libraries
copd <- merge(tcopd13, y=c(tcopd50, tcopd63, vucopd29, vucopd30, vucopd34, vucopd38, vucopd41, vucopd42, vucopd44, vucopd45, hd65, hd66, hd67, hd68, hd69, hd70, hd71, thd001, thd002, thd005a, thd005b, thd005c))

#add cellranger version metadata label
Idents(copd) <- 'orig.ident'
Idents(copd, cells = WhichCells(copd, idents = c("T_HD_001", "T_HD_002", "T_HD_005a", "T_HD_005b", "T_HD_005c", "VU_HD_65", "VU_HD_66", "VU_HD_67", "VU_HD_68", "VU_HD_69", "VU_HD_70", "VU_HD_71", "VU_COPD_30", "VU_COPD_38", "VU_COPD_41", "VU_COPD_42", "VU_COPD_44", "VU_COPD_45", 'T_COPD_13', 'T_COPD_50', 'T_COPD_63'))) <- "v5"
Idents(copd, cells = WhichCells(copd, idents = c('VU_COPD_29', 'VU_COPD_34'))) <- "v3"
copd$cellranger <- Idents(copd)

#remove low-quality cells
copd <- PercentageFeatureSet(copd, pattern = "^MT-", col.name = "percent.mt")
copd <- subset(copd, subset = percent.mt < 15 & percent.mt >0.5 & nFeature_RNA >750)

#split for SCT-integration based on cellranger version
copd.list <- SplitObject(copd, split.by = "cellranger")

for (i in 1:length(copd.list)) {
  copd.list[[i]] <- NormalizeData(copd.list[[i]], verbose = FALSE)
  copd.list[[i]] <- FindVariableFeatures(copd.list[[i]], selection.method = "vst", 
                                         nfeatures = 2000, verbose = FALSE)
}

#perform SCT on individual objects
for (i in 1:length(copd.list)) {
  copd.list[[i]] <- SCTransform(copd.list[[i]], verbose = FALSE)
}

#perform integration
copd.features <- SelectIntegrationFeatures(object.list = copd.list, nfeatures = 3000)
copd.list <- PrepSCTIntegration(object.list = copd.list, anchor.features = copd.features, 
                                verbose = FALSE)
copd.anchors <- FindIntegrationAnchors(object.list = copd.list, normalization.method = "SCT", anchor.features = copd.features, verbose = TRUE)
copd.integrated <- IntegrateData(anchorset = copd.anchors, normalization.method = "SCT", verbose = TRUE)

#Cluster and embed
copd.integrated <- RunPCA(copd.integrated, verbose = F)
copd.integrated <- RunUMAP(copd.integrated, dims = 1:35)
copd.integrated <- FindNeighbors(copd.integrated, dims = 1:35)
copd.integrated <- FindClusters(copd.integrated, resolution = 1.0)
DimPlot(copd.integrated, label=T)
DimPlot(copd.integrated, group.by='orig.ident')
DimPlot(copd.integrated, group.by='cellranger')


DefaultAssay(copd.integrated) <- 'SCT'
copd_markers <- c("PTPRC", "PECAM1", "COL1A1", 'ACTA2', "WT1", "MKI67", "EPCAM", "TP63", 'KRT5', "SCGB1A1", "SCGB3A1", "SCGB3A2", "MUC5B", "MUC5AC", "FOXJ1", "TP73", "ABCA3", "AGER", 'CLDN4', "KRT17", 'CDKN1A', 'SFN', 'CDKN2A', 'WWTR1')
DotPlot(copd.integrated, features = copd_markers)+ theme(axis.text.x = element_text(angle = 45, hjust=1))

FeaturePlot(copd.integrated, features = c('SFTPC', 'CD68', 'CD3E', 'PECAM1'))
saveRDS(copd.integrated, file = '~/Desktop/new_copd/copd_integrated_032321.rds')

#subset into epithelial, stromal and immune sub-objects
copd_epi <- subset(copd.integrated, idents = c(0,3,5,7,11,14,15,32,33,35,37))
copd_stromal <- subset(copd.integrated, idents = c(1,17,20,25,29,19,26,34))
copd_immune <- subset(copd.integrated, idents = c(2,4,6,8:10,12,13,16,18,21:23,27,28,30,31,38:40))

saveRDS(copd_epi, file = '~/Desktop/new_copd/copd_epi_032321.rds')
saveRDS(copd_stromal, file = '~/Desktop/new_copd/copd_stromal_032321.rds')
saveRDS(copd_immune, file = '~/Desktop/new_copd/copd_immune_032321.rds')



#integrate with ipf_epi for joint annotation
ipf_epi <- readRDS("~/ipf_epi.rds")

Idents(ipf_epi) <- 'Diagnosis'
Idents(ipf_epi, cells = WhichCells(ipf_epi, idents = c('IPF'))) <- "v5"
ipf_epi$cellranger <- Idents(ipf_epi)


epi <- merge(copd_epi, y=c(ipf_epi))
DefaultAssay(epi) <- 'RNA'
epi[['SCT']] <- NULL
epi[['integrated']] <- NULL

#split for SCT-integration
epi.list <- SplitObject(epi, split.by = "cellranger")

for (i in 1:length(epi.list)) {
  epi.list[[i]] <- NormalizeData(epi.list[[i]], verbose = FALSE)
  epi.list[[i]] <- FindVariableFeatures(epi.list[[i]], selection.method = "vst", 
                                        nfeatures = 2000, verbose = FALSE)
}

#perform SCT on individual objects
for (i in 1:length(epi.list)) {
  epi.list[[i]] <- SCTransform(epi.list[[i]], verbose = FALSE)
}

#perform integration
epi.features <- SelectIntegrationFeatures(object.list = epi.list, nfeatures = 3000)
epi.list <- PrepSCTIntegration(object.list = epi.list, anchor.features = epi.features, 
                               verbose = FALSE)
epi.anchors <- FindIntegrationAnchors(object.list = epi.list, normalization.method = "SCT", anchor.features = epi.features, verbose = TRUE)
epi.integrated <- IntegrateData(anchorset = epi.anchors, normalization.method = "SCT", verbose = TRUE)

#Cluster and embed
epi.integrated <- RunPCA(epi.integrated, verbose = F)
DefaultAssay(epi.integrated) <- 'integrated'
epi.integrated <- RunUMAP(epi.integrated, dims = 1:25)
epi.integrated <- FindNeighbors(epi.integrated, dims = 1:25)
epi.integrated <- FindClusters(epi.integrated, resolution = 1.0)
DimPlot(epi.integrated, label=T)
DimPlot(epi.integrated, group.by='orig.ident')

DefaultAssay(epi.integrated) <- 'SCT'
epi_markers <- c("PTPRC", "PECAM1", "COL1A1", "MKI67", "TP63", 'KRT5', "SCGB1A1", "SCGB3A1", "SCGB3A2", "MUC5B", "MUC5AC", "FOXJ1", "TP73", "ABCA3", "AGER", 'CLDN4', "KRT17", 'CDKN1A', 'SFN', 'CDKN2A', 'WWTR1')
DotPlot(epi.integrated, features = epi_markers)+ theme(axis.text.x = element_text(angle = 45, hjust=1))

#remove doublets
DefaultAssay(epi.integrated) <- 'integrated'
epi.integrated <- subset(epi.integrated, idents = c(0:19,22,23,26,32))
epi.integrated <- RunPCA(epi.integrated, verbose = F)
DefaultAssay(epi.integrated) <- 'integrated'
epi.integrated <- RunUMAP(epi.integrated, dims = 1:18)
epi.integrated <- FindNeighbors(epi.integrated, dims = 1:18)
epi.integrated <- FindClusters(epi.integrated, resolution = 1.5)
DimPlot(epi.integrated, label=T)
DimPlot(epi.integrated, group.by='orig.ident')

DefaultAssay(epi.integrated) <- 'SCT'
epi_markers <- c("PTPRC", "PECAM1", "COL1A1", "MKI67", "TP63", 'KRT5', "SCGB1A1", "SCGB3A1", "SCGB3A2", "MUC5B", "MUC5AC", "FOXJ1", "TP73", "ABCA3", "AGER", 'HOPX', 'NKX2-1', 'CLDN4', "KRT17", 'KRT8', 'KRT19', 'CDKN1A', 'SFN', 'CDKN2A', 'WWTR1', 'CALCA', 'FOXI1', 'COL4A4')
DotPlot(epi.integrated, features = epi_markers)+ theme(axis.text.x = element_text(angle = 45, hjust=1))

Idents(epi.integrated) <- 'orig.ident'
Idents(epi.integrated, cells = WhichCells(epi.integrated, idents = c("T_HD_001", "T_HD_002", "T_HD_005a", "T_HD_005b", "T_HD_005c", "VU_HD_65", "VU_HD_66", "VU_HD_67", "VU_HD_68", "VU_HD_69", "VU_HD_70", "VU_HD_71"))) <- "Control"
Idents(epi.integrated, cells = WhichCells(epi.integrated, idents = c('T_COPD_13', 'T_COPD_50', 'T_COPD_63', 'VU_COPD_29', 'VU_COPD_30', 'VU_COPD_34', 'VU_COPD_38', 'VU_COPD_41', 'VU_COPD_42', 'VU_COPD_44', 'VU_COPD_45'))) <- "COPD"
Idents(epi.integrated, cells = WhichCells(epi.integrated, idents = c('T_ILD_001', 'T_ILD_006a', 'T_ILD_006b', 'T_ILD_010', 'T_ILD_015a', 'T_ILD_015b', 'T_ILD_028a', 'T_ILD_028b', 'VU_ILD_053', 'VU_ILD_059a', 'VU_ILD_059b', 'VU_ILD_060a', 'VU_ILD_060b', 'VU_ILD_064', 'VU_ILD_087a', 'VU_ILD_087b'))) <- "IPF"
epi.integrated$diagnosis <- Idents(epi.integrated)


#remove doublets
DefaultAssay(epi.integrated) <- 'integrated'
epi.integrated <- subset(epi.integrated, idents = c(0:18, 20:31))
epi.integrated <- RunPCA(epi.integrated, verbose = F)
DefaultAssay(epi.integrated) <- 'integrated'
epi.integrated <- RunUMAP(epi.integrated, dims = 1:18)
epi.integrated <- FindNeighbors(epi.integrated, dims = 1:18)
epi.integrated <- FindClusters(epi.integrated, resolution = 1.5)
DimPlot(epi.integrated, label=T)
DimPlot(epi.integrated, group.by='orig.ident')

DefaultAssay(epi.integrated) <- 'SCT'
epi_markers <- c("PTPRC", "PECAM1", "COL1A1", "MKI67", "TP63", 'KRT5', "SCGB1A1", "SCGB3A1", "SCGB3A2", "MUC5B", "MUC5AC", "FOXJ1", "TP73", "LAMP3", "ABCA3", "CEACAM6", "AGER", 'HOPX', 'NKX2-1', 'CLDN18', 'CLDN4', "KRT17", 'KRT8', 'KRT19', 'CDKN1A', 'SFN', 'CDKN2A', 'WWTR1', 'COL4A4', 'FBLN5', 'ELN')
DotPlot(epi.integrated, features = epi_markers)+ theme(axis.text.x = element_text(angle = 45, hjust=1))




#Reintegrate after doublet removal
DefaultAssay(epi) <- 'RNA'
epi[['SCT']] <- NULL
epi[['integrated']] <- NULL

#split for SCT-integration
epi.list <- SplitObject(epi, split.by = "cellranger")

for (i in 1:length(epi.list)) {
  epi.list[[i]] <- NormalizeData(epi.list[[i]], verbose = FALSE)
  epi.list[[i]] <- FindVariableFeatures(epi.list[[i]], selection.method = "vst", 
                                        nfeatures = 2000, verbose = FALSE)
}

#perform SCT on individual objects
for (i in 1:length(epi.list)) {
  epi.list[[i]] <- SCTransform(epi.list[[i]], verbose = FALSE)
}

#perform integration
epi.features <- SelectIntegrationFeatures(object.list = epi.list, nfeatures = 3000)
epi.list <- PrepSCTIntegration(object.list = epi.list, anchor.features = epi.features, 
                               verbose = FALSE)
epi.anchors <- FindIntegrationAnchors(object.list = epi.list, normalization.method = "SCT", anchor.features = epi.features, verbose = TRUE)
epi.integrated <- IntegrateData(anchorset = epi.anchors, normalization.method = "SCT", verbose = TRUE)

#Cluster and embed
epi.integrated <- RunPCA(epi.integrated, verbose = F)
DefaultAssay(epi.integrated) <- 'integrated'
epi.integrated <- RunUMAP(epi.integrated, dims = 1:25)
epi.integrated <- FindNeighbors(epi.integrated, dims = 1:25)
epi.integrated <- FindClusters(epi.integrated, resolution = 1.0)
DimPlot(epi.integrated, label=T)
DimPlot(epi.integrated, group.by='orig.ident')

DefaultAssay(epi.integrated) <- 'SCT'
epi_markers <- c("PTPRC", "PECAM1", "COL1A1", "MKI67", "TP63", 'KRT5', "SCGB1A1", "SCGB3A1", "SCGB3A2", "MUC5B", "MUC5AC", "FOXJ1", "TP73", "LAMP3", "ABCA3", "CEACAM6", "AGER", 'HOPX', 'NKX2-1', 'CLDN18', 'CLDN4', "KRT17", 'KRT8', 'KRT19', 'CDKN1A', 'SFN', 'CDKN2A', 'WWTR1', 'COL4A4', 'FBLN5', 'ELN', 'FOXI1', 'CALCA', 'AXIN2', 'NTM', 'NCKAP5',  'RTKN2', 'SCEL', 'LAMA3','COL4A1', 'COL4A2', 'COL4A3', 'COL4A5', 'COL4A6')
DotPlot(epi.integrated, features = epi_markers)+ theme(axis.text.x = element_text(angle = 45, hjust=1))


Idents(epi.integrated) <- 'seurat_clusters'
Idents(epi.integrated, cells = WhichCells(epi.integrated, idents = c(18))) <- "Transitional"
Idents(epi.integrated, cells = WhichCells(epi.integrated, idents = c(1))) <- "Secretory - SCGB3A2+"
Idents(epi.integrated, cells = WhichCells(epi.integrated, idents = c(14))) <- "Secretory - SCGB3A1+/MUC5B+"
Idents(epi.integrated, cells = WhichCells(epi.integrated, idents = c(2))) <- "Secretory - SCGB1A1+/MUC5B+"
Idents(epi.integrated, cells = WhichCells(epi.integrated, idents = c(5))) <- "Secretory - SCGB1A1+/SCGB3A2+"
Idents(epi.integrated, cells = WhichCells(epi.integrated, idents = c(28))) <- "Proliferating AT2"
Idents(epi.integrated, cells = WhichCells(epi.integrated, idents = c(22))) <- "KRT5-/KRT17+"
Idents(epi.integrated, cells = WhichCells(epi.integrated, idents = c(0,4,6,7,10,12,19,24,26,27))) <- "Ciliated"
Idents(epi.integrated, cells = WhichCells(epi.integrated, idents = c(3))) <- "Basal"
Idents(epi.integrated, cells = WhichCells(epi.integrated, idents = c(8,9,11,15,17,20,21,23))) <- "AT2"
Idents(epi.integrated, cells = WhichCells(epi.integrated, idents = c(13,16,25))) <- "AT1"
epi.integrated$celltype_final <- Idents(epi.integrated)


meta.data <- read.csv("~/Desktop/new_copd/metadata.csv", header = T)

epi.integrated@meta.data$subject <- plyr::mapvalues(x = epi.integrated@meta.data$orig.ident,
                                                    from = meta.data$orig.ident,
                                                    to = as.character(meta.data$subject))

epi.integrated@meta.data$site <- plyr::mapvalues(x = epi.integrated@meta.data$orig.ident,
                                                 from = meta.data$orig.ident,
                                                 to = as.character(meta.data$site))

saveRDS(epi.integrated, file = '~/Desktop/new_copd/epi_integrated_032521.rds')


#subset copd and control
Idents(epi.integrated) <- 'diagnosis'
copd_control_epi <-subset(epi.integrated, idents = c('COPD', 'Control'))
saveRDS(copd_control_epi, file = '~/Desktop/new_copd/copd_control_epi.rds')


VlnPlot(copd_control_epi, features = c('PIGR'), split.by = 'diagnosis', pt.size = 0.1)
