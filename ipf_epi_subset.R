#generate ipf_epithelial subset for integration

#####ipf data
tild001.data <- Read10X_h5("~/Desktop/sa_ild/tild001lf_filtered_feature_bc_matrix.h5")
tild006a.data <- Read10X_h5("~/Desktop/sa_ild/tild006lf_filtered_feature_bc_matrix.h5")
tild006b.data <- Read10X_h5("~/Desktop/sa_ild/tild006mf_filtered_feature_bc_matrix.h5")
tild010.data <- Read10X_h5("~/Desktop/sa_ild/tild010mf_filtered_feature_bc_matrix.h5")
tild015a.data <- Read10X_h5("~/Desktop/sa_ild/tild015lf_filtered_feature_bc_matrix.h5")
tild015b.data <- Read10X_h5("~/Desktop/sa_ild/tild015mf_filtered_feature_bc_matrix.h5")
tild028a.data <- Read10X_h5("~/Desktop/sa_ild/tild028lf_filtered_feature_bc_matrix.h5")
tild028b.data <- Read10X_h5("~/Desktop/sa_ild/tild028mf_filtered_feature_bc_matrix.h5")
vuild053.data <- Read10X_h5("~/Desktop/sa_ild/vuild053mf_filtered_feature_bc_matrix.h5")
vuild059a.data <- Read10X_h5("~/Desktop/sa_ild/vuild059a_filtered_feature_bc_matrix.h5")
vuild059b.data <- Read10X_h5("~/Desktop/sa_ild/vuild059b_filtered_feature_bc_matrix.h5")
vuild060a.data <- Read10X_h5("~/Desktop/sa_ild/vuild060a_filtered_feature_bc_matrix.h5")
vuild060b.data <- Read10X_h5("~/Desktop/sa_ild/vuild060b_filtered_feature_bc_matrix.h5")
vuild064.data <- Read10X_h5("~/Desktop/sa_ild/vuild064_filtered_feature_bc_matrix.h5")
vuild087a.data <- Read10X_h5("~/Desktop/sa_ild/vuild087a_filtered_feature_bc_matrix.h5")
vuild087b.data <- Read10X_h5("~/Desktop/sa_ild/vuild087b_filtered_feature_bc_matrix.h5")

tild001 <- CreateSeuratObject(counts = tild001.data, min.features = 750, project = "T_ILD_001")
tild006a <- CreateSeuratObject(counts = tild006a.data, min.features = 750, project = "T_ILD_006a")
tild006b <- CreateSeuratObject(counts = tild006b.data, min.features = 750, project = "T_ILD_006b")
tild010 <- CreateSeuratObject(counts = tild010.data, min.features = 750, project = "T_ILD_010")
tild015a <- CreateSeuratObject(counts = tild015a.data, min.features = 750, project = "T_ILD_015a")
tild015b <- CreateSeuratObject(counts = tild015b.data, min.features = 750, project = "T_ILD_015b")
tild028a <- CreateSeuratObject(counts = tild028a.data, min.features = 750, project = "T_ILD_028a")
tild028b <- CreateSeuratObject(counts = tild028b.data, min.features = 750, project = "T_ILD_028b")
tild030a <- CreateSeuratObject(counts = tild030a.data, min.features = 750, project = "T_ILD_030a")
tild030b <- CreateSeuratObject(counts = tild030b.data, min.features = 750, project = "T_ILD_030b")
vuild053 <- CreateSeuratObject(counts = vuild053.data, min.features = 750, project = "VU_ILD_053")
vuild059a <- CreateSeuratObject(counts = vuild059a.data, min.features = 750, project = "VU_ILD_059a")
vuild059b <- CreateSeuratObject(counts = vuild059b.data, min.features = 750, project = "VU_ILD_059b")
vuild060a <- CreateSeuratObject(counts = vuild060a.data, min.features = 750, project = "VU_ILD_060a")
vuild060b <- CreateSeuratObject(counts = vuild060b.data, min.features = 750, project = "VU_ILD_060b")
vuild064 <- CreateSeuratObject(counts = vuild064.data, min.features = 750, project = "VU_ILD_064")
vuild087a <- CreateSeuratObject(counts = vuild087a.data, min.features = 750, project = "VU_ILD_087a")
vuild087b <- CreateSeuratObject(counts = vuild087b.data, min.features = 750, project = "VU_ILD_087b")

adata <- merge (tild001, y=c(tild006a, tild006b, tild010, tild015a, tild015b, tild028a, tild028b, vuild053, vuild059a, vuild059b, vuild060a, vuild060b, vuild064, vuild087a, vuild087b))

Idents(adata) <- 'orig.ident'
Idents(adata, cells = WhichCells(adata, idents = c('T_ILD_001', 'T_ILD_006a', 'T_ILD_006b', 'T_ILD_010', 'T_ILD_015a', 'T_ILD_015b', 'T_ILD_028a', 'T_ILD_028b', 'VU_ILD_053', 'VU_ILD_059a', 'VU_ILD_059b', 'VU_ILD_060a', 'VU_ILD_060b', 'VU_ILD_064', 'VU_ILD_087a', 'VU_ILD_087b'))) <- "IPF"
adata$diagnosis <- Idents(adata)

adata <- PercentageFeatureSet(adata, pattern = "^MT-", col.name = "percent.mt")
adata <- subset(adata, subset = percent.mt < 15 & percent.mt >0.5 & nFeature_RNA >750)
adata <- SCTransform(adata, vars.to.regress = 'percent.mt', conserve.memory=TRUE)
adata <- RunPCA(adata, verbose = F)
adata <- RunUMAP(adata, dims = 1:35)
adata <- FindNeighbors(adata, dims = 1:35)
adata <- FindClusters(adata, resolution = 1.0)
DimPlot(adata, label = T, cols = 'polychrome')

Idents(adata) <- 'seurat_clusters'

DotPlot(adata, features = c("PTPRC", "EPCAM", "PECAM1", "CD3E", "COL1A1", "CD19", 'ACTA2', 'WT1'))

epi <- subset(adata, idents = c(2,5,13,15,18,21:25,32,35,37))




epi <- SCTransform(epi, vars.to.regress = 'percent.mt', conserve.memory = TRUE)
epi <- RunPCA(epi, verbose = F)
epi <- RunUMAP(epi, dims = 1:25)
epi <- FindNeighbors(epi, dims = 1:25)
epi <- FindClusters(epi, resolution = 1.0)
DimPlot(epi, label = T, cols = 'polychrome')
DimPlot(epi, split.by = 'diagnosis', cols = 'polychrome')

DotPlot(epi, features = c("PTPRC", "PECAM1", "COL1A1",  "EPCAM", "SFTPC", "SCGB3A2", "KRT17", "AGER", "FOXJ1", "MUC5B", "MKI67", "CALCA", "FOXI1"))