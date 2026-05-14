library(Seurat)
library(SeuratWrappers)

load.and.create.seurat <- function(mtx_path, genes_path, cells_path, metadata_path) {
  curr.data <-ReadMtx(
    mtx=mtx_path,
    features=genes_path,
    cells=cells_path,
    feature.column = 1,
    feature.sep = "\n",
    cell.sep = ",",
    skip.cell = 1
  )
  
  metadata <- read.csv(metadata_path)
  
  cells <- read.csv(cells_path)
  row.names(cells) <- cells$cell_name
  
  seurat <- CreateSeuratObject(
    counts = curr.data,
    meta.data = metadata
  )
  seurat <- AddMetaData(seurat, cells)
  Idents(seurat) <- seurat$cell_type
  
  return(seurat)
}

seurat.pipeline <- function(pbmc, is_tpm = F) {
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern="^MT-")
  pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  
  all.genes = rownames(pbmc)
  
  if(is_tpm) {
    pbmc <- SetAssayData(pbmc, slot="data", new.data = log1p(pbmc@assays$RNA$counts))
  } else{
    pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  pbmc <- ScaleData(pbmc, features = all.genes)
  pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
  pbmc <- FindNeighbors(pbmc)
  pbmc <- FindClusters(pbmc)
  pbmc <- RunUMAP(pbmc, dims=1:20)
  
  return (pbmc)
}

###

