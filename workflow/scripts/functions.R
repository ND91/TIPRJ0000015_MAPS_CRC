#This script contains mostly accessory functions that are used in several analyses in one form or another.

# seuratDE is a function to perform pseudobulk differential expression on all celltypes in a column listed in one of the SeuratObject columns. It essentially creates a pseudobulk dataset and wraps around DESeq2 to perform DE analyses.
#   input:
#     seuratobj: The SeuratObject of interest.
#     sampleinfo: A dataframe containing the sample metadata, which will be fed to DESeq2 as the colData.  Note that this is not the same as the cell metadata (albeit a reduced form thereof).
#     cellsampleID: Column name in the SeuratObject metadata containing the sample identifier.  
#     cellclusterID: Column name in the SeuratObject by which the cells should be aggregated. If left empty, function will perform pseudobulk DE on all cells.
#     design: The design formula to be used in the analysis. Note that the columns should be present in \"sampleinfo\".
#     contrast: The comparison of interest per the DESeq2 requirements. Either \"name\" or \"contrast\" should be provided.
#     name: The name of the result you want from the DESeq2 output. Either \"name\" or \"contrast\" should be provided.
#   output: A list of celltypes each of which includes a DESeq2 object and the results per the contrasts indicated previously ranked by p-value.

seuratDE <- function(seuratobj, cellsampleID, cellclusterID = NULL, sampleinfo, design, contrast = NULL, name = NULL){
  
  if(!class(seuratobj) != "SeuratObj") stop("\"seuratobj\" must be of the class SeuratObject")
  if(!cellsampleID %in% colnames(seuratobj@meta.data)) stop("\"cellsampleID\" cannot be found among the cell metadata column names in the provided \"seuratobj\"")
  if(!cellclusterID %in% colnames(seuratobj@meta.data)) stop("\"cellclusterID\" cannot be found among the cell metadata column names in the provided \"seuratobj\"")
  if(is.null(contrast) & is.null(name)) stop("Either \"name\" or \"contrast\" should be provided, not both.")
  
  if(!is.null(cellclusterID)){
    seuratobj_list <- Seurat::SplitObject(object = seuratobj, split.by = cellclusterID)
  } else{
    seuratobj_list <- list(seuratobj)
  }
  
  de_list <- lapply(seuratobj_list, function(seuratentry){
    
    tryCatch(expr = {
      counts_mat <- GetAssayData(seuratentry, slot = "counts")
      
      sourceID <- paste0("~0+", cellsampleID)
      
      counts_sample <- counts_mat %*% model.matrix(as.formula(sourceID), data = seuratentry@meta.data)
      colnames(counts_sample) <- gsub(cellsampleID, "", colnames(counts_sample))
      
      counts_sample <- counts_sample[which(Matrix::rowSums(counts_sample) != 0),]
      
      dds <- DESeqDataSetFromMatrix(countData = counts_sample,
                                    colData = sampleinfo[colnames(counts_sample),],
                                    design = as.formula(design))
      dds <- DESeq(dds)
      
      if(is.null(name)){
        degs <- results(dds, contrast = contrast, independentFiltering = T)
      } else{
        degs <- results(dds, name = name, independentFiltering = T)
      }
      
      degs <- degs[order(degs$pvalue),]
      degs <- degs[!is.na(degs$padj),]
      degs$gene <- rownames(degs)
      return(list(degs = degs,
                  dds = dds))
    },
    error = function(cond) {
      message("Failed with message:")
      message(cond)
      
      return(NULL)
    })
  })
  return(de_list)
}

# seuratDA is a function to perform differential abundance analysis in a relative fashion (i.e. CD4Treg relative to CD4T) for all celltypes in parent population listed in one of the SeuratObject columns. It essentially creates a frequency dataset and wraps around DESeq2 or speckle::propeller to perform DA analyses.
#   input:
#     seuratobj: The SeuratObject of interest.
#     sampleinfo: sampleinfo A dataframe containing the sample metadata, which will be fed to DESeq2 as the colData.  Note that this is not the same as the cell metadata (albeit a reduced form thereof).
#     cellsampleID: Column name in the SeuratObject metadata containing the sample identifier.  
#     cellclusterID: Column name in the SeuratObject by which the cells should be aggregated. If left empty, function will perform pseudobulk DE on all cells.
#     clusterparentID: An optional column in the SeuratObject to which the clusters belong. Must be a super set of the factors in the \"cellclusterID\". If it is left empty (default), cell abundances will be relative to the total cell counts. If this argument is not empty, the result will be a list of tables, each representing a particular parent cluster. If the unique parent cluster is equal to the child cluster, a NULL will be returned.
#     ncell_threshold: The minimum number of cells in a cluster to be considered for differential abundance analysis. Akin to the minimum number of reads in bulk RNAseq. Defaults to 20.
#     design: If method = \"DESeq2\", the design argument should be the design formula to be used in the analysis. If method = \"propeller\", the design argument should be column name in the SeuratObject by which the cells are grouped for comparative analyses. Note that the columns should be present in \"sampleinfo\".
#     contrast: The comparison of interest per the DESeq2 requirements. Either \"name\" or \"contrast\" should be provided.
#     name: The name of the result you want from the DESeq2 output. Either \"name\" or \"contrast\" should be provided.
#   output: A list of celltypes each of which includes a DESeq2 object and the results per the contrasts indicated previously ranked by p-value.

if(require(speckle) == FALSE){
  devtools::install_github("phipsonlab/speckle", ref = "09e202f") # This commit version is the last commit before speckle was made for 4.2 only.
  require(speckle)
}

seuratDA <- function(seuratobj, sampleinfo, cellsampleID, cellclusterID, clusterparentID = NULL, ncell_threshold = 20, method = c("propeller", "DESeq2", "edgeR"), design, contrast = NULL, name = NULL){
  
  method <- match.arg(method)
  
  if(!cellsampleID %in% colnames(seuratobj@meta.data)) stop("\"cellsampleID\" cannot be found among the cell metadata column names in the provided \"seuratobj\"")
  if(!cellclusterID %in% colnames(seuratobj@meta.data)) stop("\"cellclusterID\" cannot be found among the cell metadata column names in the provided \"seuratobj\"")
  if(!is.null(clusterparentID)){
    if (!clusterparentID %in% colnames(seuratobj@meta.data)) stop("\"clusterparentID\" cannot be found among the cell metadata column names in the provided \"seuratobj\"")
  } 
  if(is.null(contrast) & is.null(name) & method == "DESeq2") stop("For DESeq2, either \"name\" or \"contrast\" should be provided, not both.")
  
  if(!is.null(clusterparentID)){
    #Remove parent cell clusters whose cell counts lie below the threshold.
    seuratobj <- seuratobj[,which(seuratobj@meta.data[,clusterparentID] %in% names(which(table(seuratobj@meta.data[,clusterparentID]) >= ncell_threshold)))]
    
    seuratobj_list <- Seurat::SplitObject(object = seuratobj, split.by = clusterparentID)
  } else{
    seuratobj_list <- list(seuratobj)
  }
  
  da_list <- lapply(seuratobj_list, function(seuratentry){
    
    proper_clusters <- which(seuratentry@meta.data[,cellclusterID] %in% names(which(table(seuratentry@meta.data[,cellclusterID]) >= ncell_threshold)))
    
    if(length(proper_clusters) == 0){
      return(NULL)
    }
    
    #Remove clusters whose cell counts lie below the threshold.
    seuratentry <- seuratentry[, proper_clusters]
    
    if(length(unique(seuratentry@meta.data[,cellclusterID])) == 1){
      return(NULL)
    } else{
      tryCatch(expr = {
        if(method == "propeller"){
          if(!design %in% colnames(seuratobj@meta.data)) stop("\"design\" cannot be found among the cell metadata column names in the provided \"seuratobj\"")
          
          da <- speckle::propeller(clusters = seuratentry@meta.data[,cellclusterID], 
                                   sample = seuratentry@meta.data[,cellsampleID],
                                   group = seuratentry@meta.data[,design])
          return(da)
        } else if(method == "DESeq2"){
          counts_df <- data.frame(table(seuratentry@meta.data[,cellclusterID], seuratentry@meta.data[,cellsampleID])) %>%
            dplyr::rename(cellclusterID = Var1,
                          cellsampleID = Var2) %>%
            tidyr::pivot_wider(id_cols = cellclusterID, 
                               names_from = cellsampleID, 
                               values_from = Freq)
          counts_df <- data.frame(counts_df[,-1], row.names = counts_df$cellclusterID)
          
          dds <- DESeqDataSetFromMatrix(countData = counts_df,
                                        colData = sampleinfo[colnames(counts_df),],
                                        design = as.formula(design))
          
          dds <- DESeq(dds)
          
          if(is.null(name)){
            da <- results(dds, contrast = contrast, independentFiltering = T)
          } else{
            da <- results(dds, name = name, independentFiltering = T)
          }
          da <- da[order(da$pvalue),]
          return(list(da = da,
                      dds = dds))
        } else if(method == "edgeR"){
          cat("This has not yet been implemented")
        }
      },
      error = function(cond) {
        
        message("Failed with message:")
        message(cond)
        
        return(NULL)
      })
    }
  })
  return(da_list)
}
