#' Read CNR from .cnr files
#'
#' Read all .cnr files from a dir and output a CNR object
#'
#'
#' @param CNR A dir which contain cnr files
#' @param chr.use You can specify chromosome to use
#'
#' @return Returns a CNR object from those files
#'
#' @export
#'
#' @rdname ReadCNR
#' @export ReadCNR
#'
ReadCNR<-function(dir,chr.use=paste0('chr',c(seq(1,22),'X'))){
  CNR_list <- list.files(path = dir, pattern = '.cnr$')
  if(length(CNR_list)==0){stop(paste0('No cnr file found in ',dir))}
  CNR_file_list <- lapply(CNR_list, function(x){fread(file = paste0(dir,'/',x))})
  metadata <- CNR_file_list[[1]][,c('chromosome','start','end','gene')]
  if(length(CNR_list)>1){
    read_error <- sapply(lapply(CNR_file_list,function(x){
      subset(x,select = c('chromosome','start','end','gene'))}),
      FUN = identical, metadata)
    if(!all(read_error)){
      stop(paste(CNR_list[which(!read_error)],'index is not equal to',CNR_list[1],sep=' '))
    }
  }
  CNR <- list()
  CNR$metadata <- metadata
  CNR$log2 <- as.data.frame(lapply(CNR_file_list,function(x){subset(x,select = 'log2')}))
  CNR$depth <- as.data.frame(lapply(CNR_file_list,function(x){subset(x,select = 'depth')}))
  CNR$weight <- as.data.frame(lapply(CNR_file_list,function(x){subset(x,select = 'weight')}))
  for (i in c('log2','depth','weight')){names(CNR[[i]])<-CNR_list}
  class(CNR)<-'CNR'
  return(CNR)
}

#' subset CNR object
#'
#' Remove some cells you want to remove
#'
#'
#' @param CNR_obj A CNR object
#' @param cells Cells you want to keep
#'
#' @return Returns a CNR object after subset
#'
#' @export
#'
#' @rdname SubsetCNR
#' @export SubsetCNR
#'
SubsetCNR<-function(CNR_obj, cells){
  CNR_cells <- colnames(CNR_obj$log2)
  if(length(setdiff(cells,CNR_cells))>0){stop('cells contain extra cells')}
  kept <- match(cells,CNR_cells)
  CNR_obj$log2<-CNR_obj$log2[,kept]
  CNR_obj$depth<-CNR_obj$depth[,kept]
  CNR_obj$weight<-CNR_obj$weight[,kept]
  if("segment" %in% names(CNR_obj)){CNR_obj$segment<-CNR_obj$segment[,kept]}
  if("CNV_cluster" %in% names(CNR_obj)){CNR_obj$CNV_cluster<-CNR_obj$CNV_cluster[kept]}
  if("CNV_UMAP" %in% names(CNR_obj)){CNR_obj$CNV_UMAP<-CNR_obj$CNV_UMAP[kept,]}
  if("ploidy_segment" %in% names(CNR_obj)){CNR_obj$ploidy_segment<-CNR_obj$ploidy_segment[kept,]}
  return(CNR_obj)
}

#' Filter CNR object
#'
#' Remove some low quality cell in CNR object which major reads covered in a tiny region
#'
#'
#' @param CNR_obj A CNR object
#' @param threshold threshold of low quality cells
#' @param topn topn sites' depth to compare with median depth
#'
#' @return Returns a CNR object after filter
#'
#' @export
#'
#' @rdname FilterCNR
#' @export FilterCNR
#'
FilterCNR<-function(CNR_obj, threshold=2.5,topn=10){
  regions <- dim(CNR_obj$metadata)[1]
  depth.ratio <- unlist(apply(CNR_obj$depth,2,function(x){
    log10(mean(sort(x,decreasing = T)[1:topn])/median(x))
  }))
  kept <- which(depth.ratio<threshold)
  for(slot in setdiff(names(CNR_obj),'metadata')){
    CNR_obj[[slot]]<-CNR_obj[[slot]][,kept]
  }
  return(CNR_obj)
}

#' get CNV assay slot
#'
#' get CNR, CNV, ploidyCNV and meta slot from Seurat or CNR object
#'
#'
#' @param object A CNR or Seurat object
#' @param slot CNR, CNV, ploidyCNV or meta slot
#'
#' @return Returns a matrix or data.frame of slot
#'
#' @export
#'
#' @rdname getCNVslot
#' @export getCNVslot
#'
getCNVslot<-function(object, slot){
  if(class(object)=='Seurat'){
    if(slot=='CNR'){return(as.matrix(object@assays$inferCNV@counts))}
    else if(slot=='CNV'){return(as.matrix(object@assays$inferCNV@data))}
    else if(slot=='ploidyCNV'){return(as.matrix(object@assays$inferCNV@scale.data))}
    else if(slot=='meta'){return(object@assays$inferCNV@meta.features)}
    else{stop('slot parameter should be CNR, CNV, ploidyCNV or meta')}
  }
  else if(class(object)=='CNR'){
    if(slot=='CNR'){return(as.matrix(object$log2))}
    else if(slot=='CNV'){return(as.matrix(object$segment))}
    else if(slot=='ploidyCNV'){return(as.matrix(object$ploidy_segment))}
    else if(slot=='meta'){return(object$segment_info)}
    else{stop('slot parameter should be CNR, CNV, ploidyCNV or meta')}
  }
  else{stop('object should be Seurat or CNR')}
}

InferCNR.default<-function(object, window.size, cut.off, ref.cell, ref.shift, exp.limit,verbose,
                           ref.quantile, chr.use, step.len, band, gene_list){
  #1. Remove genes not in gene list
  raw_gene_count<-dim(object)[1]
  used_ensembl <- intersect(gene_list$gene_id, rownames(object))
  used_symbol <- intersect(gene_list$Gene, rownames(object))
  if(length(used_ensembl)>=length(used_symbol)){
    object<-object[match(used_ensembl,rownames(object)),]
    gene_list<-gene_list[match(used_ensembl,gene_list$gene_id),]
  }else{
    object<-object[match(used_symbol,rownames(object)),]
    gene_list<-gene_list[match(used_symbol,gene_list$Gene),]
  }
  if(verbose){print(paste0(raw_gene_count-dim(object)[1],' genes removed due to not in gene order file.'))}
  #2. Calculate TPM matrix
  total_counts<-apply(object,2,sum)
  TPM<-1000000*sweep(object, 2, total_counts, `/`)
  #3. Remove genes which has low global expression level
  gene_reserved<-sort(which(log2(apply(TPM,1,mean)+1)>=cut.off))
  TPM<-TPM[gene_reserved,]
  gene_list<-gene_list[gene_reserved,]
  if(verbose){print(paste0(dim(object)[1]-length(gene_reserved),' genes removed due to mean TPM less than cut off'))}
  if(verbose){print(paste0('Total ',dim(TPM)[1],' genes reserved after filter.'))}
  #4. Calculate infercnv relative expression matrix
  Exp_matrix<-log2((TPM/10)+1)
  Exp_gene_mean<-apply(Exp_matrix,1,mean)
  Exp_matrix<-apply(Exp_matrix,2,function(x){x-Exp_gene_mean})
  over_limit_ratio<-length(which(Exp_matrix>exp.limit|Exp_matrix<(-exp.limit)))/(dim(Exp_matrix)[1]*dim(Exp_matrix)[2])
  if(verbose){print(paste0('Expression limit is ',exp.limit,' and ',round(100*over_limit_ratio,2),'% of all elements are over limit.'))}
  Exp_matrix <- replace(Exp_matrix,Exp_matrix>exp.limit,exp.limit)
  Exp_matrix <- replace(Exp_matrix,Exp_matrix<(-exp.limit),(-exp.limit))
  if(verbose){print("Prepare infercnv input file over.")}
  #5. Split inferCNV into chromosome and band
  used_chr<-intersect(unique(gene_list$Chr),chr.use)
  gene_pos<-split(1:length(gene_list$Chr),factor(gene_list$Chr))[used_chr]
  if(band){
    bands<-split(1:length(gene_list$Chr),factor(substring(gene_list$band,0,1)))[c('p','q')]
    gene_pos<-unlist(lapply(gene_pos,function(x){list(p=intersect(x,bands$p),q=intersect(x,bands$q))}),recursive=F)
  }
  for(pos in names(gene_pos)){
    if(length(gene_pos[[pos]])<window.size){
      warning(paste0('Genes in ', pos, ' is less than window.size: ', window.size,', so skip it.'))
      gene_pos[pos]<-NULL
    }
  }
  #6. Calculate inferCNV
  CNV_list <- lapply(gene_pos,function(x){
    Pos_exp_matrix<-Exp_matrix[sort(x),]
    Pos_gene_list<-gene_list[sort(x),]
    CNV_start_pos<-rollapply(Pos_gene_list$GeneStart, window.size, min, by=step.len)
    CNV_end_pos<-rollapply(Pos_gene_list$GeneEnd, window.size, max, by=step.len)
    Pos_CNV<-apply(Pos_exp_matrix,2,function(y){rollapply(y, window.size, mean, by=step.len)})
    rownames(Pos_CNV)<-1:dim(Pos_CNV)[1]
    CNV_result<-list(CNV=t(Pos_CNV),start_pos=CNV_start_pos,end_pos=CNV_end_pos)
    return(CNV_result)
  })
  CNV<-as.data.frame(lapply(CNV_list,function(x){x$CNV}))
  start_pos<-unlist(lapply(CNV_list,function(x){x$start_pos}))
  end_pos<-unlist(lapply(CNV_list,function(x){x$end_pos}))
  CNV_cell_mean<-apply(CNV,1,mean)
  CNV<-t(apply(CNV,2,function(x){x-CNV_cell_mean}))
  if(!is.null(ref.cell)){
    ref_cells<-intersect(ref.cell,colnames(CNV))
    if(length(ref_cells)!=length(ref.cell)){
      warning('Some ref.cell not in seurat object!')
    }
    Max_baseline<-apply(CNV[,ref_cells,drop=F],1,function(x){quantile(x,ref.quantile)})
    Min_baseline<-apply(CNV[,ref_cells,drop=F],1,function(x){quantile(x,(1-ref.quantile))})
    CNV_pos<-rownames(CNV)
    CNV<-apply(CNV,2,function(x){
      cnvi<-rep(0,length(x))
      out_max<-which(x>(Max_baseline+ref.shift))
      out_min<-which(x<(Min_baseline-ref.shift))
      cnvi[out_max]<-(x-Max_baseline)[out_max]
      cnvi[out_min]<-(x-Min_baseline)[out_min]
      return(cnvi)
    })
    rownames(CNV)<-CNV_pos
  }
  if(verbose){print("Calculate infercnv over.")}
  chromosome<-gsub('\\..*$','',rownames(CNV),perl=T)
  meta_feature<-data.frame(chromosome)
  if(band){meta_feature$band<-gsub('^(.*\\..*)\\..*$','\\1',rownames(CNV),perl=T)}
  meta_feature$start_pos<-start_pos
  meta_feature$end_pos<-end_pos
  rownames(meta_feature)<-rownames(CNV)
  return(list(metadata=meta_feature,InferCNV=CNV))
}
InferCNR.Seurat<-function(object,window.size=51,cut.off=3.5,ref.cell=NULL,ref.shift=0.2,exp.limit=3,verbose=T,
                          ref.quantile=0.95,chr.use=paste0('chr',c(seq(1,22),'X')), step.len=1, band=F,
                          gene_list=Human_hg38.bed){
  exp.mat <- as.matrix(GetAssayData(object, slot = "counts"))
  infercnv <- InferCNR.default(exp.mat, window.size, cut.off, ref.cell, ref.shift, exp.limit,verbose,
                               ref.quantile, chr.use, step.len, band, gene_list)
  if("inferCNV" %in% names(object@assays)){object[["inferCNV"]] <- NULL}
  object[["inferCNV"]] <- CreateAssayObject(counts = infercnv$InferCNV)
  object@assays$inferCNV@meta.features <- infercnv$metadata
  object$inferCNV_score<-NA
  object$inferCNV_score<-apply(infercnv$InferCNV,2,function(x){mean(x^2)})
  return(object)
}

#' Map DNA CNR to Seurat InferCNV slot
#'
#' Map DNA CNR region to RNA to jointy analysis
#'
#'
#' @param CNR A CNR object
#' @param Seurat_obj A Seurat object with InferCNV slot to be mapped
#' @param CNV.ratio Specify a CNV ratio if you know it
#'
#' @return Returns a CNR object with mapped slot
#'
#' @export
#'
#' @rdname mapCNR2Seurat
#' @export mapCNR2Seurat
#'
mapCNR2Seurat <- function(CNR, Seurat_obj, CNV.ratio=NULL){
  Seurat_meta<-Seurat_obj@assays$inferCNV@meta.features
  chr<-Seurat_meta$chromosome
  start<-Seurat_meta$start_pos
  end<-Seurat_meta$end_pos
  InferCNV<-as.matrix(Seurat_obj@assays$inferCNV@counts)
  meta<-CNR$metadata
  match_pos<-lapply(1:length(chr),function(x){
    which(meta$chromosome==chr[x]&meta$end>start[x]&meta$start<end[x])
  })
  DNA_CNV<-t(as.data.frame(lapply(match_pos,function(x){colMeans(CNR$log2[x,])})))
  if(!is.numeric(CNV.ratio)){
    CNV.ratio<-Matrix::nnzero(InferCNV)/prod(dim(InferCNV))
    print(paste0("Estimated CNV percent is ",100*CNV.ratio,'%.'))
  }
  max_thes<-quantile(DNA_CNV,1-CNV.ratio/2)
  min_thes<-quantile(DNA_CNV,CNV.ratio/2)
  DNA_CNV <- replace(DNA_CNV,DNA_CNV>min_thes&DNA_CNV<max_thes,0)
  DNA_value<-as.numeric(DNA_CNV)
  RNA_value<-as.numeric(InferCNV)
  range_diff<-mean(abs(RNA_value[which(RNA_value!=0)]))/mean(abs(DNA_value[which(DNA_value!=0)]))
  DNA_CNV<-DNA_CNV*range_diff
  DNA_CNV <- replace(DNA_CNV,DNA_CNV>max(InferCNV),max(InferCNV))
  DNA_CNV <- replace(DNA_CNV,DNA_CNV<min(InferCNV),min(InferCNV))
  rownames(DNA_CNV)<-rownames(InferCNV)
  CNR$mapped<-as.data.frame(DNA_CNV)
  return(CNR)
}
NormalizeCNR.default<-function(object, chr, loc, window.size, probe.min,
                               mad.threshold, gamma, do.smooth){
  object <- cbind.data.frame(loc,object)
  chr_index <- unique(chr)
  object <- cbind.data.frame(match(chr,chr_index),object)
  CNR_smooth <- winsorize(data=object, tau=mad.threshold, k=floor(window.size/2),
                          assembly='hg38',verbose=FALSE)
  if(do.smooth){Y<-NULL}else{Y<-object}
  CNR_segment <- multipcf(data=CNR_smooth, Y=Y, gamma=gamma, assembly='hg38',verbose=FALSE)
  # probe.min not use in this version
  seg <- CNR_segment[,-c(1:5)]
  info <- CNR_segment[,c(1:5)]
  info$chrom <- chr_index[info$chrom]
  seg_list <- list(meta=info,segments=seg)
  return(seg_list)
}
NormalizeCNR.Seurat<-function(object, window.size=30, mad.threshold=2.5, probe.min=0,
                              gamma=40, do.smooth=F){
  raw_CNV <- as.matrix(GetAssayData(object, assay='inferCNV', slot='counts'))
  raw_CNV_feature <- object@assays$inferCNV@meta.features
  raw_CNV_loc <- round((raw_CNV_feature$start_pos+raw_CNV_feature$end_pos)/2)
  raw_CNV_chr <- raw_CNV_feature$chromosome
  seg_list <- NormalizeCNR.default(object=raw_CNV, chr=raw_CNV_chr, loc=raw_CNV_loc,
                                   window.size=window.size, mad.threshold=mad.threshold,
                                   gamma=gamma, do.smooth=do.smooth)
  object@assays$inferCNV@data <- as.matrix(seg_list[['segments']])
  seg <- seg_list[['meta']]$n.probes
  segment_region <- unlist(lapply(1:length(seg), function(i){rep(i,seg[i])}))
  object@assays$inferCNV@meta.features$segment <- segment_region
  return(object)
}
NormalizeCNR.CNR<-function(object, window.size=200, mad.threshold=2.5, probe.min=20,
                           gamma=200, do.smooth=T){
  CNR_mat <- object$log2
  CNR_loc <- round((object$metadata$start+object$metadata$end)/2)
  CNR_chr <- object$metadata$chromosome
  seg_list <- NormalizeCNR.default(object=CNR_mat, chr=CNR_chr, loc=CNR_loc,
                                   window.size=window.size, mad.threshold=mad.threshold,
                                   gamma=gamma, do.smooth=do.smooth)
  CNV_Info <- seg_list[['meta']]
  chromosome <- CNV_Info$chrom
  n_probe <- CNV_Info$n.probes
  end.probe <- cumsum(n_probe)
  start.probe <- c(1,head(end.probe,-1)+1)
  start_pos <- object$metadata$start[start.probe]
  end_pos <- object$metadata$end[end.probe]
  object$segment_info <- data.frame(chromosome,n_probe,start_pos,end_pos)
  object$segment <- seg_list[['segments']]
  return(object)
}

PloidyCNV.default<-function(object,method){
  data<-as.numeric(as.matrix(object))
  if(method=='kmean'){
    kmean_result<-kmeans(data,centers = c(min(data),0,max(data)))
    cluster_rank<-order(as.numeric(kmean_result$centers))
    lower_limit<-max(data[kmean_result$cluster==cluster_rank[1]])
    upper_limit<-min(data[kmean_result$cluster==cluster_rank[3]])
  }else if(method=='sd'){
    sd_val <- sd(data)
    lower_limit<-(-sd_val)
    upper_limit<-sd_val
  }else if(method=='hard'){
    lower_limit<-(-0.1)
    upper_limit<-0.1
  }else{stop('method should be hard, sd or kmean')}
  object <- replace(object,object<upper_limit&object>lower_limit,0)
  object <- replace(object,object>=upper_limit,1)
  object <- replace(object,object<=lower_limit,-1)
  return(object)
}
PloidyCNV.Seurat<-function(object,method='sd'){
  data <- object@assays$inferCNV@data
  ploidy_data <- PloidyCNV.default(object=data,method=method)
  object@assays$inferCNV@scale.data <- ploidy_data
  return(object)
}
PloidyCNV.CNR<-function(object,method='sd'){
  data <- object$segment
  ploidy_data <- PloidyCNV.default(object=data,method=method)
  object$ploidy_segment <- ploidy_data
  return(object)
}

RunCNVPCA.default<-function(object, CNV_info, method='length', dims=50){
  if(method=='n_probe'){
    segment_weight <- CNV_info$n_probe
  }else if(method=='length'){
    segment_weight <- CNV_info$end_pos-CNV_info$start_pos
  }else if(method=='raw'){
    segment_weight <- rep(1,dim(CNV_info)[1])
  }else{
    stop('not support yet')
  }
  new_obj <- t(object*segment_weight)
  pca.result<-prcomp(new_obj, rank. = dims, center = F)
  return(pca.result)
}
RunCNVPCA.Seurat<-function(object,method='length',dims=50){
  CNV_mat <- object@assays$inferCNV@data
  CNR_info<-object@assays$inferCNV@meta.features
  if(!all(c('segment','chromosome','start_pos','end_pos')) %in% colnames(CNR_info)){stop('segment, chromosome, start_pos, end_pos must in meta.features in inferCNV assay')}
  CNV_list<-split(CNR_info,cut(CNR_info$segment,0:max(CNR_info$segment)))
  CNV_info <- data.frame(t(data.frame(lapply(CNV_list,function(x){
    if(length(unique(x$segment))>1||length(unique(x$chromosome))>1){
      stop('segment or chromosome not unique')
    }
    return(c(min(x$start_pos),max(x$end_pos),length(x$segment)))
  }))))
  colnames(CNV_info)<-c('start_pos','end_pos','n_probe')
  pca.result <- RunCNVPCA.default(CNV_mat, CNV_info=CNV_info, method=method, dims=dims)
  sdev <- pca.results$sdev
  feature.loadings <- pca.results$rotation
  cell.embeddings <- pca.results$x
  rownames(x = feature.loadings) <- rownames(x = CNV_mat)
  colnames(x = feature.loadings) <- paste0('PC_', 1:dims)
  rownames(x = cell.embeddings) <- colnames(x = CNV_mat)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = 'inferCNV',
    stdev = sdev,
    key = 'PC_'
  )
  object[['CNVpca']]<-reduction.data
  return(object)
}
RunCNVPCA.CNR<-function(object,method='length',dims=50){
  CNV_mat <- object$segment
  CNV_info <- object$segment_info
  pca.result <- RunCNVPCA.default(CNV_mat, CNV_info=CNV_info, method=method, dims=dims)
  CNVpca<-pca.result$x
  colnames(CNVpca)<-paste0('PC_', 1:dims)
  rownames(CNVpca)<-colnames(CNV_mat)
  object$CNVpca<-CNVpca
  return(object)
}

ClusterCNV.default<-function(data.use, method='euclidean', k=5){
  clusters <- cutree(hclust(dist(data.use, method=method),method = 'ward.D2'), k=k)
  clusters<-paste('CNV',clusters,sep='_')
  names(clusters)<-colnames(object)
  return(clusters)
}
ClusterCNV.Seurat<-function(object, slot='pca', method='euclidean', k=5, dims=10){
  if(slot=='pca'){
    data.use <- Embeddings(object[['CNVpca']])[,1:dims]
  }else if(slot=='umap'){
    data.use <- Embeddings(object[['CNVumap']])
  }else if(slot=='raw'){
    data.use <- t(object@assays$inferCNV@data)
  }else{
    stop('not support yet')
  }
  Clus <- ClusterCNV.default(data.use, method=method, k=k)
  object$CNV_cluster<-Clus
  return(object)
}
ClusterCNV.CNR<-function(object, slot='pca', method='euclidean', k=5, dims=10){
  if(slot=='pca'){
    data.use <- object$CNVpca[,1:dims]
  }else if(slot=='umap'){
    data.use <-object$CNVumap
  }else if(slot=='raw'){
    data.use <- t(object$segment)
  }else{
    stop('not support yet')
  }
  Clus <- ClusterCNV.default(data.use, method=method, k=k)
  object$CNV_cluster<-Clus
  return(object)
}

RunCNVUMAP.default<-function(data, use.dims=1:10){
  data <- data[, use.dims]
  UMAP <- umap(data)
  return(UMAP)
}
RunCNVUMAP.Seurat<-function(object, use.dims=1:10){
  data.use <- Embeddings(object[['CNVpca']])
  CNVumap <- RunCNVUMAP.default(data.use, use.dims=use.dims)
  colnames(CNVumap)<-c('CNVumap_1','CNVumap_2')
  rownames(CNVumap)<-Cells(object)
  object[['CNVumap']] <- CreateDimReducObject(
    embeddings = CNVumap,
    key = 'CNVumap_',
    assay = 'inferCNV',
    global = TRUE
  )
  return(object)
}
RunCNVUMAP.CNR<-function(object, use.dims=1:10){
  data.use <- object$CNVpca
  CNVumap <- RunCNVUMAP.default(data.use, use.dims=use.dims)
  colnames(CNVumap)<-c('CNVumap_1','CNVumap_2')
  rownames(CNVumap)<-colnames(object$segment)
  object$CNVumap<-CNVumap
  return(object)
}
#' Get CNV changed between 2 CNV subcluster
#'
#' Get CNV changed between 2 CNV subcluster
#'
#'
#' @param object A CNR object or Seurat object
#' @param group.1 A character of a CNV subcluster
#' @param group.2 A character of another CNV subcluster
#' @param min.cutoff Min CNV change threshold
#'
#' @return Returns a vector contain changed CNV
#'
#' @export
#'
#' @rdname CNVChanged
#' @export CNVChanged
#'
CNVChanged<-function(object,group.1,group.2,min.cutoff=0.3){
  group1<-which(object$CNV_cluster==group.1)
  group2<-which(object$CNV_cluster==group.2)
  ploidy_data<-getCNVslot(object,'ploidyCNV')
  CNV_change<-apply(ploidy_data,1,function(y){
    p.val<-wilcox.test(x=y[group1],y=y[group2])$p.value
    distance<-mean(y[group2])-mean(y[group1])
    if(distance>min.cutoff&p.val<0.05){return('gain')}
    if(distance<(-min.cutoff)&p.val<0.05){return('loss')}
    return('stable')
  })
  return(list(gain=which(CNV_change=='gain'),loss=which(CNV_change=='loss')))
}

#' Get CNV region located chromosome
#'
#' Get CNV region located chromosome
#'
#'
#' @param object Some CNR object and Seurat object
#' @param region A character contain CNV region serial number
#'
#' @return Returns a character contain chromosomes which CNV regions located
#'
#' @export
#'
#' @rdname CNVRegionLocated
#' @export CNVRegionLocated
#'
CNVRegionLocated<-function(object,region){
  sites<-which(getCNVslot(object,'meta')$segment==region)
  chr<-unique(getCNVslot(object,'meta')$chromosome[sites])
  if(is.null(chr)){stop(paste0('region ',region,' not in object'))}
  if(length(chr)>1){stop(paste0('region ',region,' in more than 1 chromosome'))}
  return(chr)
}

#' Merge CNR object and Seurat object's CNV
#'
#' Merge CNR object and Seurat object's CNV
#'
#'
#' @param ... Some CNR object and Seurat object
#' @param CNVlist A list contain CNR object and Seurat object
#' @param chr.use Used chromosome
#'
#' @return Returns a list with megered CNV and its metadata
#'
#' @export
#'
#' @rdname MergeCNV
#' @export MergeCNV
#'
MergeCNV <- function(..., CNVlist = NULL,min.interval=100000,scale=T,
                     chr.use=paste0('chr',c(seq(1,22),'X'))){
  CNV_objs <- c(list(...), CNVlist)
  MLIST<-lapply(CNV_objs,function(x){
    if(class(x)=='CNR'){
      return(x$segment_info)
    }else if(class(x)=='Seurat'){
      CNR_info<-x@assays$inferCNV@meta.features
      CNV_list<-split(CNR_info,cut(CNR_info$segment,0:max(CNR_info$segment)))
      CNV_info <- data.frame(t(data.frame(lapply(CNV_list,function(x){
        if(length(unique(x$segment))>1||length(unique(x$chromosome))>1){
          stop('segment or chromosome not unique')
        }
        return(c(unique(x$chromosome),min(x$start_pos),max(x$end_pos),length(x$segment)))
      }))))
      colnames(CNV_info)<-c('chromosome','start_pos','end_pos','n_probe')
      CNV_info$start_pos<-as.numeric(CNV_info$start_pos)
      CNV_info$end_pos<-as.numeric(CNV_info$end_pos)
      CNV_info$n_probe<-as.numeric(CNV_info$n_probe)
      return(CNV_info)
    }else(stop('unrecognized object class in list'))
  })
  DLIST<-lapply(CNV_objs,function(x){
    if(class(x)=='CNR'){
      return(x$segment)
    }else if(class(x)=='Seurat'){
      return(x@assays$inferCNV@data)
    }else(stop('unrecognized object class in list'))
  })
  MLIST<-lapply(MLIST,function(data){
    pos_fixed<-lapply(chr.use,function(chr){
      x<-data[which(data$chromosome==chr),]
      if(dim(x)[1]>=2){
        for(i in 2:dim(x)[1]){
          if(x[i,'start_pos'] < x[i-1,'end_pos']){
            mid_point<-round((x[i,'start_pos']+x[i-1,'end_pos'])/2)
            x[i,'start_pos']<-mid_point
            x[i-1,'end_pos']<-mid_point
          }
        }
      }
      return(x)
    })
    return(do.call("rbind", pos_fixed))
  })
  chr.pos<-lapply(chr.use,function(chr){
    chr_pos<-sort(unique(unlist(lapply(MLIST,function(x){
      as.numeric(as.matrix(x[which(x$chromosome==chr),c('start_pos','end_pos')]))
    }))))
    chr_region<-c()
    for(i in 1:(length(chr_pos)-1)){
      if(chr_pos[i+1]-chr_pos[i]<=min.interval){next}
      chr_region<-rbind(chr_region,c(chr_pos[i],chr_pos[i+1]))
    }
    return(list(name=chr,region=chr_region))
  })
  CNV_match_pos<-lapply(MLIST,function(x){
    unlist(lapply(chr.pos,function(chr){
      apply(chr$region,1,function(i){
        CNV_match<-which(x$chromosome==chr$name&x$start_pos<=i[1]&x$end_pos>=i[2])
        return(ifelse(length(CNV_match)==1,CNV_match,NA))
      })
    }))
  })
  Merged_DLIST<-as.data.frame(lapply(1:length(MLIST),function(x){
    new_Data<-DLIST[[x]][CNV_match_pos[[x]],]
    new_Data<-replace(new_Data,is.na(new_Data),0)
    rownames(new_Data)<-1:dim(new_Data)[1]
    if(scale){
      upper_limit<-quantile(as.numeric(as.matrix(new_Data)),0.99)
      lower_limit<-quantile(as.numeric(as.matrix(new_Data)),0.01)
      new_Data<-replace(new_Data,new_Data>0,new_Data[new_Data>0]/upper_limit)
      new_Data<-replace(new_Data,new_Data<0,new_Data[new_Data<0]/abs(lower_limit))
    }
    return(new_Data)
  }))
  Merged_meta<-do.call("rbind", lapply(chr.pos,function(chr){
    region_length<-unlist(apply(chr$region,1,function(x){x[2]-x[1]}))
    for(i in 1:(length(region_length)-1)){
      middle_val<-round((chr$region[i,2]+chr$region[i+1,1])/2)
      chr$region[i,2]<-middle_val
      chr$region[i+1,1]<-middle_val
    }
    return(cbind(chr$name,chr$region))
  }))
  Merged_meta<-as.data.frame(Merged_meta)
  colnames(Merged_meta)<-c('chromosome','start_pos','end_pos')
  Merged_meta$start_pos<-as.numeric(Merged_meta$start_pos)
  Merged_meta$end_pos<-as.numeric(Merged_meta$end_pos)
  Merged_CNR<-list(segment_info=Merged_meta,segment=Merged_DLIST)
  class(Merged_CNR)<-'CNR'
  return(Merged_CNR)
}

#' Write CNV region into RNA assay's each gene
#'
#' Write CNV region into RNA assay's each gene
#'
#'
#' @param object A Seurat object
#'
#' @return Returns Seurat object with CNV region in each gene in RNA assay
#'
#' @export
#'
#' @rdname GeneCNVRegion
#' @export GeneCNVRegion
#'
GeneCNVRegion<-function(object,gene_list=Human_hg38.bed,verbose=T){
  #1. Remove genes not in gene list
  if(!all(c('gene_symbol','ensembl_id') %in% colnames(object@assays$RNA@meta.features))){
    stop('gene_symbol and ensembl_id should in RNA meta.feature slot')
  }
  raw_gene_name<-object@assays$RNA@meta.features$gene_symbol
  raw_gene_id<-object@assays$RNA@meta.features$ensembl_id
  object@assays$RNA@meta.features$chromosome<-NA
  object@assays$RNA@meta.features$start_pos<-NA
  object@assays$RNA@meta.features$end_pos<-NA
  used_ensembl <- intersect(gene_list$gene_id, raw_gene_id)
  used_symbol <- intersect(gene_list$Gene, raw_gene_name)
  if(length(used_ensembl)>=length(used_symbol)){
    gene_list<-gene_list[match(used_ensembl,gene_list$gene_id),]
    matched<-match(used_ensembl,raw_gene_id)
    object@assays$RNA@meta.features$chromosome[matched]<-gene_list$Chr
    object@assays$RNA@meta.features$start_pos[matched]<-gene_list$GeneStart
    object@assays$RNA@meta.features$end_pos[matched]<-gene_list$GeneEnd
  }else{
    gene_list<-gene_list[match(used_symbol,gene_list$Gene),]
    matched<-match(used_symbol,raw_gene_name)
    object@assays$RNA@meta.features$chromosome[matched]<-gene_list$Chr
    object@assays$RNA@meta.features$start_pos[matched]<-gene_list$GeneStart
    object@assays$RNA@meta.features$end_pos[matched]<-gene_list$GeneEnd
  }
  if(verbose){print(paste0(length(raw_gene_name)-length(matched),
                           ' genes removed due to not in gene order file.'))}
  CNV_info<-GetSeuratCNVinfo(object)
  gene_info<-object@assays$RNA@meta.features[matched,c('chromosome','start_pos','end_pos')]
  gene_in_CNV<-unlist(apply(gene_info,1,function(x){
    gene_CNV<-which(CNV_info$chromosome==x[1]&
                      as.numeric(CNV_info$start_pos)<=as.numeric(x[2])&
                      as.numeric(CNV_info$end_pos)>=as.numeric(x[3]))
    if(length(gene_CNV)==1){return(gene_CNV)}
    return(NA)
  }))
  object@assays$RNA@meta.features$CNV_region<-NA
  object@assays$RNA@meta.features$CNV_region[matched]<-gene_in_CNV
  return(object)
}

GetSeuratCNVinfo<-function(object){
  CNR_info<-object@assays$inferCNV@meta.features
  CNV_list<-split(CNR_info,cut(CNR_info$segment,0:max(CNR_info$segment)))
  CNV_info <- data.frame(t(data.frame(lapply(CNV_list,function(x){
    if(length(unique(x$segment))>1||length(unique(x$chromosome))>1){
      stop('segment or chromosome not unique')
    }
    return(c(unique(x$chromosome),min(x$start_pos),max(x$end_pos),length(x$segment)))
  }))))
  colnames(CNV_info)<-c('chromosome','start_pos','end_pos','n_probe')
  return(CNV_info)
}

#' Calculate each gene's CV and write it into RNA assay
#'
#' Calculate each gene's CV and write it into RNA assay
#'
#'
#' @param object A Seurat object
#' @param CNV_correction Do CNV correction or not(For Cells which have CNV)
#'
#' @return Returns Seurat object with CV value in each gene in RNA assay
#'
#' @export
#'
#' @rdname CalculateGeneCV
#' @export CalculateGeneCV
#'
CalculateGeneCV<-function(object,min.cell=3,min.rvalue=-Inf,CNV_correction=F){
  if('CNV_region' %in% colnames(object@assays$RNA@meta.features)){
    object@assays$RNA@meta.features$CNV_region<-0
  }
  CV<-unlist(lapply(1:dim(object@assays$RNA@data)[1],function(i){
    CNVRegion<-object@assays$RNA@meta.features$CNV_region[i]
    if(is.na(CNVRegion)){return(NA)}
    CNVExp<-object@assays$RNA@data[i,]
    CNV_data<-data.frame(expression=exp(CNVExp),CNV=object@assays$inferCNV@data[CNVRegion,])
    if(length(which(CNV_data$expression>1))<min.cell){return(NA)}
    lm_result<-lm(expression~CNV, CNV_data)
    r_value<-sign(coef(lm_result)[2])*sqrt(summary(lm_result)$r.squared)
    if(r_value<min.rvalue){return(NA)}
    fixed_exp<-CNV_data$expression
    if(CNV_correction){
      fixed_exp<-fixed_exp-(CNV_data$CNV*coef(lm_result)[2])
    }
    return(sd(fixed_exp)/abs(mean(fixed_exp)))
  }))
  object@assays$RNA@meta.features$CV<-CV
  return(object)
}

#' A estimated parsimony tree by CNV
#'
#' Estimate a parsimony tree from CNV data
#'
#'
#' @param segment A segment
#'
#' @return Returns a estimated parsimony tree
#'
#' @export
#'
#' @rdname parsimony_tree
#' @export parsimony_tree
#'
parsimony_tree <- function(segment){
  weight <- log10(segment$n.probes)
  cal_matrix <- segment[,-c(1:5)]
  level=unique(as.numeric(as.matrix(cal_matrix)))
  CNV_phy <- phyDat(cal_matrix*weight, type = 'USER',levels=level)
  treeRatchet <- pratchet(CNV_phy, trace = 0)
  return(treeRatchet)
}

