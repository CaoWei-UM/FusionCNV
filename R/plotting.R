PlotCNR.default<-function(object, cell.anno, chr.anno, max.limit, min.limit, clust,
                          cutree=NA, fixed.scaling=NA, chr.scaling=NA, file=NA){
  if(!any(is.na(c(fixed.scaling,chr.scaling)))){
    stop('Cannot use two scaling methods in same time.')
  }else if(any(is.numeric(c(fixed.scaling,chr.scaling)))){
    chr.list <- lapply(unique(chr.anno),function(x){sort(which(chr.anno==x))})
    names(chr.list)<-unique(chr.anno)
    if(is.numeric(chr.scaling)){
      if(chr.scaling > min(unlist(lapply(chr.list,length)))){stop('Too many scaling parts in chr.scaling!')}
      if(chr.scaling > 1){
        chr.list <- lapply(chr.list,function(x){split(x,cut(1:length(x),chr.scaling))})
      }
      chr.anno<-unlist(lapply(unique(chr.anno),function(x){rep(x,chr.scaling)}))
    }else if(is.numeric(fixed.scaling)){
      chr.list <- lapply(chr.list,function(x){split(x,cut(1:length(x),length(x)/fixed.scaling))})
      chr.anno <- unlist(lapply(unique(chr.anno),function(x){rep(x,length(chr.list[[x]]))}))
    }
    object <- t(as.data.frame(lapply(unlist(chr.list,recursive=F),function(x){colMeans(object[x,,drop=F])})))
  }
  chr_gaps<-as.integer(cumsum(head(table(chr.anno)[unique(chr.anno)],-1)))
  chr.anno[-c(c(0,chr_gaps)+as.integer(table(chr.anno)[unique(chr.anno)]/2))]<-''
  if(!(clust|all(is.na(cell.anno)))){
    ranks<-apply(cell.anno,1,function(x){paste(x, collapse = "-")})
    names(ranks)<-rownames(cell.anno)
    cell_order <- order(ranks[colnames(object)])
    object <- object[,cell_order]
  }
  object <- replace(object,object>max.limit,max.limit)
  object <- replace(object,object<(-max.limit),(-max.limit))
  object <- replace(object,object>(-min.limit)&object<min.limit,0)
  mode_value<-modeest::mfv(object)[1]
  paletteLength=50
  myBreaks <- c(seq(min(object), mode_value, length.out=ceiling(paletteLength/2) + 1),
                seq(mode_value+(max(object)/paletteLength), max(object), length.out=floor(paletteLength/2)))
  if(is.na(file)){silent=F}else{silent=T}
  if(!is.na(chr.scaling)&&chr.scaling==1){chr_gaps=NULL}
  CNV_plot<-pheatmap(object,
                     color = colorRampPalette(c("blue","white","red"))(paletteLength),
                     cluster_rows=F, cluster_cols = clust, show_rownames=T, show_colnames=F,
                     gaps_row=chr_gaps, labels_row =chr.anno, cutree_cols=ifelse(clust,cutree,NA),
                     annotation_col=cell.anno, breaks=myBreaks,
                     filename=file, silent=silent)
  if(is.na(file)){return(CNV_plot)}else{return()}
}
PlotCNR.CNR<-function(object,max.limit=10,min.limit=1,fixed.scaling=1000,chr.scaling=NA, cell.anno=NA,clust=T,cell.use=NA,
                      cutree=NA,chr.use=paste0('chr',c(seq(1,22),'X')),file=NA){
  if(is.na(cell.use)){cell.use <- colnames(object$log2)}
  chr.use <- intersect(chr.use,unique(object$metadata$chromosome))
  pos_used <- sort(which(object$metadata$chromosome %in% chr.use))
  CNR_mat<-object$log2[pos_used,cell.use,drop=F]
  chr.anno <- object$metadata$chromosome[pos_used]
  fig<-PlotCNR.default(CNR_mat, cell.anno, chr.anno, max.limit, min.limit,
                       clust, cutree, fixed.scaling, chr.scaling, file)
  return(fig)
}
PlotCNR.Seurat<-function(object,annotation=NA,chr.use=paste0('chr',c(seq(1,22),'X')),cell.use=NA,max.limit=1,min.limit=0.1,clust=T,
                         band.use=T, cutree=NA,chr.scaling=5,fixed.scaling=NA,file=NA){
  if(!'inferCNV' %in% names(object@assays)){stop('Must run InferCNV first!')}
  metadata <- object@assays$inferCNV@meta.features
  chr.use <- intersect(chr.use,unique(metadata$chromosome))
  pos_used <- which(metadata$chromosome %in% chr.use)
  if(is.na(cell.use)){cell.use <- Cells(object)}
  if(band.use & 'band' %in% colnames(metadata)){chr.anno<-metadata$band[pos_used]}else{chr.anno<-metadata$chromosome[pos_used]}
  used_anno <- intersect(annotation,colnames(object@meta.data))
  if(length(used_anno)>0){
    cell.anno <- object@meta.data[cell.use,used_anno,drop=F]
  }else{cell.anno <- NA}
  CNR_mat<-as.matrix(GetAssayData(object,slot='counts',assay='inferCNV'))
  CNR_mat <- CNR_mat[pos_used,cell.use,drop=F]
  fig<-PlotCNR.default(CNR_mat, cell.anno, chr.anno, max.limit, min.limit,
                       clust, cutree, fixed.scaling, chr.scaling, file)
  return(fig)
}

PlotCNV.default <- function(CNV_meta, CNV_data, anno=NULL, max.limit=Inf, min.limit=0, cluster=F,
                            chr.length=Human_hg38.length, chr.use=paste0('chr',c(seq(1,22),'X')), file=NULL){
  Coord<-unname(c(0,cumsum(chr.length[chr.use])))
  names(Coord)<-c(chr.use,'end')
  if(dim(CNV_data)[1]!=dim(CNV_meta)[1]){stop('CNV_meta should have same rows with CNV_data!')}
  used_pos<-sort(which(CNV_meta$chromosome %in% chr.use))
  CNV_data<-CNV_data[used_pos,]
  CNV_meta<-CNV_meta[used_pos,]
  chr_label<-round(Coord[chr.use] + chr.length[chr.use]/2)
  chr_seg <- data.frame(y=Coord, xend=dim(CNV_data)[2])
  #trans chromosome position to global position
  CNV_meta$start_pos<-as.numeric(CNV_meta$start_pos)
  CNV_meta$end_pos<-as.numeric(CNV_meta$end_pos)
  CNV_meta$start.pos<-Coord[CNV_meta$chromosome] + CNV_meta$start_pos
  CNV_meta$end.pos<-Coord[CNV_meta$chromosome] + CNV_meta$end_pos
  #deal overlap region
  for(i in 2:dim(CNV_meta)[1]){
    if(CNV_meta[i,'start.pos'] < CNV_meta[i-1,'end.pos']){
      mid_point<-round((CNV_meta[i,'start.pos']+CNV_meta[i-1,'end.pos'])/2)
      CNV_meta[i,'start.pos']<-mid_point+1
      CNV_meta[i-1,'end.pos']<-mid_point
    }
  }
  if(cluster){
    CNV_data<-CNV_data[,hclust(dist(t(CNV_data)), method = "ward.D2" )$order]
  }else{
    if(!is.null(anno)){CNV_data<-CNV_data[,do.call(order, anno)]}
  }
  CNV_bind<-cbind(CNV_meta[,c('start.pos','end.pos')],CNV_data)
  CNV_plot<-melt(CNV_bind, id.vars = c('start.pos','end.pos'),
                 variable.name = 'cell_id',value.name='CNV_score')
  CNV_plot$CNV_score[which(CNV_plot$CNV_score>max.limit)]<-max.limit
  CNV_plot$CNV_score[which(CNV_plot$CNV_score<(-max.limit))]<-(-max.limit)
  CNV_plot$CNV_score[which(abs(CNV_plot$CNV_score)<min.limit)]<-0
  #start plot
  fig <- ggplot(CNV_plot) +
    geom_linerange(aes_string(ymin='start.pos', ymax='end.pos', x='cell_id', col='CNV_score'), size=5) +
    scale_y_continuous(breaks=chr_label, labels=names(chr_label)) +
    scale_color_gradient2(low = 'blue', high = 'red', midpoint = 0, mid = "white") +
    theme(panel.background=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),
          axis.title.y=element_blank(),axis.text.x = element_text(angle = -90,vjust=0.5)) +
    geom_segment(aes_string(x=1, xend='xend', y='y', yend='y'), data=chr_seg, col='black')
  if(!is.null(anno)){
    CNV_anno <- anno[colnames(CNV_data),,drop=F]
    CNV_anno$cell_id <- colnames(CNV_data)
    CNV_anno_rank<-seq(from = -1e8, to = -1e7, length.out = dim(anno)[2]+1)
    for(i in 1:dim(anno)[2]){
      fig <- fig + new_scale_color() +
        geom_linerange(aes_string(ymin=CNV_anno_rank[i], ymax=CNV_anno_rank[i+1],
                                  x='cell_id', col=colnames(anno)[i]), size=5, data=CNV_anno)
    }
  }
  fig <- fig + coord_flip()
  if (!is.null(file)) {
    ggsave(file, fig, device='png', width=20, height=8)
    return(NULL)
  } else {
    return(fig)
  }
}

PlotCNV.Seurat<-function(object, slot='CNV', anno=NULL, max.limit=Inf, min.limit=0,
                         cluster=F, chr.length=Human_hg38.length,
                         chr.use=paste0('chr',c(seq(1,22),'X')), file=NULL){
  if(slot=='CNV'){
    CNV_mat <- object@assays$inferCNV@data
  }else if(slot=='ploidyCNV'){
    CNV_mat <- object@assays$inferCNV@scale.data
  }else{stop('slot should be \\"CNV\\" or \\"p\\"')}
  CNR_info<-object@assays$inferCNV@meta.features
  CNV_list<-split(CNR_info,cut(CNR_info$segment,0:max(CNR_info$segment)))
  CNV_info <- data.frame(t(data.frame(lapply(CNV_list,function(x){
    if(length(unique(x$segment))>1||length(unique(x$chromosome))>1){
      stop('segment or chromosome not unique')
    }
    return(c(unique(x$chromosome),min(x$start_pos),max(x$end_pos),length(x$segment)))
  }))))
  colnames(CNV_info)<-c('chromosome','start_pos','end_pos','probes')
  rownames(CNV_info)<-1:dim(CNV_info)[1]
  CNV_info$start_pos<-as.numeric(CNV_info$start_pos)
  CNV_info$end_pos<-as.numeric(CNV_info$end_pos)
  if(!is.null(anno)){anno<-object@meta.data[,anno,drop=F]}
  fig<-PlotCNV.default(CNV_info, CNV_mat, anno=anno, max.limit=max.limit, cluster=cluster,
                       min.limit=min.limit, chr.length=chr.length, chr.use=chr.use, file=file)
  return(fig)
}

PlotCNV.CNR<-function(object, slot='CNV', anno=NULL, max.limit=Inf, min.limit=0,
                      cluster=F, chr.length=Human_hg38.length,
                      chr.use=paste0('chr',c(seq(1,22),'X')), file=NULL){
  if(slot=='CNV'){
    CNV_mat <- as.matrix(object$segment)
  }else if(slot=='ploidyCNV'){
    CNV_mat <- as.matrix(object$ploidy_segment)
  }else{stop('slot should be \\"CNV\\" or \\"ploidyCNV\\"')}
  CNV_info <- object$segment_info
  CNV_info$start_pos<-as.numeric(CNV_info$start_pos)
  CNV_info$end_pos<-as.numeric(CNV_info$end_pos)
  if(!is.null(anno)){anno<-as.data.frame(object[anno])}
  fig<-PlotCNV.default(CNV_info, CNV_mat, anno=anno, max.limit=max.limit, cluster=cluster,
                       min.limit=min.limit, chr.length=chr.length, chr.use=chr.use, file=file)
  return(fig)
}

PlotMST.default<-function(CNV_mat, cluster, weight, root='Normal', plot=T){
  direct_matrix<-function(mat,root){
    all_nodes<-colnames(mat)
    if(!(root %in% all_nodes)){stop('root not exist')}
    pro_list=c()
    scaning_list<-c(root)
    while(length(pro_list)<length(all_nodes)){
      temp<-c()
      for(i in scaning_list){
        pro_list<-append(pro_list,i)
        after_node<-setdiff(all_nodes[which(mat[i,]>0)],pro_list)
        mat[i,setdiff(all_nodes,after_node)]<-0
        temp<-append(temp,after_node)
      }
      scaning_list<-unique(temp)
    }
    return(mat)
  }
  cluster_data<-lapply(unique(cluster),function(x){
    rowMeans(CNV_mat[,which(cluster==x)])
  })
  cluster_data<-as.data.frame(cluster_data)
  colnames(cluster_data)<-unique(cluster)
  cluster_data<-cluster_data*weight
  obj_dist<-dist(t(as.matrix(cluster_data)))
  mst<-as.matrix(obj_dist)*ape::mst(obj_dist)
  direct_mst<-round(direct_matrix(mst,root),digits = 2)
  net = network(direct_mst, directed = T,ignore.eval = F,names.eval = "weights")
  if(plot){
    if(root=='Normal'){
      node_color<-c('grey',hue_pal()(dim(mst)[2]-1))
    }else{
      node_color<-hue_pal()(dim(mst)[2])
    }
    plot_node_dist<-sqrt(as.matrix(obj_dist))
    plot_node_dist<-0.5*(plot_node_dist/max(plot_node_dist))+0.5
    diag(plot_node_dist)<-0
    fig<-ggnet2(net,mode='kamadakawai',layout.par = list(elen=plot_node_dist),arrow.size = 12,arrow.gap = 0.1, color=node_color, size=30, label = TRUE, edge.label = "weights",edge.color = 'black')+ scale_x_continuous(breaks = NULL, limits = c(-0.1,1.1)) + scale_y_continuous(breaks = NULL, limits = c(-0.1,1.1))
    return(fig)
  }
  direct_list<-lapply(1:length(net$mel),function(x){
    return(c(colnames(direct_mst)[net$mel[[x]]$outl],colnames(direct_mst)[net$mel[[x]]$inl],net$mel[[x]]$atl$weights))
  })
  return(direct_list)
}
PlotMST.Seurat<-function(object,seg.weight=T,group.by='CNV_cluster',root='Normal',plot=T){
  CNV_mat<-object@assays$inferCNV@data
  cluster<-object@meta.data[[group.by]]
  if(seg.weight){
    weight<-as.numeric(table(object@assays$inferCNV@meta.features$segment))
  }else{weight<-rep(1,dim(CNV_mat)[1])}
  fig<-PlotMST.default(CNV_mat, cluster, weight, root=root, plot=plot)
  return(fig)
}
PlotMST.CNR<-function(object,seg.weight=T,group.by='CNV_cluster',root='Normal',plot=T){
  CNV_mat<-as.matrix(object$segment)
  cluster<-object[[group.by]]
  if(seg.weight){
    weight<-as.numric(object$n_probe)
  }else{weight<-rep(1,dim(CNV_mat)[1])}
  fig<-PlotMST.default(CNV_mat, cluster, weight, root=root, plot=plot)
  return(fig)
}
#' Plot Expression level vs CNV scatter plot
#'
#' Plot Expression level vs CNV scatter plot
#'
#'
#' @param object A Seurat object with CNV_region slot
#' @param gene gene name to plot
#'
#' @return Returns ggplot2 object
#'
#' @export
#'
#' @rdname PlotExpCNV
#' @export PlotExpCNV
#'
PlotExpCNV<-function(object,gene,region=NULL,plot=T){
  if(is.null(region)){
    gene_in_region<-object@assays$RNA@meta.features[gene,'CNV_region']
  }else{gene_in_region<-region}
  gene_in_chr<-object@assays$RNA@meta.features[gene,'chromosome']
  plot_obj<-cbind(expm1(object@assays$RNA@data[gene,]),object@assays$inferCNV@data[gene_in_region,])
  colnames(plot_obj)<-c('expression','CNV')
  plot_obj<-as.data.frame(plot_obj)
  drop_ratio<-sum(plot_obj$expression==0)/length(plot_obj$expression)
  lm_result<-lm(expression~CNV, plot_obj)
  r2_value<-round(summary(lm_result)$adj.r.squared,2)
  r_value<-as.numeric(sign(coef(lm_result)[2])*sqrt(summary(lm_result)$r.squared))
  p_value<-round(summary(lm_result)$coefficients[2,4],2)
  formula = paste0('y=',round(coef(lm_result)[2],2),'x+',round(coef(lm_result)[1],2))
  if(plot){
    plot_obj$CNV_cluster<-as.factor(object$CNV_cluster)
    fig<-ggplot(plot_obj, aes(CNV, expression,colour=CNV_cluster)) +
      geom_point() +
      geom_smooth(method = "lm",inherit.aes = F,aes(CNV, expression),formula = y~x) +
      labs(title = paste0("CNV vs expression in gene: ",gene),
           subtitle = paste0('chr: ',gene_in_chr,'   CNV_region: ',gene_in_region,'   r^2=',r2_value,'   r.val=',r_value,'   p.val=',p_value,'   ',formula))+ylab('normalized expression level')
    return(fig)
  }
  return(list(gene=gene,chromosome=gene_in_chr,CNV_region=gene_in_region,r2=r2_value,r.val=r_value,p.val=p_value,formula=formula,drop_ratio=drop_ratio))
}

#' Plot CNV region volcano plot of expression and variation
#'
#' Plot CNV region volcano plot of expression and variation
#'
#'
#' @param object A Seurat object with CNV_region slot
#'
#' @return Returns ggplot2 object
#'
#' @export
#'
#' @rdname PlotCNVRegion
#' @export PlotCNVRegion
#'
PlotCNVRegion<-function(object, topn.label=10,var.limit=0.001,mean.limit=0.01){
  CNV_stat<-data.frame(Var=unlist(apply(object@assays$inferCNV@data,1,var)),
                       Mean=unlist(apply(object@assays$inferCNV@data,1,mean)))
  CNV_stat$stat<-'noCNV'
  CNV_stat$stat[which(CNV_stat$Var>var.limit&CNV_stat$Mean>mean.limit)]<-'amplify'
  CNV_stat$stat[which(CNV_stat$Var>var.limit&CNV_stat$Mean<(-mean.limit))]<-'deletion'
  point_colors<-c("red","blue","black")
  names(point_colors)<-c('amplify','deletion', "noCNV")
  fig<-ggplot(CNV_stat, aes(x=Mean, y=Var)) +
    geom_point(aes(colour = stat), size=2) +
    scale_colour_manual(values = point_colors)+
    labs(title ="CNV region statistic")
  topn.label<-as.integer(topn.label)
  if(is.integer(topn.label)){
    top_genes1<-CNV_stat[which(CNV_stat$stat=='amplify'),]
    top_genes2<-CNV_stat[which(CNV_stat$stat=='deletion'),]
    top_genes1<-head(rownames(top_genes1)[order(top_genes1$Mean,-abs(top_genes1$Var))],topn.label)
    top_genes2<-head(rownames(top_genes2)[order(top_genes2$Mean,-abs(top_genes2$Var))],topn.label)
    fig<-LabelPoints(plot=fig, points=c(top_genes1,top_genes2),
                     repel = TRUE,xnudge = 0, ynudge = 0)
  }
  return(fig)
}

#' Plot CV line plot of 2 objects
#'
#' Plot CV line plot of 2 objects
#'
#'
#' @param object.1 A Seurat object with CNV_region slot
#' @param object.2 A Seurat object with CNV_region slot
#' @param show.exp Show expression level of each gene or not
#' @param region Specify regions or use all region
#'
#' @return Returns ggplot2 object
#'
#' @export
#'
#' @rdname PlotCVCompare
#' @export PlotCVCompare
#'
PlotCVCompare<-function(object.1,object.2,show.exp=F,region=NA){
  CV1<-object.1@assays$RNA@meta.features$CV
  CV2<-object.2@assays$RNA@meta.features$CV
  CV_ratio<-CV1/CV2
  exist_CV<-which(!is.na(CV_ratio))
  used_site<-exist_CV[order(CV_ratio[exist_CV])]
  if(!is.na(region)){
    used_site<-intersect(used_site,which(object.1@assays$RNA@meta.features$CNV_region %in% region))
  }
  CV1_exp<-apply(object.1@assays$RNA@data[used_site,],1,mean)
  CV2_exp<-apply(object.2@assays$RNA@data[used_site,],1,mean)
  gene_rank<-1:length(CV1)/length(CV1)
  plot_obj<-data.frame(logCV=c(log(CV1),log(CV2)),
                       belong=c(rep(object.1@project.name,length(CV1)),
                                rep(object.2@project.name,length(CV2))),
                       gene=c(gene_rank,gene_rank),
                       exp=c(CV1_exp,CV2_exp))
  plot_obj$logCV<-as.numeric(plot_obj$logCV)
  fig<-ggplot(plot_obj,aes(x=gene,color=belong))+
    geom_vline(aes(xintercept=gene_rank[which.min(abs((CV1/CV2)[used_site]-1))]),colour='black')+
    scale_x_continuous(labels = scales::percent)
  if(show.exp){
    fig<-fig+geom_line(aes(y=logCV))+geom_line(aes(y=exp),linetype="l")+
      scale_y_continuous(sec.axis = sec_axis(~.+max(c(CV1,CV2)),name='Exp'))
  }
  return(fig)
}

#' Plot CV violin plot of 2 objects divide by reference object
#'
#' Plot CV violin plot of 2 objects divide by reference object
#'
#'
#' @param object.1 A Seurat object with CNV_region slot
#' @param object.2 A Seurat object with CNV_region slot
#' @param object.ref A reference Seurat object with CNV_region slot
#' @param show.exp Show expression level of each gene or not
#' @param region Specify regions or use all region
#'
#' @return Returns ggplot2 object
#'
#' @export
#'
#' @rdname PlotCVTest
#' @export PlotCVTest
#'
PlotCVTest<-function(object.1,object.2,object.ref,show.exp=F,region=NA){
  if(!is.na(region)){
    used_site<-which(object.1@assays$RNA@meta.features$CNV_region %in% region)
  }
  CV1<-object.1@assays$RNA@meta.features$CV[used_site]
  CV2<-object.2@assays$RNA@meta.features$CV[used_site]
  CVref<-object.ref@assays$RNA@meta.features$CV[used_site]
  CV1_ratio<-CV1/CVref
  CV2_ratio<-CV2/CVref
  EXP1<-apply(object.1@assays$RNA@data,1,mean)
  EXP2<-apply(object.2@assays$RNA@data,1,mean)
  EXPref<-apply(object.ref@assays$RNA@data,1,mean)
  EXP1_ratio<-(EXP1/EXPref)[which(!is.na(CV1_ratio))]
  EXP2_ratio<-(EXP2/EXPref)[which(!is.na(CV2_ratio))]
  CV1_ratio<-CV1_ratio[which(!is.na(CV1_ratio))]
  CV2_ratio<-CV2_ratio[which(!is.na(CV2_ratio))]
  # I want use color to show expression level ,to be done
  if(mean(CV2_ratio)>mean(CV1_ratio)){tit='cancer cell have less variation'}else{tit='not'}
  plot_obj<-data.frame(CV=c(CV1_ratio,CV2_ratio),
                       EXP=c(EXP1_ratio,EXP2_ratio),
                       sample=c(rep('CV1/ref',length(CV1_ratio)),
                                rep('CV2/ref',length(CV2_ratio))))
  plot_obj$CV<-as.numeric(plot_obj$CV)
  plot_obj$EXP<-as.numeric(plot_obj$EXP)
  fig<-ggviolin(plot_obj, x = 'sample', y = 'CV', fill = 'sample',title=tit,
                add = "boxplot", add.params = list(fill = "white"))
  signifi_list<-combn(names(table(plot_obj[['sample']])),2,simplify = F)
  fig<-fig + stat_compare_means(comparisons = signifi_list, label = "p.format")
  fig<-fig+ geom_hline(aes(yintercept=1),colour='red')
  return(fig)
}
