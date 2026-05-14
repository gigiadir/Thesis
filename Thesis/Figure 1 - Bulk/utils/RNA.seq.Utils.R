### RNA.seq.Utiles.R

library(GSVA)
library(AUCell)
library(sva)
#library(data.table)


################################## percentile.expression.CCPM
# Using: ccpm.RNA.TPM.df, ccpm.RNA.TPM.Signatures.df (more inclusive...)
#        ccpm_TCGA.RNA.rank
percentile.expression.CCPM <- function(expression_mat = ccpm.RNA.TPM.df, signatures_mat = ccpm.RNA.TPM.Signatures.df, mGene = "TP53", mSig = NULL , mSample_ID = "401_T1", 
                                       samples.in.distribution = subcohort.list[["Sample_ID.RNA.Mets.TopPurity"]], 
                                       tcga_expression_mat = ccpm_TCGA.RNA.rank, is.rank.among.TCGA = F, Receptor.subset = NULL){
  
  if(!is.null(mGene)){
    if(!is.null(samples.in.distribution))  samples.in.distribution = colnames(expression_mat)
    ecdf.dist <- ecdf(expression_mat[mGene, samples.in.distribution ])
    mExpression = expression_mat[mGene,mSample_ID]
    mPercentile = ecdf.dist(mExpression)
    return(list("mPercentile" = mPercentile, "mExpression" = mExpression))
  }
  
}

################################## my.log.normalize
my.log.normalize <- function(expression_mat, log.base = 2, countPer = 1e6, is.center = F){
  if(log.base==2){
    return( log2(scale(expression_mat, center=is.center, scale=colSums(expression_mat)/countPer ) + 1) )
  }else if(log.base==10){
    return( log10(scale(expression_mat, center=is.center, scale=colSums(expression_mat)/countPer ) + 1) )
  }
}

################################## combat.wrapper
combat.wrapper <- function(expression_mat, sample.names.batch.list, is.first.log = T, is.forec.return.positive = F, log.base = 2){
  if(is.first.log){
    if(log.base==2){
      expression_mat = log2(expression_mat+1)
    }else if(log.base==10){
      expression_mat = log10(expression_mat+1)
    }
  }
  
  batch.names = names(sample.names.batch.list)
  
  batch = rep("NA",ncol(expression_mat)) # init batch with the size of mm_merge
  for(myset in batch.names){
    batch[colnames(expression_mat) %in% sample.names.batch.list[[myset]]] <- myset
  }
  print(table(batch))
  
  myvars <- apply( expression_mat  ,1, var,na.rm=TRUE)
  keep.idx = which(myvars>0)
  print(paste0( "remove ", nrow(expression_mat) - length(keep.idx) ," genes with no variance"))
  expression_mat = expression_mat[keep.idx,]
  
  cols = colnames(expression_mat)
  rows = row.names(expression_mat)
  expression_mat = matrix(as.numeric(unlist(expression_mat)),nrow=nrow(expression_mat))
  colnames(expression_mat) = cols
  row.names(expression_mat) = rows
  
  cleandat <- ComBat(dat=expression_mat, batch=batch, mod=NULL) # use svaseq for RNA-seq (counts/reads)
  
  frac.of.negative =  length(which(cleandat<0)) / length(cleandat)
  min.of.negative =  min(cleandat)
  print(paste0( "frac.of.negative= ", frac.of.negative ," min.of.negative= ", min.of.negative))
  
  # Negative values occasionally show up when adjusting for batch or surrogate values. These are expression values that were likely already very low. 
  if(is.forec.return.positive)  cleandat[cleandat<0] = 0 
  
  return(cleandat)
  
}





### calc_cpm
# # https://hemberg-lab.github.io/scRNA.seq.course/normalization-for-library-size.html
calc_cpm <-   function (expr_mat) {
    norm_factor <- colSums(expr_mat)
    return(t(t(expr_mat)/norm_factor)) * 10^6
}

### calc_uq - calculate upper quartile for library normalization
calc_uq <- function (expr_mat) 
{
  UQ <- function(x) {
    quantile(x[x > 0], 0.75)
  }
  uq <- unlist(apply(expr_mat, 2, UQ))
  norm_factor <- uq/median(uq)
  return(t(t(expr_mat)/norm_factor))
}

### calc_sf - calculate upper quartile for library normalization
# The size factor (SF) was proposed and popularized by DESeq (Anders and Huber (2010)).
calc_sf <-   function (expr_mat) {
  geomeans <- exp(rowMeans(log(expr_mat)))
  SF <- function(cnts) {
    median((cnts/geomeans)[(is.finite(geomeans) & geomeans > 
                              0)])
  }
  norm_factor <- apply(expr_mat, 2, SF)
  return(t(t(expr_mat)/norm_factor))
}


### Get log2 and centered versions of data
get.log2.and.centered.versions.of.data <- function(in.df, gct.out.file = NULL, tsv.out.file = NULL, is.center = T , is.scale = F){
  # in - data.frame or file
  # out - log2, centered, 
  # write gct, tsv or both
  # log transform
  print(dim(in.df))
  in.df.log = log2( (in.df) +1) 
  
  # center
  if(is.center){
    in.df.log.centered = sweep(in.df.log, 1, rowMeans(in.df.log))
    if(is.scale){
      in.df.log.centered <- t(apply( in.df.log  ,1, scale ,scale=is.scale)) # can do centering only, but takes longer compute
      colnames(in.df.log.centered) = colnames(in.df.log)
    }
    out.df = in.df.log.centered
  }else{
    out.df = in.df.log
  }
  
  if( !is.null(tsv.out.file)) write.table(out.df, file = tsv.out.file ,quote=F,sep = "\t", col.names =NA  ) 
  if( !is.null(gct.out.file)) write.gct.2(gct.data.frame = as.data.frame(out.df), descs = row.names(out.df),  filename=gct.out.file)
  
  return(out.df)
}




## update colnames by SampleSheet
update.colnames.by.SampleSheet <- function(df,common.prefix, in.SampleSheet.tsv){
  SampleSheet.df = read.delim(in.SampleSheet.tsv, stringsAsFactors = F, row.names = 1)
  colnames.org = colnames(df)
  colnames.new = gsub(common.prefix,"",colnames.org)
  colnames.new = make.names(paste0(SampleSheet.df[ gsub("S","A",colnames.new) ,"Sample_Name"],".", colnames.new))
  colnames(df) = colnames.new
  return(df)
}

## produce matrix - encapsulate as function
produce.TMP.matrix.from.RSEM.agg <- function(file.tpm, inputListFile, out.tpm.mat, is.inputListFile.with.header = F){
  tpm = read.delim(file.tpm, stringsAsFactors = F)
  print(dim(tpm))
  # find the number of genes, loop-size
  first.gene.name = tpm[1,1]
  
  first.gene.name.idx = which(tpm$gene_id == first.gene.name)
  number.of.genes = first.gene.name.idx[2]-1
  gene.names = tpm$gene_id[1:number.of.genes]
  
  
  df_inputListFile = read.delim(inputListFile, stringsAsFactors = F, header = is.inputListFile.with.header)
  print(dim(df_inputListFile))
  num.of.samples = dim(df_inputListFile)[1]
  
  stopifnot(dim(tpm)[1] == number.of.genes*num.of.samples)
  
  tpm.df.as.matrix = data.frame(matrix(nrow = number.of.genes, ncol=num.of.samples))
  row.names(tpm.df.as.matrix) = gene.names
  colnames(tpm.df.as.matrix) = df_inputListFile[,1]
  
  tpm = tpm[,2] # keep only the tmp values !!!
  
  for(i in 1:num.of.samples){
    sample.name = df_inputListFile[i,1]
    tpm.of.genes.per.sample = tpm[ ((i-1)*number.of.genes+1) : (i*number.of.genes) ]
    tpm.df.as.matrix[,sample.name] = tpm.of.genes.per.sample
  }
  print(dim(tpm.df.as.matrix))
  write.table(tpm.df.as.matrix, file = out.tpm.mat ,quote=F,sep = "\t", col.names =NA  ) ### NEW!
}




#source(file.path(code_source,"DISSECTOR_lib.v13.R")) # Note: can't insall packages on Broad
######## Reads a gene expression dataset in GCT format and converts it into an R data frame
MSIG.Gct2Frame <- function(filename = "NULL") {
  ds <- read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T, na.strings = "")
  descs <- ds[,1]
  ds <- ds[-1]
  row.names <- row.names(ds)
  names <- names(ds)
  return(list(ds = ds, row.names = row.names, descs = descs, names = names))
}


##################################  top.Var.Genes
top.Var.Genes <- function(df=NULL, input_file_gct=NA, number_of_top_var_genes = 1000, num.of.isoexpressed.bins = 50, is.topvar.per.isoexpressed.bin = T, output_file_gct = NA, is.already.centered = F){
  
  stopifnot(!is.null(df) | !is.na(input_file_gct));
  
  if(is.null(input_file_gct)){
    stopifnot(file.exists(input_file_gct))
    mat.obj <- MSIG.Gct2Frame(filename = input_file_gct)
    df_ <- data.matrix(mat.obj$ds)
  }else{
    df_ = df
  }
  num.of.genes = dim(df_)[1]   
  
  if(is.topvar.per.isoexpressed.bin){
    top_var_names_all_bins = NULL
    score_mat_tt = df_
    mymedian <- apply(score_mat_tt,1, mean,na.rm=TRUE) # change to mean ??? Works... 
    mymedian <- sort(mymedian,decreasing=TRUE)
    eq_expression_bins <- split(mymedian, ceiling(seq_along(mymedian)/ (num.of.genes/num.of.isoexpressed.bins) ))
    if(!is.already.centered)   stopifnot(max(eq_expression_bins[[length(eq_expression_bins)]])>0) # Error
    for( bin_name in  names(eq_expression_bins) ){
      myvars <- apply( score_mat_tt[ names(eq_expression_bins[[bin_name]]) ,]  ,1, var,na.rm=TRUE) 
      myvars <- sort(myvars,decreasing=TRUE) 
      myvars <- myvars[1: floor(number_of_top_var_genes/( length(eq_expression_bins) )  ) ]
      top_var_names_all_bins <- c(top_var_names_all_bins, names(myvars))
    }    
  }else{
    myvars <- apply( df_  ,1, var,na.rm=TRUE)
    myvars <- sort(myvars,decreasing=TRUE)
    top_var_names_all_bins <- names(myvars)[1:number_of_top_var_genes]
  }
  df_var =df_[top_var_names_all_bins,]
  df_var =df_var[which(rowSums(is.na(df_var))==0),] # patch, otherwise return NA ...
  
  if(!is.na(output_file_gct))     write.gct.2(gct.data.frame = as.data.frame(df_var), descs = row.names(df_var),  filename=output_file_gct)

  return(df_var)
  
}


##################################  ```{run.NMF}
run.NMF <- function( input_file_gct, number_of_top_var_genes = 1000, num.of.isoexpressed.bins = 20, is.topvar.per.isoexpressed.bin = T, is.NMF.from.rank = F,
                     is.run.specified.num.components =T, number_of_components=4,
                     min_k = 3, max_k = 8, rep.num = 2, nrun = 50, res.folder = NA, num.of.top.genes.per.component.for.geneset = 200 ){
  
  gene.set.databases= c("~/Dropbox/Cancer_DB/GSEA/breast_sigs.symbols.gmt",
                        "~/Dropbox/Cancer_DB/GSEA/my_pathways.symbols.gmt",
                        "~/Dropbox/Cancer_DB/GSEA/Breast.Stover.FinalSignatures_Symbol.gmt")
  
  if(is.na(res.folder)){
    #res.folder =  dirname(input_file_gct)
    res.folder = gsub(".gct", "", input_file_gct)
    dir.create(res.folder)
  } 
  mat.obj <- MSIG.Gct2Frame(filename = input_file_gct)
  df_ <- data.matrix(mat.obj$ds)
  num.of.genes = dim(df_)[1]
  
  if(is.topvar.per.isoexpressed.bin){
    top_var_names_all_bins = NULL
    score_mat_tt = df_
    mymedian <- apply(score_mat_tt,1, median,na.rm=TRUE) 
    mymedian <- sort(mymedian,decreasing=TRUE)
    eq_expression_bins <- split(mymedian, ceiling(seq_along(mymedian)/ (num.of.genes/num.of.isoexpressed.bins) ))
    for( bin_name in  names(eq_expression_bins) ){
      myvars <- apply( score_mat_tt[ names(eq_expression_bins[[bin_name]]) ,]  ,1, var,na.rm=TRUE) 
      myvars <- sort(myvars,decreasing=TRUE) 
      myvars <- myvars[1: floor(number_of_top_var_genes/( length(eq_expression_bins) )  ) ]
      top_var_names_all_bins <- c(top_var_names_all_bins, names(myvars))
    }    
  }else{
    myvars <- apply( df_  ,1, var,na.rm=TRUE)
    myvars <- sort(myvars,decreasing=TRUE)
    top_var_names_all_bins <- names(myvars)[1:number_of_top_var_genes]
  }
  df_var =df_[top_var_names_all_bins,]

  
  if(is.NMF.from.rank){
    df_rank <- apply(df_var,2, rank)
    df.2.NMF = df_rank
  }else{
    df.2.NMF = df_var
    
  }
  colnames(df.2.NMF) = gsub("^X", "", colnames(df.2.NMF))
  cat(dim(df.2.NMF))
  
  if(is.run.specified.num.components){
    nmf_ <- nmf(df.2.NMF, rank=number_of_components) 
    cat("Residuals", residuals(nmf_), "\n")
    
    W<-nmf_ @fit@W # LM genes accross k component, 1:k matix of components over the 978 landmark genes (m) [features] (use this to tag the components with pathway/states)
    H<-nmf_ @fit@H # ORFs across k components, 1:n(samples) all peturbations [samples] matrix over the k components (use this as input for t-SNE, heatmap, etc.)
    
    row.names(H) <- paste0("C",1:nrow(H))
    colnames(W) <- paste0("C",1:ncol(W))
    key = paste0("in.",dim(df.2.NMF)[1],".",dim(df.2.NMF)[2],".out.C", number_of_components)
    
    write.csv(H, file=paste(res.folder,"/", paste("H.Samples.coeff.",key,".csv",sep=""),sep=""))
    write.csv(W, file=paste(res.folder,"/", paste("W.basis",key,".csv",sep=""),sep=""),row.names=T)
    
    write.gct.2(gct.data.frame = as.data.frame(H), descs = row.names(H),  filename=paste(res.folder,"/", paste("H.Samples.coeff.",key,".gct",sep=""),sep=""))
    write.gct.2(gct.data.frame = as.data.frame(W), descs = row.names(W),  filename=paste(res.folder,"/", paste("W.basis.",key,".gct",sep=""),sep=""))
    
    ### scale the W matrix...(aming to get gene sets)
    W.zscore <- t(apply(W,1, scale)) # scale by row (compare to other components)
    colnames(W.zscore) = colnames(W)
    
    pdf(file=paste(res.folder,"/","NNMF",".",number_of_components,".pdf",sep=""))
    #basismap(nmf_ , main=paste("basismap",paste("#comp=", number_of_components))) 
    #coefmap(nmf_ ,fontsize = 8,cexCol=5, main=paste("coefmap",paste("#comp=", number_of_components))); # heatmap resulting in cluster the samples (overexpressed ORFs) with their weights across the basis (components)
    pheatmap(H, fontsize_row =6, main=paste("heatmap H", paste("#comp=", number_of_components)) ,fontsize = 8)
    pheatmap(W, fontsize_row =4, main=paste("heatmap W", paste("#comp=", number_of_components)) ,fontsize = 4)   
    
    tsne_nmf <- tsne(t(H), perplexity = 30) # default=30 (more == more "cohesive", "ball-like")
    tsne_H = data.frame(tsne_nmf)
    colnames(tsne_H) = c("x","y")
    tsne_H[,"lab"] = colnames(H);  #colnames(rank_mat)
    my_ggplot <- ggplot(data = tsne_H, aes(x = x, y = y, label = lab ) ) + theme_bw() + 
      #geom_point(aes(colour = factor(state), shape = factor(cond)) , size = 1  ) +
      geom_point() +
      geom_text(data = within(tsne_H, c(y <- y+.01, x <- x-.01)), hjust = 0, vjust = 0, size=2  ) + 
      labs(title="") 
    ggsave(filename = paste(res.folder,"/","NNMF t-SNE ",number_of_components,key,".pdf",sep=""));
    
    
    while(dev.cur()>1){dev.off();}
    
    
    ### write gmt file with gene sets for each of the components
    sink(file=paste0(res.folder,"/","Top.genes.in.components.txt"))
    for(i in 1:ncol(W)){
      cat(i,"\n")
      C.order = W.zscore[ order(W.zscore[,i],decreasing = T) ,i]
      cat( paste0( names( C.order[1:num.of.top.genes.per.component.for.geneset] ), sep="\t"), "\n"  )
    }
    sink()

    ### ssGSEA over the components
    W.file = paste(res.folder,"/", paste("W.basis.",key,".gct",sep=""),sep="")
    W.ssGSEA.file = paste0(gsub(".gct","",W.file), paste0(".Features.ssGSEA.","db",".gct" ))
    run.ssGSEA(file2read=W.file,  gene.set.databases=gene.set.databases)
    z.score.expression.tgt.by.ref(gct_tgt=W.ssGSEA.file, gct_ref=W.ssGSEA.file)
        

  }else{
    if(!exists("estimate_H_k")) estimate_H_k = list();
    first_iter = 1 #update after running
    ## Run 
    #nrun = 100 # "Each point on the graph was obtained from 50 runs of the Brunet et"
    while(dev.cur()>1){dev.off();}
    pdf(file=paste(res.folder,"/","NMF",".dim.",dim(df.2.NMF)[1],".",dim(df.2.NMF)[2],".k.",min_k,".",max_k,".",nrun,".pdf",sep=""))
    
    for(i in seq(first_iter, (rep.num+first_iter-1)) ){
      cat(i,"\n")
      estimate_H_k[[i]] <- nmfEstimateRank(df.2.NMF, seq(min_k,max_k), method='brunet', nrun=nrun) # , seed=123456
      plot(estimate_H_k[[i]])
      first_iter = first_iter+1
    }
    while(dev.cur()>1){dev.off();}    
  }
}

##################################  ```produce.weights.H.matrix.from.z.scores.gct
produce.weights.H.matrix.from.z.scores.gct <- function( input_file_gct, row.names.2.keep, outfile=NA  ){
  input_file_gct = "/Users/ofirc/Dropbox/Broad/Research_projects/Breast.ER+.Wagle/Firehose/RNA_expression_level/ER_Pos_RNASeq_ind_set.rpkm.48644.101.filtered_samples.hugo.tracker.names.intersected.Features.ssGSEA.breast_sigs.my.h.Stover.ref_self.z.scored.gct"
  if(is.na(outfile)) outfile = gsub(".gct", ".weights.gct", input_file_gct)
  row.names.2.keep = qw("EGUCHI_CELL_CYCLE_RB1_TARGETS  ERBB2.DESMEDT.18698033 DOANE_BREAST_CANCER_ESR1_UP")
  
  
  mat.obj <- MSIG.Gct2Frame(filename = input_file_gct)
  df_ <- data.matrix(mat.obj$ds)
  df_ = df_[ row.names(df_) %in%  row.names.2.keep, ]
  df_.exp = 2^df_

  df_.weight = sweep(df_.exp,2,colSums(df_.exp),`/`)
  
  write.gct.2(gct.data.frame = df_.weight , descs = row.names(df_.weight),
              filename=outfile )
  
  
}





##################################  ```{r cast columns as Factor Numeric function}
cast.columns.as.double <- function( dt, num_cols ){
  
  library(data.table)
  dt = data.table(dt)
  
  for (col in num_cols){
    e = substitute(X := as.numeric(X), list(X = as.symbol(col)))
    dt[ , eval(e)]
  }
  return(dt)
}  



# RNA-seq RPKM or Counts - perp matrix
################################## RNA-seq RPKM or Counts - perp matrix
RNAseq.names.filterQC.prep.matrix <- function(input_file_gct, imput_file.QC, type="rpkm", input_file.Firehose.samples.R.object,
                                              min.num.of.Genes.Detected = 10000, max.allowed.Cont=5 ,
                                              agg.fun.to.use = "median"){
  #input_file.Firehose.samples.R.object="/Users/ofirc/Dropbox/Broad/Research_projects/Breast.ER+.Wagle/Trackers.Summary/Tracker.df.list.05062016.RData"
  #type="reads"
  df.qc = read.delim(imput_file.QC, stringsAsFactors = F)
  dataset.2 <- MSIG.Gct2Frame(filename = input_file_gct)
  m2.breast.matrix <- data.matrix(dataset.2$ds)
  print(dim(m2.breast.matrix))
  
  ### filter out poor quality RNA-seq
  cols2keep = colnames(m2.breast.matrix) %in% make.names( df.qc$Sample[which(df.qc$Genes.Detected>=min.num.of.Genes.Detected)] )
  m2.breast.matrix = m2.breast.matrix[ ,  cols2keep ]
  print(dim(m2.breast.matrix))
  
  m2.breast.matrix = cast.columns.as.double(m2.breast.matrix, colnames(m2.breast.matrix)) # cast into double if reads
  m2.breast = data.frame(m2.breast.matrix, "Gene"=dataset.2$descs)
  
  # merge duplicate rows, Ave.
  m2_agg_duplicates <- aggregate.duplicated.row.names(m2.breast, dulicated.names.cols="Gene", agg.fun.to.use = agg.fun.to.use)
  
  # simply skip duplicated. Good after sorting 
  #df_merge.omit.dup = df_merge_[! duplicated(df_merge_$hgnc_symbol) ,]
  print(dim(m2_agg_duplicates))
  head(m2_agg_duplicates)
  row.names(m2_agg_duplicates) = m2_agg_duplicates$Gene
  m2_agg_duplicates = m2_agg_duplicates[,!grepl("Gene",colnames(m2_agg_duplicates))]
  
  # {r write gct}
  suffix.ext = gsub( paste0(".+",type,"."),"",input_file_gct)
  file_hugo_gct = paste0(gsub(suffix.ext, paste0(dim(m2_agg_duplicates), collapse = ".") ,input_file_gct),".hugo.gct") # out
  write.gct.2(gct.data.frame=m2_agg_duplicates ,descs =row.names(m2_agg_duplicates),filename= file_hugo_gct ) 
  
  # {r sample names map convert}
  load(file = input_file.Firehose.samples.R.object)
  df.t = Tracker.df.list[["Firehose.samples"]]
  
  my.colnames.2.convert = make.names(colnames(m2_agg_duplicates))
  Tracker.names = my.colnames.2.convert #init

  Tracker.names  = llply(my.colnames.2.convert  ,function(x)
    if(length(grep( x, make.names(df.t$sample_id_Firehose.samples)) )>0){
      df.t[which( make.names( df.t$sample_id_Firehose.samples) ==x ),"External.ID.fix_Firehose.samples"]  
    }else{
      x
    }
  ) %>%  unlist()
  
  BRCA.05246.mat = m2_agg_duplicates
  if(length(which(Tracker.names==""))==0){
    colnames(BRCA.05246.mat) = Tracker.names
    file_hugo_gct.tracker.names = paste0(gsub(".gct","",file_hugo_gct),".hugo.tracker.names.gct") # out
    write.gct.2(gct.data.frame=BRCA.05246.mat ,descs =row.names(BRCA.05246.mat),filename= file_hugo_gct.tracker.names ) 
  }  
  
  
  # {r filter for contamination}
  RNA.samples = colnames(BRCA.05246.mat)
  
  # which RNA samples are not in tracker
  RNA.samples[! RNA.samples %in% df_tracker$Sample_ID ] # "remove_from_cohort"]] = qw("298 302 308 312 023 041 317")
  
  cols2keep = colnames(BRCA.05246.mat) %in% df_tracker$Sample_ID [which(!df_tracker$RNA.Cont>=max.allowed.Cont)]
  BRCA.05246.mat = BRCA.05246.mat[ ,  cols2keep ]
  print(dim(BRCA.05246.mat))
  
  suffix.ext = gsub( paste0(".+",type,"."),"",file_hugo_gct.tracker.names)
  file_hugo_gct.tracker.names = paste0(gsub(suffix.ext, paste0(dim(BRCA.05246.mat), collapse = ".") ,file_hugo_gct.tracker.names),".filtered_samples.hugo.tracker.names.gct") # out
  write.gct.2(gct.data.frame=BRCA.05246.mat ,descs =row.names(BRCA.05246.mat),filename= file_hugo_gct.tracker.names ) 
  
}












################################## require NextSeq.PlateMap.Barcodes.df_
ERK2.get.Sample.ID <- function(Samples, is.return.rep = T, is.only.Allele = F){
  mbs = stri_sub(Samples,-6,-1) # get last 6 char, the barcode
  Lane = stri_sub(Samples,-8,-8) # get the 8th from last char, the Lane
  Pool =  llply(Samples  ,function(x) if(stri_sub(x,-9,-9)=="s"){ return("A")} else {return("B")} ) %>%   unlist()
  
  Sample.ID = llply(paste(Pool,Lane, mbs,sep = ".") ,function(x) NextSeq.PlateMap.Barcodes.df_[ NextSeq.PlateMap.Barcodes.df_$Key == x, "Sample.ID"  ] ) %>%   unlist()
  Sample.ID = gsub("\\+",".TRE", Sample.ID )
  Sample.ID = gsub("\\-","", Sample.ID )
  
  Allele = gsub("\\.\\S+","", Sample.ID )  
  Sample.ID.rep = paste(Sample.ID, Lane, Pool,sep = ".")
  if(is.only.Allele){
    return(Allele)
  }
  
  if(is.return.rep) return(Sample.ID.rep)
  else return(Sample.ID)
}



################################## intersect
intersect.genes.write.gct <- function(gct_tgt, gct_ref, gene.names.from.desc.tgt = T, gene.names.from.desc.ref = T) {
  
  tgt <- MSIG.Gct2Frame(filename = gct_tgt)
  tgt.mat <- data.matrix(tgt$ds) 
  if(gene.names.from.desc.tgt) row.names(tgt.mat) = tgt$descs
  print(dim(tgt.mat))
  
  ref <- MSIG.Gct2Frame(filename = gct_ref)
  ref.mat <- data.matrix(ref$ds) 
  if(gene.names.from.desc.ref) row.names(ref.mat) = ref$descs
  print(dim(ref.mat))
  
  overlap <- intersect(row.names(tgt.mat), row.names(ref.mat)) # do all that are in data_set_names_info
  
  tgt.mat.intersected = tgt.mat[overlap,]
  ref.mat.intersected = ref.mat[overlap,]
  
  write.gct.2(gct.data.frame=tgt.mat.intersected ,descs =row.names(tgt.mat.intersected),
              filename= gsub(".gct",".intersected.gct",gct_tgt) )
  write.gct.2(gct.data.frame=ref.mat.intersected ,descs =row.names(ref.mat.intersected),
              filename= gsub(".gct",".intersected.gct",gct_ref) ) 
  
}

################################## runner.fgsea
runner.fgsea <- function(fgsea.Ranks, 
                       pathways=GSDB.list.of.lists[["msigdb.c2.augmented"]],
                       nperm=10000,
                       minSize=3, maxSize = 1000,
                       pval.return = 0.01){
  
  library(fgsea)
  fgsea.Ranks = fgsea.Ranks[!duplicated(names(fgsea.Ranks))]
  
  fgseaRes <- fgsea(pathways = pathways,
                    stats = fgsea.Ranks,
                    minSize=minSize,
                    maxSize=maxSize,
                    nperm=nperm,
                    nproc=4) 
  dd = data.frame(subset(fgseaRes, pval <= pval.return))
  if(nrow(dd)>0){
    for(i in 1:nrow(dd))    dd[i,"leadingEdgeString"] = paste0( dd$leadingEdge[[i]], collapse = ", ")
    dd = dd[order(dd$NES, decreasing = T),]
    dd.writeable = dd[,qw("pathway pval padj ES NES nMoreExtreme size leadingEdgeString")]
    return(dd.writeable)
  }else{
    warning("No results")
  }
}

################################## run.ssGSEA
run.ssGSEA <- function(file2read, gene.set.databases, gsea_type = "db", sample.norm.type= "rank", weight = 0.75, output_feature_file=NA) {
  
  print(file2read)
  if(is.na(output_feature_file)){
    output_feature_file = paste0(gsub(".gct","",file2read), paste0(".Features.ssGSEA.",gsea_type,".gct" ))
  }
  
  OPAM.project.dataset.6(  
    input.ds                 = file2read, # input, normalized gct
    output.ds                = output_feature_file,
    gene.set.databases       = gene.set.databases,
    gene.set.selection       = "ALL",
    sample.norm.type         = sample.norm.type,  # "rank", "log" or "log.rank"
    weight                   = weight,
    statistic                = "area.under.RES",
    output.score.type        = "ES",
    combine.mode             = "combine.add",  # "combine.off", "combine.replace", "combine.add"
    nperm                    =  1,
    min.overlap              =  1,
    correl.type              =  "z.score")             # "rank", "z.score", "symm.rank"
  
}

################################## run.ssGSEA
run.ssGSEA.df.wrapper <- function(in.df, gene.set.databases, gsea_type = "db", sample.norm.type= "rank", weight = 0.75, output_feature_file=NA) {
  file2read = "~/in.df.tmp.gct"
  write.gct.2(gct.data.frame = in.df , descs = row.names(in.df),
              filename=file2read )
  if(is.na(output_feature_file)){
    output_feature_file = paste0(gsub(".gct","",file2read), paste0(".Features.ssGSEA.",gsea_type,".gct" ))
  }
  run.ssGSEA (file2read, gene.set.databases, gsea_type =gsea_type, sample.norm.type= sample.norm.type, weight =weight, output_feature_file=output_feature_file)
  ssGSEA.mat <- data.matrix(MSIG.Gct2Frame(filename = output_feature_file)$ds) 
  return(ssGSEA.mat)
  
}


################################## z.score.expression.tgt.by.ref
z.score.expression.tgt.by.ref <- function(gct_tgt, gct_ref, gene.names.from.desc.tgt=F, gene.names.from.desc.ref=F, ref.code = "", outfile = NA){
  
  tgt <- MSIG.Gct2Frame(filename = gct_tgt)
  tgt.mat <- data.matrix(tgt$ds) 
  if(gene.names.from.desc.tgt) row.names(tgt.mat) = tgt$descs
  print(dim(tgt.mat))
  
  ref <- MSIG.Gct2Frame(filename = gct_ref)
  ref.mat <- data.matrix(ref$ds) 
  if(gene.names.from.desc.ref) row.names(ref.mat) = ref$descs
  print(dim(ref.mat))
  
  mean_ref = data.frame(Means=rowMeans(ref.mat) )
  sd_ref = data.frame(Means=rowSds(ref.mat) )
  
  tgt.mat.z= sweep(tgt.mat,1,t(mean_ref), FUN = "-")
  tgt.mat.zz= sweep(tgt.mat.z,1,t(sd_ref), FUN = "/")

  ref.code.name = paste0(".ref",ref.code,".z.scored.gct")
  if(is.na(outfile)){
    outfile = gsub(".gct",ref.code.name,gct_tgt)
  }
  write.gct.2(gct.data.frame = tgt.mat.zz , descs = row.names(tgt.mat.zz),
              filename=outfile )
}

################################## calc.Signatures.Sum
calc.Signatures.Sum <- function(df, GSDB.list=NULL,  gene.sets.file = "~/Dropbox/Wagle.CCPM.ERpos/Data.sources/Gene.sets.Cell.Types/breast_sigs.inclusive.Stover.symbols.gmt"
                                ,min.gene.in.geneset = 3, mock.genes = qw("ESR1 ERBB2"), is.run.zscore = T ){
  data.info = data.frame( t(df[mock.genes,]) )
  expression.mat = as.matrix(df) # range from 0 to max log normalized(?) TPM values ~5, summary(rowMeans(countData@raw.data))
  
  if(is.null(GSDB.list)){
    GSDB.list = list()
    GSDB <- Read.GeneSets.db(gene.sets.file, thres.min = 0, thres.max = 10000, gene.names = NULL)
    row.names(GSDB$gs) = GSDB$gs.names
    for(sig in row.names(GSDB$gs)){
      Genes.in.set = GSDB$gs[sig,] %>% unique()
      Genes.in.set = Genes.in.set[!Genes.in.set %in% c("","null")]
      N.genes = length(Genes.in.set)
      if(N.genes>=min.gene.in.geneset & N.genes<=1000 & ! sig %in% names(GSDB.list))  GSDB.list[[sig]] = Genes.in.set
    } 
  }
  for(sig in names(GSDB.list)){
    Genes.in.set = GSDB.list[[sig]]
    if(length(which(row.names(expression.mat) %in% Genes.in.set))>=min.gene.in.geneset){
      data.info[  , sig] = as.numeric(t(((colMeans(expression.mat[ row.names(expression.mat) %in% Genes.in.set , ]) ))))
    }
  }
  
  df.signatures = t(data.info)
  for(mock.gene in mock.genes) df.signatures = df.signatures[-which(row.names(df.signatures) == mock.gene),]
  if(is.run.zscore){
    df.signatures.z = z.score.expression.tgt.by.ref.data.frames(df.signatures, df.signatures)
    return(df.signatures.z)    
  }else{
    return(df.signatures)
  }
}




################################## correlated.Gene.Set
correlated.Gene.Set <- function(df_ = ccpm.RNA.TPM.df, 
                                 gene.set.A, gene.set.B,
                                gene.set.Ctrl=NULL,
                                 is.fast.Sum = F,
                                 cor.method = "spearman"){

  # data - TCGA, CCPM(mets), CCLE, SC
  # TCGA.BRCA.mat[,grep_multi_patterns(TCGA.ERpos.Pts, colnames(TCGA.BRCA.mat))], 
  # METABRIC.RNA, 
  # ccpm.RNA.TPM.df[,]
  library(GSVA)
  library(matrixStats) # needed if using fast proxy
  
  if(is.fast.Sum){
    warning("Using calc.Signatures.Sum is not normalized... expect inflated correlation")
    gene.set.A.OE = calc.Signatures.Sum (df_ , GSDB.list=list("A"=gene.set.A), min.gene.in.geneset = 1, is.run.zscore = F)
    gene.set.B.OE = calc.Signatures.Sum (df_ , GSDB.list=list("B"=gene.set.B), min.gene.in.geneset = 1, is.run.zscore = F)
    if(!is.null(gene.set.Ctrl)){
      gene.set.Ctrl.OE = calc.Signatures.Sum (df_ , GSDB.list=list("Ctrl"=gene.set.Ctrl), min.gene.in.geneset = 1, is.run.zscore = F)
      my.df = data.frame(my.B=gene.set.B.OE, "my.Ctrl" = gene.set.Ctrl.OE)
      m1 <- lm(my.df$my.B ~ my.df$my.Ctrl)  #Create a linear model
      gene.set.B.OE = resid(m1)      
    }
  }else{
    gene.set.A.OE = calc.Signatures.GSVA (df_ , GSDB.list=list("A"=gene.set.A), min.gene.in.geneset = 1)
    gene.set.B.OE = calc.Signatures.GSVA (df_ , GSDB.list=list("B"=gene.set.B), min.gene.in.geneset = 1)
    if(!is.null(gene.set.Ctrl)){
      gene.set.Ctrl.OE = calc.Signatures.GSVA (df_ , GSDB.list=list("Ctrl"=gene.set.Ctrl), min.gene.in.geneset = 1)
      my.df = data.frame(gene.set.B.OE["B",], gene.set.Ctrl.OE["Ctrl",]  )
      m1 <- lm(my.df[,1] ~ my.df[,2])  #Create a linear model
      gene.set.B.OE = t( resid(m1) )      
    }    
  }
  
  my.cor = cor.test(t(gene.set.A.OE), t(gene.set.B.OE), method = cor.method)
  return(my.cor)

}




################################## calc.Signatures.GSVA
calc.Signatures.GSVA <- function(df, GSDB.list=NULL,  gene.sets.file = "~/Dropbox/Wagle.CCPM.ERpos/Data.sources/Gene.sets.Cell.Types/breast_sigs.inclusive.Stover.symbols.gmt", 
                                 min.gene.in.geneset = 3, max.gene.in.geneset=1000,  is.counts = F, method = "gsva",parallel.sz=4){
  expression.mat = as.matrix(df) # range from 0 to max log normalized(?) TPM values ~5, summary(rowMeans(countData@raw.data))
  
  if(is.null(GSDB.list)){
    GSDB.list = list()
    GSDB <- Read.GeneSets.db(gene.sets.file, thres.min = 0, thres.max = 10000, gene.names = NULL)
    row.names(GSDB$gs) = GSDB$gs.names
    for(sig in row.names(GSDB$gs)){
      Genes.in.set = GSDB$gs[sig,] %>% unique()
      Genes.in.set = Genes.in.set[!Genes.in.set %in% c("","null")]
      N.genes = length(Genes.in.set)
      if(N.genes>=min.gene.in.geneset & N.genes<=1000 & ! sig %in% names(GSDB.list))  GSDB.list[[sig]] = Genes.in.set
    } 
  }
  
  #library(GSVAdata) # c2BroadSets
  #library(GSEABase)
  #library(org.Hs.eg.db)
  #library(GSVA)
  kcdf="Gaussian"
  if(is.counts) kcdf="Poisson"
  enrichment.scores <- gsva((expression.mat), GSDB.list, kcdf=kcdf, mx.diff=TRUE, verbose=TRUE, parallel.sz=parallel.sz, min.sz=min.gene.in.geneset, max.sz=max.gene.in.geneset, method = method) #
  
  return(enrichment.scores)
}





################################## calc.Signatures.AUC
# https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html
# Default - top.n.percent.of.genes.in.ranking.used= 5%
calc.Signatures.AUC <- function(df, GSDB.list=NULL , 
                                 top.n.percent.of.genes.in.ranking.used = 20 , nCores =4 , is.new.version = T, is.return.z.score = F, is.setdiff.GSDB.names=NULL){
  expression.mat = as.matrix(df)
  if(!is.null(is.setdiff.GSDB.names)) GSDB.list = GSDB.list[ setdiff(names(GSDB.list), is.setdiff.GSDB.names) ]
  
  cells_rankings <- AUCell_buildRankings(expression.mat, plotStats =F)
  cells_AUC <- AUCell_calcAUC(GSDB.list, cells_rankings, aucMaxRank=nrow(cells_rankings)* (top.n.percent.of.genes.in.ranking.used/100), nCores=nCores )
  
  if(is.new.version){
    enrichment.scores =  cells_AUC@assays@data$AUC
  }else{
    enrichment.scores =  cells_AUC@assays$data$AUC
  }
  if(is.return.z.score){
    enrichment.scores = z.score.expression.tgt.by.ref.data.frames(enrichment.scores, enrichment.scores)
  }
  return(enrichment.scores)
}

################################## z.score.expression.tgt.by.ref.data.frames
z.score.expression.tgt.by.ref.data.frames <- function(tgt.mat, ref.mat){
  
  mean_ref = data.frame(Means=rowMeans(ref.mat) )
  sd_ref = data.frame(Means=rowSds(ref.mat) )
  
  tgt.mat.z= sweep(tgt.mat,1,t(mean_ref), FUN = "-")
  tgt.mat.zz= sweep(tgt.mat.z,1,t(sd_ref), FUN = "/")
  
  return(tgt.mat.zz)
}






##################################################  
### Gene Symbol conversion - Ensemble_id & HUGO
# *** NOTE Ensembl.version = "85", in gene_ens_map is GRCh38 (hg38) reference genome 
#   load(file.path(Cancer_Reference.folder,"gene_ens_map.Add.desc.RData")) # 57820 genes {Ensemble_id.no.version Ensemble_id HUGO Desc Description}
#   load(file.path(Cancer_Reference.folder,"Ensemble_id.HUGO.df.RData")) # 52576 genes with a version {HUGO Ensemble_id}
#   load(file.path(Cancer_Reference.folder, paste0("gene_ens_map.v",Ensembl.version,".RData")) ) # gene_ens_map,  58069 genes

######################### convert.Genes.Expression.Matrix.between.Ensemble_id.and.HUGO
convert.Genes.Expression.Matrix.between.Ensemble_id.and.HUGO <- function(df , agg.fun.to.use = "max", gene_ens_map.Obj.name = "gene_ens_map",
                        Ensembl.version = "85", convert.from = "Ensemble_id", convert.to = "HUGO"){
  # df with Genes as row names
  # agg.fun.to.use = {median max}
  # gene_ens_map.Obj.name = {gene_ens_map, gene_ens_map.Add.desc, Ensemble_id.HUGO.df}
  # convert.from / convert.tog  = {Ensemble_id Ensemble_id.no.version HUGO}
  #convert.to = "Ensemble_id.no.version"
  #convert.from = "HUGO"
  # *** is.ignore.Ensemble_id.version remove - not needed
  
  genes.with.no.match = setdiff( row.names(df), eval(parse(text=gene_ens_map.Obj.name))[, convert.from] )
  if(length(genes.with.no.match > 0)){
    warning(paste0("genes.with.no.match ", length(genes.with.no.match), " ", paste0(genes.with.no.match, collapse = ",")  ))
  }

  df.gene.names.map = merge(df, eval(parse(text=gene_ens_map.Obj.name)) [,c(convert.from, convert.to)], by.x = "row.names", by.y = convert.from ,
                            sort = F, all = F) # all.x = T, all.y = F
  
  genes.with.non.unique.match.idx = which( duplicated(df.gene.names.map[,convert.to]) )
  non.unique = unique(df.gene.names.map[,convert.to][genes.with.non.unique.match.idx])
  non.unique.from = unique(df.gene.names.map$Row.names[genes.with.non.unique.match.idx])
  #df.gene.names.map[ df.gene.names.map[,convert.to] %in% non.unique  ,c("Row.names", convert.to, names(df.gene.names.map)[1:10]) ]
  df.gene.names.map = df.gene.names.map[,-grep("Row.names", names(df.gene.names.map))]
  
  if(length(genes.with.non.unique.match.idx) > 0){
    warning(paste("genes.with.non.unique.match", length( non.unique ), paste0(non.unique, collapse = ",") ))
    df.gene.names.map.dedup <- aggregate.duplicated.row.names(df.gene.names.map, dulicated.names.cols=convert.to, agg.fun.to.use = agg.fun.to.use) # was median
  }
  row.names(df.gene.names.map.dedup) = df.gene.names.map.dedup[,convert.to]
  df.gene.names.map.dedup = df.gene.names.map.dedup[,!names(df.gene.names.map.dedup) %in% convert.to ]
  
  return(df.gene.names.map.dedup)
}

######################### read.Ensemble_id.return.HUGO
read.Ensemble_id.return.HUGO <- function(df , agg.fun.to.use = "max", gene_ens_map.Obj.name = "gene_ens_map",
                                         Ensembl.version = "85"){
  # Call convert.Genes.Expression.Matrix.between.Ensemble_id.and.HUGO with:
  convert.from = "Ensemble_id"
  convert.to = "HUGO"
  return( convert.Genes.Expression.Matrix.between.Ensemble_id.and.HUGO(df , convert.from = convert.from, convert.to = convert.to,
                  agg.fun.to.use = agg.fun.to.use, gene_ens_map.Obj.name = gene_ens_map.Obj.name,
                  Ensembl.version = Ensembl.version) )
  
}
######################### read.HUGO.return.Ensemble_id
read.HUGO.return.Ensemble_id <- function(df , agg.fun.to.use = "max", gene_ens_map.Obj.name = "gene_ens_map",
                                         Ensembl.version = "85"){
  # Call convert.Genes.Expression.Matrix.between.Ensemble_id.and.HUGO with:
  convert.from = "HUGO"
  convert.to = "Ensemble_id.no.version"
  return( convert.Genes.Expression.Matrix.between.Ensemble_id.and.HUGO(df , convert.from = convert.from, convert.to = convert.to,
                                                                       agg.fun.to.use = agg.fun.to.use, gene_ens_map.Obj.name = gene_ens_map.Obj.name,
                                                                       Ensembl.version = Ensembl.version) )  
}
## Sanity check
# df = POG570.cohort.list[["TPM"]]
# df.Ens = read.HUGO.return.Ensemble_id ( df )
# df.re = read.Ensemble_id.return.HUGO ( df.Ens )
# all_equal(df, df.re) # TRUE

########### get.Ensemble_id.return.HUGO: synonym to read.Ensemble_id.return.HUGO
get.Ensemble_id.return.HUGO <- function(df) {
  # synonym to read.Ensemble_id.return.HUGO
  # df with Ensemble_id as row names
  return( read.Ensemble_id.return.HUGO(df) ) 
}

######################### return.Ensemble_id.no.version
return.Ensemble_id.no.version <- function(Ensemble_ids){
  return( gsub("\\..+","",Ensemble_ids )  )
}


################################## get.HUGO.return.Ensemble_id
get.HUGO.return.Ensemble_id <- function(df, agg.fun.to.use = "max", gene_ens_map.Obj.name = "gene_ens_map.Add.desc") {
  # *** Old version - try using read.HUGO.return.Ensemble_id
  # agg.fun.to.use = {median max}
  # gene_ens_map.Obj.name = {gene_ens_map, gene_ens_map.Add.desc, Ensemble_id.HUGO.df}  
  
  if(!exists(gene_ens_map.Obj.name)) load(file.path(Cancer_Reference.folder, paste0(gene_ens_map.Obj.name,".RData")))
  
  print(dim(df))
  if(is.ignore.Ensemble_id.version){
    df = merge(df, unique(gene_ens_map.Add.desc[,qw("Ensemble_id.no.version HUGO")] ), by.x = "row.names", by.y = "HUGO", all.x = T, all.y = F, sort = F)
  }else{
    df = merge(df, unique(gene_ens_map.Add.desc[,qw("Ensemble_id HUGO")] ), by.x = "row.names", by.y = "HUGO", all.x = T, all.y = F, sort = F)
  }
  df <- aggregate.duplicated.row.names(df, dulicated.names.cols="Row.names", agg.fun.to.use = agg.fun.to.use)
  
  genes.with.NA.Ensemble_id = is.na( df$Ensemble_id ) # ***
  df = df[! is.na( df$Ensemble_id ),]
  
  row.names(df) = df$Ensemble_id
  df = df[,!grepl("Ensemble_id",colnames(df))]
  df = df[,!grepl("Row.names",colnames(df))]
  print(dim(df))
  return(df)
}



################################## read.gct.aggregate.by.HUGO.return.df
read.gct.aggregate.by.HUGO.return.df <- function(input_file_gct , agg.fun.to.use = "max") {
  mat.obj <- MSIG.Gct2Frame(filename = input_file_gct)
  df <- data.matrix(mat.obj$ds)
  df = data.frame(df, "Gene"=mat.obj$descs)   # Get HUGO
  print(paste("original dim", paste0(dim(df), collapse = ", ") ))
  
  df <- aggregate.duplicated.row.names(df, dulicated.names.cols="Gene", agg.fun.to.use = agg.fun.to.use) # was median
  print(paste("After aggregate.duplicated.row.names by HUGO", paste0(dim(df), collapse = ", ") ))
  
  row.names(df) = df$Gene
  df = df[,!grepl("Gene",colnames(df))]
  return(df)
}


####################################################################  
################################## read.gct.aggregate.by.gene.add.desc
read.gct.aggregate.by.gene.add.desc <- function(input_file_rpkm_gct, outfile=NULL, agg.fun.to.use = "max",
                                                is.remove.source.from.Description = T)
{
  
  # WARN using annotEnsembl63 (old) *** use read.gct.aggregate.by.HUGO.return.df()
  if(is.null(outfile))   outfile = gsub(".gct", ".agg.add.desc.gct", input_file_rpkm_gct)

  dataset.2 <- MSIG.Gct2Frame(filename = input_file_rpkm_gct)
  m2.breast.matrix <- data.matrix(dataset.2$ds)
  m2.breast = data.frame(m2.breast.matrix, "Gene"=dataset.2$descs)
  print(dim(m2.breast))
  
  # merge duplicate rows, Ave.
  m2_agg_duplicates <- aggregate.duplicated.row.names(m2.breast, dulicated.names.cols="Gene", agg.fun.to.use = agg.fun.to.use )
  
  #row.names(m2_agg_duplicates) = m2_agg_duplicates$Gene
  #m2_agg_duplicates = m2_agg_duplicates[,-grep("Gene",colnames(m2_agg_duplicates))]
  
  ### Add full gene's descriptions
  if(!exists("annotEnsembl63")) load(file.path(Cancer_Reference.folder, paste0("annotEnsembl63",".RData")))
  annotEnsembl63 = annotEnsembl63[,c("Symbol","Description")]
  annotEnsembl63 = unique(annotEnsembl63) # from ~50K to ~44K
  if(is.remove.source.from.Description) annotEnsembl63$Description = gsub("\\[.+\\]","", annotEnsembl63$Description) 
  annotEnsembl63 = annotEnsembl63[order(annotEnsembl63$Description,decreasing = T),] # order to get "" values down
  
  annotEnsembl63 = annotEnsembl63[-which(duplicated(annotEnsembl63$Symbol)),] # remove duplicated Symbols, needed for merge
  row.names(annotEnsembl63) = annotEnsembl63$Symbol
  
  deGenes <- merge(m2_agg_duplicates, annotEnsembl63, by.x="Gene", by.y="Symbol",all.x = T, all.y= F, sort=FALSE)  
  deGenes = reoder_cols_x_after_y(deGenes, x="Description", "Gene")
  
  row.names(deGenes) = deGenes$Gene # Reqquired for gct
  
  write.gct.2(gct.data.frame = deGenes[,-c(1,2)] , descs = deGenes$Description,
              filename=outfile )
  
}




################### add.purity.from.RNA
add.purity.from.RNA <- function(df, Cancer.cells.col = "Malignant.ER",
                                main.TME.cels.cols = qw("CAFs T.cells Macrophages Endothelial B.cells"), purity.from.RNA.new.col = "purity.from.RNA",
                                Expression.Mat = NULL,df.ID = "Sample_ID",merge.with.TPM.df.ID="Sample_ID",
                                is.print.summary = T, is.range01 = T){
  
  used.cols = union(main.TME.cels.cols,Cancer.cells.col )
  
  if(!is.null(Expression.Mat) & length( setdiff(used.cols, colnames(df)) )>0 ){
    tmp.Signatures.AUC =  calc.Signatures.AUC(Expression.Mat, GSDB.list.of.lists[["cell.types.main"]][ setdiff(used.cols, colnames(df)) ],  
                                              top.n.percent.of.genes.in.ranking.used = 20 , nCores =4, is.new.version = T)
    tmp.Signatures.AUC.df = data.frame( t(tmp.Signatures.AUC))

    df <- merge(df, tmp.Signatures.AUC.df, by.x = merge.with.TPM.df.ID, by.y="row.names", all.x=T, all.y=F, sort=F)
    row.names(df) = df[,df.ID]
    df = df[, grep("row.names", names(df), invert = T, ignore.case = T )]    
  }

  df[,purity.from.RNA.new.col] = df[,Cancer.cells.col] / rowSums(df[, used.cols ])
  if(is.range01) df[,purity.from.RNA.new.col] = range01(df[,purity.from.RNA.new.col])
  if(is.print.summary) print( summary(df[,purity.from.RNA.new.col]))
  return(df)
}


################## add. Gene.and.Signature.Expression.to.Samples.table
add.Gene.and.Signature.Expression.to.Samples.table <- function (mat.expr, Samples.annotation, 
                                                                Genes =NULL, Gene.sets = NULL, GSDB.list=GSDB.list.inclusive, sample_ID.col="Sample_ID"){
  stopifnot( length(which(duplicated(Samples.annotation[,sample_ID.col])))==0 )
  
  if(!is.null(Gene.sets)){
    Signatures.AUC =  calc.Signatures.AUC(mat.expr, GSDB.list=GSDB.list[Gene.sets] )
    Signatures.AUC.df = data.frame( t(Signatures.AUC))
    Samples.annotation = my.merge.new.cols(df1=Samples.annotation, df2=Signatures.AUC.df)
  }
  if(!is.null(Genes)){
    Genes.df = data.frame( t(mat.expr[Genes,]))
    Samples.annotation = my.merge.new.cols(df1=Samples.annotation, df2=Genes.df)
  }
  return(Samples.annotation)
}




################## add. Feature.annotation.to.Samples.table
add.Feature.annotation.to.Samples.table <- function (Samples.annotation, Samples.annotation.sample.col="Sample_ID",
                                                     Reference.annotation.df, Reference.annotation.sample.col="Sample_ID",
                                                     annotate.keep.by.cols.vals.list,
                                                     annotate.filter.by.cols.vals.list=NULL, is.verbos=F ){
  
  sample.IDs.TRUE.annotation = unique( Reference.annotation.df[,Reference.annotation.sample.col] )
  tmp.df = Reference.annotation.df
  
  for(mTest.col in intersect(names(annotate.keep.by.cols.vals.list), names(tmp.df)) ){
    if(is.verbos) print(mTest.col)
    tmp.df <- tmp.df %>% dplyr::filter( !!as.symbol(mTest.col) %in% annotate.keep.by.cols.vals.list[[mTest.col]]  )
    sample.IDs.TRUE.annotation <- intersect(sample.IDs.TRUE.annotation,  tmp.df[, Reference.annotation.sample.col ]   )
    if(is.verbos) print(length(sample.IDs.TRUE.annotation))
  }
  if(!is.null(annotate.filter.by.cols.vals.list)){
    for(mTest.col in intersect(names(annotate.filter.by.cols.vals.list), names(tmp.df)) ){
      if(is.verbos) print(mTest.col)
      tmp.df <- tmp.df %>% dplyr::filter( !!as.symbol( mTest.col) %in% annotate.filter.by.cols.vals.list[[mTest.col]])
      sample.IDs.TRUE.annotation <- setdiff(sample.IDs.TRUE.annotation,  tmp.df[, Reference.annotation.sample.col ]   )
      if(is.verbos) print(length(sample.IDs.TRUE.annotation))
    }
  }
  Samples.annotation[ , annotate.keep.by.cols.vals.list[["make.col.name"]] ] = Samples.annotation[,Samples.annotation.sample.col] %in% sample.IDs.TRUE.annotation
  
  rm(tmp.df)
  return(Samples.annotation)
}







