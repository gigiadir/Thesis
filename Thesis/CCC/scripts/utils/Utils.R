### Utils.R
## General purpose functions

library(data.table)
library(plyr)
library(dplyr)
library(magrittr)
library(ggrepel)
require(limma)
require(DESeq2)
require(edgeR)
library(matrixStats)
library(gplots) # heatmap.2
library(cluster)
library(edgeR)
library(NMF)
#library(tsne)
library(ggplot2)
library(pheatmap)
library(stringi)
library(stringr)
library(clipr)

library(RColorBrewer)


#library(Rcmdr)
#library(mygene)
#require(EDASeq)
#library(tweeDEseqCountData)

################################ perl-like qw()
qw <- function(x) { unlist(strsplit(x, "[[:space:]]+")) }
qwr <- function(x, r=",") { qw( gsub(r,"",x) ) }


#######################  convert.to.numeric.matrix
convert.to.numeric.matrix <- function(dat){
  scols = colnames(dat)
  srows = row.names(dat)
  mat = matrix(as.numeric(unlist(dat)),nrow=nrow(dat))
  colnames(mat) = scols
  row.names(mat) = srows
  return(mat)
}




#######################   Perform clustering heatmap with NAs using heatmap.2
# heatmap.2(as.matrix(cor.Genes),dendrogram="row",trace="none",
#           hclust=  hclustfunc,   distfun=distfunc, col=redgreen(75));
distfunc <- function(x) daisy(x,metric="gower")
hclustfunc <- function(x) hclust(x, method="complete")


#######################  dedup.df
dedup.df <- function(df, value, id, is.decreasing = T){
  df = df[ order(df[,value], decreasing = is.decreasing) , ]
  df = df[ which( !duplicated (df[,id]) ) , ]
  return(df)
}

####################### plot.simple.text
plot.simple.text <- function(text.to.plot, add.table=NULL, cex.table=0.5){
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, text.to.plot,   cex = 1.6, col = "black")
  if(!is.null(add.table))   addtable2plot(cex=cex.table, 0 , 0, add.table, bty="o",display.rownames=TRUE,hlines=TRUE,
                                                vlines=TRUE,title="")   
}

####################### my.merge
my.merge <- function(x, y, by.x, by.y, all.x, all.y, sort  ){
  stopifnot( by.x %in% names(x) | by.y %in% names(y))
  y.2merge = y[,c(by.y, setdiff(names(y), names(x))) ]
  return( merge(x, y.2merge, by.x = by.x, by.y = by.y, all.x=all.x, all.y=all.y,sort=sort) )
}


####################### my.setdiff
my.setdiff <- function(a, b ){
  print("a vs. b")
  print( setdiff(a, b) )
  print("b vs. a")
  print( setdiff(b, a) )  
}

################## print.dim.or.length
print.dim.or.length <- function(mvar){
  mdim = dim(mvar)
  if(is.null(mdim)){
    print(length(mvar))
  }else{
    print(mdim)
  }
}


################## all.exists
all.exists <- function(myvars){
  all.myvars.exists = T
  for(myvar in  myvars  ) all.myvars.exists = all.myvars.exists & exists(myvar)
  return(all.myvars.exists)
}



################################################################
################## loadRData
# loads an RData file, and returns it
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}


# If you're just saving a single object, use an .RDS file:
#saveRDS(x, "x.rds")
#x <- readRDS("x.rds")

####################### my.save.RData
my.save.RData <- function(myvars, myfolder = "~/Dropbox/RData.Robj", 
                          is.save.obj = T, is.print.dim = T, is.print.dim.foreach = F,
                          is.rename.file.keep.old.copy = T, 
                          is.load.obj = F, is.load.only.if.not.exists = F,
                          is.gc = F, time.stamp = "mtime" ){
  
  if(is.load.obj | is.load.only.if.not.exists) {is.save.obj=F; is.rename.file.keep.old.copy=F}
  
  for(myvar in  myvars  ){ 
    print(myvar)
    myvar.date = tail(file.info( file.path(myfolder, paste0(myvar,".RData")) )[, time.stamp]) 
    print(myvar.date)
    file.name = file.path(myfolder, paste0(myvar,".RData"))
    new.file.name.with.date = file.path(myfolder ,paste0(myvar,".", str_split_fixed(myvar.date, " ", 3)[,1]  ,".RData"))
    
    if(is.load.obj) {print("load.Rdata"); miceadds::load.Rdata( file.name , myvar )}
    #if(is.load.obj) load( file.name  ) # not loaded to the Gloval Env
    if(is.load.only.if.not.exists) {
      if(! exists(myvar)){
        print(paste(myvar, "not exists => load.Rdata..."))
        miceadds::load.Rdata( file.name , myvar )
      }else{
        print(paste(myvar, "found"))
      }
    }     
    
    if(is.rename.file.keep.old.copy & file.exists(file.name)) file.rename( file.name,  new.file.name.with.date)
    if(is.save.obj) save( list=myvar, file=file.name )
 
    # if(is.rm.obj){
    #   warning("is.rm.obj not working with Gloval Env")
    #   #assign(myvar, NULL)
    # }
    if(is.print.dim) print.dim.or.length( eval(parse(text=myvar)) )
    if(is.print.dim.foreach) for(mname in names(eval(parse(text=myvar)))){print(mname); print.dim.or.length( eval(parse(text=myvar))[[mname]] )}
    if(is.gc) gc(); # Clear memory
  }
}


####################### my.saveRDS
my.saveRDS <- function(mObject, file.name, is.print.dim = T, is.print.dim.foreach = F,
                          is.rename.file.keep.old.copy = FALSE, time.stamp = "mtime",
                          is.rm.obj = FALSE, is.gc = FALSE , file_extension = ".rds"){
  # time.stamp = {mtime ctime}, Modified or Touched
  
    myvar.date = tail(file.info( file.name )[, time.stamp])
    my.date = str_split_fixed(myvar.date, " ", 3)[,1]
    print(myvar.date)
    new.file.name.with.date = paste0( file_path_sans_ext(file.name),".", my.date , file_extension)
    if(is.rename.file.keep.old.copy) file.rename( file.name,  new.file.name.with.date)
    saveRDS( mObject, file=file.name )

    if(is.print.dim) print.dim.or.length( eval(parse(text=myvar)) )
    if(is.print.dim.foreach) for(mname in names(eval(parse(text=myvar)))){print(mname); print.dim.or.length( eval(parse(text=myvar))[[mname]] )}
    if(is.gc) gc(); # Clear memory
}


####################### find.and.load
## Go over main RData folders (MS.SF.2017.Robj.folder, ...) - if the file is found - load it
## list.of.RData.folders defined as global variable in Args.source.R
find.and.load <- function(myvar,
                          file.extention = ".RData",
                          list.of.RData.folders = "RData.Robj.Misc MS.SF.2017.Robj.folder 
                          Cancer_Reference.folder Gene.sets.Cell.Types.folder  Data.sources.folder  MS.2019.folder
                          HTAPP.fresh.folder  HTAPP.frozen.folder SC.All_QC_10x.folder"){
  
  for(myfolder in qw(list.of.RData.folders)){
    m.path = eval(parse(text=myfolder))
    file.name <- list.files(path = m.path ,pattern = paste0(paste0(myvar,file.extention)) , recursive=F, include.dirs=F)
    if(length(file.name)>0){
      print(paste("load.Rdata from ", myfolder));
      miceadds::load.Rdata( file.path(m.path, file.name) , myvar )
      return( print.dim.or.length( eval(parse(text=myvar)) ) )
    }
  }
}


########################## residual_sig
residual_sig <-function(dat, sig, regress_out_sig){
  m1 <- lm(dat[,sig] ~ dat[,regress_out_sig])  #Create a linear model
  return(resid(m1))
}

########################## my.eval
my.eval <- function(varname){
  return( eval(parse(text =varname )) )
}

########################## my.open.system
my.open.system <- function(file.or.folder.2open){
  system(paste0("open ", file.or.folder.2open ))
}

########################## write.xlsx.append.sheet
write.xlsx.append.sheet <- function(df, mFile, sheetName){
  mFile = gsub("~/",path.expand("~/"), mFile) # Replacing "~/" with "Users//" works for Mac 
  if(file.exists(mFile)) m.append = T else m.append =F
  if(is.numeric(nrow(df)))  xlsx::write.xlsx(df, file=mFile, sheetName=sheetName, append=m.append, row.names=FALSE) # Niice! - write multiple sheets  
}

########################## read_excel_allsheets
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame) # tidyverse tibbles (the default with read_excel)
  names(x) <- sheets
  x
}
####################### my.read.table.by.extension.path
my.read.table.by.extension.path <- function(path,filename, is.all.sheets = T,  sheet = NULL, na = "", agg.fun.to.use = "max"){
  my.read.table.by.extension(file.path(path,filename), 
                             is.all.sheets = is.all.sheets,  sheet = sheet, na = na, agg.fun.to.use = agg.fun.to.use)
}
####################### read.table.by.extension / my.read.table.by.extension
my.read.table.by.extension <- function(filename, is.all.sheets = T,  sheet = NULL, na = "", agg.fun.to.use = "max"){
  # Invoke functions to read tables in {xlsx, csv, tsv/txt, gct}
  file_extension = file_ext(filename)
  print(paste("Read file=",filename," by file_extension ", file_extension))
  
  if(file_extension == "xlsx"){
    if(is.all.sheets){
      df = read_excel_allsheets( filename )
      if(length(df)==1) df = data.frame(df[[1]])
    }else{
      df = data.frame( read_xlsx( filename , sheet = sheet, na = na) )
    }
  }else if(file_extension == "csv"){
    df = read.csv( filename )
  }else if(file_extension %in% qw("tsv txt")){
    df = read.delim( filename )
  }else if(file_extension %in% qw("gct")){
    warning(paste("Read gct file assumes RNA.seq matrix. for HUGO - agg.fun.to.use=", agg.fun.to.use))
    df = read.gct.aggregate.by.HUGO.return.df ( filename, agg.fun.to.use = "max" )     
  }else{
    warning(paste0("file_extension ", file_extension, " not a valid option. Attempt Tab delim"))
    df = read.delim( filename ) 
  }
  return(df)
}
####################### 
read.table.by.extension <- function(filename, is.all.sheets = T,  sheet = NULL, na = "", agg.fun.to.use = "max"){
  my.read.table.by.extension(filename, is.all.sheets ,  sheet , na , agg.fun.to.use)
}


###################### convert.list.to.df
convert.list.to.df <- function(m.list, by.col = T){
  df = data.frame(matrix())
  for(mItem in names(m.list)){
    n.df = data.frame(mItem = m.list[[mItem]]); names(n.df) = mItem;
    df = rowr::cbind.fill(df, n.df ,fill="")
  }
  df = df[,-1]
  
  if(!by.col) df = t(df)
  return(df)
}

###################### 
col.rename <- function(df, old.names, new.names){
  names(df)[grep_multi_patterns_index(qw(old.names), names(df) )] <- qw(new.names)
  return(df)
}

####################### my.merge.new.cols
my.merge.new.cols <- function(df1, df2, row.names.col=NULL, merge.by = "row.names",all.x=T,all.y=F,sort=F){
  # Add cols in df1, Assume row.names are unique IDs, merge.by = "row.names"
  new.cols = setdiff(names(df2), names(df1))
  if(!is.null(row.names.col)) row.names(df1) = df1[,row.names.col]
  if(merge.by == "row.names"){
    df <- merge(df1, df2[,new.cols], by=merge.by, all.x=all.x, all.y=all.y, sort=sort)
  }else{
    df <- merge(df1, df2[,c(new.cols,merge.by)], by=merge.by, all.x=all.x, all.y=all.y, sort=sort)
  }
  if(!is.null(row.names.col)) row.names(df) = df[,row.names.col]
  if(merge.by == "row.names") df = df[, grep("Row.names", names(df), invert = T )]
  #row.names(df) = df[,row.names.col]
  return(df)
}


####################### my.write.table
#  col.names=NA => a blank column name is added, which is the convention used for CSV files to be read by spreadsheets
my.write.table <- function(df, file, sep="\t", row.names = F, quote = F,  col.names=T){
  if(row.names == T) col.names=NA;
  write.table(df, file ,sep=sep, row.names = row.names, quote = quote, col.names=col.names)
}

#######################  cat with new line
catnl <- function(x) cat("\n",paste0(x, "\n"),"\n")

#######################  my.order
my.order <- function(X, decreasing = F){
  X[order(X, decreasing = decreasing)]
}
  


#######################  table.sorted
table.sorted <- function(vector, decreasing = T, is.return.df = F){
  # table + sort
  my.sorted.table = table(vector)
  my.sorted.table = my.sorted.table[order(my.sorted.table, decreasing = decreasing)]
  if(is.return.df){
    return(data.frame(my.sorted.table))
  }else{
    return (my.sorted.table)
  }
}
my.table.sorted <- function(vector, decreasing = T,is.return.df = F){
  return ( table.sorted(vector, decreasing, is.return.df) )
}
my.table <- function(vector, decreasing = T,is.return.df = F){
  return ( table.sorted(vector, decreasing, is.return.df) )
}
length.unique <- function(X){
  return ( length(unique(X)) )
}


#######################  my.cbind.intersect.rownames
my.cbind.intersect.rownames <- function(df1, df2){
  # useful to merge Expression Matrix (keep joint Genes)
  rownames.keep = intersect(row.names(df1), row.names(df2))
  return(cbind(df1[rownames.keep,], df2[rownames.keep,]))
}

#######################  my.write_clip
my.write_clip <- function(dat) write_clip( data.frame(dat) )



#######################  Dates
library(lubridate) # year
myyear <- function(x, year=1968){
  m <- year(x) %% 100
  year(x) <- ifelse(m > year %% 100, 1900+m, 2000+m)
  x
}
foo <- function(x, year=1968){
  return(myyear(x,year))
}
### my.Date {%m/%d/%Y, %Y-%m-%d}
my.Date <- function(x, format="%m/%d/%Y"){
  myyear(as.Date(x, format))
}




#######################  get.myvars
# e.g.,  with Expression Matrix - return Genes vector ranked by the var
get.myvars <- function(mat, row.or.col = 1, na.rm=T){
  myvars <- apply( mat  , row.or.col , var , na.rm=na.rm)
  myvars = myvars[order(myvars, decreasing = T)]
  return( myvars )
}

#######################  get.last.n
# get last n items in vector
get.last.n <- function(x, n){
  return( x[(length(x)-n) : length(x)]  )
}

#######################  substrRight
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}



################################################################
#######################  fisher.test.vec.input
fisher.test.vec.input <- function(AB, AnB, BnA, nAnB , alt="two.sided"){
  mat <- matrix(as.numeric(c( AB,
                              AnB,
                              BnA,
                              nAnB  )), ncol=2)
  f <- fisher.test(as.table(mat), alt=alt)
  return(f)
}
#######################  fisher.test.vec.set.input
fisher.test.vec.set.input <- function(AB, AnB, BnA, nAnB , alt="two.sided", is.print.mat = T){
  mat <- matrix(as.numeric(c( length(AB),
                              length(AnB),
                              length(BnA),
                              length(nAnB)  )), ncol=2)
  if(is.print.mat) print(mat)
  f <- fisher.test(as.table(mat), alt=alt)
  return(f)
}


########################## prep.fgsea.Ranks
prep.fgsea.Ranks <- function(my_table, sort.by = "logFC", name.by = "Gene"){
  my_table = my_table[order(my_table[,sort.by], decreasing = T),]
  fgsea.Ranks = my_table[,sort.by]
  names(fgsea.Ranks) = my_table[,name.by]
  return(fgsea.Ranks)
}

#######################  fisher.test.Gene.Sets.from.fgsea.Ranks
fisher.test.Gene.Sets.from.fgsea.Ranks <- function(fgsea.Ranks, mysigName = "", pathways, U=NULL, mfile = NULL, ntop.genes = 200, q.value.write = 0.25, is.intersect.set.with.U = T, 
                                                   is.write.row.names = T, is.already.names = F, is.return.fisher.test.Gene.Sets.df = F){
  
  if(is.already.names) names(fgsea.Ranks) = fgsea.Ranks
  if(is.null(U)) U = names(fgsea.Ranks)
  
  A = names(fgsea.Ranks[1:ntop.genes])
  if(is.intersect.set.with.U) A = intersect(A, U) # not needed if U is based on fgsea.Ranks
  
  fisher.test.Gene.Sets.df = data.frame()
  for(mySigMain in names(pathways)){
    B = pathways[[mySigMain]]
    if(is.intersect.set.with.U) B = intersect(B, U)
    key = paste(mysigName,".vs.", mySigMain,".top",ntop.genes)
    fisher.test.Gene.Sets.df[key, "mySigMain" ] = mySigMain
    
    mU = unique(c(A,B,U)) # not needed if is.intersect.set.with.U, otherwise - change the universe
    tmp = fisher.test.groups.input(A, B, mU, is.return.list=T, is.print = F)
    fisher.test.Gene.Sets.df[key, "estimate" ] = signif( tmp$estimate, digits = 3)  
    fisher.test.Gene.Sets.df[key, "p.value" ] =  signif( tmp$p.value , digits = 3)  
    fisher.test.Gene.Sets.df[key, "mat" ] = paste0( as.vector(tmp$mat) , collapse = ", ") 
    fisher.test.Gene.Sets.df[key, "AB" ] = paste0(tmp$AB, collapse = ", ")
    fisher.test.Gene.Sets.df[key, "n.AB" ] = length(tmp$AB)
  }
  fisher.test.Gene.Sets.df = fisher.test.Gene.Sets.df[order(fisher.test.Gene.Sets.df$p.value),]
  fisher.test.Gene.Sets.df$q.value = p.adjust( fisher.test.Gene.Sets.df$p.value, method = "BH")
  mSummary.estimate = summary(fisher.test.Gene.Sets.df$estimate)
  print(mSummary.estimate)
  fisher.test.Gene.Sets.df = subset(fisher.test.Gene.Sets.df, q.value <= q.value.write)
  
  if(!is.null(mfile)){
    write.csv(fisher.test.Gene.Sets.df, mfile, row.names = is.write.row.names)
  }else if(!is.return.fisher.test.Gene.Sets.df){
    write_clip(fisher.test.Gene.Sets.df)
    return(mSummary.estimate)
  }else{
    return(fisher.test.Gene.Sets.df)
  }
}

################################ fisher.test.Gene.Sets.from.fgsea.Ranks
gsea.test.Gene.Sets.from.fgsea.Ranks <- function(fgsea.Ranks, 
                                       genes.sets.2test = "msigdb.c2.augmented.SABCS19",
                                       mfile = NULL,
                                       nperm = 10000, pval.return = 0.1 
                                       ){
 
  dd.writeable = runner.fgsea(fgsea.Ranks,
                                pathways = genes.sets.2test,
                                nperm=nperm,
                                pval.return=pval.return)
    
  if(!is.null(mfile)){
    write.csv(dd.writeable, mfile  )
  }else{
    write_clip(dd.writeable)
  }
}





#######################  fisher.test.vec.set.raw.input
fisher.test.vec.set.raw.input <- function(A, nA, B, nB , alt="two.sided", is.print.sets = T, is.print.mat = T){
  AB = intersect(A,B)
  AnB = intersect(A,nB)
  BnA = intersect(nA,B)
  nAnB = intersect(nA,nB)
  
  mat <- matrix(as.numeric(c( length(AB),
                              length(AnB),
                              length(BnA),
                              length(nAnB)  )), ncol=2)
  
  if(is.print.mat) print(mat)
  if(is.print.sets){
    print("AB"); print(AB)
    print("AnB"); print(AnB)    
  }

  f <- fisher.test(as.table(mat), alt=alt)
  return(f)
}

#######################  fisher.test.groups.input
fisher.test.groups.input <- function(A, B, U , alt="two.sided", is.return.list = F, is.print = T, is.print.nonoverlap = F){
  AB = intersect(A,B)
  BnA = base::setdiff(B, A)
  AnB = base::setdiff(A, B)
  nAnB = base::setdiff(U,c(B,A))
  
  mat <- matrix(as.numeric(c( length(AB),
                              length(AnB),
                              length(BnA),
                              length(nAnB) )), ncol=2)
  f <- fisher.test(as.table(mat), alt=alt)
  if(is.print){
    print(mat)
    print("AB"); print(AB)
  }
  if(is.print.nonoverlap){ print("AnB"); print(AnB); print("BnA"); print(BnA);}
  


  if(is.return.list){
    return(list(f=f,mat=mat,  AB=AB, BnA=BnA, AnB=AnB, estimate = f$estimate, p.value=f$p.value))
  }else{
    return(f)
  }
}



#######################  inverse.minus.log (power)
inverse.minus.log <- function(x, base = 10){
  # -log10 (10 ^ - as.numeric(DE.set.2.write$minus.log.q.val), inverse function for -log10
  return(base ^ (-as.numeric( x )) )
}


#######################  setdiff.twosides
setdiff.twosides <- function(x, y, is.print.intersection = T){
  setdiff.twosides.list = list()
  if(is.print.intersection) setdiff.twosides.list[["xy"]] = intersect(x,y)
  setdiff.twosides.list[["x.not.y"]] = setdiff(x,y)
  setdiff.twosides.list[["y.not.x"]] = setdiff(y,x)
  return(setdiff.twosides.list)
}



######################### add.GO.annotation.to.df
add.GO.annotation.to.df <- function(df, HUGO.col.name){
  Genes.to.annotate = df[,HUGO.col.name]
  res <- mygene::queryMany(Genes.to.annotate, scopes='symbol', fields=c("symbol",'entrezgene', 'go'), species='human',return.as="DataFrame")
  res = subset(res, !is.na(symbol))
  res = unique(res)
  dim(res)
  for(go.type in qw("go.BP go.CC go.MF")){ #  BP = Biological Process, MF = Molecular Function,  CC = Cellular Component.
    print(go.type)
    res [,paste0(go.type,".","id")] = llply(1:nrow(res)  ,function(x) {paste(res[x, go.type][[1]]$id, collapse = ", ")   } ) %>%  unlist()    
    res [,paste0(go.type,".","term")] = llply(1:nrow(res)  ,function(x) {paste(res[x, go.type][[1]]$term, collapse = ", ")   } ) %>%  unlist()    
  }
  res = res[order(res$go.BP.id, decreasing = T),]
  res = res[!duplicated(res$symbol),]
  row.names(res) = res$symbol
  res.final = res[,qw("symbol go.BP.id go.BP.term go.CC.id go.CC.term go.MF.id go.MF.term")]
  
  df.merged = merge(df, res.final, by.x = HUGO.col.name, by.y = "symbol", all.x = T, all.y = F, sort = F )
  print(dim(df.merged))
  return( data.frame( df.merged) )
}


######################### replace.col.name.in.data.frame
replace.col.name.in.data.frame <- function(df, old.name, new.name){
  names(df)[names(df) == old.name] <- new.name
  return(df)
}






########## reorder_cols_movetolast
reorder_cols_movetolast <- function(data, move) {
  data[c(setdiff(names(data), move), move)]
}


########## add col x after col y
reoder_cols_x_after_y <- function(df, x, y=NULL){
  col_idx <- which(colnames(df)==x)
  if(is.null(y)){
    keep_order_idx = 0
    col_idx_after = ((keep_order_idx+1):ncol(df))[ ! ((keep_order_idx+1):ncol(df)) %in%   col_idx  ]    
    df <- df[, c(col_idx, col_idx_after) ]
  }else{
    keep_order_idx = which(colnames(df)==y)
    col_idx_after = ((keep_order_idx+1):ncol(df))[ ! ((keep_order_idx+1):ncol(df)) %in%   col_idx  ]
    df <- df[, c(1:keep_order_idx , col_idx, col_idx_after) ]
  }
  return(df)
}
reorder_cols_x_after_y <- function(df, x, y=NULL){
  return(reoder_cols_x_after_y(df,x,y))
}


################################################################
#######################  gg.point.XY.plot
gg.point.XY.plot <- function(dat, X, Y,  lab,
                             point.size=NULL, 
                             point.size.number = NULL,
                             point.shape=NULL,
                             point.color =NULL,
                             label.size = 2,
                             is.print.label=T, is.repel = T, max.overlaps=50,
                             subset.lables = NULL, label.color = NULL, # NULL?
                             is.scale.color = F, midpoint.col =NULL,
                             is.trend = F,trend.line.method = "glm",is.rug=F,
                             filename = NULL, is.auto.filename=F,
                             scale=1, width=10, height=8,
                             panel.grid.element_blank=T,
                             scale_color_manual.values=NULL, use.distinctColorPalette=F,
                             title="", subtitle=""){
  


  myaes <- aes_string(x = X, y = Y)
  myaes.full <- aes_string(x = X, y = Y, label = lab)
  
  if(!is.null(point.size.number)) myaes$size <- point.size.number
  if(!is.null(point.size)) myaes$size <- as.name(point.size)  # as category to use
  if(!is.null(point.color)) myaes$colour <-  as.name(point.color)
  if(!is.null(point.shape)) myaes$shape <- as.name(point.shape)

  p = ggplot(data = dat) + geom_point(myaes)
  
  if(is.rug) p = p + geom_rug(myaes)

  if(panel.grid.element_blank){
    p = p +  theme_bw() + theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(),
                                                panel.background = element_blank()) #+ geom_point(size = point.size)
  }else{
    p = p + theme_bw() #+ geom_point(size = point.size)
  }
  
  
  if(!is.null(scale_color_manual.values)){
    print(scale_color_manual.values)
    p = p + scale_color_manual(values=scale_color_manual.values)
  }else if(use.distinctColorPalette){
    my_color_palette <- distinctColorPalette( length(unique(data[,point.color])) )
    p = p + scale_color_manual(values=my_color_palette)
  }else if(is.scale.color){
    if(is.null(midpoint.col)) midpoint.col= my.median( dat[,point.color] )
    p = p + scale_color_gradient2(low = "darkcyan", midpoint = midpoint.col, mid = "grey", high = "red4") #+ theme_bw()
  }
  
  if(!is.null(subset.lables)){
    dat.for.label = subset(dat, dat[,label] %in% subset.lables)
  }else{
    dat.for.label = dat
  }
  if(is.trend){
    p = p + geom_smooth(myaes, se = T, method = trend.line.method) # aes_string(),aes_string("RB.relative.to.cycle"),
    m.cor = cor.test(dat[,X], dat[,Y], use = "complete.obs")
    subtitle= paste0(subtitle, " R=" , signif(m.cor$estimate, digits=3), ", p=", signif(m.cor$p.value, digits=3)  )
  }
  if(is.print.label){
    myaes = myaes.full
    if(is.repel){
      p = p + geom_text_repel(data=dat.for.label, myaes, size=label.size, max.overlaps=max.overlaps, color="black") # , color = label.color
    }else{
      p = p + geom_text(data=dat.for.label, myaes, size=label.size, max.overlaps=max.overlaps, color="black") # ,color =label.color
    }
  }
  
  
  p = p + labs(title=title, subtitle= subtitle)
  print(p)
  
  if(is.auto.filename) filename = paste("Figures/", title, Y, " vs ", X,".pdf")
  if(!is.null(filename)) ggsave(filename = filename, device="pdf", scale = scale, width = width, height = height)   
  
}




library(scales) 
#######################  gg.histogram
gg.histogram <- function(dat, 
                         hist.by, 
                         color.by = NULL, color = "blue", fill = "grey50",
                         is.percent = T,
                         is.add.density = T, density.factor = 10,
                         title.size=26, subtitle.size = 24, text.size =18, axis.text.size=22,axis.title.size=22,
                         title = "",subtitle="", y = NULL, x=NULL ) {
  # aes(y=..count../sum(..count..)), aes(y=..density..)
  if(!is.null(color.by) ){
    g <- ggplot(dat, aes_string(x=hist.by,  color=color.by, fill=color.by)) + geom_histogram(aes(y=..count../sum(..count..)), position="dodge", bins=nrow(dat)/3)
  }else{
    g <- ggplot(dat, aes_string(x=hist.by)) + geom_histogram(aes(y=..count../sum(..count..)), color = "blue", fill="white", position="dodge", bins=nrow(dat)/3)
  }

  if(is.add.density) g <- g + geom_density(aes(y=..count../(sum(..count..)/density.factor)), alpha=.2)
  if(is.percent) g <- g + scale_y_continuous(labels=percent)
  g <- g +  labs(title=title, subtitle= subtitle)
  if(!is.null(y)) g <- g +  labs(y=y) # Note! - this call of labs() is not overriding the pervious one - incremental! (works with basics? geom_bar?)
  if(!is.null(x)) g <- g +  labs(x=x)
  g <- g +  theme_classic()
  g <- g + theme(text = element_text(size=text.size), 
                 plot.title=element_text(size=title.size),
                 plot.subtitle=element_text(size=subtitle.size, hjust=0.2, face="italic"),
                 axis.text=element_text(size=axis.text.size),axis.title=element_text(size=axis.title.size,face="bold") )
  return(g)
}

#######################  gg.bar
gg.bar <- function(dat, bar.by,
                   bar.color = "black", bar.fill.color = "lightblue",
                   title = "",   subtitle="", y=NULL,x=NULL, 
                   title.size=26, subtitle.size = 24, text.size =18, axis.text.size=22,axis.title.size=22,
                   add.label.count = F,hjust= -0.05, vjust= 0,
                   is.coord_flip = T, is.theme_minimal = T, 
                   is.order.bars = T, is.max.on.bottom = F) {
  
  if(is.order.bars){
    x.order = table.sorted(dat[,bar.by], decreasing = is.max.on.bottom)
    dat[,bar.by]<-factor(dat[,bar.by], levels=names(x.order) )    
  }

  g <- ggplot(dat, aes_string(x=bar.by)) # , y=PATIENT_ID
  g <- g + geom_bar( color=bar.color, fill=bar.fill.color) # stat="identity", only if already counted.
  if(add.label.count) g <- g + geom_text(stat='count', aes(label=..count..), hjust=hjust, vjust=vjust)
  g <- g +  labs(title=title, subtitle= subtitle )
  if(!is.null(y)) g <- g +  labs(y=y) # Note! - this call of labs() is not overriding the pervious one - incremental! (works with basics? geom_bar?)
  if(!is.null(x)) g <- g +  labs(x=x)
  
  if(is.theme_minimal) g <- g + theme_minimal()
  if(is.coord_flip) g <- g + coord_flip()
  
  g <- g + theme(text = element_text(size=text.size), 
                 plot.title=element_text(size=title.size),
                  plot.subtitle=element_text(size=subtitle.size, hjust=0.2, face="italic"),
                  axis.text=element_text(size=axis.text.size),axis.title=element_text(size=axis.title.size,face="bold") )
  
  return(g)
}



#######################  gg.box.plot.with.names
gg.box.plot.with.names <- function(dat, category.x.col=NULL, quantitative.y.col=NULL, lab.col=NULL,
                                   y.values.name="", txt.size = 3, point.size =1, jitter.factor = 1.5, is.repel = F, max.overlaps=30,
                                   is.anova = F, is.kruskal=F, is.levene = F, is.Welch = T,
                                   title = "",  title.size=16, subtitle.size = 12,
                                   ylim = NA, fill=NA,is.print.label=T,
                                   scale_color_manual.values=NULL,colors_palette_list=NULL,
                                   geom_line.group=NULL, notch=F, outlier.shape=19,lwd=0.5,geom_boxplot.linecol="black", 
                                   outlier.size=0,
                                   ylim1=0,ylim2=NULL,
                                   colour = "x",
                                   is.reorder.groups.by.median = F,order.decreasing=T,
                                   panel.grid.element_blank=F, x.angle = NA, hjust=1, vjust=0.2){
  
  if(!is.null(category.x.col))     dat$x = dat[,category.x.col]
  if(!is.null(quantitative.y.col)) dat$y = dat[,quantitative.y.col]
  if(!is.null(lab.col)) dat$lab = dat[,lab.col]
  
  dat$colour = dat[,colour]
  
  ### order by Median of: y.axis.measured, for category.x.col
  x.order = as.character(unique(dat$x))
  if(is.reorder.groups.by.median){
    formula=as.formula(paste0(quantitative.y.col," ~ ", category.x.col))
    x.order.df = aggregate( formula,   my.median, data=dat)
    x.order <- x.order.df[order(x.order.df[,quantitative.y.col], decreasing = order.decreasing), category.x.col  ]
    dat$x<-factor(dat$x, levels=x.order)    
  }
  if(!is.null(colors_palette_list)){
    scale_color_manual.values.df.list = list()
    categories = as.character( unique(dat$x) )
    #scale_color_manual.values = c()
    for(Cols in names(colors_palette_list)){
      #scale_color_manual.values = c(scale_color_manual.values,  colorRampPalette(brewer.pal(7,  Cols ))(length(colors_palette_list[[Cols]])+1)[-1] )
      scale_color_manual.values.df.list[[Cols]] = data.frame(cat=colors_palette_list[[Cols]], color=colorRampPalette(brewer.pal(7,  Cols ))(length(colors_palette_list[[Cols]])+1)[-1])
      categories = setdiff(categories, colors_palette_list[[Cols]])
    }
    if(length(categories)>0) {
      scale_color_manual.values.df.list[["Other"]] =  data.frame(cat=categories, color=brewer.pal(length(categories), "BrBG"))
    }
    scale_color_manual.values.df = as.data.frame(bind_rows(scale_color_manual.values.df.list))
    scale_color_manual.values = scale_color_manual.values.df$color
    names(scale_color_manual.values) = scale_color_manual.values.df$cat
    scale_color_manual.values = scale_color_manual.values[as.character( x.order)]
  }
  
  datJit <- dat
  if("fill" %in% colnames(dat)){
    fill.factor = 1/length(unique(dat$fill))-0.1
    mymin = -fill.factor
    mymax = fill.factor
    val = as.numeric(factor(dat$fill))
    #val.normalized.shifted = (val-min(val))/(max(val)-min(val)) 
    val.normalized.shifted= (mymax-mymin)/(max(val)-min(val))*(val-max(val) )+mymax
    datJit$xj <- jitter(as.numeric(factor(dat$x))+val.normalized.shifted, factor=jitter.factor)
  }else{
    datJit$xj <- jitter(as.numeric(factor(dat$x)), factor=jitter.factor)
  }
  

  


  test.type = ""
  my.p.value = ""
  if("z" %in% colnames(dat)){
    Welch = oneway.test(y~z , data=dat, na.action=na.omit, var.equal=FALSE); test.type = "Welch"
    my.p.value = my.Welch.p.value = my.signif(Welch$p.value)
  }else if(is.anova){
    my.aov <- aov(y ~ x, data = dat);    test.type = "ANOVA"
    my.p.value = my.aov.p.value = my.signif(summary(my.aov)[[1]][["Pr(>F)"]][1])    
  }else if(is.kruskal){
    my.kruskal.test = kruskal.test(dat$y~factor(dat$x))
    my.p.value = my.kruskal.test.pvalue = my.signif(my.kruskal.test$p.value)  
  }else if(is.levene){
    my.levene.test = levene.test(dat$y~factor(dat$x))
    my.p.value = my.levene.p.value = my.signif((my.levene.test)[["Pr(>F)"]][1])
  }else if(is.Welch){
    Welch = oneway.test(y~x, data=dat, na.action=na.omit, var.equal=FALSE); test.type = "Welch"
    my.p.value = my.Welch.p.value = my.signif(Welch$p.value)     
  }
  

  pp = ggplot(dat,aes(x=x,y=y))
  if("fill" %in% colnames(dat)) pp = ggplot(dat,aes(x=x,y=y,fill=fill))
  
  pp = pp + ylab(y.values.name) +
    geom_boxplot(notch = notch, outlier.shape = outlier.shape,lwd=lwd,colour = geom_boxplot.linecol, outlier.size=outlier.size) # geom_violin
  if(!is.null(ylim2)) pp = pp + coord_cartesian(ylim=c(ylim1,ylim2))

  if("fill" %in% colnames(dat)){
    if("shape" %in% colnames(dat)){
      pp = pp + geom_point(data=datJit,aes(x=xj, colour=colour, fill=fill, shape=shape ),size=point.size)
    }else{
      pp = pp + geom_point(data=datJit,aes(x=xj, colour=colour, fill=fill),size=point.size)
    }
  }else{
    if("shape" %in% colnames(dat)){
      pp = pp + geom_point(data=datJit,aes(x=xj, colour=colour, shape=shape ),size=point.size)
    }else{
      pp = pp + geom_point(data=datJit,aes(x=xj, colour=colour),size=point.size)
    }    
  }
  if(!is.null(scale_color_manual.values) ){ 
    pp = pp + scale_color_manual(values=scale_color_manual.values)
  }
  

  pp = pp +  #geom_text(data=datJit,aes(x=xj,label=lab),size=txt.size ) +
    theme(axis.text=element_text(size=12) ,axis.title=element_text(size=12) )
  
  if(is.numeric(my.p.value))   pp = pp +   labs(title=title,
                   subtitle= paste0(test.type, " " ,signif(my.p.value, digits = 3))  ) # x, y, title, subtitle  
  
  if(!is.na(ylim)) pp = pp + ylim(ylim)
  if(is.print.label){
    if(is.repel){
      pp = pp + geom_text_repel(data=datJit,aes(x=xj,label=lab),size=txt.size ,max.overlaps=max.overlaps)
    }else{
      pp = pp + geom_text(data=datJit,aes(x=xj,label=lab),size=txt.size )
  }
  }
  if(!is.null(geom_line.group))   pp <- pp + geom_line(data=datJit, aes_string(x="xj",group = geom_line.group))
  
  pp <- pp + theme_bw()
  if(panel.grid.element_blank){
    pp <- pp  + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      #panel.border = element_blank(),
                      #,panel.grid.minor = element_blank(), panel.background= element_blank()) # axis.ticks = element_blank(), panel.grid.major = element_blank()
                      plot.title=element_text(size=title.size),
                      plot.subtitle=element_text(size=12, hjust=0.2, face="italic") )    
  }else if(!is.na(x.angle)){
    pp <- pp + theme(axis.text.x=element_text(angle=x.angle,hjust=1,vjust=0.2),
                     plot.title=element_text(size=title.size),
                     plot.subtitle=element_text(size=subtitle.size, hjust=0.2, face="italic") ) 
  }else{
    pp <- pp + theme(plot.title=element_text(size=title.size),
                     plot.subtitle=element_text(size=subtitle.size, hjust=0.2, face="italic") )
  }
  
  return(pp)
} 






####################### gg.box.plot.with.names.2facors.fill
# https://cmdlinetips.com/2019/02/how-to-make-grouped-boxplots-with-ggplot2/
# dat %>% 
#   ggplot(aes(x=x, y=y, fill=factor(fill) ) ) + # , colour = (Perturbation)
#   geom_boxplot(outlier.size=0) + 
#   labs(fill = "Perturbation") + 
#   geom_point(position=position_jitterdodge(seed=1),alpha=0.5) +
#   theme_bw(base_size = 10)
gg.box.plot.with.names.2facors.fill <- function(dat, y.values.name, txt.size = 3, point.size =1, jitter.factor = 1.5, 
                                                is.print.label = T, is.repel = F,
                                                scale_color_manual.values=NULL,
                                                is.anova = F, is.kruskal=F, is.levene = F, title = title,  ylim = NA, fill=NA){
  # dat with x==label of groups, y==values to measure, lab==names of datapoints
  pp = ggplot(dat,aes(x = x, y = y, fill=fill)) + 
    geom_boxplot() + 
    geom_point( position=position_jitterdodge(jitter.width=0.1) ) ##+ geom_text(aes(label=lab),size=txt.size )
  
  if(!is.null(scale_color_manual.values)) pp = pp + scale_fill_manual(values=scale_color_manual.values)
  if(!is.na(ylim)) pp = pp + ylim(ylim)
  
  if(is.print.label){
    if(is.repel){
      pp = pp + geom_text_repel(data=datJit,aes(x=xj,label=lab),size=txt.size )
    }else{
      pp = pp + geom_text(data=datJit,aes(x=xj,label=lab),size=txt.size )
    }
  }
  
  pp <- pp + theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black"))
  
  pp <- pp + theme_bw()
  pp
  
}


#######################  Volcano.plot.DE
Volcano.plot.DE <- function( DE.set,
                             Genes.to.plot = NULL,
                             mytitle="",
                             all.genes.quantile.to.add = 0.999,
                             selected.gene.labels = NULL,
                             text.size.specify = NULL,
                             point.size = 1,text.size=4,
                             max.overlaps=10, # getOption("ggrepel.max.overlaps", default = 10),
                             keyplot = "",
                             Significance.by = 2,
                             element_text_size = 12,
                             panel.grid.element_blank = T,
                             Significance.col = "minus.log.p.val",
                             logFC.col = "logFC",
                             Gene.col = "Gene",
                             max.number.of.geom_text_repel.efficient.print = 500,
                             mCell.type="", 
                             DE.Volcano.plot.folder = "~/Dropbox",
                             is.lower.pval.for.label = F,
                             genes.autoadd=F,genes.autoadd.LFC=F,
                             width = 10, height = 8, is.open.plot.folder = F){
  
  
  library(ggrepel)
  #is.lower.pval.for.label = is.lower.pval.for.label | is.gene.of.interest.Signature | is.gene.of.interest.specific | is.gene.of.interest.all.cancer
  
  
  if(grepl("minus.log", Significance.col) ){
    DE.set$Significance = as.numeric( DE.set[,Significance.col])
  }else{
    DE.set$Significance = as.numeric( -log10(DE.set[,Significance.col])  )
  }
  DE.set$minus.log.p.val = DE.set$Significance
  DE.set$logFC = as.numeric( DE.set[,logFC.col])
  DE.set$Gene = DE.set[,Gene.col]
  
  min.logFC.label.low = 1
  min.minus.log.p.val.low.quantile = 0.75 # 0.9
  
  
  DE.set.sub = DE.set # init
  genes.autoadd = intersect(genes.autoadd, subset(DE.set.sub, (minus.log.p.val>10 & abs(logFC)>0.5) | abs(logFC) >= genes.autoadd.LFC)$Gene)
  if(!is.null(Genes.to.plot)){
    DE.set.sub= subset(DE.set.sub, DE.set.sub[,Gene.col] %in% Genes.to.plot )
    point.size = 3 # 2
    text.size = 5 # 4
  }
  if(!is.null(text.size.specify)) text.size = text.size.specify
  
  
  min.logFC.label =   quantile(abs(DE.set.sub$logFC), all.genes.quantile.to.add)
  min.minus.log.p.val =  quantile(DE.set.sub$minus.log.p.val, all.genes.quantile.to.add)
  min.minus.log.p.val.low = min( quantile(DE.set.sub$minus.log.p.val, min.minus.log.p.val.low.quantile), min.minus.log.p.val) # 4
  
  
  DE.set.sub$Significance = DE.set.sub$minus.log.p.val
  
  p = ggplot(DE.set.sub, aes(logFC,  Significance, colour = logFC) ) + geom_point(size=point.size)  
  
  p = p + scale_x_continuous(breaks = round(seq(min(DE.set.sub$logFC), max(DE.set.sub$logFC), by = 1),0))
  p = p + scale_y_continuous(breaks = round(seq(min(DE.set.sub$Significance), max(DE.set.sub$Significance), by = Significance.by),0))
  
  
  if(panel.grid.element_blank){
    p = p + scale_color_gradient2(low = "darkcyan", midpoint = 0, mid = "grey", high = "red4") + theme_bw() + theme(panel.grid.major = element_blank(),
                                                                                                                    panel.grid.minor = element_blank(),
                                                                                                                    panel.background = element_blank(),
                                                                                                                    axis.text.x = element_text(size=element_text_size), axis.text.y = element_text(size=element_text_size)) + 
      geom_point(size = point.size) # was c("blue","violet","red")
  }else{
    p = p + scale_color_gradient2(low = "darkcyan", midpoint = 0, mid = "grey", high = "red4") + theme_bw() + 
      theme( axis.text.x = element_text(size=element_text_size),axis.text.y = element_text(size=element_text_size))
  }
  p = p + ggtitle(mytitle)
  
  
  if(is.lower.pval.for.label){
    dat.2.plot.labels = subset(DE.set.sub, abs(logFC) >= min.logFC.label.low | minus.log.p.val> min.minus.log.p.val.low | Gene %in% genes.autoadd )
    keyplot = paste(keyplot,".p",signif(min.minus.log.p.val.low,digits=3),sep=".") 
  }else{
    dat.2.plot.labels = subset(DE.set.sub, abs(logFC) >= min.logFC.label | minus.log.p.val > min.minus.log.p.val | Gene %in% genes.autoadd)
    keyplot = paste(keyplot,".p",signif(min.minus.log.p.val,digits=3),sep=".")
  }
  
  
  if(!is.null(selected.gene.labels) && length(selected.gene.labels)>0 ){
    dat.2.plot.labels = subset(dat.2.plot.labels, Gene %in% selected.gene.labels)
    keyplot = paste(keyplot,".selected.genes",length(selected.gene.labels),sep=".")
  }
  
  p=p+geom_text_repel(data=dat.2.plot.labels, aes(label=Gene), size = text.size, max.overlaps=max.overlaps) 
  
  
  if(nrow(dat.2.plot.labels) < max.number.of.geom_text_repel.efficient.print){
    p
    ggsave(filename = file.path(DE.Volcano.plot.folder, paste0( "Volcano.",keyplot,".PDF") ), device="pdf",
           scale = 1, width = width, height = height)
  }else{
    warning("Too many labels to print...")
  }
  
  if(is.open.plot.folder) system(paste0("open ", DE.Volcano.plot.folder) )
  
}












#######################  Merge Two Lists
appendList <- function (x, val) 
{
  stopifnot(is.list(x), is.list(val))
  xnames <- names(x)
  for (v in names(val)) {
    x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]])) 
      appendList(x[[v]], val[[v]])
    else c(x[[v]], val[[v]])
  }
  return(x)
}

#######################  geometric mean
gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}


########################## t.test.single.vec
t.test.single.vec <- function(vec, min.vec.length = 2){
  # test if mean is not equal to 0 for vector of numbers
  if(length(vec) >= min.vec.length){
    my.t.test = t.test(vec.logFC)
  }else{
    my.t.test$estimate = NA;   my.t.test$p.value = NA  
  }
  return(my.t.test)
}
#######################  my.t.test
my.t.test <- function(x,y){
  if(length(which(!is.na(x)))>1 & length(which(!is.na(y)))>1){
    return(t.test(x,y))
  }else{
    return("NA")
  }
}

#######################  aggregate.duplicated.row.names - New
# *** agg.fun.to.use = median or max? (changed to max - verify consistency/reproducibility)
# alternative - sort - to get top and remove duplicates (edgeR) - sames as Max!
aggregate.duplicated.row.names <- function(m.2, dulicated.names.cols="hgnc_symbol", agg.fun.to.use = "max", is.test.hugo.symbols = F, is.make.row.names=F){
  # Aim: aggregated df with duplicated HUGO symbols
  
 
  X <- as.data.table(m.2)

  if(agg.fun.to.use == "mean"){
    m2_agg_duplicates <- X[,lapply(.SD,my.mean),dulicated.names.cols]
  }else if(agg.fun.to.use == "median"){
    m2_agg_duplicates <- X[,lapply(.SD,my.median),dulicated.names.cols]
  }else if(agg.fun.to.use == "max"){
    m2_agg_duplicates <- X[,lapply(.SD,max),dulicated.names.cols]
  }else if(agg.fun.to.use == "min"){
    m2_agg_duplicates <- X[,lapply(.SD,min),dulicated.names.cols]
  }else if(agg.fun.to.use == "which.max"){
    m2_agg_duplicates <- X[,lapply(.SD,which.max),dulicated.names.cols]
  }else if(agg.fun.to.use == "which.min"){
    m2_agg_duplicates <- X[,lapply(.SD,which.min),dulicated.names.cols]    
  }else if(agg.fun.to.use == "sum"){
    m2_agg_duplicates <- X[,lapply(.SD,sum),dulicated.names.cols]
  }else if(agg.fun.to.use == "sd"){
    m2_agg_duplicates <- X[,lapply(.SD,sd),dulicated.names.cols]
  }else if(agg.fun.to.use == "prod"){
    m2_agg_duplicates <- X[,lapply(.SD,prod),dulicated.names.cols]    
  }else{
    cat("Error - agg.fun.to.use:",agg.fun.to.use,", was not found")
  }    
  print(dim(m2_agg_duplicates))
  if(is.test.hugo.symbols) row.names(m2_agg_duplicates) = fix_HUGO_excel_symbol(row.names(m2_agg_duplicates))
  m2_agg_duplicates = as.data.frame.matrix( m2_agg_duplicates)
  
  if(is.make.row.names){
    NA.gene.num = which(is.na(m2_agg_duplicates[, 1 ]))
    if(length(NA.gene.num)>0)  m2_agg_duplicates = m2_agg_duplicates[-NA.gene.num, ]
    row.names(m2_agg_duplicates) = m2_agg_duplicates[,1]
    m2_agg_duplicates = m2_agg_duplicates[,-1]
  }
  
  return(  m2_agg_duplicates )
}


#######################  aggregate.duplicated.row.names - old
# aggregate.duplicated.row.names.old <- function(m.2, dulicated.names.col = NULL, agg.fun.to.use = "mean", is.test.hugo.symbols = T){
#   
#   if(is.null( dulicated.names)){
#     dulicated.names = list(row.names(m.2)) # can't assign duplicated row.names, thus might be missing in m.2
#   }else{
#     cat("Field (colname) for the duplicated row naumes is given as input",dulicated.names.col,"\n")
#     dulicated.names = list(m.2[,dulicated.names.col])
#     m.2 = m.2[ ,which(colnames(m.2) != dulicated.names.col ) ] # remove the non-numeric col
#   }
#   
#   print(dim(m.2))
#   if(agg.fun.to.use == "mean"){
#     m2_agg_duplicates = as.data.frame(apply(m.2, 2, function(x) aggregate(  as.numeric(x), dulicated.names, my.mean))) 
#   }else if(agg.fun.to.use == "median"){
#     m2_agg_duplicates = as.data.frame(apply(m.2, 2, function(x) aggregate(  as.numeric(x), dulicated.names, my.median))) 
#   }else if(agg.fun.to.use == "max"){
#     m2_agg_duplicates = as.data.frame(apply(m.2, 2, function(x) aggregate(  as.numeric(x), dulicated.names, my.max ))) 
#   }else if(agg.fun.to.use == "sum"){
#     m2_agg_duplicates = as.data.frame(apply(m.2, 2, function(x) aggregate(  as.numeric(x), dulicated.names, my.sum ))) 
#   }else{
#     cat("Error - agg.fun.to.use:",agg.fun.to.use,", was not found")
#   }
#   print(dim(m2_agg_duplicates))
#   row.names(m2_agg_duplicates) = m2_agg_duplicates[,1] # first col, Description consists of gene names
#   m2_agg_duplicates = m2_agg_duplicates[, grepl(".x$", colnames(m2_agg_duplicates))] # keep data in ".x" col
#   if(is.test.hugo.symbols) row.names(m2_agg_duplicates) = fix_HUGO_excel_symbol(row.names(m2_agg_duplicates))
#   colnames(m2_agg_duplicates) = gsub(".x$", "", colnames(m2_agg_duplicates))
#   print(dim(m2_agg_duplicates))
#   return(m2_agg_duplicates)
# }










### *** DEGs via limma (move to RNA.seq.Utils)
#######################  diff.expression.limma
diff.expression.limma <- function(mat.expr, 
                                  is_log_transform = F, is_normalize = T, is_voom = F, test.cols.numbers=NULL, ctrl.cols.numbers=NULL,test.cols.names=NULL, ctrl.cols.names=NULL, 
                                  design.batch.vector = NULL, design.numeric.vector = NULL,
                                  batch.cols.chars=NULL, is.batch.last = T, numeric.covariate=NULL, numeric.covariate.name = "purity",
                                  trend = T, robust = T, proportion = 0.01){

  # If mat.expr TPM, process with my.log.normalize() first and run with is_log_transform = is_voom = F, & is_normalize = T
  stopifnot(!is.null(test.cols.numbers) | !is.null(test.cols.names))
  stopifnot(!is.null(ctrl.cols.numbers) | !is.null(ctrl.cols.names))
  
  if(is.null(test.cols.numbers)) test.cols.numbers = which(colnames(mat.expr) %in% c(test.cols.names) )
  if(is.null(ctrl.cols.numbers)) ctrl.cols.numbers = which(colnames(mat.expr) %in% c(ctrl.cols.names) )
  
  mat_data = mat.expr[ , c(test.cols.numbers,ctrl.cols.numbers)]
  mat_data = mat_data[rowSums(abs(mat_data))>0  ,] # remove non-expressed genes
  
  
  #### design
  if(is.null(design.batch.vector) & is.null(design.numeric.vector)){
    design = return.design.based.on.colnames (test.cols.numbers=test.cols.numbers, ctrl.cols.numbers=ctrl.cols.numbers,  colnames(mat_data),
                                              batch.cols.chars=batch.cols.chars, is.batch.last = is.batch.last,
                                              numeric.covariate=numeric.covariate, numeric.covariate.name = numeric.covariate.name)
  }else{
    print("design based on vectors")
    design = return.design.based.on.vectors (test.cols.numbers, ctrl.cols.numbers, design.batch.vector , design.numeric.vector, numeric.covariate.name)
  }
  # if(is.null(batch.cols.chars)){
  #   if(is.null(numeric.covariate)){
  #     design=cbind(test=c(rep(1,length(test.cols.numbers)),rep(0,length(ctrl.cols.numbers))),
  #                  ctrl=c(rep(0,length(test.cols.numbers)),rep(1,length(ctrl.cols.numbers))))
  #   }else{
  #     group <- as.factor(c(rep("test",length(test.cols.numbers)), rep("ctrl",length(ctrl.cols.numbers))  ))
  #     samples.names = paste0( str_split_fixed(colnames(mat_data), "_", 3)[,1],"_",str_split_fixed(colnames(mat_data), "_", 3)[,2]) # assume batch in name e.g., "065_Prospective1_Liver"...
  #     design <- model.matrix(~0 + group + numeric.covariate [ samples.names ] )
  #     colnames(design) <- gsub("group", "", colnames(design))
  #     colnames(design)[ncol(design)] = numeric.covariate.name
  #     write_clip(design)
  #   }
  # }else{
  #   group <- as.factor(c(rep("test",length(test.cols.numbers)), rep("ctrl",length(ctrl.cols.numbers))  ))
  #   batch = rep("", ncol(mat_data) )
  #   for(batch.char in batch.cols.chars){
  #     if(is.batch.last){
  #       batch[grep(paste0(batch.char,"$"), colnames(mat_data))] = batch.char
  #     }else{
  #       batch[grep(batch.char, colnames(mat_data))] = batch.char
  #     }
  #   }
  #   batch = as.factor(batch)
  # 
  #   if(is.null(numeric.covariate)){
  #     design <- model.matrix(~0 + group + batch)
  #     colnames(design) <- gsub("group", "", colnames(design))
  #   }else{
  #     samples.names = paste0( str_split_fixed(colnames(mat_data), "_", 3)[,1],"_",str_split_fixed(colnames(mat_data), "_", 3)[,2]) # assume batch in name e.g., "065_Prospective1_Liver"...
  #     design <- model.matrix(~0 + group + batch + numeric.covariate [ samples.names ] )
  #     colnames(design) <- gsub("group", "", colnames(design))
  #     colnames(design)[ncol(design)] = numeric.covariate.name
  #     write_clip(design)
  #   }
  # }
  ### END of design
      

  # voom: Mean-variance trend. Also do log2. But what if already in log space - don't run voom...
  if(is_voom){
    print("voom - correcting Mean-variance trend. Also do log2 (i.e., assuming it wasn't already in log...)")
    y <- DGEList(mat_data)
    y <- calcNormFactors(y)  # not clear it's needed
    v <- voom(y,design,plot = F)
    ## duplicateCorrelation "using duplicateCorrelation with limma+voom for RNA-seq data"
    #dupcor <- duplicateCorrelation(v, design = v$design, block = df$Subject) # map df$Subject into code...
    
  }else{
    if(is_log_transform) mat_data<- log2(mat_data+1)
    if(is_normalize){
      #mat_data = normalize( mat_data, norm.method = "quantile", within = FALSE, data.type = c("rseq") )
      #mat_data <- withinLaneNormalization(mat_data)
      mat_data <- normalizeBetweenArrays(mat_data) # Not vetted.. (seems OK)
    } 
    v = mat_data
  }

  fit <- lmFit(v, design)
  cont.matrix <- makeContrasts(STvsCO=test-ctrl, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2,robust=robust, trend=trend, proportion=proportion) # logical, should an intensity-trend be allowed for the prior variance? Default is that the prior variance is constant.

  return(fit2)
}
####################### process.topTable
process.topTable <- function(fit2, number=Inf, p.value = 1, adjust.method="BH",sort.by="p", round.digits = 5, is.round = F, is.sig=T,
                             cols.2.keep.qw = "Gene logFC  minus.log.p.val minus.log.q.val" ){
  
  my_table = topTable(fit2, number=number, p.value=p.value, adjust.method=adjust.method, sort.by=sort.by)
  my_table$minus.log.q.val = -log10(my_table$adj.P.Val)
  my_table$minus.log.p.val = -log10(my_table$P.Value)
  
  if(is.round)   my_table = round(my_table, round.digits)
  if(is.sig)   my_table = signif(my_table, digits = round.digits)
  my_table$Gene = row.names(my_table)
  my_table = my_table[,qw(cols.2.keep.qw)] # AveExpr
}

#######################  return.design.based.on.colnames
return.design.based.on.colnames <- function(test.cols.numbers, ctrl.cols.numbers,  colnames,
                                            batch.cols.chars=NULL, is.batch.last = T,
                                            numeric.covariate=NULL, numeric.covariate.name = "purity"){
  if(is.null(batch.cols.chars)){
    if(is.null(numeric.covariate)){
      design=cbind(test=c(rep(1,length(test.cols.numbers)),rep(0,length(ctrl.cols.numbers))),
                   ctrl=c(rep(0,length(test.cols.numbers)),rep(1,length(ctrl.cols.numbers))))
    }else{
      group <- as.factor(c(rep("test",length(test.cols.numbers)), rep("ctrl",length(ctrl.cols.numbers))  ))
      samples.names = paste0( str_split_fixed(colnames, "_", 3)[,1],"_",str_split_fixed(colnames, "_", 3)[,2]) # assume batch in name e.g., "065_Prospective1_Liver"...   
      design <- model.matrix(~0 + group + numeric.covariate [ samples.names ] )
      colnames(design) <- gsub("group", "", colnames(design))
      colnames(design)[ncol(design)] = numeric.covariate.name
      write_clip(design)
    }
  }else{
    group <- as.factor(c(rep("test",length(test.cols.numbers)), rep("ctrl",length(ctrl.cols.numbers))  ))
    batch = rep("", length(ctrl.cols.numbers)+length(test.cols.numbers) )
    for(batch.char in batch.cols.chars){
      if(is.batch.last){
        batch[grep(paste0(batch.char,"$"), colnames)] = batch.char
      }else{
        batch[grep(batch.char, colnames)] = batch.char
      }
    }
    batch = as.factor(batch)
    
    if(is.null(numeric.covariate)){
      design <- model.matrix(~0 + group + batch)
      colnames(design) <- gsub("group", "", colnames(design))
    }else{
      samples.names = paste0( str_split_fixed(colnames, "_", 3)[,1],"_",str_split_fixed(colnames, "_", 3)[,2]) # assume batch in name e.g., "065_Prospective1_Liver"...   
      design <- model.matrix(~0 + group + batch + numeric.covariate [ samples.names ] )
      colnames(design) <- gsub("group", "", colnames(design))
      colnames(design)[ncol(design)] = numeric.covariate.name
      write_clip(design)
    }
  }
  return(design)
}
#######################  return.design.based.on.vectors
return.design.based.on.vectors <- function(test.cols.numbers, ctrl.cols.numbers, design.batch.vector , design.numeric.vector, numeric.covariate.name = "purity"){
  
  group <- as.factor(c(rep("test",length(test.cols.numbers)), rep("ctrl",length(ctrl.cols.numbers))  ))
  if(is.null(design.batch.vector)){
    if( is.null(design.numeric.vector)){
      print("No covariates in design")
      design=cbind(test=c(rep(1,length(test.cols.numbers)),rep(0,length(ctrl.cols.numbers))),
                   ctrl=c(rep(0,length(test.cols.numbers)),rep(1,length(ctrl.cols.numbers))))           
    }else{
      print("Only numeric covariate in design")
      design <- model.matrix(~0 + group + design.numeric.vector )
      colnames(design) <- gsub("group", "", colnames(design))
      colnames(design)[ncol(design)] = numeric.covariate.name
      write_clip(design)        
    }
    
  }else{
    batch = as.factor(design.batch.vector)
    if(is.null(design.numeric.vector)){
      print("Only batch covariate in design")
      design <- model.matrix(~0 + group + batch)
      colnames(design) <- gsub("group", "", colnames(design))
    }else{
      print("Dual batch and numric covariates in design")
      design <- model.matrix(~0 + group + batch + design.numeric.vector )
      colnames(design) <- gsub("group", "", colnames(design))
      colnames(design)[ncol(design)] = numeric.covariate.name
      write_clip(design)
    }
  }
  return(design)
}





#######################  my.association.test
my.association.test <- function(X, Y, 
                                sigX =NULL, sigY =NULL, 
                                top_X_names.in=NULL, top_Y_names.in =NULL, 
                                num.top.X=100, num.top.Y=100, 
                                is.top.vs.bottom = T ) # if is.top.vs.bottom=F, hard to detect depletion ration
{
  if(is.top.vs.bottom){ # compare top vs. bottom third
    num.top.X = floor(ncol(X)/3)
    num.top.Y = floor(ncol(Y)/3)
  }else{ # compare top quartile to other
    num.top.X = floor(ncol(X)/4)
    num.top.Y = floor(ncol(Y)/4)
  }
  
  
  set.of.sig.X = sigX;
  if(is.null(sigX)) set.of.sig.X= row.names(X)
  set.of.sig.Y = sigY;
  if(is.null(sigY)) set.of.sig.Y = row.names(Y)
  
  Fisher.test.df = data.frame(matrix(ncol = 1, nrow = 1))
  
  for(sig.X in set.of.sig.X){
    if(is.null(top_X_names.in) | is.null(sigX)){
      vec.X = X[sig.X, ]
      vec.X = vec.X[order(vec.X, decreasing = T)]
      top_X = vec.X[1:num.top.X]
      top_X_names = names(top_X)
      if(is.top.vs.bottom){
        vec.X = vec.X[order(vec.X, decreasing = F)]
        top_X = vec.X[1:num.top.X]
        bottome_X_names = names(top_X)        
      }
    }else{
      top_X_names = top_X_names.in
    }
    for(sig.Y in set.of.sig.Y){
      if(sig.Y == sig.X) next;
      if(any(paste(sig.Y,sig.X,sep = ".") %in% row.names(Fisher.test.df))) next;
      if(is.null(top_Y_names.in) | is.null(sigX)){
        vec.Y = Y[sig.Y, ]
        vec.Y = vec.Y[order(vec.Y, decreasing = T)]
        top_Y = vec.Y[1:num.top.Y]
        top_Y_names = names(top_Y)
        if(is.top.vs.bottom){
          vec.Y = vec.Y[order(vec.Y, decreasing = F)]
          top_Y = vec.Y[1:num.top.Y]
          bottom_Y_names = names(top_Y)            
        }
      }else{
        top_Y_names = top_Y_names.in
      }
      X.and.Y = length( which( top_X_names %in% top_Y_names) )
      if(is.top.vs.bottom){
        X.and.nY = length( which( top_X_names %in% bottom_Y_names) )
      }else{
        nY = colnames(Y)[!colnames(Y) %in% top_Y_names]
        X.and.nY = length( which( top_X_names %in% nY) )        
      }
      if(is.top.vs.bottom){
        nX.and.Y = length(which(bottome_X_names %in% top_Y_names))
        nX.and.nY = length(which(bottome_X_names %in% bottom_Y_names))
      }else{
        nX = colnames(X)[!colnames(X) %in% top_X_names]
        nX.and.Y = length(which(nX %in% top_Y_names))
        nX.and.nY = length(which(nX %in% nY))        
      }
      my.fisher.test = fisher.test(matrix(c(X.and.Y, nX.and.Y,X.and.nY, nX.and.nY),ncol = 2))
      Fisher.test.df[paste(sig.X,sig.Y,sep = "."),"p.value"] = my.fisher.test$p.value
      Fisher.test.df[paste(sig.X,sig.Y,sep = "."),"estimate"] = my.fisher.test$estimate
      
      mutual.cols = colnames(X)[colnames(X) %in% colnames(Y)]
      my.cor.test = cor.test(as.numeric(Y[sig.Y, mutual.cols]), as.numeric(X[sig.X, mutual.cols]))
      Fisher.test.df[paste(sig.X,sig.Y,sep = "."),"cor.p.value"] = my.cor.test$p.value
      Fisher.test.df[paste(sig.X,sig.Y,sep = "."),"cor.estimate"] = my.cor.test$estimate
    }
  }
  dim(Fisher.test.df)
  Fisher.test.df= Fisher.test.df[order(Fisher.test.df$p.value),]
  Fisher.test.df = Fisher.test.df[,-1]
  
  return(Fisher.test.df)
}









#######################  aggregate_values_by_gene_sets
aggregate_values_by_gene_sets <-function(expr_, Signature_gene_sets, is_median_for_gene_sets = F, verbos_level = 1){

  gene_sets_counts <- data.frame(matrix(ncol = ncol(expr_), nrow = 1))
  colnames(gene_sets_counts) = colnames(expr_)
  
  signature_values = list();
  Signatures_found = ""
  
  for(gene_set_name in names(Signature_gene_sets)){
    print(paste(gene_set_name));
    genes_found = Signature_gene_sets[[gene_set_name]] [ Signature_gene_sets[[gene_set_name]] %in% row.names(expr_)  ]
    if(length(genes_found)<1) next;
    Signatures_found = c(Signatures_found, gene_set_name)
    if(verbos_level>1) print(paste("length=",length(genes_found), paste(genes_found, collapse="," )))
    signature_values[[gene_set_name]] <- expr_[which(row.names(expr_) %in% Signature_gene_sets[[gene_set_name]]),]
    if(is_median_for_gene_sets){
      gene_sets_counts[gene_set_name,] <- apply(signature_values[[gene_set_name]], 2, median)
    }else{
      gene_sets_counts[gene_set_name,] <- apply(signature_values[[gene_set_name]], 2, sum)
    }
  }
  gene_sets_counts = gene_sets_counts[row.names(gene_sets_counts) %in% names(Signature_gene_sets),]  
  return(gene_sets_counts)
}



#######################  grep_multi_patterns
grep_multi_patterns <- function(pattern2match, x, ignore.case = F, invert = F, also.pattern2match = NULL){
  matched.pattern =   unique (grep(paste(pattern2match,collapse="|"), 
                      x, value=TRUE,ignore.case=ignore.case,invert=invert))
  if(!is.null(also.pattern2match)) matched.pattern = grep_multi_patterns(also.pattern2match, matched.pattern)
  
  return(matched.pattern)
}
#######################  grepl_multi_patterns
grepl_multi_patterns <- function(pattern2match, x, ignore.case = F, invert = F, also.pattern2match = NULL){
  matched.pattern =   unique (grep(paste(pattern2match,collapse="|"), 
                                   x, value=TRUE,ignore.case=ignore.case,invert=invert))
  if(!is.null(also.pattern2match)) matched.pattern = grep_multi_patterns(also.pattern2match, matched.pattern)
  
  return(x %in% matched.pattern)
}
#######################  grep_multi_patterns_index
grep_multi_patterns_index <- function(pattern2match, x, ignore.case = F, invert = F, exact = F){
  if(exact) {
    pattern2match <- paste0("^",pattern2match,"$")
  }
  grep(paste(pattern2match,collapse="|"), 
       x, value=F,ignore.case=ignore.case,invert=invert)
}



#######################  check for a package, install it otherwise and load again
pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}


####################### fix dates
fix_HUGO_excel_symbol <- function(x){
  # Replace date-like format (issue when saved in excel) with gene names (e.g., DEC1 FEB1 MARC1 MARCH1 SEP15 etc.)
  gene.names = x
  
  names_fixed = qw("DEC1 FEB1 FEB2 FEB4 FEB5 FEB6 FEB7 MARC1 MARC2 MARCH1 MARCH10 MARCH11 MARCH2 MARCH3 MARCH4 MARCH5 MARCH6 MARCH7 MARCH8 MARCH9 SEP15 SEPT1 SEPT10 SEPT11 SEPT12 SEPT14 SEPT2 SEPT3 SEPT4 SEPT5 SEPT6 SEPT7 SEPT8 SEPT9")
  names_as_dates = qw("1-Dec  1-Feb   2-Feb   4-Feb   5-Feb   6-Feb   7-Feb   1-Mar   2-Mar   1-Mar   10-Mar  11-Mar  2-Mar   3-Mar   4-Mar   5-Mar   6-Mar   7-Mar   8-Mar   9-Mar   15-Sep  1-Sep   10-Sep  11-Sep  12-Sep  14-Sep  2-Sep   3-Sep   4-Sep   5-Sep   6-Sep   7-Sep   8-Sep   9-Sep")
  names_of_genes_ambig = qw("1-Mar   2-Mar")

  for(i in 1:length(names_of_genes_ambig)){
    if(any(grepl(names_of_genes_ambig[i],gene.names))) print(paste("WARN !!! can't fix ambig gene symbol name", gene.names[grepl(names_of_genes_ambig[i],gene.names)]))
  }
  for(i in 1:length(names_as_dates)){
    print( gene.names[grepl(names_as_dates[i],gene.names)] )
    gene.names = gsub( names_as_dates[i], names_fixed[i] , gene.names )
  }
  return(gene.names)
}



####################### auto fontsize for pheatmap
my_pheatmap <- function(x){
  fontsize_row = ifelse(dim(x)[1] < 50, 10, max(3,(10 - dim(x)[1]/50)))
  fontsize_col = ifelse(dim(x)[2] < 50, 10, max(3,(10 - dim(x)[2]/50)))
  pheatmap(x,fontsize_col=fontsize_col,fontsize_row=fontsize_row)
}

####################### get_number_from_str
get_number_from_str <-function(x) {  
  return( as.numeric( gsub("[^0-9]", "", x)))  
}











######################## query.COR.mat.values
query.COR.mat.values <- function(mat, query.A, query.B, is.print = F){
  # return values either mat[query.A, query.B] or mat[query.B, query.A], whichever is not NA

  All.gene.in.mat = unique(c(row.names(mat), colnames(mat)))
  if(is.print){
    print(paste0("No available query for", paste0(setdiff(query.A, All.gene.in.mat), collapse=",") ))
    print(paste0("No available query for", paste0(setdiff(query.B, All.gene.in.mat), collapse=",") ))    
  }

  mat.fill.NA = data.frame(matrix())
  mat.fill.NA [query.A, query.B] = NA
  mat.fill.NA = mat.fill.NA[-1,-1, drop=F]
  
  query.rows = intersect(query.A,intersect(query.A, row.names(mat)))
  query.cols = intersect(query.B,intersect(query.B, colnames(mat)))    
  mat.fill.NA [query.rows, query.cols] = mat[query.rows, query.cols]
  
  # for(qA in intersect(query.A, colnames(mat)) ){
  #   for(qB in intersect(query.B, row.names(mat)) ){
  #     ##if(!is.na(mat[qB, qA]) ) mat.fill.NA [qB, qA] = mat[qB, qA]
  #     if(!is.na(mat[qB, qA]) ) mat.fill.NA [qA, qB] = mat[qB, qA]
  #   }
  # }
  t.mat = t(mat)
  myList <- list();
  query.rows =  intersect(query.A,row.names(t.mat))
  query.cols = intersect(query.B,colnames(t.mat))
  myList[[1]] <- mat.fill.NA[ query.rows ,  query.cols ]
  myList[[2]] <- t.mat      [ query.rows,   query.cols ]
  mat.fill.NA[query.rows,  query.cols] = Reduce(`+`, lapply(myList,function(x) {x[is.na(x)] <-0;x}))
  return(mat.fill.NA)
}


######################## producing a full, symmetrical matrix with minimal NA
produce.full.symetrical.mat <- function(mat){
  # fill the NA pairs in the matrix [x,y] = [y,x]
  print(dim(mat))
  mat = as.data.frame(mat)
  col.row.names = unique(c(row.names(mat), colnames(mat)))
  mat.fix = query.COR.mat.values(mat, col.row.names, col.row.names)
  
}
######################## produce_LDA_projection
produce_LDA_projection <-function(Exp_Mat, Annot_table, Annot_col_name, number_of_top_var_genes=2000){
  Top.Var.genes = top_var_genes(Exp_Mat, number_of_top_var_genes, number_of_expression_bins=1)
  lda.data = merge( t(Exp_Mat[ Top.Var.genes ,]), Annot_table[,Annot_col_name, drop=F], by = "row.names", all=F  )
  row.names(lda.data) = lda.data$Row.names
  lda.data = lda.data[,-grep("Row.names", colnames(lda.data))]
  mformula=as.formula(paste0( Annot_col_name," ~ ", " ."))
  
  mylda <- lda(formula = mformula, data = lda.data) #  method = {t, moment, mle, mve
  plda <- predict(object = mylda,
                  newdata = lda.data)
  df = data.frame(LD1= plda$x[,1], LD2=plda$x[,2], lab = row.names(plda$x))
  df = merge(df, Annot_table, by = "row.names")
  row.names(df) = df$Row.names
  return(list(df=df, mylda=mylda, plda=plda))
}










################################################################
### "my" short Arithmetic functions
my.ecdf <- function(x,perc) { ecdf(x)(perc) }
my.signif <- function(x, digits=3){ signif(x, digits=digits) } 
range01 <- function(x, na.rm = T){(x-min(x,na.rm=na.rm))/(max(x,na.rm=na.rm)-min(x,na.rm=na.rm))} # *** range01 Not obvious it's working.. use pnorm() instead to get percentile for Z-score 
my.max <- function(x) {max(x, na.rm=TRUE)}
my.min <- function(x) {min(x, na.rm=TRUE)}
my.mean <- function(x) {mean(x, na.rm=TRUE)}
my.sd <- function(x) {sd(x, na.rm=TRUE)}
my.median <- function(x) {median(x, na.rm=TRUE)}
my.quantile <- function(x) {quantile(x, probs =0.9, na.rm=TRUE)}
my.Q1 <- function(x) {quantile(x, probs =0.25, na.rm=TRUE)}
my.Q3 <- function(x) {quantile(x, probs =0.75, na.rm=TRUE)}
my.sum <- function(x) {sum(x, na.rm=TRUE)}

if.T.keep.element <- function(X) return(X[X])
if.Pos.keep.element <- function(X) return(X[X>0])

######################## top_var_genes
top_var_genes <- function(expression_matrix, number_of_top_var_genes=1000, number_of_expression_bins = 4) {
  if(number_of_expression_bins>1){
    Top.Var.genes = c()
    mymeans <- apply( expression_matrix  ,1, mean,na.rm=TRUE)
    prev_mean_quantile = 0
    for(nbin in 1:number_of_expression_bins){
      print(nbin)
      mean_quantile = quantile(mymeans, prob = nbin / number_of_expression_bins)
      genes_in_bin = names(mymeans)[mymeans >= prev_mean_quantile & mymeans <= mean_quantile]
      prev_mean_quantile = mean_quantile
      
      myvars <- apply( expression_matrix[genes_in_bin,]  ,1, var,na.rm=TRUE) 
      myvars <- sort(myvars,decreasing=TRUE) 
      Top.Var.genes <- c(Top.Var.genes, names(myvars)[1:number_of_top_var_genes/number_of_expression_bins])
    }
    return( unique(Top.Var.genes) )
  }else{
    myvars <- apply( expression_matrix  ,1, var,na.rm=TRUE) 
    myvars <- sort(myvars,decreasing=TRUE) 
    Top.Var.genes <- names(myvars)[1:number_of_top_var_genes]
    return( Top.Var.genes )
  }
}  


### my "+" operator with na.rm = TRUE
# https://stackoverflow.com/questions/45311490/is-it-possible-to-skip-na-values-in-operator
`%+%` <- function(x, y)  mapply(sum, x, y, MoreArgs = list(na.rm = TRUE))
### element-wise summation of matrices
# https://stackoverflow.com/questions/14147655/is-there-an-r-function-for-the-element-wise-summation-of-the-matrices-stored-as
my.add <- function(x) Reduce("%+%", x)


######################## LRT.simple - Get p-values for 2 logL values of nested models with n df:
LRT.simple <-function(logL0, logL1, degFree=1)
{
  if((logL1 >0) || (logL0 >0)){
    cat("The log likelihood values need to be negative", "\n", sep="");
    return()   
  }
  LLog  <- 2*(logL1 - (logL0));
  pVal <- pchisq(LLog, df=degFree, lower.tail=FALSE);
  cat("p-value of models with logL of: M0= ", logL0," M1= ", logL1, " with df= ",degFree,"\n",pVal,"\n", sep="");  
  return(pVal)
}

#######################  my.scale.zscore
my.scale.zscore <- function(Mat, scale.by = "rows", center = TRUE, scale = TRUE){
  Zscored = Mat
  if(scale.by == "rows"){
    Zscored <- t(apply(Mat,1, scale, center = center, scale = scale)) # scale by row (compare to other components)
    colnames(Zscored) = colnames(Mat)    
  }else if(scale.by == "cols"){
    Zscored <- apply(Mat,2, scale, center = center, scale = scale) # scale by row (compare to other components)
    row.names(Zscored) = row.names(Mat)    
  }else{
    print('scale.by = {rows, cols}')
  }
  return(Zscored)
}

######################### 3 new cols - median mean        sd
produce.z.scored.matrix <- function(mat, is.by.row=T, is.keep.only.orginal.cols=T){
  
  if(is.by.row){
    mat_ = cbind( mat, 
                  median= apply(mat, 1, median),
                  mean= apply(mat, 1, mean),
                  sd= apply(mat, 1, sd)
    )
    mat.z.scored= sweep(mat_,1,t(mat_[,"mean"]), FUN = "-")
    mat.z.scored= sweep(mat.z.scored,1,t(mat_[,"sd"]), FUN = "/")
    if(is.keep.only.orginal.cols){
      mat.z.scored = mat.z.scored[,-which(colnames(mat.z.scored) %in% c("mean", "sd", "median"))]
    }
  }
  
  print(dim(mat.z.scored))
  return(mat.z.scored)
}


#######################  my.summary
my.summary <- function(df){
  df <- df[,which(unlist(lapply(df, is.numeric))),drop=F ]
  my.summary.df =   do.call(data.frame, 
                            list(mean = apply(df, 2, my.mean),
                                 sd = apply(df, 2, my.sd),
                                 Q1 = apply(df, 2,  my.Q1),
                                 median = apply(df, 2, my.median),
                                 Q3 = apply(df, 2,  my.Q3),
                                 min = apply(df, 2, my.min),
                                 max = apply(df, 2, my.max),
                                 n = apply(df, 2, length)))
  
  return(my.summary.df)
}

################################ my.head
my.head <- function(df, nrows = 5, ncols = 5) print(df[1:nrows, 1:ncols])

################################ my.table.as.data.frame
my.table.as.data.frame <- function(df, by.col, count = "Freq"){
  my.table =  data.frame( table( df[,by.col] ) )
  colnames(my.table) = c(by.col, count)
  return(my.table)
}

################################ fill.missing.values.in.df
fill.missing.values.in.df <- function(df, df.fill, match.by = "Sample_ID", cols.2.fill.str = "BX_LOCATION	BX_ER	BX_PR	BX_HER2OVERALL"){
  # Make match.by ("Sample_ID") as Row.names (assume unique)
  stopifnot( match.by %in% names(df) &&  match.by %in% names(df.fill))
  row.names(df) = df [,match.by]
  row.names(df.fill) = df.fill [,match.by]
  cols.2.fill = qw(cols.2.fill.str)
  for(mCol in cols.2.fill){
    keys.missing = df[ which( is.na(df[,mCol]) | df[,mCol]=="" ) ,match.by]
    keys.missing = df[ which( ! is.valid(df[,mCol]) ) ,match.by]
    for(mKey in keys.missing){ if( is.valid(df.fill[mKey,mCol]) ) df[mKey,mCol] = df.fill[mKey,mCol]   }
  }
  # df     [keys.missing,cols.2.fill]
  # df.fill[keys.missing,cols.2.fill]
  return(df)
}






################################ 
################################ generic.quantitative.associations.test
generic.quantitative.associations.test <- function(df, feat1, feat2.set=NULL, feat2.set.2ignore = c()){
  # make quantitative.associations.df - quantify association of numeric.cols (correlation), character.cols/logical.cols (Welch's t-test)
  all.cols = names(df)
  numeric.cols  <- unlist(lapply(df, is.numeric))  %>% if.T.keep.element %>% names
  character.cols  <- unlist(lapply(df, is.character))  %>% if.T.keep.element %>% names
  logical.cols  <- unlist(lapply(df, is.logical))  %>% if.T.keep.element %>% names
  categorical.cols <- c(logical.cols, character.cols)
  
  if(!is.null(feat2.set)){
    numeric.cols = intersect(numeric.cols, feat2.set)
    categorical.cols = intersect(categorical.cols, feat2.set)
  }
  
  quantitative.associations.df = data.frame(matrix())
  for(mC in numeric.cols  ){
    if(mC %in% feat2.set.2ignore) next;
    key = paste(mC, feat1, sep=".by.")
    quantitative.associations.df[key, "feat1"] = feat1
    quantitative.associations.df[key, "Feature"] = mC
    dat = subset(df, !is.na(df[,feat1]))
    dat$x = dat[,mC]
    dat$y = as.numeric( dat[,feat1] )
    if(length(which(!is.na(dat$x)))<20 ) next;
    Cor = cor.test(dat$y, dat$x, use = "pairwise.complete.obs")
    quantitative.associations.df[key, "p.value"] = Cor$p.value
    quantitative.associations.df[key, "R"] = Cor$estimate    
  }    
  for(mC in categorical.cols ){
    if(mC %in% feat2.set.2ignore) next;
    key = paste(mC, feat1, sep=".by.")
    quantitative.associations.df[key, "feat1"] = feat1
    quantitative.associations.df[key, "Feature"] = mC
    dat = subset(df, !is.na(df[,feat1]))
    dat$x = dat[,mC]
    dat[which(dat$x=="" | is.na(dat$x)),"x"] = "_"
    dat$y = as.numeric( dat[,feat1] )
    table.x = table(dat[,mC])
    dat = subset(dat, x %in% names(table.x[table.x>5]))
    if(nrow(dat) < 50) next;
    if(length(table(dat[,mC])) < 2) next;
    Welch = oneway.test(y~x, data=dat, na.action=na.omit, var.equal=FALSE)
    quantitative.associations.df[key, "p.value"] = Welch$p.value
    summ = tapply(dat$y, dat$x , summary)
    for(mG in names(summ)){
      quantitative.associations.df[key, paste0(mG, ".Median")] = summ[[mG]][ "Median" ]
    }
  }
  quantitative.associations.df = quantitative.associations.df[-1,-1]
  return(quantitative.associations.df)
} 



############################ annotate.Genes.with.Cell.Type
annotate.Genes.with.Cell.Type <- function(my_table, Gene.col = "Gene",
                                          Cell.type.Top.marker.Genes.list.key="inclusive1000+canonical.Cell.type",
                                          is.q_val_select_genes = T,
                                          mCell.type.FDR = "Epithelial.ER.breast", p.val.col = "P.Value",is.keep.only.Cell.Type.Genes = F,
                                          is.annotate.all.Cell.Types = T){
  
  # Add cols to table - is.Cell.Type?
  # Compute q_val_select_genes (require: mCell.type.FDR, p.val.col)
  
  stopifnot(exists("Cell.type.Top.marker.Genes.list.of.lists"))
  
  Cell.Type.Genes = Cell.type.Top.marker.Genes.list.of.lists [[Cell.type.Top.marker.Genes.list.key]] [[mCell.type.FDR]]
  print(paste("Cell.Type", mCell.type.FDR, "with ", length(Cell.Type.Genes), "Cell.Type.Genes"))
  if(is.annotate.all.Cell.Types){
    for(mC in names(Cell.type.Top.marker.Genes.list.of.lists [[Cell.type.Top.marker.Genes.list.key]]) ){
      mC.Genes = Cell.type.Top.marker.Genes.list.of.lists [[Cell.type.Top.marker.Genes.list.key]] [[mC]]
      my_table[,mC] = my_table[,Gene.col] %in% mC.Genes 
    }
  }else{
    my_table[, mCell.type.FDR] = my_table[,Gene.col] %in% Cell.Type.Genes 
  }
  my_table.cell.type = subset(my_table, my_table[,Gene.col] %in% Cell.Type.Genes)
  if(is.q_val_select_genes){
    my_table.cell.type$q_val_select_genes = p.adjust(my_table.cell.type$P.Value, method = "BH")
    if(is.keep.only.Cell.Type.Genes){
      return(my_table.cell.type)
    }else{
      return(merge(my_table, my_table.cell.type[,c(Gene.col, "q_val_select_genes")],by=Gene.col, all=T, sort=F ) )
    }    
  }else{
    return(my_table)
  }
  
}





################################ is.valid - not NULL/NA/""?
is.valid.val <- function(x){
  return( !( is.null(x) | is.na(x) | x == "")  )
}
is.valid.var <- function(x){
  return( exists(deparse(substitute(x))) && is.valid.val(x)  )
}
is.valid <- function(x){  return(is.valid.var(x)) }

################################ my.str_split_fixed.and.paste
my.str_split_fixed.and.paste <- function(splt.string, splt.pattern, splt.n, keep.n, paste.pattern = "_"){
  str_split = str_split_fixed(splt.string, splt.pattern ,splt.n)
  if(length(keep.n) > 1){
    return( apply( str_split[,keep.n] , 1 , paste , collapse = paste.pattern ) )
  }else{
    return( str_split[,keep.n]  )
  }
}



################################################################ 
################################ write.gct.2 from Pablo's lib
write.gct.2 <- function(gct.data.frame, descs = "", filename) 
{
  f <- file(filename, "w")
  cat("#1.2", "\n", file = f, append = TRUE, sep = "")
  cat(dim(gct.data.frame)[1], "\t", dim(gct.data.frame)[2], "\n", file = f, append = TRUE, sep = "")
  cat("Name", "\t", file = f, append = TRUE, sep = "")
  cat("Description", file = f, append = TRUE, sep = "")
  
  colnames <- colnames(gct.data.frame)
  cat("\t", colnames[1], file = f, append = TRUE, sep = "")
  
  if (length(colnames) > 1) {
    for (j in 2:length(colnames)) {
      cat("\t", colnames[j], file = f, append = TRUE, sep = "")
    }
  }
  cat("\n", file = f, append = TRUE, sep = "\t")
  
  oldWarn <- options(warn = -1)
  m <- matrix(nrow = dim(gct.data.frame)[1], ncol = dim(gct.data.frame)[2] +  2)
  m[, 1] <- row.names(gct.data.frame)
  if (length(descs) > 1) {
    m[, 2] <- descs
  } else {
    m[, 2] <- row.names(gct.data.frame)
  }
  index <- 3
  for (i in 1:dim(gct.data.frame)[2]) {
    m[, index] <- gct.data.frame[, i]
    index <- index + 1
  }
  write.table(m, file = f, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
  close(f)
  options(warn = 0)
}


################################ read.gct from Pablo's lib
MSIG.Gct2Frame <- function(filename = "NULL") { 
  # Reads a gene expression dataset in GCT format and converts it into an R data frame
  ds <- read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T, na.strings = "")
  descs <- ds[,1]
  ds <- ds[-1]
  row.names <- row.names(ds)
  names <- names(ds)
  return(list(ds = ds, row.names = row.names, descs = descs, names = names))
}



#### https://stackoverflow.com/questions/6979917/how-to-unload-a-package-without-restarting-r
# multiple versions of a package loaded at once
detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}










################################################################ 
# http://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
#######################  qqunif.plot
library(lattice)
qqunif.plot<-function(pvalues, 
                      names=NULL,cex.label=0.8,
                      should.thin=T, thin.obs.places=2, thin.exp.places=2, 
                      xlab=expression(paste("Expected (",-log[10], " p-value)")),
                      ylab=expression(paste("Observed (",-log[10], " p-value)")), 
                      draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
                      already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
                      par.settings=list(superpose.symbol=list(pch=pch)), ...) {
  
  
  #error checking
  if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
  if(!(class(pvalues)=="numeric" || 
       (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
    stop("pvalue vector is not numeric, can't draw plot")
  if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
  if (already.transformed==FALSE) {
    if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
  } else {
    if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
  }
  
  
  grp<-NULL
  n<-1
  exp.x<-c()
  if(is.list(pvalues)) {
    nn<-sapply(pvalues, length)
    rs<-cumsum(nn)
    re<-rs-nn+1
    n<-min(nn)
    if (!is.null(names(pvalues))) {
      grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
      names(pvalues)<-NULL
    } else {
      grp=factor(rep(1:length(pvalues), nn))
    }
    pvo<-pvalues
    pvalues<-numeric(sum(nn))
    exp.x<-numeric(sum(nn))
    for(i in 1:length(pvo)) {
      if (!already.transformed) {
        pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
        exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
      } else {
        pvalues[rs[i]:re[i]] <- pvo[[i]]
        exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
      }
    }
  } else {
    n <- length(pvalues)+1
    if (!already.transformed) {
      exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
      pvalues <- -log10(pvalues)
    } else {
      exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
    }
  }
  
  
  #this is a helper function to draw the confidence interval
  panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
    require(grid)
    conf.points = min(conf.points, n-1);
    mpts<-matrix(nrow=conf.points*2, ncol=2)
    for(i in seq(from=1, to=conf.points)) {
      mpts[i,1]<- -log10((i-.5)/n)
      mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
      mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
      mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
    }
    grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
  }
  
  #reduce number of points to plot
  if (should.thin==T) {
    if (!is.null(grp)) {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places),
                                grp=grp))
      grp = thin$grp
    } else {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places)))
    }
    pvalues <- thin$pvalues
    exp.x <- thin$exp.x
  }
  gc()
  
  prepanel.qqunif= function(x,y,...) {
    A = list()
    A$xlim = range(x, y)*1.02
    A$xlim[1]=0
    A$ylim = A$xlim
    return(A)
  }
  
  #draw the plot
  #draw the plot
  if(is.null(names)){
    xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
           prepanel=prepanel, scales=list(axs="i"), pch=pch,
           panel = function(x, y, ...) {
             if (draw.conf) {
               panel.qqconf(n, conf.points=conf.points, 
                            conf.col=conf.col, conf.alpha=conf.alpha)
             };
             panel.xyplot(x,y, ...);
             panel.abline(0,1);
           }, par.settings=par.settings, ...
    )
  }else{
    xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
           prepanel=prepanel, scales=list(axs="i"), pch=pch,
           panel = function(x, y, ...) {
             if (draw.conf) {
               panel.qqconf(n, conf.points=conf.points, 
                            conf.col=conf.col, conf.alpha=conf.alpha)
             };
             panel.xyplot(x,y, ...);
             ltext(x=x, y=y, labels=names, offset=0, cex=cex.label); # , pos=1
             panel.abline(0,1);
           }, par.settings=par.settings, ...
    )
  }
}


######################## rmd2rscript
require("roxygen2")
rmd2rscript <- function(infile){
  # read the file
  flIn <- readLines(infile)
  # identify the start of code blocks
  cdStrt <- which(grepl(flIn, pattern = "```{r*", perl = TRUE))
  # identify the end of code blocks
  cdEnd <- sapply(cdStrt, function(x){
    preidx <- which(grepl(flIn[-(1:x)], pattern = "```", perl = TRUE))[1]
    return(preidx + x)
  })
  # define an expansion function
  # strip code block indacators
  flIn[c(cdStrt, cdEnd)] <- ""
  expFun <- function(strt, End){
    strt <- strt+1
    End <- End-1
    return(strt:End)
  }
  idx <- unlist(mapply(FUN = expFun, strt = cdStrt, End = cdEnd, 
                       SIMPLIFY = FALSE))
  # add comments to all lines except code blocks
  comIdx <- 1:length(flIn)
  comIdx <- comIdx[-idx]
  for(i in comIdx){
    flIn[i] <- paste("#' ", flIn[i], sep = "")
  }
  # create an output file
  nm <- paste(strsplit(infile, split = "\\.")[[1]][-length(strsplit(infile, split = "\\.")[[1]])], collapse = ".") #strsplit(infile, split = "\\.")[[1]][1]
  flOut <- file(paste(nm, ".R", sep = ""), "w")
  for(i in 1:length(flIn)){
    cat(flIn[i], "\n", file = flOut, sep = "\t")
  }
  close(flOut)
}






##################################################################### Ad-hoc functions
# #######################  subset.of.CMJ.ORFs
# subset.of.CMJ.ORFs <- function(set, only.Drug.subset){
#   if(only.Drug.subset == "ERK"){
#     subset = set[grepl("_E$", set )]
#     length(subset)
#   }else if(only.Drug.subset == "RAF_MEK"){
#     subset = set[grepl("_RM$", set )]
#     length(subset)
#   }else if(only.Drug.subset == "DMSO"){
#     subset = set[!grepl("_RM$", set ) & !grepl("_E$", set )]
#     length(subset)
#   }else if(only.Drug.subset == "non-DMSO"){
#     subset = set[grepl("_RM$", set ) | grepl("_E$", set )]
#     length(subset)  
#   }else{
#     return(set)
#   }
#   return(subset)
# }



# #######################  gg.box.plot.with.names - OLD version
# gg.box.plot.with.names <- function(dat, y.values.name, txt.size = 3, point.size =1, jitter.factor = 1.5, is.repel = F, is.anova = F, is.kruskal=F, is.levene = F, title = title,  
#                                    ylim = NA, fill=NA,is.print.label=T,scale_color_manual.values=NULL){
#   datJit <- dat
#   if("fill" %in% colnames(dat)){
#     fill.factor = 1/length(unique(dat$fill))-0.1
#     mymin = -fill.factor
#     mymax = fill.factor
#     val = as.numeric(factor(dat$fill))
#     #val.normalized.shifted = (val-min(val))/(max(val)-min(val)) 
#     val.normalized.shifted= (mymax-mymin)/(max(val)-min(val))*(val-max(val) )+mymax
#     datJit$xj <- jitter(as.numeric(factor(dat$x))+val.normalized.shifted, factor=jitter.factor)
#   }else{
#     datJit$xj <- jitter(as.numeric(factor(dat$x)), factor=jitter.factor)
#   }
#   
#   if(is.anova){
#     my.aov <- aov(y ~ x, data = dat)
#     my.aov.p.value =   signif(summary(my.aov)[[1]][["Pr(>F)"]][1], digits = 3)
#   }
#   if(is.kruskal){
#     my.kruskal.test = kruskal.test(dat$y~factor(dat$x))
#     my.kruskal.test.pvalue = signif(my.kruskal.test$p.value, digits = 3)  
#   }
#   if(is.levene){
#     my.levene.test = levene.test(dat$y~factor(dat$x))
#     my.levene.p.value =   signif((my.levene.test)[["Pr(>F)"]][1], digits = 3)
#   }
#   
#   pp = ggplot(dat,aes(x=x,y=y))
#   if("fill" %in% colnames(dat)) pp = ggplot(dat,aes(x=x,y=y,fill=fill))
#   
#   pp = pp + ylab(y.values.name) +
#     geom_boxplot()
#   
#   if("fill" %in% colnames(dat)){
#     if("shape" %in% colnames(dat)){
#       pp = pp + geom_point(data=datJit,aes(x=xj, colour=x, fill=fill, shape=shape ),size=point.size)
#     }else{
#       pp = pp + geom_point(data=datJit,aes(x=xj, colour=x, fill=fill),size=point.size)
#     }
#   }else{
#     if("shape" %in% colnames(dat)){
#       pp = pp + geom_point(data=datJit,aes(x=xj, colour=x, shape=shape ),size=point.size)
#     }else{
#       pp = pp + geom_point(data=datJit,aes(x=xj, colour=x),size=point.size)
#     }    
#   }
#   
#   if(!is.null(scale_color_manual.values)) pp = pp + scale_color_manual(values=scale_color_manual.values)
#   
#   pp = pp +  #geom_text(data=datJit,aes(x=xj,label=lab),size=txt.size ) +
#     theme(axis.text=element_text(size=12) ,axis.title=element_text(size=12) ) 
#   if(is.anova) pp = pp +   labs(title=title,
#                                 subtitle= paste0("ANOVA ", my.aov.p.value)  ) # x, y, title, subtitle
#   if(is.kruskal) pp = pp +   labs(title=title,
#                                   subtitle= paste0("kruskal ", my.kruskal.test.pvalue)  ) # x, y, title, subtitle  
#   if(is.levene) pp = pp +   labs(title=title,
#                                  subtitle= paste0("levene ", my.levene.p.value)  ) # x, y, title, subtitle    
#   
#   if(!is.na(ylim)) pp = pp + ylim(ylim)
#   if(is.print.label){
#     if(is.repel){
#       pp = pp + geom_text_repel(data=datJit,aes(x=xj,label=lab),size=txt.size )
#     }else{
#       pp = pp + geom_text(data=datJit,aes(x=xj,label=lab),size=txt.size )
#     }
#   }
#   pp <- pp + theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black")) + theme_bw()  
#   pp
# } 


