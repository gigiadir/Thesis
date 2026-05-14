### Utiles.FH.R

# library(plyr)
# library(dplyr)
# library(magrittr)
# library(GenVisR)

################################## unique.Pts
unique.Pts <- function(Sample_IDs, ncharPt = 3){
  return(unique(substr(Sample_IDs,1, ncharPt)))
}
  


################################## add.cols.to.MAF
# in cols: Genome_Change Protein_Change Hugo_Symbol
add.cols.to.MAF <- function(mMAF)
{
  # Protein_Change.fill - Fill Protein_Change with Genome_Change - Splice site
  idx.missing.Protein_Change = which(mMAF$Protein_Change == "")
  mMAF[,"Protein_Change.fill"] = mMAF[,"Protein_Change"]
  mMAF[idx.missing.Protein_Change,"Protein_Change.fill"] = mMAF[idx.missing.Protein_Change,"Genome_Change"]
  
  # Pt.Allele
  mMAF =  within(mMAF, Allele <- paste(Hugo_Symbol,Protein_Change.fill, sep="."))
  mMAF =  within(mMAF, Sample.Allele <- paste(Sample_ID,Allele, sep="."))
  mMAF =  within(mMAF, Sample.Genome_Change <- paste(Sample_ID,Genome_Change, sep="."))
  
  mMAF =  within(mMAF, Sample.Gene <- paste(Sample_ID,Hugo_Symbol, sep="."))
  mMAF =  within(mMAF, Pt.Allele <- paste(Pt,Allele, sep="."))
  mMAF =  within(mMAF, Allele.Pt <- paste(Allele,Pt, sep="."))
  
  return(mMAF)
}


################################## add.Mut.Type.to.MAF with BRCA and MBC for hotspots
# "GENIE.Mut.MBC.Alleles.table")) # *** get latest GENIE (freeze/record)
# CosmicMutant.Breast.Alleles.table")) # get latest Cosmic (freeze/record)
# allAnnotatedVariants.oncokb")) # *** new version - v1.19_patch_1
# catalog_of_validated_oncogenic_mutations")) # *** check version 
# *** replace with oncokb annotation CGC.TSG.genes.aug, CGC.oncogene.genes.aug

add.Mut.Type.to.MAF <- function(mMAF, hotfix.Alleles = NULL)
{
  print("Mut.Type annotation for MAF, Categories (later overwrite earlier): Other, NonFunc, missense, NonCoding, CodingBounds, In_Frame_Indel, warmspot.Onco, warmspot.TSG, hotspot, hotspot.Onco, hotspot.TSG, truncating, ValidAllele, OncoAllele, LOF, GOF")
  
  requisite.Objects = qw("GENIE.Mut.MBC.Alleles.table CosmicMutant.Breast.Alleles.table allAnnotatedVariants.oncokb catalog_of_validated_oncogenic_mutations CGC.oncogene.genes.aug CGC.TSG.genes.aug
                         Variant_Classification.NonFunc Variant_Classification.NonCoding Variant_Classification.In_Frame_Indel")
  
  if(! all.exists(requisite.Objects) ) my.save.RData(requisite.Objects, is.load.only.if.not.exists = T)
    
  mMAF[,"Mut.Type"] = "Other" # Default, init
  idx.ignore = which(mMAF$Variant_Classification %in% Variant_Classification.NonFunc) # Variant_Classification.2ignore, "Silent" "Intron" "IGR"
  mMAF[idx.ignore,"Mut.Type"] = "NonFunc" 
  
  idx.Missense = which(mMAF$Variant_Classification %in% qw("Missense_Mutation"))
  mMAF[idx.Missense,"Mut.Type"] = "missense"
  
  idx.NonCoding = which(mMAF$Variant_Classification %in% Variant_Classification.NonCoding ) #  qw("3'UTR 5'Flank 5'UTR lincRNA RNA")
  mMAF[idx.NonCoding,"Mut.Type"] = "NonCoding"  
  
  idx.CodingBounds = which(mMAF$Variant_Classification %in% Variant_Classification.CodingBounds ) # qw("De_novo_Start_InFrame Nonstop_Mutation Start_Codon_SNP")
  mMAF[idx.CodingBounds,"Mut.Type"] = "CodingBounds"  
  
  idx.In_Frame_Indel = which(mMAF$Variant_Classification %in% Variant_Classification.In_Frame_Indel ) # qw("In_Frame_Del In_Frame_Ins")
  mMAF[idx.In_Frame_Indel,"Mut.Type"] = "In_Frame_Indel"  
 
  idx.warmspot = which(mMAF$Allele %in% names(GENIE.Mut.MBC.Alleles.table[GENIE.Mut.MBC.Alleles.table>=1]) 
                       | mMAF$Allele %in% names(CosmicMutant.Breast.Alleles.table[CosmicMutant.Breast.Alleles.table>=1])
                       | mMAF$Allele %in% names(CosmicMutant.Breast.Alleles.table[CosmicMutant.all.Alleles.table>=1])  )
  idx.warmspot.Onco = intersect(idx.warmspot , which( mMAF$Hugo_Symbol %in%  CGC.oncogene.genes.aug  )) # c(CGC.oncogene.genes,Oncogenes.Vogel.2013)
  idx.warmspot.TSG = intersect(idx.warmspot , which( mMAF$Hugo_Symbol %in%  CGC.TSG.genes.aug  )) # c(CGC.TSG.genes,TSGs.Vogel.2013)
  mMAF[idx.warmspot.Onco,"Mut.Type"] = "warmspot.Onco" 
  mMAF[idx.warmspot.TSG,"Mut.Type"] = "warmspot.TSG"  
  
  idx.hotspot = which(mMAF$Allele %in% names(GENIE.Mut.MBC.Alleles.table[GENIE.Mut.MBC.Alleles.table>=3]) 
                      | mMAF$Allele %in% names(CosmicMutant.Breast.Alleles.table[CosmicMutant.Breast.Alleles.table>=3])
                      | mMAF$Allele %in% names(CosmicMutant.Breast.Alleles.table[CosmicMutant.all.Alleles.table>=5]))
  idx.hotspot.Onco = intersect(idx.hotspot , which( mMAF$Hugo_Symbol %in%  CGC.oncogene.genes.aug  )) # c(CGC.oncogene.genes,Oncogenes.Vogel.2013)
  idx.hotspot.TSG = intersect(idx.hotspot , which( mMAF$Hugo_Symbol %in%  CGC.TSG.genes.aug  )) # c(CGC.TSG.genes,TSGs.Vogel.2013)
  
  mMAF[idx.hotspot,"Mut.Type"] = "hotspot" 
  mMAF[idx.hotspot.Onco,"Mut.Type"] = "hotspot.Onco" 
  mMAF[idx.hotspot.TSG,"Mut.Type"] = "hotspot.TSG" 
  
  idx.truncating = which(mMAF$Variant_Classification %in% qw("Splice_Site Frame_Shift_Ins Frame_Shift_Del Nonsense_Mutation De_novo_Start_OutOfFrame") ) # Nonstop_Mutation
  mMAF[idx.truncating,"Mut.Type"] = "truncating" 
  
  idx.GOF = which(mMAF$Allele %in% subset(allAnnotatedVariants.oncokb, Mutation.Effect=="Gain-of-function" )$Allele   )
  mMAF[idx.GOF,"Mut.Type"] = "GOF" 
  
  idx.LOF = which(mMAF$Allele %in% subset(allAnnotatedVariants.oncokb, grepl("Loss", Mutation.Effect))$Allele   )
  mMAF[idx.LOF,"Mut.Type"] = "LOF" 
  
  mMAF$is.annot.cancergenomeinterpreter = mMAF$Allele %in% catalog_of_validated_oncogenic_mutations$Allele
  idx.add.OncoAllele = which(mMAF$is.annot.cancergenomeinterpreter & mMAF$Hugo_Symbol %in%  CGC.oncogene.genes.aug & !  mMAF$Mut.Type %in% qw("GOF LOF") )
  mMAF[idx.add.OncoAllele,"Mut.Type"] = "OncoAllele"

  idx.add.ValidAllele = which(mMAF$is.annot.cancergenomeinterpreter & !  mMAF$Mut.Type %in% qw("GOF LOF OncoAllele") )
  mMAF[idx.add.ValidAllele,"Mut.Type"] = "ValidAllele"
  if(!is.null(hotfix.Alleles)){
    idx.add.ValidAllele = which(mMAF$Allele %in%  hotfix.Alleles)
    mMAF[idx.add.ValidAllele,"Mut.Type"] = "ValidAllele"
  }
  
  return(mMAF)
}



################################## add.Mut.Type.to.CNA
add.Mut.Type.to.CNA <- function(mMAF, col.read="rescaled_total_cn.ploidyDiff", 
                                filter.for.FocalAmp = "Num.genes",
                                col.focal.genes ="Num.genes",max.ngenes.for.FocalAmp = 100,
                                col.focal.length ="length", max.length.for.FocalAmp = 3000000,
                                min.for.Gain = 1.5,min.for.Amp = 3,min.for.HighAmp = 6,min.for.FocalAmp = 9,HOMDEL.col = "Bi_Allelic.inc")
{

  
  mMAF[,"Mut.Type"] = ""
  
  CNA.events.idx.list = list()
  GAIN.idx = which( mMAF[,col.read] >= min.for.Gain )
  mMAF[GAIN.idx,"Mut.Type"] = "GAIN"

  AMP.idx = which( mMAF[,col.read] >= min.for.Amp )
  mMAF[AMP.idx,"Mut.Type"] = "AMP"
  
  HAMP.idx = which( mMAF[,col.read] >= min.for.HighAmp )
  mMAF[HAMP.idx,"Mut.Type"] = "HighAMP"

  
  FAMP.idx = which( mMAF[,col.read] >= min.for.FocalAmp)
  
  FAMP.idx.length = which(mMAF[,col.focal.length] <= max.length.for.FocalAmp )
  CNA.events.idx.list[[col.focal.length]] = FAMP.idx.length
  
  FAMP.idx.genes = which( mMAF[,col.read] >= min.for.FocalAmp & mMAF[,col.focal.genes] <= max.ngenes.for.FocalAmp )
  CNA.events.idx.list[[col.focal.genes]] = FAMP.idx.genes
  
  if(! is.null(filter.for.FocalAmp))  FAMP.idx = intersect(FAMP.idx, CNA.events.idx.list[[filter.for.FocalAmp]])
  mMAF[FAMP.idx,"Mut.Type"] = "FocalAMP"
  
  HOMDEL.idx.naive = which( mMAF[,HOMDEL.col] )
  # *** if the SNV is truncatin, foce HOMDEL
  # *** HOMDEL only if not AMP
  mMAF[ setdiff(HOMDEL.idx.naive, AMP.idx) , "Mut.Type"] = "HOMDEL" 
  return(mMAF)
  
}


################################## plot.Mut.Type.freq.PerGene
plot.Mut.Type.freq.PerGene <- function(MAF,
                                       Mut.Type.order.str = "GOF hotspot.Onco missense LOF hotspot.TSG truncating  Other",
                                       Mut.Type.order.cols.str = "red tomato3 black blue4 darkcyan blue3  grey40",
                                       Mut.Type.2.omit.str = "NonFunc"
                                       ){
  Mut.Type.order = qw(Mut.Type.order.str)
  Mut.Type.order.cols = qw(Mut.Type.order.cols.str)
  Mut.Type.2.omit = qw(Mut.Type.2.omit.str)
  
  PerGene.Type.Stats = data.frame(matrix())
  for(mmgene in MAF$Hugo_Symbol){
    aaa = subset(MAF, Hugo_Symbol==mmgene)
    for(myMut.Type in unique(MAF$Mut.Type) ){
      if(myMut.Type %in% Mut.Type.2.omit ) next;
      bbb = subset(aaa, myMut.Type == Mut.Type)
      PerGene.Type.Stats[mmgene, myMut.Type ] = length( unique( bbb$Pt) )
    }
    PerGene.Type.Stats[mmgene, "Total" ] = length(unique(aaa$Pt))
  }
  PerGene.Type.Stats =   PerGene.Type.Stats[-1,-1]
  PerGene.Type.Stats$Gene = row.names(PerGene.Type.Stats)
  
  
  xxx = data.table::melt( PerGene.Type.Stats[,colnames(PerGene.Type.Stats)] , id.vars="Gene" )
  xxx = subset(xxx, variable %in% Mut.Type.order)
  xxx$frac = xxx$value / length(unique(MAF$Pt))
  
  xxx$Gene <- factor(xxx$Gene, levels =  PerGene.Type.Stats[order(PerGene.Type.Stats$Total, decreasing = T),"Gene"]    )
  xxx$variable <- factor(xxx$variable, levels =  Mut.Type.order)
  
  #xxx = subset(xxx, ! variable %in% Mut.Type.2.omit)
  gg =ggplot(xxx, aes(x = Gene, y = frac, fill = variable)) + geom_bar(stat = "identity") # Already counted in xxx
  gg = gg + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=10)) + scale_fill_manual( values= Mut.Type.order.cols )
  return(gg)
  
  
}









#######################  assign.High.Low.quantile
assign.High.Low.quantile <- function(df, param.2.use, new.param.name, qtile=0.25, index.2.use =NULL){
  df[,new.param.name] = ""
  if(is.null(index.2.use)){
    df[ which( df[,param.2.use] >= quantile(df[,param.2.use], 1-qtile, na.rm=T)) ,new.param.name] = "High"
    df[ which( df[,param.2.use] < quantile(df[,param.2.use], qtile, na.rm=T)) ,new.param.name] = "Low"    
  }else{
    df[ which( df[index.2.use,param.2.use] >= quantile(df[index.2.use,param.2.use], 1-qtile, na.rm=T)) ,new.param.name] = "High"
    df[ which( df[index.2.use,param.2.use] < quantile(df[index.2.use,param.2.use], qtile, na.rm=T)) ,new.param.name] = "Low"    
  }
  return(df)
}

############## run.GenVisR func
run.GenVisR <- function( MAF , sample_id="Sample_ID" , mainLabelCol = "Protein_Change", geneOrder = NULL, PatientOrder.watefall = T, out.file = NULL, 
                         mainLabelSize=2, main_geneLabSize=10, height=8, width=12, rmvSilent = T, mainDropMut = T, plotMutBurden=F,
                         clinData=NULL, clinLegCol=1){
  mutationset = MAF
  mutationset$Variant_Classification = with(mapping.FH.MAF, MAF[match(mutationset$Variant_Classification, FH)])
  
  print(dim(mutationset))
  #catnl(unique(mutationset$Pt))
  
  mutationset$Tumor_Sample_Barcode = mutationset[,sample_id]
  
  if(PatientOrder.watefall){ #order by Patiet
    sampOrder = unique(mutationset$Tumor_Sample_Barcode)
    sampOrder = sampOrder[order(sampOrder)]
  }else{
    sampOrder = NULL
  }
  
  if(!is.null(out.file))  pdf(file= out.file, height=height, width=width)
  
  waterfall(mutationset, fileType="MAF", mainLabelCol=mainLabelCol, mainLabelSize=mainLabelSize, main_geneLabSize = main_geneLabSize, mainXlabel=TRUE,
            geneOrder=geneOrder, sampOrder=sampOrder,  rmvSilent = rmvSilent, mainDropMut=mainDropMut, plotMutBurden=plotMutBurden, plot_proportions=F,
            clinData=clinData, clinLegCol=clinLegCol)
  
  if(!is.null(out.file)) dev.off()  
  
}




#######################  add.Sample_ID.Pt.2.FCset
# entity.sample_id
add.Sample_ID.Pt.2.FCset <- function(mFCset, settype = "sample"){
  if(settype == "sample" ){
    if( any( grepl("^BRCA", mFCset$entity.sample_id) ))   warning("old FH format");
    
    sample_id.split = str_split_fixed(mFCset$entity.sample_id, "_", 7)
    mFCset$Pt = substrRight( sample_id.split[,4],3)
    mFCset$Sample_ID = paste0(  mFCset$Pt, "_", sample_id.split[,5])
    
  }else{
    warning("settype ", settype); return(NULL)
  }
  
  return(data.frame(mFCset))
}

#######################  add.Sample_ID.Pt.FHformat
# external_id_capture, external_id_rna
add.Sample_ID.Pt.FHformat <- function(mFCset, settype = "sample"){
  if(settype == "sample" ){
    
    external_id_capture = mFCset$external_id_capture
    idx.valid.external_id_capture = grep("_CCPM", external_id_capture)
    external_id_capture.split = str_split_fixed(external_id_capture, "_", 4)
    mFCset[ idx.valid.external_id_capture , "Pt" ] = substrRight( external_id_capture.split[ idx.valid.external_id_capture ,3], 3)
    mFCset[ idx.valid.external_id_capture , "Sample_ID" ] = paste0(mFCset[ idx.valid.external_id_capture , "Pt" ], "_", external_id_capture.split[ idx.valid.external_id_capture ,4])
    
    # RNA only
    external_id_rna = mFCset$external_id_rna
    idx.valid.external_id_rna = grep("_CCPM", external_id_rna)
    idx = setdiff(idx.valid.external_id_rna, idx.valid.external_id_capture)
    external_id.split = str_split_fixed(external_id_rna, "_", 4)
    mFCset[ idx , "Pt" ] = substrRight( external_id.split[ idx ,3], 3)
    mFCset[ idx , "Sample_ID" ] = paste0(mFCset[ idx , "Pt" ], "_", external_id.split[ idx ,4])
    
    # Blood biopsies
    entity.sample_id = mFCset$entity.sample_id
    entity.sample_id.split = str_split_fixed(entity.sample_id, "-", 3)
    idx = which(is.na(mFCset$Sample_ID))
    mFCset[ idx , "Pt" ] = substrRight( entity.sample_id.split[ idx ,2], 3)
    mFCset[ idx , "Sample_ID" ] = paste0(mFCset[ idx , "Pt" ], "_", entity.sample_id.split[ idx ,3])
    
    mFCset$Sample_ID = gsub("-", "_", mFCset$Sample_ID)
    
    # # Blood biopsies - short format
    # entity.sample_id.split = str_split_fixed(entity.sample_id, "-", 4)
    # idx = grep("_$", mFCset$Sample_ID)
    # mFCset[ idx , "Sample_ID" ] = paste0(mFCset[ idx , "Pt" ], "_", entity.sample_id.split[ idx ,4])
    
  }else{
    warning("settype ", settype); return(NULL)
  }
  
  return(data.frame(mFCset))
}



############# 459_BB2_noT_FC19310108 => 459_FC19310108
getSimplified.BB.ID <- function(input_ID ){
  output_ID.simplified = input_ID # initialize
  BB.idx = grep("_BB", input_ID)
  output_ID.simplified[BB.idx] = paste0(str_split_fixed(input_ID[BB.idx], "_", 4)[,1],  "_",  str_split_fixed(input_ID[BB.idx], "_", 4)[,4])
  return(output_ID.simplified)
}


### function - add "Sample_ID" based on FH-sample (later)
# use FH Sample (latest)
# make FH RData. Check if exists, load
getID.from.FH_Sample_ID_dictionary <- function(input_ID ,intype = "sample_id", outtype = "Sample_ID", is.make.names = T, re.make.Obj = F){
  ### def - outtype={Sample_ID, individual_id, }
  samples_FH.filename = paste0("~/Dropbox/Wagle.CCPM.ERpos/MS.2017/Samples.Clinical.Tracker/Samples.Tracker.2017.freeze/samples.BRCA_05246_inclusive.Tumor.vs.Normal.one.per.case.tsv")
  pairs_FH.filename = paste0("~/Dropbox/Wagle.CCPM.ERpos/MS.2017/Samples.Clinical.Tracker/Samples.Tracker.2017.freeze/pairs.BRCA_05246_inclusive.Tumor.vs.Normal.one.per.case.tsv")
  
  ### Not yet tested/used
  # Samples.Clinical.Tracker.folder = "~/Dropbox/Wagle.CCPM.ERpos/MS.2017/Samples.Clinical.Tracker/Samples.Tracker.2017.freeze"
  # mypset="BRCA_05246_inclusive.Tumor.vs.Normal.one.per.case.SampleUnique"
  # 
  # samples_FH.filename =  file.path(Samples.Clinical.Tracker.folder, paste0("samples.",mypset,".tsv"))
  # pairs_FH.filename =file.path(Samples.Clinical.Tracker.folder, paste0("pairs.",mypset,".tsv"))
  
  
  Obj.file = file.path(DataHub.Shiny.ERpos,"samples_FH_Sample_ID_dictionary.RData")
  
  if(intype == "sample_id" & outtype == "Sample_ID"){
    BB.elements = grep_multi_patterns (qw("_BB -BB bloodbiopsy cfDNA ^FC _FC"), input_ID)
    BB.index = which(input_ID %in% BB.elements)
    output_ID = rep("", length(input_ID))
    if(length(BB.index)>0){
      output_ID[-BB.index] = getID.from.FH_Sample_ID_dictionary (input_ID[-BB.index] ,intype , outtype , is.make.names , re.make.Obj)
      output_ID[BB.index] = gsub("BRCA-05246_CCPM_", "",input_ID[BB.index])
      output_ID[BB.index] = gsub("-","_",output_ID[BB.index])
      Pt = substrRight( str_split_fixed( output_ID[BB.index] , "_", 2)[,1], 3)
      output_ID[BB.index] = paste0( Pt ,"_", str_split_fixed( output_ID[BB.index] , "_", 2)[,2]  )      
      return(output_ID)
    }
  }

  if(file.exists(Obj.file) & !re.make.Obj){
    load(Obj.file)
  }else{
    stopifnot(file.exists(samples_FH.filename))
    samples_FH <- read.delim(paste0(samples_FH.filename), stringsAsFactors=F)
    #samples_FH$external_id_capture

    samples_FH$Sample_ID_capture = get.Sample_IDs.from.full.names(samples_FH$external_id_capture)
    samples_FH$Sample_ID_rna = get.Sample_IDs.from.full.names(samples_FH$external_id_rna)
    samples_FH$Sample_ID = samples_FH$Sample_ID_capture
    index_missing = which(!grepl("\\d", samples_FH$Sample_ID)) 
    samples_FH[index_missing,"Sample_ID"] = samples_FH$Sample_ID_rna[index_missing]
    samples_FH_Sample_ID_dictionary = samples_FH[,qw("Sample_ID sample_id  clean_bam_file_capture")] # individual_id
    
    stopifnot(file.exists(pairs_FH.filename))
    pairs_FH <- read.delim(paste0(pairs_FH.filename), stringsAsFactors=F)
    pairs_FH_samples_dictionary = pairs_FH[,qw("pair_id case_sample control_sample")] # individual_id
    save(samples_FH_Sample_ID_dictionary,pairs_FH_samples_dictionary, file=Obj.file)    
  }
  if(is.make.names){
    if(intype %in% colnames(samples_FH_Sample_ID_dictionary)){
      output_ID = samples_FH_Sample_ID_dictionary[ match( make.names(input_ID), make.names(samples_FH_Sample_ID_dictionary[,intype]) ) ,  grep(outtype, colnames(samples_FH_Sample_ID_dictionary)) ]
    }else if(intype %in% colnames(pairs_FH_samples_dictionary)){
      output_ID = pairs_FH_samples_dictionary[ match( make.names(input_ID), make.names(pairs_FH_samples_dictionary[,intype]) ) ,  grep(outtype, colnames(pairs_FH_samples_dictionary)) ]
    }
  }else{
    if(intype %in% colnames(samples_FH_Sample_ID_dictionary)){
      output_ID = samples_FH_Sample_ID_dictionary[ match(input_ID, samples_FH_Sample_ID_dictionary[,intype]) ,  grep(outtype, colnames(samples_FH_Sample_ID_dictionary)) ]
    }else if(intype %in% colnames(pairs_FH_samples_dictionary)){
      output_ID = pairs_FH_samples_dictionary[ match(input_ID, pairs_FH_samples_dictionary[,intype]) ,  grep(outtype, colnames(pairs_FH_samples_dictionary)) ]
    }    
  }
  missing.idx = which(is.na(output_ID))
  #output_ID[missing.idx] = input_ID[missing.idx]
  output_ID[missing.idx] = get.Sample_ID.from.samples.FC ( input_ID[missing.idx] )
  
  return(output_ID)
}


################################ getID.from.FC.MBCProject.Sample_ID
getID.from.FC.MBCProject.Sample_ID <- function(samples_vector, is.return.full.names.table = F, MBCproject.clinical.samples = MBCproject.clinical.samples){
  # A-la getID.from.FH_Sample_ID_dictionary
  df = merge(  data.frame(samples_vector ) , MBCproject.clinical.samples, 
               by.x = "samples_vector", by.y = "entity.sample_id", all.x=T, all.y=F, sort=F )
  
  row.names(df) = df$samples_vector
  
  matched = NA
  for(mS in subset(df, is.na(Sample_ID.MBC) & grepl("SM-", samples_vector) )$samples_vector ){
    SM.non.matched = my.str_split_fixed.and.paste( mS, "-", 4, 4)
    matched = subset( MBCproject.clinical.samples, grepl(SM.non.matched, MBCproject.clinical.samples ))$Sample_ID.MBC
    if(length(matched)==1) df[ mS , "Sample_ID.MBC" ] = matched
  }
  
  df$Sample_ID.format1 = my.str_split_fixed.and.paste( df$samples_vector, "_", 7, c(3,4), "_")
  df$Sample_ID.format2 = gsub("MBCProject_", "", my.str_split_fixed.and.paste( df$samples_vector, "-", 6, c(2,5), "_"), ignore.case = T )
  df$Sample_ID.format3 = gsub("MBCProject_", "", my.str_split_fixed.and.paste( df$samples_vector, "-", 4, c(1,4), "_"), ignore.case = T )
  
  df$Sample_ID = gsub("MBCPROJECT_", "",df$Sample_ID.MBC)
  df[is.na(df$Sample_ID) ,"Sample_ID"] = df[is.na(df$Sample_ID) ,"Sample_ID.format1"]
  df[is.na(as.numeric( substr(df$Sample_ID, 1, 4))) ,"Sample_ID"] = df[is.na(as.numeric( substr(df$Sample_ID, 1, 4))) ,"Sample_ID.format2"]
  df[is.na(as.numeric( substr(df$Sample_ID, 1, 4))) ,"Sample_ID"] = df[is.na(as.numeric( substr(df$Sample_ID, 1, 4))) ,"Sample_ID.format3"]
  
  if(is.return.full.names.table){
    return(df)
  }else{
    return(df$Sample_ID)
  }
}



### Get Individual from full names
get.Individuals.from.full.names <- function(sample_id_full.vec){
  
  Individuals = 
    llply(
      sample_id_full.vec, 
      function(x) strsplit(x, '_')[[1]][3] %>% 
        (function(x) substrRight(x,3))  %>% 
        unlist()
    ) %>% unlist()
  
  return(Individuals)
  
}



### Get Sample_ID from full names
get.Sample_IDs.from.full.names <- function(sample_id_full.vec){
  
  # BB.elements = grep_multi_patterns (qw("_BB -BB bloodbiopsy cfDNA ^FC"), sample_id_full.vec)
  # BB.index = which(sample_id_full.vec %in% BB.elements)
  # Sample_IDs = rep("", length(sample_id_full.vec))
  # if(length(BB.index)>0){
  #   Sample_IDs[-BB.index] = get.Sample_IDs.from.full.names(sample_id_full.vec[-BB.index])
  #   Sample_IDs[BB.index] = str_extract(sample_id_full.vec[BB.index], "FC[0-9]+")
  #   return(Sample_IDs)
  # }
    
  example.sample = sample_id_full.vec[which(sample_id_full.vec!="")][1] 
  num.of._.split = length( strsplit(example.sample, '_')[[1]] )
  PtNum.idx = num.of._.split-1
  Sample.suffix.idx = num.of._.split
  
  Individuals = 
    llply(
      sample_id_full.vec, 
      function(x) strsplit(x, '_')[[1]][PtNum.idx] %>% 
        (function(x) substrRight(x,3))  %>% 
        unlist()
    ) %>% unlist()
  
  suffix = 
    llply(
      sample_id_full.vec, 
      function(x) strsplit(x, '_')[[1]][Sample.suffix.idx] %>% 
        unlist()
    ) %>% unlist()
  
  Sample_IDs = paste0(Individuals,"_",suffix)
  Sample_IDs = gsub("[[:space:]]", "",Sample_IDs)
  
  return(Sample_IDs)
  
}

## Get Sample_ID from bam_file_rna
get.Sample_IDs.from.bam_file_rna <- function(bam_file_rna.vec){
  res = basename(bam_file_rna.vec)
  res = gsub("\\.bam$","",res)
  res = get.Sample_IDs.from.full.names(res)
  return(res)
  
}

#######################  get.Sample_ID.from.samples.FC
get.Sample_ID.from.samples.FC <- function(entity.sample_ids, samples.ref.file = "~/Dropbox/Wagle.CCPM.ERpos/MS.2017/Samples.Clinical.Tracker/sample.FC.An_ER_Pos_WAGLE-ms_2017.tsv" ){
  if(!exists("sample.FC.RNA.ws")) sample.FC.RNA.ws =  read.delim(samples.ref.file)
  
  entity.sample_ids = make.names(entity.sample_ids)
  sample.FC.RNA.ws.sub = subset(sample.FC.RNA.ws, make.names(entity.sample_id) %in%  entity.sample_ids )
  row.names(sample.FC.RNA.ws.sub) = make.names( sample.FC.RNA.ws.sub$entity.sample_id )
  
  bam_file_rnas = sample.FC.RNA.ws.sub[entity.sample_ids,"bam_file_rna"]
  get.Sample_IDs.from.bam_file_rna( bam_file_rnas )
  
}






