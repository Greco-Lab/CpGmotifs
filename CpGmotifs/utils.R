library(tibble)
library(gplots)
library(dendextend)
library(foreach)
library(doParallel)
library(XML)
library(methods)
library(Biostrings)

#change this to your dreme execution path

# this script uses tomtom utility and meme motif databeses, make sure you have them installed 
# on your system and provide the corresponding paths here
meme_bin_path  <-  "/meme/bin/"
motif_dbs <- "/motif_databases/"
tf_db <- "JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme"

#tf_db <- "HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"

#retrives methylation information from the support set of the motif
get_meth_sign <- function(mot_cpg,idm_type) {
  print(idm_type)
  
  
  idm_type$hyper <- NA
  idm_type$hypo  <- NA
  #define hyper- of hypo- methylation status
  if ("Status" %in% colnames(idm_type)){
    idm_type$hyper <-  (idm_type$Status == "hyper")
    idm_type$hypo <-  (idm_type$Status == "hypo")
  } else {
    idm_type$hyper <-  NA
    idm_type$hypo <- NA
  }
  
  if (sum(!is.na(idm_type$hyper)) > 0 ){
    idm_type <- idm_type[!is.na(idm_type$hyper),]
    
    hyper_type <- rownames(idm_type)[idm_type$hyper]
    n_hyper  <-  length(intersect(hyper_type,mot_cpg$id))
    
    hypo_type <- rownames(idm_type)[idm_type$hypo]
    n_hypo <- length(intersect(hypo_type,mot_cpg$id))
    
    tot_hyper <- length(hyper_type)
    tot_hypo <- length(hypo_type)
    
    #test for enrichment of hyper- or hypo-methylated CpGs in the support set
    contingency_mat <- matrix(c(n_hyper, n_hypo, tot_hyper-n_hyper, tot_hypo-n_hypo) ,nrow=2,ncol=2)
    fisher <- fisher.test(contingency_mat)
    
    #save collected info in a data.frame
    motif_out_data <- data.frame(recovered_pos= nrow(mot_cpg), 
                                 n_hyper=n_hyper,n_hypo=n_hypo, 
                                 tot_hyper=tot_hyper, 
                                 tot_hypo = tot_hypo, 
                                 fisher.p = fisher$p.value)
  } else {
    motif_out_data <- data.frame(recovered_pos= nrow(mot_cpg), 
                                 n_hyper=0,n_hypo=0, 
                                 tot_hyper=0, 
                                 tot_hypo = 0, 
                                 fisher.p = NA)
  }
  
  motif_out_data
}


get_mot_sum <- function(iDMCs, meme_out) {
  types_mot <- tibble::tibble(type = character(), motif = character(),length = numeric(),
                             positives = numeric(), tot_pos = numeric(),
                             negatives = numeric(), tot_neg = numeric(),
                             pos_ratio = numeric(), neg_ratio = numeric(),
                             pvalue = numeric(), evalue = numeric(),
                             recovered_pos=numeric(),
                             n_hyper= numeric(), n_hypo = numeric(),
                             tot_hyper= numeric(),tot_hypo = numeric(),
                             fisher.p = numeric()
  )
  
  for (type in names(iDMCs)){
    type_cpg_of_int <- iDMCs[[type]]
    
    types_tab <- read.table(paste0(meme_out,"/",type,"/summary.txt"), header = T,stringsAsFactors = F)
    if (nrow(types_tab) > 0){
      types_tab <- cbind(types_tab, 
                      do.call(rbind, 
                              apply(types_tab, 1 ,
                                    function(x,idmcs){
                                      motif_file <- paste0(meme_out,"/",x["type"],"/",x["motif"],".txt")
                                      mot_cpg = read.table(motif_file , quote="\"", comment.char="",header = T, stringsAsFactors=FALSE)
                                      get_meth_sign(mot_cpg,idmcs)
                                    }, 
                                    idmcs=type_cpg_of_int)
                      )
      )
      types_mot <- dplyr::bind_rows(types_mot, types_tab)
    }
  }
  
  types_mot
}


assign_main_dir <- function(types_mot, unbalance_tresh, evalue_tresh, 
                            pvalue_tresh, fisher_thesh){
  
  #code methylation status of the support set as hyper if > unbalance_tresh % are hyper, 
  # hypo if > unbalance_tresh % are hypo, unbal otherwise
  types_mot$main_dir <- "neut"
  
  if(any(types_mot$n_hyper/(types_mot$recovered_pos) > unbalance_tresh)){
    types_mot[types_mot$n_hyper/(types_mot$recovered_pos) > unbalance_tresh, ]$main_dir <- "hyper"
  }
  if(any(types_mot$n_hypo/(types_mot$recovered_pos) > unbalance_tresh)){
    types_mot[types_mot$n_hypo/(types_mot$recovered_pos) > unbalance_tresh, ]$main_dir <- "hypo"  
  }
  
  types_mot <- types_mot %>% dplyr::filter(evalue <= evalue_tresh, fisher.p <= fisher_thesh, pvalue <= pvalue_tresh)
  
  types_mot
}





# motif_similiarity -------------------------------------------------------

motif_similiarity <- function(types_mot) {
  types_mot$name  <-  types_mot$type
  if(nrow(types_mot)>1){
    cat(paste0(">",types_mot$motif,"-->",types_mot$type,"\n",types_mot$motif,"\n",collapse = ""))
    names= paste0(types_mot$name,"-->",types_mot$motif,"-->",types_mot$main_dir)
    
    #create similiarity matrix
    mat_aln_all <- matrix(NA,nrow = nrow(types_mot), ncol = nrow(types_mot),dimnames = list(names,names))
    
    for(i in 1:nrow(mat_aln_all))
      for(j in 1:nrow(mat_aln_all))
        #here you can provide another function to compute the similarity between i-th and j-th motif
        #for example you can provide a function that akes two PWMs in input
        mat_aln_all[i,j] <- pairwiseAlignment(types_mot$motif[i], types_mot$motif[j], 
                                              substitutionMatrix=nucleotideSubstitutionMatrix(),
                                              type="global",scoreOnly=T,gapOpening=1, gapExtension=0.1)
    
    # min max normalization in the range 0-1
    range01 <- function(x){
      (x - min(x))/(max(x) - min(x))
    }
    
    mat_aln_all_0_1 <- apply(mat_aln_all, 2, range01)
    #mat_aln_d=apply(mat_aln_all,2,function(x) 1 - x/max(x))
    mat_aln_d_all <- 1 - mat_aln_all_0_1
  }
  else{
    mat_aln_d_all <- NULL
  }
  mat_aln_d_all
}

annotate_tfbs <- function(types,types_mot, meme_out_fold, meme_bin_path, motif_dbs, tf_db,progress_ind = NULL) {
  
  #setup parallel backend to use 8 processors
  #if you don't want to to run the analisys on multiple cores comment the following twolines
  # cl<-makeCluster(8)
  # registerDoParallel(cl)
  
  res_all <- foreach(type_can = types,.combine = rbind) %do% {
    require(XML)
    library("methods")
    
    progress_ind$inc(1/length(types), detail = paste("Analyzing ", type_can))
    
    res=NULL
    tomtom_comm <- paste0(meme_bin_path,"/tomtom -no-ssc -oc tmp -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 10.0 ",
                          meme_out_fold,"/",type_can,"/dreme.txt ",motif_dbs,"/",tf_db, " -oc ",meme_out_fold,"/",type_can,"/tomtom")
    #print(tomtom_comm)
    system(tomtom_comm)
    
    #res <- read.table(paste0(meme_out_fold,"/",type_can,"/tomtom/tomtom.tsv"),
    res <- tryCatch(read.table(paste0(meme_out_fold,"/",type_can,"/tomtom/tomtom.tsv"),
                      stringsAsFactors = F,
                      comment.char = "#"), error= function(e){
                        message("No motifs for ", type_can)
                      })
    if (!is.null(res)){
        colnames(res) <- c("Query_ID",	"Target_ID"	,"Optimal_offset",	"p.value"	,
                           "E.value",	"q.value",	"Overlap"	,"Query_consensus"	,
                           "Target_consensus",	"Orientation")
        res$type <- type_can
        res$tfbs <- sapply(res$Target_ID,function(x) strsplit(x,"_")[[1]][1])
    
        #parse dreme xml output
        xml_out <- xmlParse(paste0(meme_out_fold,"/",type_can,"/tomtom/tomtom.xml"))
        rootnode <- xmlRoot(xml_out)
        motifs <- xmlChildren(rootnode[["targets"]])
        tab_names <- do.call(rbind,
                            lapply(as.list(motifs),
                                   function(x) {
                                     c(as.character(xmlAttrs(x)["id"]),as.character(xmlAttrs(x)["alt"]))
                                   }
                            )
        )
        matches <- match(res$Target_ID,tab_names[,1])
        res$tfbs2 <- tab_names[matches,2]
    #Sys.sleep(1)
    }
    return(res)
  }
  
  #comment this if you dont use multi-core
  # stopCluster(cl)
  
  ann_tfbs <- types_mot %>%
              dplyr::select(type,motif,main_dir) %>%
              dplyr::inner_join(res_all, by = c( "type"="type","motif"="Query_ID"))
  
  print(ann_tfbs)
  ann_tfbs
  
}


filter_tfbs <- function(tfbs_tab, tf_pvalue_tresh, tf_evalue_tresh, tf_qvalue_thresh, tf_overlap_tresh){
  tfbs_tab %>% dplyr::filter(p.value <= tf_pvalue_tresh, E.value <= tf_evalue_tresh, 
                      q.value <= tf_qvalue_thresh, Overlap >= tf_overlap_tresh)
}

