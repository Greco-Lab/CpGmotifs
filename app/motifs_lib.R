# given a target bed file and a background bed file 
#finds enriched motifs for a given flanking region size around target sites

#this script calls dreme tool, make sure you have it installed on your system 
# and set the path to the folder containing its executable in dreme_path

library(BSgenome.Hsapiens.UCSC.hg19)
require(XML)
library(methods)
library(tibble)
library(gtools)
library(Biostrings)



#change this to your dreme execution path
dreme_path = "/meme/bin/dreme-py3"

hg19 = BSgenome.Hsapiens.UCSC.hg19

#returns a sequence object containting hg19 sequences flanking each target site
get_fasta = function(ids, chrs, pos, flankSize,bsgen){
  #create dataframe       
  query = data.frame(id=as.character(ids),chr=as.character(chrs),start=pos-flankSize,end=pos+flankSize+1)
  seqs  = getSeq(bsgen,query$chr,start=query$start,end=query$end)
  names(seqs) = query$id
  return(seqs)
}

args <- commandArgs(trailingOnly = TRUE)

targ_name    = args[1] ### here goes the id of target 
back_name    = args[2] ### here goes the id of the background
bed_folder   = args[3] ### here goes the path of the bed folder
fasta_folder = args[4] ### here goes the path of the fasta folder
memeOutFold  = args[5] ### here goes the path of the meme_out folder
flankSize    = as.numeric(args[6]) ### here goes the size of the flanking region

#read target and background regions files 
back_pos  = read.table(paste0(bed_folder,"/",back_name,".bed"),sep = "\t", header = T,stringsAsFactors = F)
targ_pos  = read.table(paste0(bed_folder,"/",targ_name,".bed"),sep = "\t", header = T,stringsAsFactors = F)

back_fasta =  paste0("./fasta/",targ_name,"_back.fasta")
targ_fasta = paste0("./fasta/",targ_name,".fasta")

# get backgorund sequences and write them to a fasta file
back_seq = get_fasta(ids = back_pos$id, chrs = back_pos$chr, pos = back_pos$pos, flankSize = flankSize, bsgen = hg19)
writeXStringSet(back_seq, filepath=back_fasta, append=FALSE, format="fasta")


# get target sequences and write each one to a fasta file
targ_seq = get_fasta(targ_pos$id, targ_pos$chr, targ_pos$pos, flankSize, hg19)
writeXStringSet(targ_seq, filepath=targ_fasta, append=FALSE, format="fasta")


# call dreme with default parameters on the fasta target  
system(paste(dreme_path,"-oc",paste0(memeOutFold,"/",targ_name),"-p", 
              targ_fasta,"-n", back_fasta,">",paste(memeOutFold,"/",targ_name,"_log.txt",sep="")))

## after the dreme is finished 
# the system call is blocking by default 




###### retrive imput sites associted with each motif using the XML output file

# eventual filter on the nominal pvalue or the evalue 
pThresh = 0.05


types_mot = tibble::tibble(type = character(), motif = character(),length = numeric(), 
                           positives = numeric(), tot_pos = numeric(),
                           negatives = numeric(), tot_neg = numeric(),
                           pos_ratio = numeric(), neg_ratio = numeric(),
                           pvalue = numeric(), evalue = numeric())


imput_seqs = as.data.frame(targ_seq)


#read the dreme xml output to get enriched motifs
xml_out = xmlParse(paste0(memeOutFold,"/",targ_name,"/dreme.xml"))

rootnode = xmlRoot(xml_out)

model = rootnode[["model"]]
motifs = rootnode[["motifs"]]

allPositives = as.numeric(xmlAttrs(model[["positives"]])["count"])
allNegatives = as.numeric(xmlAttrs(model[["negatives"]])["count"])

for (i in 1:xmlSize(motifs)) {
  cgMot = c()
  motif  = motifs[[i]]
  mot_RE = as.character(xmlAttrs(motif)["seq"])
  length = as.numeric(xmlAttrs(motif)["length"])
  positives = as.numeric(xmlAttrs(motif)["p"])
  negatives = as.numeric(xmlAttrs(motif)["n"])
  pvalue = as.numeric(xmlAttrs(motif)["pvalue"])
  evalue = as.numeric(xmlAttrs(motif)["evalue"])

  #use evalue here if you want to filter by e-values 
  if ( length(pvalue) > 0 && evalue < pThresh){
    
    types_mot = add_row(types_mot, type = targ_name, motif = mot_RE ,length = length, 
                        positives = positives, negatives = negatives,
                        tot_pos = allPositives, tot_neg = allNegatives,
                        pos_ratio = positives/allPositives, neg_ratio = negatives/allNegatives,
                        pvalue = pvalue, evalue = evalue)
    
    matches = motif["match"]
    for (j in 1:xmlSize(matches)){
      match  = matches[[j]]
      seq    = as.character(xmlAttrs(match)["seq"])
      revseq = as.character(reverseComplement(DNAString(seq)))

      # search the cpgIds containg seq AND the reverse complement of seq. in the fasta file
      cgMot = union(cgMot,rownames(imput_seqs)[ grepl(seq,imput_seqs$x,fixed = TRUE) | grepl(revseq,imput_seqs$x,fixed = TRUE) ])
    }
  } 

  targ_mot = targ_pos[targ_pos$id %in% cgMot,]
  write.table(targ_mot,paste0(memeOutFold,"/",targ_name,"/",mot_RE,".txt"), quote = F, row.names = F, col.names = T)
}

write.table(types_mot,paste0(memeOutFold,"/",targ_name,"/summary.txt"), quote = F, row.names = F, col.names = T)
