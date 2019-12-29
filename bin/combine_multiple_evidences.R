library("tidyverse")
library("optparse")

################################################################################
## sub functions

add_rna_var_evidence=function(x,y){
  if("neoflow_var_ID" %in% names(x)){
    stop("add_rna_var_evidence: neoflow_var_ID was found in x!\n")
  }
  cat("add variant evidence from RNA-Seq data ...\n")
  x$neoflow_var_ID <- paste(x$Chr,x$Start,x$End,x$Ref,x$Alt,sep="_")
  y$neoflow_var_ID <- paste(y$Chr,y$Start,y$End,y$Ref,y$Alt,sep="_")
  x$rna_var_evidence <- ifelse(x$neoflow_var_ID %in% y$neoflow_var_ID,1,0)
  x$neoflow_var_ID <- NULL
  return(x)
}


add_gene_exp_evidence=function(x,y){
  
  cat("add gene expression evidence from RNA-Seq data ...\n")
  ## currently, we use gene (not transcript) expression data
  dat <- y %>% select(gene_id,TPM) %>% 
    distinct() %>%
    rename(Gene=gene_id,gene_expression=TPM)
  x <- left_join(x,dat,by="Gene")    
  return(x)
}


add_pro_exp_evidence=function(x,y){
  
  cat("add protein expression from proteomics data ...\n")
  ## currently, we use gene (not transcript) expression data
  dat <- y %>% select(protein,iBAQ) %>% 
    distinct() %>%
    rename(mRNA=protein,protein_expression=iBAQ)
  x <- left_join(x,dat,by="mRNA")    
  return(x)
}

is_variant_protein=function(x,prefix="VAR"){
  ids <- str_split(x$Protein,pattern = ";") %>% unlist()
  isvar <- FALSE
  if(all(str_detect(ids,pattern = "^VAR"))){
    isvar <- TRUE  
  }
  x$is_var <- isvar
  return(x)
}

split_proteins=function(x){
  ids <- str_split(x$Protein,pattern = ";") %>% unlist()
  ids_var <- str_replace_all(ids,pattern = "\\(.*\\)",replacement = "")
  pres <- data.frame(Variant_ID=ids_var,stringsAsFactors = FALSE)
  pres$Protein <- x$Protein
  if("peptide" %in% names(x)){
    pres$protein_var_evidence_pep <- x$peptide
  }
  return(pres)
  
}

add_pro_var_evidence=function(x,y,fdr=0.01){
  y_res <- y %>% filter(QValue<=fdr) %>% select(Protein,Peptide) %>%
    distinct() %>% 
    group_by(Protein,Peptide) %>%
    do(is_variant_protein(.)) %>% 
    filter(is_var==TRUE) %>%
    ungroup() %>%
    group_by(Protein) %>%
    summarise(peptide=paste(Peptide,collapse = ";")) %>%
    ungroup() %>%
    group_by(Protein) %>%
    do(split_proteins(.)) %>%
    ungroup() %>%
    select(Variant_ID,protein_var_evidence_pep) %>%
    distinct()
  
  vres <- left_join(x,y_res,by="Variant_ID")
  
  vres$protein_var_evidence <- ifelse(is.na(vres$protein_var_evidence_pep),0,1)
  return(vres)
  
}


################################################################################


option_list = list(
  make_option(c("-i","--sample_id"), type="character", default="test", 
              help="Sample ID [default= %default]", metavar="character"),
  make_option(c("-m","--mhc"), type="character", default=NULL, 
              help="MHC-binding prediction file [default= %default]", metavar="character"),
  make_option(c("-f","--mhc_filter"), type="integer", default=150, 
              help="The cutoff of MHC-binding affinity [default= %default]", metavar="150"),
  make_option("--pro_exp", type="character", default=NULL, 
              help="Protein expression data file [default= %default]", metavar="character"),
  make_option("--pro_var", type="character", default=NULL,
              help="Variant expression file (protein) [default= %default]", metavar="character"),
  make_option("--rna_exp", type="character", default=NULL,
              help="Gene expression data file [default= %default]", metavar="character"),
  make_option("--rna_var", type="character", default=NULL,
              help="Variant expression evidence file (mRNA)[default= %default]", metavar="character"),
  make_option("--db", type="character", default=NULL,
              help="protein reference database[default= %default]", metavar="character"),
  make_option(c("-o","--out_dir"), type="character", default="./", 
              help="Output folder [default= %default]", metavar="./")
  
  
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

## 
for(para in names(opt)){
  cat(" --",para,": ",opt[[para]],"\n",sep="")
}


if (is.null(opt[['mhc']])){
  print_help(opt_parser)
  stop("No valid input for --mhc\n")
}

## read MHC prediction result
res <- read_csv(opt[['mhc']])

if(!is.null(opt[['db']])){
  ref_db = opt[['db']]
  res %>% select(Neoepitope) %>% write_tsv(path = "pep.txt",col_names = FALSE)
  out_pep_map = "pep2pro_ref.txt";
  system(paste("java -jar /opt/pepmap.jar -i pep.txt -d ",ref_db," -o ",out_pep_map,sep=" "))
  pep_map_res <- read_tsv(out_pep_map)
  cat("The number of peptides mapped to reference proteins:",nrow(pep_map_res),"\n")
  res <- res %>% filter(!(Neoepitope %in% pep_map_res$peptide))
  
}

## read variant evidence from RNA-Seq data
if(!is.null(opt[['rna_var']])){
  rna_var_evidence <- read_tsv(opt[['rna_var']],na = "#_#_#")
  res <- add_rna_var_evidence(res,rna_var_evidence)
}else{
  cat("No variant evidence information from RNA-Seq data!\n")
}

## read gene expression data
if(!is.null(opt[['rna_exp']])){
  rna_exp_data <- read_tsv(opt[['rna_exp']])
  res <- add_gene_exp_evidence(res,rna_exp_data)
}else{
  cat("No gene expression evidence information from RNA-Seq data!\n")
}

## read variant evidence from proteomics data
if(!is.null(opt[['pro_exp']])){
  pro_var_data <- read_tsv(opt[['pro_var']])
  res <- add_pro_var_evidence(res,pro_var_data)
}else{
  cat("No protein expression evidence information from RNA-Seq data!\n")
}

## read protein expression data
if(!is.null(opt[['pro_exp']])){
  pro_exp_evidence <- read_tsv(opt[['pro_exp']])
  res <- add_pro_exp_evidence(res,pro_exp_evidence)
}else{
  cat("No protein expression evidence information from RNA-Seq data!\n")
}

outdir <- opt[["out_dir"]]
sample_id <- opt[["sample_id"]]
all_res <- paste(outdir,"/",sample_id,"-all_neoantigens.tsv",sep="")
write_tsv(res,path = all_res)

################################################################################
## output final result
## filter neoantigen based on binding affinity
mhc_filter <- opt[["mhc_filter"]]
res_filtered <- res %>% filter(netMHCpan_binding_affinity_nM <= mhc_filter)
res_filtered_file <- paste(outdir,"/",sample_id,"-filtered_neoantigens.tsv",sep="")
write_tsv(res_filtered,path = res_filtered_file)


################################################################################
## output summary information
cat("Summary information:\n")
output_stat <- list()
output_stat$neoantigen_somatic_mutations <- res_filtered %>% 
  select(Chr,Start,End,Ref,Alt) %>% 
  distinct() %>% nrow()
output_stat$neoantigen_somatic_mutations_with_rna_var <- res_filtered %>% 
  filter(rna_var_evidence==1) %>% 
  select(Chr,Start,End,Ref,Alt) %>% 
  distinct() %>% nrow()
output_stat$neoantigen_somatic_mutations_with_pro_var <- res_filtered %>% 
  filter(protein_var_evidence==1) %>% 
  select(Chr,Start,End,Ref,Alt) %>% 
  distinct() %>% nrow()
output_stat$neoantigen_somatic_mutations_with_rna_pro_var <- res_filtered %>% 
  filter(rna_var_evidence==1,protein_var_evidence==1) %>% 
  select(Chr,Start,End,Ref,Alt) %>% 
  distinct() %>% nrow()
output_stat$neoantigen <- res_filtered %>% 
  select(Neoepitope) %>% 
  distinct() %>% nrow()
output_stat$neoantigen_with_rna_var <- res_filtered %>% 
  filter(rna_var_evidence==1) %>% 
  select(Neoepitope) %>% 
  distinct() %>% nrow()
output_stat$neoantigen_with_pro_var <- res_filtered %>% 
  filter(protein_var_evidence==1) %>% 
  select(Neoepitope) %>% 
  distinct() %>% nrow()
output_stat$neoantigen_with_rna_pro_var <- res_filtered %>% 
  filter(rna_var_evidence==1,protein_var_evidence==1) %>% 
  select(Neoepitope) %>% 
  distinct() %>% nrow()

res_stat <- data.frame(item=names(output_stat),value=unlist(output_stat))
row.names(res_stat) <- NULL
print(res_stat)
res_stat_file <- paste(outdir,"/",sample_id,"-stat.tsv",sep="")
write_tsv(res_stat,path = res_stat_file)


