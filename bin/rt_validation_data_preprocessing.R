library(tidyverse)


is_novel_peptides=function(protein_ids){
    ids <- str_split(protein_ids,pattern=";") %>% unlist()
    is_novel <- all(str_detect(ids,pattern="^VAR"))
    return(is_novel)
}

format_modifications=function(peptide, modification, search_engine="comet"){
    
    ## Only consider Oxidation on M. Update this function to support 
    ## more modifications.
    
    ## X!Tandem or MS-GF+
    ## GSCYPCPETVDVK
    ## example: 57.02147@3C;57.02147@6C
    rt_valid_pep = TRUE
    aa <- str_split(peptide,pattern = "") %>% unlist()
    if(search_engine == "xtandem" || search_engine == "msgf"){
        if(str_detect(modification,pattern = "@")){
            mods = str_split(modification,pattern = ";") %>% unlist()
            for(mod in mods){
                ## 
                mod_info <- str_split(mod,pattern = "@") %>% unlist
                if(str_detect(mod_info[1],pattern = "^57\\.") && str_detect(mod_info[2],pattern = "C$")){
                    ## 57 on C        
                }else if(str_detect(mod_info[1],pattern = "^15\\.9") && str_detect(mod_info[2],pattern = "M$")){
                    pos <- str_replace_all(mod_info[2],pattern = "M",replacement = "") %>% as.numeric()
                    aa[pos] <- "1"
                }else{
                    rt_valid_pep = FALSE
                    cat(paste("Modification:", mod, " is not supported!\n"))
                }
            }
        }
        
    }else if(search_engine == "comet"){
        ## http://comet-ms.sourceforge.net/parameters/parameters_201801/output_txtfile.php
        ## 11_S_57.021464;16_S_57.021464
        ## 6_V_15.994900
        ## Comet
        ## VCLLGCGISTGYGAAVNTAK
        ## example: 2_S_57.021464;6_S_57.021464
        
        if(str_detect(modification,pattern = "_")){
            mods = str_split(modification,pattern = ";") %>% unlist()
            for(mod in mods){
                ## 
                mod_info <- str_split(mod,pattern = "_") %>% unlist
                if(str_detect(mod_info[3],pattern = "^57\\.") && str_detect(aa[as.numeric(mod_info[1])],pattern = "C")){
                    ## 57 on C        
                }else if(str_detect(mod_info[3],pattern = "^15\\.9") && str_detect(aa[as.numeric(mod_info[1])],pattern = "M")){
                    pos <- mod_info[1] %>% as.numeric()
                    aa[pos] <- "1"
                }else{
                    rt_valid_pep = FALSE
                    cat(paste("Modification:", mod, " is not supported!\n"))
                }
            }
        }
	}
    
    
    format_pep = paste(aa,collapse = "")
    res <- data.frame(peptide = peptide,
                      mods=modification,
                      format_peptide=format_pep,
                      rt_valid=rt_valid_pep,
                      stringsAsFactors = FALSE)
    return(res)
}

load_rt_data=function(rt_dir="./"){
    rt_files <- list.files(path = rt_dir, pattern = ".txt",full.names = TRUE,include.dirs = TRUE)
    rt_dat <- lapply(rt_files, function(x){
        rt <- read.delim(x,stringsAsFactors = FALSE) %>% mutate(fraction=basename(x) %>% str_replace_all(pattern = ".mgf_index.txt$",replacement = ""))
        return(rt)
    })
    
    res <- dplyr::bind_rows(rt_dat) %>% mutate(rt=rt/60.0) %>%
        rename(msms_index=index,msms_title=Title)
    return(res)
}


prepare_training_testing_data=function(psm_file, pep_file, rt_file_dir = "./",
                                       novel_prefix="VAR",
                                       out_dir="./",search_engine="comet"){
    
    ## 
    rt_fractions_name = NULL
    
    ## peptide level FDR result
    peptide_level_data <- read_tsv(pep_file)
    all_peps <- peptide_level_data %>% separate(index, c("fraction", "msms_index"),sep=":",remove=FALSE) %>%
        mutate(fraction=basename(fraction)) %>%
        filter(Qvalue<=0.01) %>% 
        mutate(is_novel=sapply(protein,is_novel_peptides)) %>%
        filter(isdecoy==FALSE) %>%
        dplyr::select(-rt)
    
    ## PSM level FDR result
    psm_level_data <- read_tsv(psm_file)
    all_psms <- psm_level_data %>% separate(index, c("fraction", "msms_index"),sep=":",remove=FALSE) %>%
        mutate(fraction=basename(fraction)) %>%
        filter(Qvalue<=0.01) %>% 
        mutate(is_novel=sapply(protein,is_novel_peptides)) %>%
        filter(isdecoy==FALSE) %>%
        dplyr::select(-rt)
    
    all_psms <- all_psms %>% filter(peptide %in% all_peps$peptide)
    
    ref_psms   <- all_psms %>% filter(is_novel!=TRUE)
    novel_psms <- all_psms %>% filter(is_novel==TRUE)
    
    if(nrow(novel_psms) <= 0){
        cat("No novel peptide identification!\n")
    }else{
        
        rt_data <- load_rt_data(rt_file_dir)
        all_psms_rt <- merge(all_psms,rt_data, by=c("fraction","msms_index"))
        if(nrow(all_psms_rt) != nrow(all_psms)){
            stop("RT data is not complete!")
        }
        ## format modification for all peptides
        format_pep_data <- lapply(1:nrow(all_psms_rt), function(i){
            res <- format_modifications(all_psms_rt$peptide[i],all_psms_rt$mods[i],search_engine = search_engine)
            return(res)
        })
        format_pep_data <- dplyr::bind_rows(format_pep_data) %>% distinct()
        all_psm_rt_tmp <- merge(all_psms_rt,format_pep_data,by=c("peptide","mods")) 
        
        train_data <- all_psm_rt_tmp %>% filter(is_novel==FALSE) %>% 
            dplyr::filter(rt_valid == TRUE) %>% 
            dplyr::select(fraction,format_peptide,rt) %>%
            group_by(fraction,format_peptide) %>% 
            dplyr::summarise(rt_range=max(rt)-min(rt),rt_mean=mean(rt)) %>%
            ungroup()
        
        ## remove PSMs with rt range > 3 minutes
        peptides_remove <- train_data %>% filter(rt_range > 3)
        cat("Remove peptides with large rt range:", nrow(peptides_remove),"\n")
        test_data <- all_psm_rt_tmp %>% filter(is_novel==TRUE) %>% 
            dplyr::filter(rt_valid == TRUE) 
        ## Only need to perform training and prediction for those runs who have novel identifications 
        rt_fractions_name <- unique(test_data$fraction)
        cat("The number of runs with identified novel peptides:",length(rt_fractions_name),"\n")
        ## save training data
        train_data_dir <- paste(out_dir,"/train_data",sep="")
        dir.create(train_data_dir,showWarnings = FALSE)
        for(fn in rt_fractions_name){
            fd <- train_data %>% filter(fraction == fn) %>%
                dplyr::select(format_peptide,rt_mean) %>% 
                dplyr::rename(x=format_peptide,y=rt_mean)
            write_tsv(fd,path = paste(train_data_dir,"/",fn,"_train.tsv",sep=""))
            
            fd_test <- test_data %>% filter(fraction == fn) %>%
                dplyr::select(format_peptide,rt) %>% 
                dplyr::rename(x=format_peptide,y=rt)
            write_tsv(fd_test,path = paste(train_data_dir,"/",fn,"_test.tsv",sep=""))
        }
    }
}

cmdargs <- commandArgs(TRUE)
psm_file <- cmdargs[1]
pep_file <- cmdargs[2]
rt_file_dir <- cmdargs[3]
search_engine <- cmdargs[4]
out_dir <- cmdargs[5]

if(!dir.exists(out_dir)){
    dir.create(out_dir,showWarnings = FALSE,recursive = TRUE) 
}

prepare_training_testing_data(psm_file = psm_file, 
                              pep_file = pep_file, 
                              rt_file_dir = rt_file_dir, 
                              search_engine = search_engine,
                              out_dir = out_dir)

