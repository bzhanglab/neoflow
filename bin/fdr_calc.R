library(PGA)
library(dplyr)
library(stringr)
library(readr)

# search_engine: comet, xtandem or msgf
pre_processing=function(id_file_dir, search_engine="xtandem"){
    all_data <- NULL
    
    if(search_engine == "xtandem" || search_engine == "msgf"){
        input_files <- list.files(id_file_dir,pattern = ".mzid",full.names = TRUE,include.dirs = TRUE)
        
        dat_res <- lapply(input_files, function(file){
            sample_name <- sub(x=file,pattern=".mzid$",replacement="")
            parserGear(file=file,
                       db=database,
                       decoyPrefix="XXX_",
                       novelPrefix="VAR",
                       prefix=sample_name,
                       xmx=10,
                       thread=8,
                       alignment=0,
                       outdir=output_dir)
            data <- read.delim(paste(output_dir,"/", sample_name, "-rawPSMs.txt", sep=""),stringsAsFactors = FALSE)
            data <- data %>% mutate(index = as.numeric(index) + 1)
            data$index <- paste(sample_name, data$index, sep=":")
            colnames(data)[2] <- "score"
            return(data)
        })
        
        all_data <- dplyr::bind_rows(dat_res)

    }else if(search_engine == "comet"){
        input_files <- list.files(id_file_dir, pattern = "_rawResults.txt$", full.names = TRUE, include.dirs = TRUE) 
        dat_res <- lapply(input_files, function(file){
            sample_name <- sub(x=file,pattern="_rawResults.txt$",replacement="")
            comet_result <- read.delim(file, stringsAsFactors = FALSE,check.names = FALSE) %>%
                dplyr::select(scan, xcorr, charge, exp_neutral_mass, plain_peptide, protein, modifications)
            names(comet_result) <- c("index", "score", "charge", "mass", "peptide", "protein", "mods")
            comet_result$protein <- str_replace_all(comet_result$protein, ",", ";")
            comet_result$mods <- str_replace_all(comet_result$mods, ",", ";")
            comet_result <- comet_result %>% mutate(index = as.numeric(index))
            comet_result$index <- paste(sample_name, comet_result$index, sep=":")
            return(comet_result)
        })
        
        all_data <- dplyr::bind_rows(dat_res)
    }else{
        stop(paste("Search engine: ", search_engine, " is not supported!\n",sep = ""))
    }
    return(all_data)
}



args <- commandArgs(TRUE)

id_file_dir   <- args[1]
database      <- args[2]
batchID       <- args[3]
search_engine <- args[4]
output_dir    <- args[5]


all_data <- pre_processing(id_file_dir, search_engine = search_engine)
psm_file <- paste(output_dir, "/", batchID, "-rawPSMs.txt", sep="")
write_tsv(all_data, psm_file)

calculateFDR(psmfile=psm_file,
        db=database,
        fdr=0.01,
        decoyPrefix="XXX_",
        better_score_lower=FALSE,
        remap=FALSE,
        peptide_level=TRUE,
        score_t = 0,
        protein_inference=FALSE,
        out_dir=paste(output_dir, "/", "peptide_level/global_fdr", sep=""),
        xmx=10)

calculateFDR(psmfile=paste(output_dir, "/", batchID, "-rawPSMs.txt", sep=""),
        db=database,
        fdr=0.01,
        decoyPrefix="XXX_",
        better_score_lower=FALSE,
        remap=FALSE,
        peptide_level=FALSE,
        score_t = 0,
        protein_inference=FALSE,
        out_dir=paste(output_dir, "/", "psm_level/global_fdr", sep=""),
        xmx=10)
