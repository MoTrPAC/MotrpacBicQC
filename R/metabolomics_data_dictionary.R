
# All about the RefMet data dictionary


#' @title Get and validate the entire RefMet database from Metabolomics Workbench
#'
#' @description Get and validate Metabolomics Data Dictionary from
#' Metabolomics Workbench
#' @param remove_duplications (logical) if `TRUE``, removes duplications.
#' @return (vector) PHASE code
#' @export

# get_and_validate_mdd <- function(remove_duplications = FALSE){
#   
#   refmet <- MotrpacBicQC::metabolomics_data_dictionary
# 
#   if(remove_duplications){
#     # REMOVE DUPLICATIONS
#     if( any(duplicated(refmet$refmet_name)) ){
#       duplirefmets <- length(refmet$refmet_name[(duplicated(refmet$refmet_name))])
#       # message("WARNING: [ ",duplirefmets, " ] DUPLICATION(s) found in data dictionary!")
#       refmet <- refmet[!(duplicated(refmet$refmet_name)),]
#     }
#   }
#   return(refmet)
# }

get_and_validate_mdd <- function(remove_duplications = FALSE){

  .id = name = NULL

  # REST metabolomics workbench data dictionary
  # Previous REST version (motrpac only)
  # refmetjson <- jsonlite::fromJSON("https://www.metabolomicsworkbench.org/rest/refmet/motrpac")
  refmetjson <- jsonlite::fromJSON("https://www.metabolomicsworkbench.org/rest/refmet/all/")

  dt_list <- purrr::map(refmetjson, as.data.table)
  dt <- data.table::rbindlist(dt_list, fill = TRUE, idcol = T)
  df <- as.data.frame(dt)

  colnames(df) <- tolower(colnames(df))

  if( !("name" %in% colnames(df)) ){
    stop("<refmet_name> column not found in the Metabolomics Workbench data dictionary")
  }else{
    df <- rename(df, refmet_name = name)
  }

  # CHECK DUPLICATIONS
  if(any(duplicated(df$refmet_name))){
    duplications <- df[duplicated(df$refmet_name),]
    message("Duplicated ids: ", length(duplications$refmet_name))
    message("IDS: ", paste(duplications$refmet_name, collapse = ", "))
    message("DUPLICATIONS IN REFMET ONLINE (REST VERSION)")
  }

  refmet <- subset(df, select = -c(.id))

  if(remove_duplications){
    # REMOVE DUPLICATIONS
    if( any(duplicated(refmet$refmet_name)) ){
      duplirefmets <- length(refmet$refmet_name[(duplicated(refmet$refmet_name))])
      # message("WARNING: [ ",duplirefmets, " ] DUPLICATION(s) found in data dictionary!")
      refmet <- refmet[!(duplicated(refmet$refmet_name)),]
    }
  }
  return(refmet)
}

#' @title Validate refmet_name
#'
#' @description Validate the refment_name using the Metabolomics Workbench API
#' @param dataf (data.frame) Data frame with the refmet_name column
#' @param verbose (logical) `TRUE` (default) shows messages)
#' @return (numeric) number of refmet_name ids not available in RefMet
#' @export
validate_refmetname <- function(dataf, verbose){
  
  # dataf <- data.frame(refmet_name=c("Valine-d8 [iSTD]", "Phenylalanine-d8 [iSTD]", "2-amino-6-hydroxyhexanoic acid", "O-Ethylhomoserine", "isoquinoline-1,5-diol", "Tridecylamine", "N-Acetylisoputreanine", "Piperine", "p-Methylhippuric acid", "4-Hydroxyhippuric acid", "Palmitoleoyl-EA", "Stearamide", "Tetradecylamine", "N-Acetylglycine", "N-Acetylvaline", "N-Acetylglutamic acid", "Lys-Gly", "Glu-Gly", "Asp-Pro", "Ile-Glu", "Diethanolamine", "Triethanolamine", "Tributylamine", "N,N-Dicyclohexylamine", "Myo-inositol", "Heme", "Protoporphyrin", "Pyridoxal", "13-cis-Retinoic acid", "CAR 2:0", "CAR 3:0", "CAR 4:0", "CAR 5:1", "CAR 5:0", "CAR 4:0;OH", "CAR 6:0", "CAR DC3:0;2Me", "CAR 7:0", "CAR 6:0;OH", "CAR 8:0", "CAR DC6:0", "CAR 8:0;OH", "CAR 10:1", "CAR 10:0", "CAR 12:1", "CAR 12:0", "CAR 12:0;OH", "CAR 14:2", "CAR 14:1", "CAR 14:0", "CAR 16:2", "CAR 16:0", "CAR 16:0;OH", "CAR 18:2", "CAR 18:1", "CAR 18:0", "CAR 18:1;OH", "CAR 18:0;OH", "CAR 20:4", "CAR 26:0", "LPC 14:0", "LPC 15:1", "LPC 15:0", "LPC 16:1", "LPC 16:0", "LPC 18:3", "LPC 18:2", "LPC 18:1", "LPC 18:0", "LPC 19:0", "LPC 20:5", "LPC 20:4", "LPC 20:3_hp_a", "LPC 20:3_hp_b", "LPC 20:1", "LPC 20:0", "LPC 22:6_hp_a", "LPC 22:6_hp_b", "LPC 22:5", "LPC 22:4_hp_a", "LPC 22:4_hp_b", "LPC 22:0", "LPC 24:0", "LPC P-16:0 or LPC O-16:1", "LPC O-16:0", "LPC P-18:0 or LPC O-18:1_hp_a", "LPC P-18:0 or LPC O-18:1_hp_b", "LPC O-18:0", "PC O-32:0", "PC P-34:2 or PC O-34:3", "PC P-34:0 or PC O-34:1", "PC P-36:4 or PC O-36:5", "PC P-36:3 or PC O-36:4", "PC P-35:1 or PC O-35:2", "PC P-38:6 or PC O-38:7", "PC P-38:5 or PC O-38:6", "PC P-38:4 or PC O-38:5", "PC 30:0", "PC 31:1", "PC 32:2", "PC 32:0", "PC 33:1", "PC 33:0", "PC 34:4", "PC 34:3", "PC 35:4", "PC 35:2", "PC 36:5", "PC 36:2", "LPE 16:0", "LPE 17:0", "LPE 18:3", "LPE 18:2_hp_a", "LPE 18:2_hp_b", "LPE 18:1", "LPE 20:4_hp_a", "LPE 20:4_hp_b", "LPE 20:1", "LPE 20:0", "LPE 22:6", "LPE P-18:0 or LPE O-18:1", "PE P-34:2 or PE O-34:3", "PE P-34:1 or PE O-34:2", "PE P-34:0 or PE O-34:1", "PE P-36:4 or PE O-36:5", "PE P-36:2 or PE O-36:3", "PE P-36:1 or PE O-36:2", "PE P-38:6 or PE O-38:7", "PE P-38:5 or PE O-38:6", "PE P-38:4 or PE O-38:5", "PE P-40:6 or PE O-40:7", "PE P-40:4 or PE O-40:5", "PE 34:2", "PE 34:0", "PE 36:4", "PE 36:2", "PE 37:4", "PE 38:6", "PE 38:4", "PE 40:6", "PS P-36:1 or PS O-36:2", "Cer 18:1;O2/15:0", "Cer 18:1;O2/16:1", "Cer 18:1;O2/16:0", "Cer 18:1;O2/20:0", "Cer 18:1;O2/22:0", "Cer 18:1;O2/24:1", "GlcCer 18:2;O2/22:0", "C20 Sphinganine", "C16 Sphingosine", "C20 Sphingosine", "C22 Sphingosine", "SM 18:1;O2/14:0", "SM 18:1;O2/16:1", "SM 18:1;O2/16:0", "SM 18:1;O2/15:0", "SM 18:1;O2/18:3", "SM 18:1;O2/18:1", "SM 36:1;O2", "SM 18:1;O2/19:0", "SM 18:1;O2/20:1", "SM 18:1;O2/20:0", "SM 18:1;O2/22:1", "SM 18:1;O2/22:0", "SM 18:1;O2/23:2", "DG 36:4", "DG 36:2", "DG 38:5", "TG 46:2" ))
  irm <- 0
  if((dim(dataf)[1] > 20) & (verbose)){
    pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent",
                                     total = dim(dataf)[1],
                                     complete = "=",   # Completion bar character
                                     incomplete = "-", # Incomplete bar character
                                     current = ">",    # Current bar character
                                     clear = FALSE,    # If TRUE, clears the bar when finish
                                     width = 100)      # Width of the progress bar  
  }
  
  for(i in 1:dim(dataf)[1]){
    rn <- dataf$refmet_name[i]
    
    if((dim(dataf)[1] > 20) & (verbose)){
      pb$tick()
    }
    # Remove the isotope label
    rn <- "LPC 20:3_hp_a"
    tosearch <- "_hp_||_rp_||_rn_||_in_||_lp_||_ln_"
    if(grepl(tosearch, rn)){
      rn <- gsub("(.*)(_\\w{2}_\\w{1})", "\\1", rn)
    }
    
    search_api <- paste0("https://www.metabolomicsworkbench.org/rest/refmet/match/",URLencode(rn),"/name/")
    here <- jsonlite::fromJSON(search_api)
    if(length(here) == 7){
      value <- here$refmet_name[1]
      if(value != rn){
        if(verbose) message(paste0("      (-) refmet_name [", rn, "] not available in RefMet. Please, contact MW/BIC (Error RN2)"))
        irm <- irm + 1
      }
    }else{
      if(verbose) message(paste0("      (-) refmet_name [", rn, "] not available in RefMet. Please, contact MW/BIC (Error RN1)"))
      irm <- irm + 1
    }
  }
  return(irm)
}
