
# Datasets
#______________________________________________________________________________


#' Motrpac assay codes
#'
#' @format A data frame with assay codes
#' \describe{
#'   \item{omics_code}{Main OMIC group}
#'   \item{submission_code}{Assay submission code}
#'   \item{assay_code}{Assay release code}
#'   \item{assay_name}{Assay long name}
#'   \item{cas_code}{MoTrPAC chemical analysis site}
#' }
#' @examples
#' \dontrun{
#'  assay_codes
#' }
"assay_codes"


#' BIC Tissue Code and Name
#'
#' @format A data frame with tissue codes
#' \describe{
#'   \item{bic_tissue_code}{`T` + integer, as encoded in vial_label by DMAQC}
#'   \item{bic_tissue_name}{Tissue name (as provided by DMAQC)}
#'   \item{motrpac_tissue_code}{Both code and tissue name}
#'   \item{tissue_name_release}{Code + tissue name as used in releases}
#'   \item{abbreviation}{Tissue abbreviation}
#'   \item{tissue_hex_colour}{Tissue Hex Colour code}
#'   
#' }
#' @examples
#' \dontrun{
#'  bic_animal_tissue_code
#' }
"bic_animal_tissue_code"

#' Results test dataset for NAMED metabolites
#'
#' @format A data frame with 58 rows and 98 columns:
#' \describe{
#'   \item{metabolite_name}{metabolite names}
#'   \item{90001013104}{Sample id. Unit: Peak Area}
#'   \item{90005013104}{sample id. Unit: Peak Area}
#'   \item{90007013104}{sample id. Unit: Peak Area}
#'   \item{90008013104}{sample id. Unit: Peak Area}
#'   \item{90009013104}{sample id. Unit: Peak Area}
#'   \item{90010013104}{sample id. Unit: Peak Area}
#'   \item{90011013104}{sample id. Unit: Peak Area}
#'   \item{90012013104}{sample id. Unit: Peak Area}
#'   \item{90013013104}{sample id. Unit: Peak Area}
#'   \item{90014013104}{sample id. Unit: Peak Area}
#'   \item{90015013104}{sample id. Unit: Peak Area}
#'   \item{90017013104}{sample id. Unit: Peak Area}
#'   \item{90018013104}{sample id. Unit: Peak Area}
#'   \item{90020013104}{sample id. Unit: Peak Area}
#'   \item{90023013104}{sample id. Unit: Peak Area}
#'   \item{90025013104}{sample id. Unit: Peak Area}
#'   \item{90026013104}{sample id. Unit: Peak Area}
#'   \item{90027013104}{sample id. Unit: Peak Area}
#'   \item{90028013104}{sample id. Unit: Peak Area}
#'   \item{90029013104}{sample id. Unit: Peak Area}
#'   \item{90031013104}{sample id. Unit: Peak Area}
#'   \item{90032013104}{sample id. Unit: Peak Area}
#'   \item{90033013104}{sample id. Unit: Peak Area}
#'   \item{90034013104}{sample id. Unit: Peak Area}
#'   \item{90037013104}{sample id. Unit: Peak Area}
#'   \item{90038013104}{sample id. Unit: Peak Area}
#'   \item{90039013104}{sample id. Unit: Peak Area}
#'   \item{90040013104}{sample id. Unit: Peak Area}
#'   \item{90042013104}{sample id. Unit: Peak Area}
#'   \item{90043013104}{sample id. Unit: Peak Area}
#'   \item{90044013104}{sample id. Unit: Peak Area}
#'   \item{90045013104}{sample id. Unit: Peak Area}
#'   \item{90046013104}{sample id. Unit: Peak Area}
#'   \item{90047013104}{sample id. Unit: Peak Area}
#'   \item{90048013104}{sample id. Unit: Peak Area}
#'   \item{90050013104}{sample id. Unit: Peak Area}
#'   \item{90052013104}{sample id. Unit: Peak Area}
#'   \item{90053013104}{sample id. Unit: Peak Area}
#'   \item{90109013104}{sample id. Unit: Peak Area}
#'   \item{90110013104}{sample id. Unit: Peak Area}
#'   \item{90112013104}{sample id. Unit: Peak Area}
#'   \item{90114013104}{sample id. Unit: Peak Area}
#'   \item{90115013104}{sample id. Unit: Peak Area}
#'   \item{90117013104}{sample id. Unit: Peak Area}
#'   \item{90118013104}{sample id. Unit: Peak Area}
#'   \item{90119013104}{sample id. Unit: Peak Area}
#'   \item{90120013104}{sample id. Unit: Peak Area}
#'   \item{90121013104}{sample id. Unit: Peak Area}
#'   \item{90122013104}{sample id. Unit: Peak Area}
#'   \item{90123013104}{sample id. Unit: Peak Area}
#'   \item{90124013104}{sample id. Unit: Peak Area}
#'   \item{90126013104}{sample id. Unit: Peak Area}
#'   \item{90127013104}{sample id. Unit: Peak Area}
#'   \item{90128013104}{sample id. Unit: Peak Area}
#'   \item{90129013104}{sample id. Unit: Peak Area}
#'   \item{90133013104}{sample id. Unit: Peak Area}
#'   \item{90134013104}{sample id. Unit: Peak Area}
#'   \item{90135013104}{sample id. Unit: Peak Area}
#'   \item{90136013104}{sample id. Unit: Peak Area}
#'   \item{90138013104}{sample id. Unit: Peak Area}
#'   \item{90139013104}{sample id. Unit: Peak Area}
#'   \item{90140013104}{sample id. Unit: Peak Area}
#'   \item{90141013104}{sample id. Unit: Peak Area}
#'   \item{90142013104}{sample id. Unit: Peak Area}
#'   \item{90143013104}{sample id. Unit: Peak Area}
#'   \item{90144013104}{sample id. Unit: Peak Area}
#'   \item{90145013104}{sample id. Unit: Peak Area}
#'   \item{90146013104}{sample id. Unit: Peak Area}
#'   \item{90147013104}{sample id. Unit: Peak Area}
#'   \item{90148013104}{sample id. Unit: Peak Area}
#'   \item{90150013104}{sample id. Unit: Peak Area}
#'   \item{90152013104}{sample id. Unit: Peak Area}
#'   \item{90155013104}{sample id. Unit: Peak Area}
#'   \item{90156013104}{sample id. Unit: Peak Area}
#'   \item{90157013104}{sample id. Unit: Peak Area}
#'   \item{90159013104}{sample id. Unit: Peak Area}
#'   \item{90160013104}{sample id. Unit: Peak Area}
#'   \item{CS00000MP-10}{sample id. Unit: Peak Area}
#'   \item{CS00000MP-11}{sample id. Unit: Peak Area}
#'   \item{CS00000MP-12}{sample id. Unit: Peak Area}
#'   \item{CS00000MP-2}{sample id. Unit: Peak Area}
#'   \item{CS00000MP-3}{sample id. Unit: Peak Area}
#'   \item{CS00000MP-4}{sample id. Unit: Peak Area}
#'   \item{CS00000MP-5}{sample id. Unit: Peak Area}
#'   \item{CS00000MP-6}{sample id. Unit: Peak Area}
#'   \item{CS00000MP-7}{sample id. Unit: Peak Area}
#'   \item{CS00000MP-8}{sample id. Unit: Peak Area}
#'   \item{CS00000MP-9}{sample id. Unit: Peak Area}
#'   \item{CS000QCMP-1}{sample id. Unit: Peak Area}
#'   \item{CS000QCMP-2}{sample id. Unit: Peak Area}
#'   \item{CS000QCMP-4}{sample id. Unit: Peak Area}
#'   \item{CS000QCMP-5}{sample id. Unit: Peak Area}
#'   \item{CS000QCMP-6}{sample id. Unit: Peak Area}
#'   \item{CSMR80014-1}{sample id. Unit: Peak Area}
#'   \item{CSMR80014-2}{sample id. Unit: Peak Area}
#'   \item{CSMR80015-1}{sample id. Unit: Peak Area}
#'   \item{CSMR80015-2}{sample id. Unit: Peak Area}
#' }
#' @examples
#' \dontrun{
#'  results_named
#' }
"results_named"

#' Results test dataset for UNNAMED metabolites
#'
#' @format A data frame with 58 rows and 98 columns:
#' \describe{
#'   \item{metabolite_name}{metabolite names based on mz_rt}
#'   \item{90001013104}{Sample id. Unit: Peak Area}
#'   \item{90005013104}{Sample id. Unit: Peak Area}
#'   \item{90007013104}{Sample id. Unit: Peak Area}
#'   \item{90008013104}{Sample id. Unit: Peak Area}
#'   \item{90009013104}{Sample id. Unit: Peak Area}
#'   \item{90010013104}{Sample id. Unit: Peak Area}
#'   \item{90011013104}{Sample id. Unit: Peak Area}
#'   \item{90012013104}{Sample id. Unit: Peak Area}
#'   \item{90013013104}{Sample id. Unit: Peak Area}
#'   \item{90014013104}{Sample id. Unit: Peak Area}
#'   \item{90015013104}{Sample id. Unit: Peak Area}
#'   \item{90017013104}{Sample id. Unit: Peak Area}
#'   \item{90018013104}{Sample id. Unit: Peak Area}
#'   \item{90020013104}{Sample id. Unit: Peak Area}
#'   \item{90023013104}{Sample id. Unit: Peak Area}
#'   \item{90025013104}{Sample id. Unit: Peak Area}
#'   \item{90026013104}{Sample id. Unit: Peak Area}
#'   \item{90027013104}{Sample id. Unit: Peak Area}
#'   \item{90028013104}{Sample id. Unit: Peak Area}
#'   \item{90029013104}{Sample id. Unit: Peak Area}
#'   \item{90031013104}{Sample id. Unit: Peak Area}
#'   \item{90032013104}{Sample id. Unit: Peak Area}
#'   \item{90033013104}{Sample id. Unit: Peak Area}
#'   \item{90034013104}{Sample id. Unit: Peak Area}
#'   \item{90037013104}{Sample id. Unit: Peak Area}
#'   \item{90038013104}{Sample id. Unit: Peak Area}
#'   \item{90039013104}{Sample id. Unit: Peak Area}
#'   \item{90040013104}{Sample id. Unit: Peak Area}
#'   \item{90042013104}{Sample id. Unit: Peak Area}
#'   \item{90043013104}{Sample id. Unit: Peak Area}
#'   \item{90044013104}{Sample id. Unit: Peak Area}
#'   \item{90045013104}{Sample id. Unit: Peak Area}
#'   \item{90046013104}{Sample id. Unit: Peak Area}
#'   \item{90047013104}{Sample id. Unit: Peak Area}
#'   \item{90048013104}{Sample id. Unit: Peak Area}
#'   \item{90050013104}{Sample id. Unit: Peak Area}
#'   \item{90052013104}{Sample id. Unit: Peak Area}
#'   \item{90053013104}{Sample id. Unit: Peak Area}
#'   \item{90109013104}{Sample id. Unit: Peak Area}
#'   \item{90110013104}{Sample id. Unit: Peak Area}
#'   \item{90112013104}{Sample id. Unit: Peak Area}
#'   \item{90114013104}{Sample id. Unit: Peak Area}
#'   \item{90115013104}{Sample id. Unit: Peak Area}
#'   \item{90117013104}{Sample id. Unit: Peak Area}
#'   \item{90118013104}{Sample id. Unit: Peak Area}
#'   \item{90119013104}{Sample id. Unit: Peak Area}
#'   \item{90120013104}{Sample id. Unit: Peak Area}
#'   \item{90121013104}{Sample id. Unit: Peak Area}
#'   \item{90122013104}{Sample id. Unit: Peak Area}
#'   \item{90123013104}{Sample id. Unit: Peak Area}
#'   \item{90124013104}{Sample id. Unit: Peak Area}
#'   \item{90126013104}{Sample id. Unit: Peak Area}
#'   \item{90127013104}{Sample id. Unit: Peak Area}
#'   \item{90128013104}{Sample id. Unit: Peak Area}
#'   \item{90129013104}{Sample id. Unit: Peak Area}
#'   \item{90133013104}{Sample id. Unit: Peak Area}
#'   \item{90134013104}{Sample id. Unit: Peak Area}
#'   \item{90135013104}{Sample id. Unit: Peak Area}
#'   \item{90136013104}{Sample id. Unit: Peak Area}
#'   \item{90138013104}{Sample id. Unit: Peak Area}
#'   \item{90139013104}{Sample id. Unit: Peak Area}
#'   \item{90140013104}{Sample id. Unit: Peak Area}
#'   \item{90141013104}{Sample id. Unit: Peak Area}
#'   \item{90142013104}{Sample id. Unit: Peak Area}
#'   \item{90143013104}{Sample id. Unit: Peak Area}
#'   \item{90144013104}{Sample id. Unit: Peak Area}
#'   \item{90145013104}{Sample id. Unit: Peak Area}
#'   \item{90146013104}{Sample id. Unit: Peak Area}
#'   \item{90147013104}{Sample id. Unit: Peak Area}
#'   \item{90148013104}{Sample id. Unit: Peak Area}
#'   \item{90150013104}{Sample id. Unit: Peak Area}
#'   \item{90152013104}{Sample id. Unit: Peak Area}
#'   \item{90155013104}{Sample id. Unit: Peak Area}
#'   \item{90156013104}{Sample id. Unit: Peak Area}
#'   \item{90157013104}{Sample id. Unit: Peak Area}
#'   \item{90159013104}{Sample id. Unit: Peak Area}
#'   \item{90160013104}{Sample id. Unit: Peak Area}
#'   \item{CS00000MP-10}{Sample id. Unit: Peak Area}
#'   \item{CS00000MP-11}{Sample id. Unit: Peak Area}
#'   \item{CS00000MP-12}{Sample id. Unit: Peak Area}
#'   \item{CS00000MP-2}{Sample id. Unit: Peak Area}
#'   \item{CS00000MP-3}{Sample id. Unit: Peak Area}
#'   \item{CS00000MP-4}{Sample id. Unit: Peak Area}
#'   \item{CS00000MP-5}{Sample id. Unit: Peak Area}
#'   \item{CS00000MP-6}{Sample id. Unit: Peak Area}
#'   \item{CS00000MP-7}{Sample id. Unit: Peak Area}
#'   \item{CS00000MP-8}{Sample id. Unit: Peak Area}
#'   \item{CS00000MP-9}{Sample id. Unit: Peak Area}
#'   \item{CS000QCMP-1}{Sample id. Unit: Peak Area}
#'   \item{CS000QCMP-2}{Sample id. Unit: Peak Area}
#'   \item{CS000QCMP-4}{Sample id. Unit: Peak Area}
#'   \item{CS000QCMP-5}{Sample id. Unit: Peak Area}
#'   \item{CS000QCMP-6}{Sample id. Unit: Peak Area}
#'   \item{CSMR80014-1}{Sample id. Unit: Peak Area}
#'   \item{CSMR80014-2}{Sample id. Unit: Peak Area}
#'   \item{CSMR80015-1}{Sample id. Unit: Peak Area}
#'   \item{CSMR80015-2}{Sample id. Unit: Peak Area}
#' }
#' @examples
#' \dontrun{
#'  results_unnamed
#' }
"results_unnamed"

#' Metabolomics data dictionary
#'
#' @format A data frame with > 2000 rows and 15 columns (variables):
#' \describe{
#'   \item{CURRENT_REFMET_NAME}{Updated version of RefMet (after December 2020)}
#'   \item{refmet_name}{Old version of refmet_name}
#'   \item{metabolite_name}{Metabolite name provided by each laboratory}
#'   \item{is_standard}{is a reference standard? 1: yes, 0: no}
#'   \item{super_class}{Metabolite super class}
#'   \item{main_class}{Metabolite main class}
#'   \item{sub_class}{Metabolite sub class}
#'   \item{formula}{Metabolite formula}
#'   \item{exactmass}{Metaboliteexact mass}
#'   \item{pubchem_cid}{Metabolite pubchem id}
#'   \item{kegg_id}{Metabolite kegg id}
#'   \item{inchi_key}{Metabolite inchi key}
#'   \item{lm_id}{Metabolite lm id}
#'   \item{hmdb_id}{Metabolite hmdb id}
#'   \item{chebi_id}{Metabolite chebi id}
#' }
#' @examples
#' \dontrun{
#'  metabolomics_data_dictionary
#' }
"metabolomics_data_dictionary"

#' Metadata Metabolites test dataset for NAMED metabolites
#'
#' @format A data frame with 53 rows and 6 columns:
#' \describe{
#'   \item{metabolite_name}{metabolite names}
#'   \item{refmet_name}{RefMet id (metabolomics workbench)}
#'   \item{mz}{mass over charge}
#'   \item{rt}{retention time}
#'   \item{formula}{chemical formula}
#'   \item{neutral_mass}{neutral mass}
#' }
#' @examples
#' \dontrun{
#'  data(metadata_metabolites_named)
#' }
"metadata_metabolites_named"

#' Metadata Metabolites test dataset for UNNAMED metabolites
#'
#' @format A data frame with 53 rows and 4 columns:
#' \describe{
#'   \item{metabolite_name}{metabolite names based on mz_rt}
#'   \item{mz}{mass over charge}
#'   \item{rt}{retention time}
#'   \item{neutral_mass}{neutral mass}
#' }
#' @examples
#' \dontrun{
#'  data(metadata_metabolites_unnamed)
#' }
"metadata_metabolites_unnamed"


#' Metadata Sample test dataset for NAMED metabolites
#'
#' @format A data frame with 97 rows and 4 columns:
#' \describe{
#'   \item{sample_id}{metabolite names based on mz_rt}
#'   \item{sample_type}{mass over charge}
#'   \item{sample_order}{retention time}
#'   \item{raw_file}{neutral mass}
#' }
#' @examples
#' \dontrun{
#'  data(metadata_sample_named)
#' }
"metadata_sample_named"

#' Metadata Sample test dataset for UNNAMED metabolites
#'
#' @format A data frame with 97 rows and 4 columns:
#' \describe{
#'   \item{sample_id}{metabolite names based on mz_rt}
#'   \item{sample_type}{mass over charge}
#'   \item{sample_order}{retention time}
#'   \item{raw_file}{neutral mass}
#' }
#' @examples
#' \dontrun{
#'  metadata_sample_unnamed
#' }
"metadata_sample_unnamed"

#' Motrpac Phenotypes PASS1A 6 Months
#'
#' @format A data frame with 8616 rows and 6 columns:
#' \describe{
#'   \item{tissue_code}{MoTrPAC Tissue code}
#'   \item{vial_label}{Vial label id (only available for MoTrPAC Samples)}
#'   \item{tissue_name}{Tissue name}
#'   \item{group_time_point}{Group time point}
#'   \item{sex}{Animal sex}
#'   \item{site_code}{cas}
#' }
#' @examples
#' \dontrun{
#'  phenotypes_pass1a06_short
#' }
"phenotypes_pass1a06_short"

#' assay abbreviations
"assay_abbr"

#' assay order
"assay_order"

#' group abbreviations
"group_abbr"

#' group colors
"group_cols"

#' sex abbreviations
"sex_abbr"

#' sex colors
"sex_cols"

#' tissue abbreviations
"tissue_abbr"

#' tissue colors
"tissue_cols"

#' tissue order
"tissue_order"


