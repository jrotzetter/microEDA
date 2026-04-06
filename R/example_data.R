#' SRS014476-Supragingival_plaque data
#'
#' Example data of a single MetaPhlAn4 output file in `metaphlanProfile`
#' format, containing the computed taxon abundances.
#'
#' @format ## `single_abundance_profile`
#' An object of class \linkS4class{metaphlanProfile} with slots `mpa_version` and `mpa_taxProfile`:
#' \describe{
#'   \item{mpa_version}{A `Character` vector indicating the MetaPhlAn database version used.}
#'   \item{mpa_taxProfile}{A data frame with 20 rows and 2 columns:
#'
#'   `clade_name`: The taxonomic lineage of the taxon reported on this row.
#'   Taxon names are prefixed with one-letter codes to help indicate their rank.
#'
#'   `single_sample_profile`: Each taxon's relative abundance in %.}
#' }
#' @source <https://github.com/biobakery/MetaPhlAn/releases/tag/4.1.0>
"single_metaphlan_profile"


#' Merged MetaPhlAn4 profiles data
#'
#' Example data of multiple MetaPhlAn4 profiles that were merged into one file
#' in `metaphlanProfile` format.
#'
#' @format ## `merged_abundance_profiles`
#' An object of class \linkS4class{metaphlanProfile} with slots `mpa_version` and `mpa_taxProfile`:
#' \describe{
#'   \item{mpa_version}{A `Character` vector indicating the MetaPhlAn database version used.}
#'   \item{mpa_taxProfile}{A data frame with 59 rows and 7 columns:
#'
#'   `clade_name`: The taxonomic lineage of the taxon reported on this row.
#'   Taxon names are prefixed with one-letter codes to help indicate their rank.
#'
#'   `SRS0144...`: The sample names for which sample-specific abundance
#'    profiles are found in a column.}
#' }
#' @source <https://github.com/biobakery/MetaPhlAn/releases/tag/4.1.0>
"merged_metaphlan_profiles"
