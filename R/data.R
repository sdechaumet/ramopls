#' liver.toxicity sample data from R
#'
#' Data available from R mixOmics
#'
#' @format A list with three levels:
#' \describe{
#'   \item{gene}{Observation x variables}
#'   \item{clinic}{Observation x variables}
#'   \item{treatment}{Observation x variables}
#'   \item{gene.ID}{Observation x variables}
#' }
"liver.toxicity"

#' Dataset used in Ruiz et al. (2017)
#'
#' Data used in the article from Ruiz et al. 2017 for the calcul of VIP²
#' There is 127 identified compounds with two duplicated names (L-ARGININE and O-ACETYL-L-CARNITINE) whereas there are 126
#' compounds selected in the article. Since no evidence were available to choose the 126 compounds from the 127,
#' the duplicated ones were renamed as CPD_NAME_1 and CPD_NAME_2 and included in this dataset?
#'
#' @format a list of two level with
#' \describe{
#'   \item{datamatrix}{a data.table with measurment results:  Observation x variables (18 x 127)}
#'   \item{samplemetadata}{a data.table with sample groups:  Observation x metadata (18 x 3)}
#' }
"data_Ruiz2017"

#' Minimal dataset for testing
#'
#' Minimal dataset created from Ruiz et al. 2017 to run test and exemples.
#'
#' @format a list of two level with
#' \describe{
#'   \item{datamatrix}{a data.table with measurment results:  Observation x variables}
#'   \item{samplemetadata}{a data.table with sample groups:  Observation x metadata}
#' }
"data_exemple"

#' Dataset used in Boccard et al. (2016)
#'
#' Data used in the article from Boccard et al. 2016 for the implementation of AMOPLS
#' and available at http://www.metaboanalyst.ca/resources/data/cress_time.csv
#'
#' \describe{
#'   \item{datamatric}{sample x variables: 48 observation described by 1468 variables}
#'   \item{samplemetadata}{sample x metadata: class }
#'   \item{variablemetadata}{variable x metadata: information about variables}
#' }
"data_Boccard2016"

#' Dataset used in Boccard et al. (2019)
#'
#' Data used in the article from Boccard et al. 2019 for the implementation of
#' stratified subsampling.
#'
#' \describe{
#'   \item{datamatric}{sample x variables: 49 observation described by 131 variables}
#'   \item{samplemetadata}{sample x metadata: Chemical (7 levels) and Batch (1 to 3)}
#' }
"data_Boccard2019"

