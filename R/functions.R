utils::globalVariables(c("."))


.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
      paste0(
        "To use this package, kopls functionalities must be installed.\n",
        "Please install `kopls` by one of the following method:\n",
        "  - Use the internal function: `rAMOPLS::install_kopls()`\n",
        "  - Download the original source code from  `http://kopls.sourceforge.net/download.shtml` and compile it manually with `devtools::install()`"
      )
    )
}

#' Install the kopls package included in rAMOPLS
#'
#' This function try to install the kopls package from rAMOPLS.
#' It needs Rtools and devtools to run
#'
#' @param ... Argument passed to install
#'
#' @examples
#' install_kopls()
#'
#' @export
install_kopls <- function(...) {
  ## Check Rtools env
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Package \"devtools\" needed for this function to work, please install it using install.packages('devtools').",
         call. = FALSE)
  }

  Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
  Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

  ## Unzip package in temporary folder
  temp_dir_path <- tempdir()
  zip::unzip(file.path(system.file("package", package = "rAMOPLS"), "kopls.zip"), exdir = temp_dir_path)
  ## Install from unzipped temporary folder
  devtools::install(temp_dir_path, quick = T, ...)
  ## Remove temporary folder
  unlink(temp_dir_path)


  if (!requireNamespace("kopls")) {
    warning("kopls installation from local file has failed. You can install it directly from the original authors: \n
            http://kopls.sourceforge.net/download.shtml")
  } else {message("kopls was successfully installed")}
}


#' Return factor names and interaction
#'
#' Return the names of all the studied factors (+ interactions) of the given dataframe.
#'
#' @param data_factors  Dataframe to study
#' @param studied_factors String of the studied factors indexes according to their column number
#'
#' @examples
#' M <- data.frame(matrix(nrow = 3, ncol = 3, 1:9))
#' colnames(M) <- c('Dose','Time','Age')
#' fun_factor_names(M,'1,2,23')
#'
#' @export
fun_factor_names <- function(data_factors, studied_factors) {
  factor_names <- c()
  s <- studied_factors %>% {
    strsplit(as.character(.), ",")
  } %>% unlist() %>% as.numeric() %>% as.list()
  for (l in s) {
    if (nchar(l) == 1) {
      factor_names <- c(factor_names, colnames(data_factors)[l])
    }
    if (nchar(l) == 2) {
      factor1_name <- colnames(data_factors)[as.numeric(substr(l, 1, 1))]
      factor2_name <-
        colnames(data_factors)[as.numeric(substr(l, 2, 2))]
      factor_int_name <- paste(factor1_name, 'x', factor2_name)
      factor_names <- c(factor_names, factor_int_name)
    }
    if (nchar(l) > 2) {
      stop("only 2 factors interaction")
    }
  }
  return(factor_names)
}


#' Load and format data
#'
#' @param datamatrix Datamatrix as a matrix, data.frame or data.table
#' @param samplemetadata Metadata on samples
#' @param factor_names Column(s) name(s) from samplemetadata to use
#'
#' @return a list of two object : dataset and factors
#' @export
#'
fun_load_data <- function(datamatrix, samplemetadata, factor_names) {
  Data <- list()
  Data$dataset <- switch(class(datamatrix)[[1]],
                         "matrix" = {datamatrix},
                         "data.frame" = {
                           ## Check if first column is numeric
                           if (is.numeric(datamatrix[, 1][[1]])) {
                             # Return as matrix
                             temp <- as.matrix(datamatrix)
                             temp
                           } else {
                             temp <- as.matrix(datamatrix[, -1])
                             ## add rownames
                             rownames(temp) <- datamatrix[, 1][[1]]
                             temp
                           }
                         },
                         "data.table" = {
                           ## Check if first column is numeric
                           if (is.numeric(datamatrix[, 1][[1]])) {
                             # Return as matrix
                             temp <- as.matrix(datamatrix)
                             temp
                           } else {
                             temp <- as.matrix(datamatrix[, -1])
                             ## add rownames
                             rownames(temp) <- datamatrix[, 1][[1]]
                             temp
                           }
                         },
                         stop("datamatrix must be a matrix, a data.frame or a data.table")
  )

  if (data.table(Data$dataset) %>% .[, lapply(.SD, function(x) {any(is.na(x))})] %>% {any(. == T)}) {
    stop("NA are not allowed in the datamatrix, check and treat the NA before AMOPLS analysis")
  }

  ## Format samplemetadata
  ### Search sampleid in samplemetadata
  samplemetadata <- as.data.table(samplemetadata, keep.rownames = T)
  sampleid_col <- samplemetadata %>% .[, lapply(.SD, function(x) {all(x %in% rownames(Data$dataset))})] %>% {which(. == T)}
  if (length(sampleid_col) == 0) {
    stop("No column in samplemetadata corresponds to the sample IDs.")
  }

  ## Check factor_names
  if (!any(factor_names %in% colnames(samplemetadata))) {
    stop("Some factor_names are not present in samplemetadata, check samplemetadata column names")
  }

  Data$factors <- samplemetadata[, factor_names, with = F] %>% as.matrix(., rownames = samplemetadata[, sampleid_col, with = F][[1]])

  return(Data)
}



#' Get Row Repeats
#'
#' Return the unique row patterns from a given matrix and the indices of the corresponding repeated rows.
#'
#' @param mat Matrix to study
#'
#' @return \item{result}{The single row patterns and the lists of the corresponding indices for each pattern}
#'
#' @examples
#' M <- matrix(nrow = 3, ncol = 3, 1:9)
#' colnames(M) <- c('Dose','Time','Age')
#' rAMOPLS:::fun_GetRowRepeats(M)
#'
fun_GetRowRepeats <- function(mat) {
  if (!is.matrix(mat)) {
    mat <- matrix(mat)
  }
  result <- list()
  no.rows <- dim(mat)[1]
  no.cols <- dim(mat)[2]
  result$row.patterns <- matrix(nrow = 0, ncol = no.cols)
  no.patterns.found <- 0
  result$indices.per.pattern <- list()

  fun_IsIdenticalIgnoreNames <- function(x, y) {
    # a new function to check for identical matrices was needed, one that ignores column and row names
    x.nameless <- c(x)
    y.nameless <- c(y)
    if (length(x.nameless) != length(y.nameless)) {
      return(FALSE)
    }
    for (i in 1:length(x.nameless)) {
      if (x.nameless[i] != y.nameless[i]) {
        return(FALSE)
      }
    }
    return(TRUE)
  }

  for (r in 1:no.rows) {
    # go through all the rows in the input matrix, and check whether that row-pattern was already discovered before
    pattern.number <- which(apply(result$row.patterns, 1, fun_IsIdenticalIgnoreNames, y = mat[r, ]) ==
                              TRUE) # which pattern does this row match, if any
    if (length(pattern.number) == 0) {
      # if the row does not match a previous pattern, then add it to the list of patterns
      result$row.patterns <- rbind(result$row.patterns, mat[r, ])
      no.patterns.found <- no.patterns.found + 1
      result$indices.per.pattern[[no.patterns.found]] <- r
    } else {
      # the row does match a previous pattern, therefore remember the index of this row as an occurence of that pattern in the matrix
      result$indices.per.pattern[[pattern.number]] <-
        c(result$indices.per.pattern[[pattern.number]], r)
    }
  }
  factors.order <- sort(as.numeric(apply(result$row.patterns, 1, paste, collapse = "")), index.return = TRUE)$ix # sort the patterns by numerical order
  result$indices.per.pattern <- result$indices.per.pattern[factors.order]
  result$row.patterns <- result$row.patterns[factors.order, , drop = FALSE]

  return(result)
}


#' Get Equation Element
#'
#' Internal function : Return the mean average of the dataset on every level-combination of the given factor.
#'
#' @param model Dataset matrix to study
#' @param evaluation String : index name of the factors to study - ex:'1,2,12'
#' @param previous.model Use of a previous model independant from the data - NULL by default
#'
#' @return \code{s} List of results containing :
#' @return \item{level.combinations}{List containing all the levels combinations of the factor and the corresponding row indices}
#' @return \item{means.matrix}{Means matrix of the dataset on all the level combinations related to the selected factor (or interactions)}
#'
fun_GetEquationElement <- function(model, evaluation, previous.model) {
  s <- list()
  if (!is.null(previous.model)) {
    s$level.combinations <- previous.model[[paste(evaluation, collapse = "")]]$level.combinations
  } else {
    s$level.combinations <- fun_GetRowRepeats(model$general$factors[, evaluation, drop = FALSE])
  }
  s$means.matrix <- matrix(nrow = dim(model$general$data)[1],
                           ncol = dim(model$general$data)[2])
  for (p in 1:dim(s$level.combinations$row.patterns)[1]) {
    mean.for.this.level.combination <-
      colMeans(model$general$data[s$level.combinations$indices.per.pattern[[p]], , drop = FALSE])
    for (i in s$level.combinations$indices.per.pattern[[p]]) {
      s$means.matrix[i, ] <- mean.for.this.level.combination
    }
  }
  return(s)
}



#' Function to undersample a dataset
#'
#' This function will randomly take n observation in each unique groups of Data$factors
#' with n equal to the smallest subgroup.
#'
#' @inheritParams fun_AMOPLS
#'
#' @return a list with subsampled datamatrix and their corresponding factors
fun_balance_data <- function(Data) {
  Factors <- as.matrix(Data$factors[, colnames(Data$factors)])
  Factors_id <- as.matrix(unique(Factors)) # name of the associated factors

  Pattern_indices <- list()
  Length_indices  <- list()

  for (l in 1:dim(Factors_id)[1]) {
    # l <- 2
    line_indices <- apply(Factors, 1, identical, Factors_id[l, ]) %>% which()
    Pattern_indices[[l]] <- line_indices
    Length_indices[[l]] <- length(line_indices)
  }

  m <- min(unlist(Length_indices))
  Data_balanced <- matrix(m, nrow = 0, ncol = dim(Data$data)[2])
  Factors_balanced <- matrix(m, nrow = 0, ncol = dim(Data$factors)[2])

  if (all(Length_indices == m)) {
    message('Data already balanced')
  } else {
    for (l in 1:dim(Factors_id)[1]) {
      Pattern_indices[[l]] <- sample(Pattern_indices[[l]], m)
      Data_balanced <- rbind(Data_balanced, Data$data[Pattern_indices[[l]], ])
      Factors_balanced <- rbind(Factors_balanced, Data$factors[Pattern_indices[[l]], ])
    }
    Data$dataset <- Data_balanced
    Data$factors <- Factors_balanced
  }

  return(Data)
}

#' Check if data are balanced in selected factors
#'
#' @inheritParams run_AMOPLS
#' @inheritParams fun_AMOPLS
#' @return TRUE if the data are balanced, FALSE otherwise.
#' @export
fun_is_balanced <- function(Data, factor_names, interaction_level){
  # Data <- Data
  # factor_names <- factor_names
  # interaction_level <- 1
  factor_index <- which(factor_names %in% colnames(Data$factors))

  ## Check for each factors and their interaction if asked
  ## The number of samples by subgroups
  names_fac <- list()
  result <- sapply(0:interaction_level, function(int) {
    utils::combn(factor_index, int + 1, simplify = F)
    # Data$factors[, as.vector(temp)]
  }) %>%
    unlist(., recursive = F) %>% {
      lapply(1:length(.), function(y) {
        x <- .[[y]]
        col_sel <- colnames(Data$factors)[x]
        names_fac[[y]] <<- paste(col_sel, collapse = " x ")
        temp <- as.data.table(Data$factors)[, .(sple_nb = .N), by = col_sel]
        length(unique(temp$sple_nb)) == 1
      })
    }

  if (any(result == F)) {
    message("Data are unbalanced in:")
    message(paste(names_fac[!unlist(result)], collapse = "\n"))
    return(F)
  } else {
    return(T)
  }
}



#' Data pre-processing
#'
#' Data pre-processing step before the ANOVA decomposition.
#'
#' @inheritParams run_AMOPLS
#' @inheritParams fun_AMOPLS
#'
#' @return \code{s$general} List of results containing general information about the pre-processed dataset :
#' @return \item{Nb_compo_ortho}{Number of orthogonal components}
#' @return \item{studied_factors}{String of the indices of the studied factors}
#' @return \item{equation.elements}{List of numeric indices of the studied factors}
#' @return \item{order.to.evaluate.ee}{List of numeric indices corresponding to the order of evaluation of the factors under study}
#' @return \item{data}{Pre-processed dataset to use in the model}
#' @return \item{ssq.mean}{Mean sum of square of the dataset}
#' @return \item{ssq}{Sum of square of the dataset}
#' @return \item{factors}{Pure factors dataset}
#' @return \item{factor_names}{List of the names of the factors under study (interactions included)}
#' @return \item{factors_data}{Pure factors + interactions dataset}
#' @return \item{PCA}{Principal component analysis on the main dataset}
fun_pre_processing <- function(Data,
                               Nb_compo_ortho,
                               equation.elements,
                               scaling,
                               only.means.matrix = FALSE) {
  s <- list()
  s$general <- list("Nb_compo_ortho" = Nb_compo_ortho,
                    "studied_factors" = equation.elements)
  # Data processing : reduced - centered :
  dataAdjusted <- MetStaT.ScalePip(Data$dataset,
                                   center = TRUE,
                                   scale = scaling,
                                   quietly = TRUE)
  s$general$data <- dataAdjusted$data

  # Format factors
  if (!is.numeric(Data$factors)) {
    # message("The supplied factors are not numeric. Converting levels to numeric values")
    temp <- apply(Data$factors, 2, function(x) {as.numeric(as.factor(x))})
    rownames(temp) <- rownames(Data$factors)
    colnames(temp) <- colnames(Data$factors)
    Data$factors <- temp
  }

  s$general$factors <- Data$factors

  if (is.character(equation.elements)) {
    equation.elements <- lapply(strsplit(strsplit(equation.elements, split = ",")[[1]], split =
                        ""), as.numeric)
  }

  for (ee in equation.elements) {
    for (f in ee)
      if (f > dim(Data$factors)[2] ||
          f < 1) {
        stop(paste("Factor ", f, " is beyond scope of study-design", sep = ""))
      }
  }
  if (nrow(Data$dataset) != nrow(Data$factors)) {
    stop(
      paste(
        "Number of rows in data (",
        dim(Data$dataset)[1],
        ") and study design (",
        dim(Data$factors)[1],
        ") do not match",
        sep = ""
      )
    )
  }

  ## Establisment of the order of study of the factor, by decreasing complexity : Main effects then interactions
  order.to.evaluate.ee <- sort(
    as.numeric(
      unlist(
        lapply(equation.elements, paste, collapse = "")
      )), index.return = TRUE)$ix

  s$general$equation.elements <- equation.elements
  s$general$order.to.evaluate.ee <- order.to.evaluate.ee



  s$general$ssq.mean <- sum(rep(dataAdjusted$center.vector / dataAdjusted$scale.vector ,nrow(Data$dataset)) ^ 2)
  s$general$ssq <- sum(Data$dataset ^ 2)


  s$general$factor_names <- fun_factor_names(s$general$factors, s$general$studied_factors)

  ## Clustering of data in one matrix for score plots
  # pure effect factors :
  s$general$factors_data <- fun_factors_data(s$general$factors, s$general$studied_factors)

  ## PCA on the dataset :
  if (!only.means.matrix) {
    s$general$PCA <- MetStat.PCA.Calculate(s$general$data)
  }

  return(s)
}


#' Return the factors_data slot
#'
#' @param factors The samplemetadata with observations groups
#' @param studied_factors The factors under study
#'
#' @return Return a data.table with observations groups for each factors under study and their interactions
fun_factors_data <- function(factors, studied_factors) {
  factor_names <- fun_factor_names(factors, studied_factors)
  ## Clustering of data in one matrix for score plots
  # pure effect factors :
  temp_data <- factors %>% data.table()
  temp_fac_list <- studied_factors %>% {strsplit(as.character(.), ",")[[1]]}

  temp_data <- mapply(function(x, z) {
    temp <- lapply(1:nchar(x), function(y){
      col_ind <- as.numeric(substr(x, y, y))
      factors[, col_ind]
    })
    output <- data.table(interaction(temp, sep = "."))
    setnames(output, z)

  }, temp_fac_list , factor_names, SIMPLIFY = F) %>%
    {Reduce(cbind, .)}
  return(temp_data)
}

#' @title Permutation setting function (optional)
#' @description Local permutation of the selected matrix, according to the selected element - cf J.Boccard et al. (2016)
#'
#' @param s List corresponding to the dataset to study
#' @param ee Internal parameter corresponding to the factor element under study
#' @param perm_t Element to permute
#'
#' @return \item{s$general$dataset}{Permuted dataset}
#' @export
fun_perm_settings <- function(s, ee, perm_t) {
  ee.name <- paste(s$general$equation.elements[[ee]],collapse="")
  ## If ee corresponds to a pure effect ('1', or '2' - size  <2 <-> pure effect) and to perm :
  if (ee.name == perm_t & nchar(perm_t)<2){
    Permterm <- colnames(s$general$factors)[as.integer(perm_t)] # Permterm : studied factor name - ex: Dose.Group
    Perm_factors <- s$general$factors[, -as.integer(perm_t), drop = F]
    # Perm_factors <- as.matrix(s$general$factors) %>% {.[grepl(Permterm, colnames(.)), drop = F]} # perte du nom de colonne
    factors_id <- unique(Perm_factors) # name of the associated factors

    Permuted_factor <- matrix(0, 1,dim(s$general$factors)[1])
    s$general$permuted_factors <- matrix(0, dim(s$general$factors)[1],dim(s$general$factors)[2])

    ## For each level non-associated to Permterm :
    for (l in 1:dim(factors_id)[1]){
      # l <- 1
      # line_indices <- apply(Perm_factors, 1, identical, factors_id[l, , drop = F]) %>% which()
      line_indices <- which(interaction(data.frame(Perm_factors)) %in% interaction(data.frame(factors_id)[l,]))
      randperm_indices <- sample(line_indices)
      # For each permuted indices dataset : permutation on the column corresponding to ee
      for (k in 1:length(line_indices)){
        # k <- 1
        Permuted_factor[line_indices[k]] <- s$general$factors[randperm_indices[k],Permterm]
        s$general$permuted_factors[line_indices[k],] <- s$general$factors[randperm_indices[k],]
      }
    }
    Permuted_factor <- t(Permuted_factor)
    s$general$factors[,Permterm] <- Permuted_factor
  }
  return(s)
}


#' @title ANOVA decomposition
#' @description ANOVA decomposition of the dataset according to the selected factors and their interactions.
#'
#' @param s List containing general information about the dataset and the factors - output of fun_pre_processing
#' @param only.means.matrix 'FALSE' by default
#' @param perm_t Term to permute if Perm = 'TRUE'
#'
#' @return \code{s$decompo_ANOVA} List of results of the ANOVA decomposition for each studied factor (including interactions + residuals) :
#' @return \item{residuals}{Residuals matrix}
#' @return \item{factor_index}{List of results related to the selected factor (interactions) :}
#' \itemize{
#'    \item \code{level.combinations} Main results on the level
#'    \item \code{means.matrix} Mean average matrix on every level combination of the selected factor
#'    \item \code{svd} Result of the Singular Value Decomposition on the means matrix - cf fun_svd_extraction
#'    \item \code{means.matrix_res} Mean average + residuals matrix on every level combination of the selected factor}
fun_decompo_ANOVA <- function(s,
                              only.means.matrix = FALSE,
                              perm_t = NULL) {
  # only.means.matrix <- F
  # perm_t <- "1"
  # Data <- list("dataset" = as.matrix(liver.toxicity$gene[, 1:20]),
  #              "factors" = as.matrix(data.table(liver.toxicity$treatment)[, .(Dose = Dose.Group, Time = Time.Group)]))
  # Nb_compo_ortho <- 1
  # equation.elements <- "1,2,12"
  # scaling <- F
  # s <-fun_pre_processing(Data = Data,
  #                        s = list(),
  #                        Nb_compo_ortho = Nb_compo_ortho,
  #                        equation.elements = equation.elements,
  #                        scaling = scaling,
  #                        only.means.matrix = only.means.matrix
  # )

  s$decompo_ANOVA <- list()
  s$decompo_ANOVA$residuals <- s$general$data # Residuals matrix

  ## For each factor (or interactions) i under study :
  temp_remainder <- s$general$data

  for (ee in s$general$order.to.evaluate.ee) {
    # ee <- s$general$order.to.evaluate.ee[[1]]
    ## Permutation on a pure effect factor or interaction - depending on perm_t
    if (!is.null(perm_t)) {s <- fun_perm_settings(s, ee, perm_t)}

    ## Calculation of the mean submatrix <X_i>, related to the average on the samples of every level of the i factor
    # <X_i> : in new.equation.element
    # reductions : for an interaction (ij) -> selection of the associated pure effects factors (i,j)
    new.equation.element <- fun_GetEquationElement(s, s$general$equation.elements[[ee]], previous.model = NULL) # $means.matrix in new.equation.element

    if (length(s$general$equation.elements[[ee]]) > 1) {
      for (r in s$general$equation.elements[[ee]]) {
        new.equation.element$means.matrix <- new.equation.element$means.matrix - s$decompo_ANOVA[[c(paste(r, collapse = ""))]]$means.matrix
      }
    }
    ## For an interaction (ij) : the pure effect submatrices are removed from <X_ij>
    # <X_int(i,j)> = <X_ij> - <X_i> - <X_j>
    # for (r in reductions) {
    #   new.equation.element$means.matrix <- new.equation.element$means.matrix - s$decompo_ANOVA[[c(paste(r, collapse = ""))]]$means.matrix
    # }

    if (nchar(s$general$equation.elements[ee]) > 1 & !is.null(perm_t)) {
      # print('interaction')
      new.equation.element$means.matrix <- new.equation.element$means.matrix[sample(1:nrow(new.equation.element$means.matrix)[1]), ]
    }

    ## Residual matrix :
    # For each factor i : <X_i> is removed from the dataset
    # s$decompo_ANOVA$residuals : residuals matrix
    #
    # For each factor i under study : Singular Value Decomposition (PCA) on <X_i> - result in s$new.equation.element$svd
    # Selection of the non-zero eigen values and vectors selected : fix or factor-dependent threshold

    if (!only.means.matrix) {
      temp_remainder <- temp_remainder - new.equation.element$means.matrix
      new.equation.element <- fun_svd_extraction(new.equation.element, threshold = 0.0001)
    }

    ## Results in s$'ee.name' - ee.name : name of the selected factor

    ee.name <- paste(s$general$equation.elements[[ee]], collapse = "")
    s$general$ee.names <- c(s$general$ee.names, ee.name)
    s$decompo_ANOVA[[ee.name]] <- new.equation.element
  }
  s$decompo_ANOVA$residuals <- temp_remainder
  s$general$ee.names <- c(s$general$ee.names, 'residuals')

  for (ee.name in s$general$ee.names[s$general$ee.names != 'residuals']) {
    s$decompo_ANOVA[[ee.name]]$means.matrix_res <- s$decompo_ANOVA[[ee.name]]$means.matrix + s$decompo_ANOVA$residuals
  }
  return(s)
}


#' @title SVD extraction
#' @description Singular Value Decomposition (PCA) extraction for the given matrix and selection of the non-zero eigen values and vectors.
#'
#' @param new.equation.element Matrix to study
#' @param threshold Threshold for the non-zero eigen values selection
#'
#' @return \item{svd}{List of results corresponding to :}
#' \itemize{
#'    \item \code{d, v, var.explained, t} Results from the SVD : cf PCA.Calculate
#'    \item \code{non_zero_eigen_vect} Non-zero singular eigen vectors from the SVD results
#'    \item \code{non_zero_eigen_val}  Non-zero singular eigen values from the SVD results}
#' @references From MetStaT : PCA.Calculate
#' @export
fun_svd_extraction <- function(new.equation.element, threshold) {
  new.equation.element$svd <- MetStat.PCA.Calculate(new.equation.element$means.matrix) # SVD sur la matrice moyenn?e correspondante (interaction ou effet pur)
  new.equation.element$svd$non_zero_eigen_vect <- as.matrix(new.equation.element$svd$t[, new.equation.element$svd$d > threshold])
  new.equation.element$svd$non_zero_eigen_val <- as.matrix(new.equation.element$svd$d[new.equation.element$svd$d > threshold])

  return(new.equation.element)
}


#' ANOVA PCA
#'
#' Principal Component Analysis on every residual-augmented experimental submatrix from ANOVA decomposition (interactions included).
#'
#' @inheritParams fun_outputs
#'
#' @return \item{s$ANOVA_PCA}{List of results including for each factor (+ interactions and residuals) the main results of the PCA}
#'
#' @references From MetStaT : PCA.Calculate
fun_ANOVA_PCA <- function(s) {
  s$ANOVA_PCA <- list()
  for (ee.name in s$general$ee.names[s$general$ee.names != 'residuals']) {
    s$ANOVA_PCA[[ee.name]]$pca <-
      MetStat.PCA.Calculate(s$decompo_ANOVA[[ee.name]]$means.matrix_res)
  }
  return(s)
}


#' Multiblock clustering
#'
#' Clustering of the ANOVA-decomposed experimental submatrices and the non-zero eigen vectors from the SVD analysis.
#'
#' @details The number of predictive components is imposed by the previous SVD - it corresponds to the number of non-zero eigen vectors.
#'
#' @inheritParams fun_outputs
#' @return \code{general$Nb_compo_pred} Number of predictive components
#' @return \code{Multiblock} List including all the results from the Multiblock clustering :
#' \itemize{
#'    \item \code{Wmat} Clustering of all the outcome mean matrices from ANOVA decomposition
#'    \item \code{AMat} AMat matrix for each factor (+ interactions and residuals) - AMat_i = t<X_i+res>*<X_i+res> (normalized - Frobenius norm)
#'    \item \code{Y} Clustering of all the non-zero eigenvectors from the SVD of the ANOVA-decomposed matrices}
fun_Multiblock <- function(s) {
  s$Multiblock$W_mat <- dim(s$general$data)[1] %>% {matrix(0, ncol = ., nrow = .)}
  s$Multiblock$Y <- NULL

  ## For each factor i under study :
  # Calculation of AMat_i = t<X_i+res>*<X_i+res> normalised - Forbenius norm ('F')
  # W_Mat : Addition of all the AMat_i matrices

  for (ee.name in s$general$ee.names[s$general$ee.names != 'residuals']) {
    # ee.name <- s$general$ee.names[s$general$ee.names != 'residuals'][[1]]
    s$Multiblock$AMat[[ee.name]] <- s$decompo_ANOVA[[ee.name]]$means.matrix_res %*% t(s$decompo_ANOVA[[ee.name]]$means.matrix_res) /
      norm(s$decompo_ANOVA[[ee.name]]$means.matrix_res %*% t(s$decompo_ANOVA[[ee.name]]$means.matrix_res), 'F')
    s$Multiblock$W_mat <- s$Multiblock$W_mat + s$Multiblock$AMat[[ee.name]]
  }

  s$Multiblock$AMat[['residuals']] <- s$decompo_ANOVA$residuals %*% t(s$decompo_ANOVA$residuals) /
    norm(s$decompo_ANOVA$residuals %*% t(s$decompo_ANOVA$residuals), 'F')
  s$Multiblock$W_mat <- s$Multiblock$W_mat + s$Multiblock$AMat[['residuals']]

  # Clustering of the non-zero eigen vectors into the Y matrix :
  for (ee.name in s$general$ee.names[s$general$ee.names != 'residuals']) {
    s$Multiblock$Y <- Reduce("cbind", list(s$Multiblock$Y, s$decompo_ANOVA[[ee.name]]$svd$non_zero_eigen_vect)
    ) # $svd$non_zero_eigen_vects already weighted in t
  }

  s$general$Nb_compo_pred <- ncol(s$Multiblock$Y)

  return(s)
}


#' K-OPLS training function
#'
#' Application of the K-OPLS model training function to the multiblock dataset from kopls package.
#'
#' @inheritParams koplsModel_custom
#'
#' @return Results from the K-OPLS model function
fun_kopls <- function(K, Y, A, nox) {
  for (j in 1:100) {
    # j <- 1
    tol <- 10^-(5*j)

    result <- tryCatch(
      koplsModel_custom(K = K,
                        Y = Y,
                        A = A,
                        nox = nox,
                        preProcK = "mc",
                        preProcY = "mc",
                        tol = tol),
      error = function (e) {e})

    if (inherits(result, "error")) {
      message("Collinerarity problem in solve function, setting tolerance to: ", tol)
      j <- j+1
    } else {
      j <- 100
    }
  }

  if (inherits(result, "error")) {
    message(result)
    return(NULL)
  } else {
    return(result)
  }
}



#' Plot R2Y and p-value for each orthogonal models
#'
#' @inheritParams fun_outputs
#'
#' @import ggplot2
#' @import data.table
#' @import magrittr
#'
#' @export
fun_plot_ortho <- function(s) {
  # s <- result
  ## check if there is multiple result
  `R2Y p-value` <- Ortho_nb <- R2Y <- R2Y_pval <- Iteration <- x <- y <- NULL
  if (!is.list(s)) {stop("Perform run_AMOPLS with multiple nb_compo_orthos")}
  lapply(s, function(x) {
    # x <- s[[1]]
    data.table("R2Y" = x$kOPLS$R2Yhat %>% {.[length(.)]},
               "R2Y_pval" = x$output$Permutation_result %>% {.[, sum(`R2Y p-value`) / .N]},
               "Ortho_nb" = x$general$Nb_compo_ortho,
               "Iteration" = x$output$Permutation_result$PermNb %>% unique())
  }) %>%
    rbindlist(use.names = T, idcol = "rn") %>% {
      ggplot(., aes(Ortho_nb, R2Y)) +
        geom_bar(stat = "identity", color = "black", fill = "grey") +
        theme_bw() +
        ylim(0,1) +
        geom_text(aes(label = formatC(R2Y_pval, digits = 3, format = "f")), vjust = -0.5) +
        labs(title = "R2Y and p-value for orthogonal component selection",
             subtitle = paste0("p-value calculated from ", .[, unique(Iteration)], " iterations"),
             x = "Number of orthogonal component",
             y = "R2Y")
    }
}

#' RSS score
#' Internal function : Calculation of the Relative Sum of Squares (RSS) + Sum of Squares (SSQ) scores.
#'
#' @inheritParams fun_outputs
#'
#' @return \code{SSQ} SSQ score table for each experimental submatrix (+ interactions and residuals)
#' @return \code{RSS} RSS score table for each factor (+ interactions and residuals)
#' @export
fun_Sum_of_Squares <- function(s) {
  s[['residuals']]$ssq <- sum(s$decompo_ANOVA$residuals ^ 2)
  s$ssq_tot <- sum(s$decompo_ANOVA$residuals ^ 2)

  for (ee.name in s$general$ee.names[s$general$ee.names != 'residuals']) {
    s[[ee.name]]$ssq <- sum(s$decompo_ANOVA[[ee.name]]$means.matrix ^ 2)
    s$ssq_tot <- s$ssq_tot + s[[ee.name]]$ssq
  }

  for (ee.name in s$general$ee.names) {
    s[[ee.name]]$RSS <- s[[ee.name]]$ssq / s$ssq_tot
  }

  RSS <- c()
  SSQ <- c()
  for (ee.name in s$general$ee.names) {
    RSS <- c(RSS, s[[ee.name]]$RSS)
    SSQ <- c(SSQ, s[[ee.name]]$ssq)
  }

  RSS <- data.table(t(RSS))
  names(RSS) <- s$general$ee.names
  # RSS <- t(matrix(RSS))

  SSQ <- data.table(t(SSQ))
  names(SSQ) <- s$general$ee.names

  return(list(SSQ, RSS))
}


#' Block saliences
#' Internal function : Calculation of the contribution of each factor on every component (predictive + orthogonal).
#'
#' @inheritParams fun_outputs
#'
#' @return  \code{block_saliences} Table of block saliences for each factor (row) and component (column) - raw
#' @return  \code{block_saliences_norm} Table of block saliences for each factor (row) and component (column) - normalized
#' @export
fun_block_saliences <- function(s) {
  i <- 0
  block_saliences <- matrix(0, nrow = length(s$general$ee.names), ncol = s$general$Nb_compo_pred + s$general$Nb_compo_ortho)

  for (ee.name in s$general$ee.names) {
    i <- i + 1
    for (d in 1:s$general$Nb_compo_pred) {
      block_saliences[i, d] <- t(as.matrix(s$kOPLS$T[, d])) %*% s$Multiblock$AMat[[ee.name]] %*% as.matrix(s$kOPLS$T[, d])
    }
    for (o in 1:s$general$Nb_compo_ortho) {
      block_saliences[i, s$general$Nb_compo_pred + o] <- t(as.matrix(s$kOPLS$To[, o])) %*% s$Multiblock$AMat[[ee.name]] %*% as.matrix(s$kOPLS$To[, o])
    }
  }
  block_saliences_norm <- block_saliences

  # Normalisation of the block_saliences :
  block_saliences_norm <- sapply(1:(s$general$Nb_compo_pred + s$general$Nb_compo_ortho), function(x) {
    block_saliences[, x] / sum(block_saliences[, x])
  })
  return(list(block_saliences, block_saliences_norm))
}



#' Most influent factor per component
#'
#' Internal function : Return for each component the index of most influent factor (corresponding to the maximum of block_saliences among all the factors + residuals).
#'
#' @inheritParams fun_outputs
#'
#' @return \code{Most_influent_factor} Table of indexes of the most influent factor for each component (predictive + orthogonal)
#'
#' @export
fun_Most_influent_factor <- function(s) {
  Most_influent_factor <- matrix(0, nrow = 1, ncol = s$general$Nb_compo_pred + s$general$Nb_compo_ortho)

  for (k in 1:(s$general$Nb_compo_pred + s$general$Nb_compo_ortho)) {
    Most_influent_factor[, k] <- which.max(s$outputs$block_saliences_norm[, k])
  }
  return(Most_influent_factor)
}


#' RSR score
#'
#' Internal function : Calculation of the Residual Structure Ratio (RSR) score.
#'
#' @inheritParams fun_outputs
#'
#' @return \code{RSR} Table of the RSR score for each factor (+ interactions and residuals)
#' @export
fun_RSR <- function(s) {
  s$outputs$block_saliences <- fun_block_saliences(s)[[1]]
  RSR <- sapply(1:length(s$general$ee.names), function(fact) {
    s$outputs$block_saliences[dim(s$outputs$block_saliences)[1], s$general$Nb_compo_pred + 1] / s$outputs$block_saliences[fact, s$general$Nb_compo_pred + 1]
  }) %>%
    t() %>%
    data.table()
  setnames(RSR, s$general$ee.names)
  return(RSR)
}


#' X-score calculation
#'
#' Internal function : Clustering of the x-scores from the kOPLS::kOPLSModel function (predictive and orthogonal).
#'
#' @inheritParams fun_outputs
#'
#' @return \code{x-scores} Matrix of the X-scores from the kOPLS model for every predictive and orthogonal component
#' @references From kopls : koplsModel
#' @export
fun_xscores <- function(s) {
  x_scores <- cbind(s$kOPLS$T, s$kOPLS$To)
  return(x_scores)
}

#' X-loadings
#'
#' Internal function : Calculation of the x-loadings corresponding to the contribution of each variable on the components.
#'
#' @inheritParams fun_outputs
#'
#' @return \code{x_loadings} Table corresponding to the X-loadings for each variable (row) and component (column) - predictive + orthogonal
#' @export
fun_xloadings <- function(s) {

  x_loadings <- matrix(0, nrow = dim(s$general$data)[2], ncol = s$general$Nb_compo_pred + s$general$Nb_compo_ortho)

  for (d in 1:s$general$Nb_compo_pred) {
    if (s$general$ee.names[s$outputs$Most_influent_factor[, d]] != 'residuals') {
      x_loadings[, d] <- t(s$decompo_ANOVA[[s$general$ee.names[s$outputs$Most_influent_factor[, d]]]]$means.matrix_res) %*% as.matrix(s$kOPLS$T[, d]) / as.numeric((t(s$kOPLS$T[, d]) %*% s$kOPLS$T[, d]))
    } else {
      x_loadings[, d] <- t(s$decompo_ANOVA$residuals) %*% as.matrix(s$kOPLS$T[, d]) / as.numeric((t(s$kOPLS$T[, d]) %*% s$kOPLS$T[, d]))
    }
  }

  for (o in 1:s$general$Nb_compo_ortho) {
    if (s$general$ee.names[s$outputs$Most_influent_factor[, o]] != 'residuals') {
      x_loadings[, s$general$Nb_compo_pred + o] <- t(s$decompo_ANOVA[[s$general$ee.names[s$outputs$Most_influent_factor[, o]]]]$means.matrix_res) %*% as.matrix(s$kOPLS$To[, o]) / as.numeric((t(s$kOPLS$T[, o]) %*% s$kOPLS$T[, o]))
    } else {
      x_loadings[, o] <- t(s$decompo_ANOVA$residuals) %*% as.matrix(s$kOPLS$T[, o]) / as.numeric((t(s$kOPLS$T[, o]) %*% s$kOPLS$T[, o]))
    }
  }
  return(x_loadings)
}


#' Y-loadings
#'
#' Internal function : Calculation of the loadings related to the Y matrix.
#'
#' @inheritParams fun_outputs
#'
#' @return \code{y_loadings} Y-loadings matrix defined as described in T. Mehmood et al.
#' @references T. Mehmood et al. (2012)
#' @export
fun_yloadings <- function(s) {

  y_loadings <- matrix(0, nrow = 1, ncol = dim(s$kOPLS$T)[2])

  for (i in 1:dim(s$kOPLS$T)[2]) {
    y_loadings[, i] <-
      t(s$Multiblock$Y[, i]) %*% s$kOPLS$T[, i] / as.numeric(t(s$kOPLS$T[, i]) %*% s$kOPLS$T[, i])

  }

  return(y_loadings)
}



#' SSa score
#'
#' Internal function : Calculation of the SSa score. SSa : variance of Y explained by the a-th component - in the calculation of the VIP formula
#'
#' @inheritParams fun_outputs
#'
#' @return \code{SSa} score for each component a
#' @return \code{var_explained} Variation explained by each component
#' @export
fun_SSa <- function(s) {
  SSa <- diag((t(s$kOPLS$T) %*% s$kOPLS$T) %*% (t(s$outputs$y_loadings) %*% s$outputs$y_loadings)) #%>% abs()
  return(list(SSa))
}


#' VIP score
#' Internal function : Calculation of the Variable
#'
#' @inheritParams fun_outputs
#'
#' @return \code{VIP} VIP Table scores for each variable (row) and factor (column) - interactions included
#' @references  T. Mehmood et al. (2012) - DOI : 188 (2012) 62-69
#' @export
fun_VIP <- function(s) {
  # ### Selection of the related components with Most_influent_factor - ex: factor Dose <-> tp2, tp4
  # ### Calculation of the SSa score for each of those components - ex: SS2, SS4
  #     For each variable j :
  #       Calculation on every a-component of the contribution of the variable of j on component a : ex: W2j, W4j
  #        formula : Waj = t(Xi_res) * Y (cf Mehmood et al. 2012)
  #        with
  #           Xi_res : i factor related matrix  (cf block_saliences calculation)
  #           Y : clustered eigen vectors matrix
  # For a factor a, for the variable j : VIP_a(j) = sum[a-components]{SSa*Waj}/sum[a-components]{SSa} - ex: VIP_a(j) = (SS2*W2j + SS4*W4j) / (SS2 +SS4)
  # s <- result$orthoNb_1
  names <- data.table(colnames(s$general$data))
  p_var <- nrow(names)

  X <- s$general$data

  ## new
  VIP_n <- lapply(1:(length(s$general$ee.names) - 1), function(i) {
    # i <- 2
    ## [Update 07/01/2020] Calculate only on predictive component
    ## Removed
    # pred_compos <- which(s$outputs$Most_influent_factor[, 1:s$general$Nb_compo_pred] %>% {.[-length(.)]} == i)
    pred_compos <- which(s$outputs$Most_influent_factor %>% {.[-length(.)]} == i)
    if (length(pred_compos) == 0) {
      VIP_term <- cbind(names, NA)
      setnames(VIP_term, c("id", s$general$factor_names[i]))
      return(VIP_term)
    }
    Term <- c()
    W <- c()
    Q <- c()

    for (p in pred_compos) {
      q <- t(s$kOPLS$Up) %*% s$kOPLS$T[, p] / as.numeric(t(s$kOPLS$T[, p]) %*% s$kOPLS$T[, p])
      u <- s$kOPLS$Up %*% q / as.numeric(t(q) %*% q)
      w <- t(X) %*% u / as.numeric(t(u) %*% u)
      w <- w / norm(w, '2')

      Term <- cbind(Term, s$kOPLS$T[, p])
      W <- cbind(W, w)
      Q <- cbind(Q, q)
    }

    Q <- t(as.matrix(Q))
    SS <- diag(t(Term) %*% Term %*% Q %*% t(Q))

    VIP <- lapply(1:p_var, function(j) {
      # j <- 1
      weight <- c()
      for (p in 1:length(pred_compos)) {
        weight <- cbind(weight, (W[j, p] / norm(as.matrix(W[, p]), '2')) ^ 2)
      }
      weight <- t(weight)
      q <- SS %*% weight
      # VIP_term <- sqrt(p_var * q / sum(SS)) ## Not needed since erased in the next step ?
      VIP_term <- p_var * q / sum(SS)
      # VIP_term <- data.table(names[j, 1], VIP_term)
      # setnames(VIP_term, c("id", s$general$factor_names[i]))
      return(VIP_term)
    }) %>%
      unlist() %>%
      {data.table("id" = names, .)}
    setnames(VIP, c("id", s$general$factor_names[i]))
    return(VIP)
  }) %>%
    {Reduce(function(z, w) {merge(z, w, by = "id", all = T)}, .)}
  VIP_n <- VIP_n[names$V1]
  return(VIP_n)
}


#' @title Outputs wrapper
#' @description Wrapper function to cluster all the outputs of the AMOPLS model.
#' @param s List containing all the information from the AMOPLS model
#' @return \code{s$outputs} Main outputs of AMOPLS corresponding to :
#' \itemize{
#'    \item \code{SSQ} Sum of squares scores
#'    \item \code{RSS} Relative sum of squares scores
#'    \item \code{block_saliences} Block saliences scores
#'    \item \code{block_saliences_norm} Normalized block saliences scores
#'    \item \code{Most_influent_factor} Most influent factor per component
#'    \item \code{RSR} Residual Structure Ratio score
#'    \item \code{x_scores} x-scores from K-OPLS
#'    \item \code{x_loadings} x-loadings from K-OPLS
#'    \item \code{y_loadings} y-loadings from K-OPLS
#'    \item \code{SSa} SSa scores
#'    \item \code{var_explained} Table of the variance explained by every component
#'    \item \code{VIP} VIP score}
fun_outputs <- function(s) {
  s$outputs$SSQ <- fun_Sum_of_Squares(s)[[1]]
  s$outputs$RSS <- fun_Sum_of_Squares(s)[[2]]
  s$outputs$block_saliences <- fun_block_saliences(s)[[1]]
  s$outputs$block_saliences_norm <- fun_block_saliences(s)[[2]]
  s$outputs$Most_influent_factor <- fun_Most_influent_factor(s)
  s$outputs$RSR <- fun_RSR(s)
  s$outputs$x_scores <- fun_xscores(s)
  s$outputs$x_loadings <- fun_xloadings(s)
  s$outputs$y_loadings <- fun_yloadings(s)
  s$outputs$SSa <- fun_SSa(s)
  s$outputs$VIP <- fun_VIP(s)
  s$outputs$R2Y <- s$kOPLS$R2Yhat[2]
  return(s)
}


#' @title AMOPLS wrapper
#' @description Wrapper function to process all the steps for the AMOPLS model.
#'
#' @param Data List of 2 numeric matrices - Data$dataset : raw data; Data$factors : factors matrix
#' @param equation.elements String with column indices containing factors and interactions to study; ex: "1,2,12"
#' @param scaling Should scaling be performed : 'TRUE' or 'FALSE'
#' @param only.means.matrix Should the means matrix only be returned : 'TRUE' or 'FALSE'
#' @param use.previous.model Should a previous model be used :'TRUE' or 'FALSE'
#' @param Nb_compo_ortho Number of orthogonal component
#' @param perm_t ... to permute
#'
#' @return \code{s} List containing all the information about the AMOPLS model, organized in 6 groups :
#' \itemize{
#'     \item \code{general} General information about the parameters
#'     \item \code{decompo_ANOVA} Outcomes from the ANOVA decomposition : experimental submatrices (+ residuals) and svd
#'     \item \code{ANOVA_PCA} Outcomes from the ANOVA-PCA (for ANOVA-PCA model)
#'     \item \code{Multiblock} Outcomes from the Multiblock_clustering function : multiblock X and Y-matrices
#'     \item \code{kOPLS} Outcomes from the kOPLS Model - cf kopls::koplsModels
#'     \item \code{outcomes} Main outcomes from AMOPLS - cf fun_outputs for details}
fun_AMOPLS <- function(Data,
                       equation.elements = "1",
                       scaling = FALSE,
                       only.means.matrix = FALSE,
                       use.previous.model = NULL,
                       Nb_compo_ortho = 1,
                       perm_t = NULL) {
  t_pre <- fun_pre_processing(
    Data = Data,
    Nb_compo_ortho = Nb_compo_ortho,
    equation.elements = equation.elements,
    scaling = scaling,
    only.means.matrix = only.means.matrix
  )
  t_ANOVA <- fun_decompo_ANOVA(
    s = t_pre,
    only.means.matrix = only.means.matrix,
    perm_t = perm_t
  )
  t_ANOVA_PCA <- fun_ANOVA_PCA(s = t_ANOVA)
  t_multiblock <- fun_Multiblock(s = t_ANOVA_PCA)

  t_multiblock$kOPLS <- t_multiblock %>% {
    fun_kopls(.$Multiblock$W_mat,
              .$Multiblock$Y,
              .$general$Nb_compo_pred,
              .$general$Nb_compo_ortho)
  }

  # t <- fun_outputs(t) most long task, perform only on demand
  return(t_multiblock)
}


#' Function to calculate permutation
#'
#' @param iter Number of iterations to compute
#' @inheritParams run_AMOPLS
#' @inheritParams fun_AMOPLS
#'
fun_temp_perm <- function(Data, equation.elements, scaling, Nb_compo_ortho, perm_t, iter) {
  P_Results <- fun_AMOPLS(Data = Data,
                          equation.elements = equation.elements,
                          scaling = scaling,
                          only.means.matrix = FALSE,
                          use.previous.model = NULL,
                          Nb_compo_ortho = Nb_compo_ortho,
                          perm_t = perm_t)

  if (is.null(P_Results)) {return(NULL)}

  factors_element <- strsplit(P_Results$general$studied_factors, ",") %>% unlist()
  output <- data.table(Iter = iter,
                       RSR = fun_RSR(P_Results)[, which(perm_t %in% factors_element), with = F],
                       RSS = fun_Sum_of_Squares(P_Results)[[2]][, which(perm_t %in% factors_element), with = F],
                       R2Y = P_Results$kOPLS$R2Yhat %>% {.[length(.)]})
  setnames(output, c("Iter", "RSR", "RSS", "R2Y"))
  return(output)
}


#' Wrapper to run AMOPLS models
#'
#' This function is a wrapper to perform AMOPLS model with permutation and
#' subsampling if data are unbalanced.
#'
#' @param scaling Logical for unit variance scaling of the data before running the model
#' @param nb_perm Number of permutation for each effect to compute p-values
#' @param nb_compo_orthos Number of orthogonal component to model
#' @param parallel Number of process to run in parallel using future and furrr
#' @param debug Logical to run a logger with debug messages
#' @param datamatrix The datamatrix with observations id in the first column (observations x variables)
#' @param samplemetadata The observations metadata with groups and levels (the first column must be the observations id)
#' @param factor_names Name of the column in samplemetadata to use for effect decomposition
#' @param interaction_level Order of interaction to consider (0 = pure effect only, 1 first order interaction between each effect)
#' @param subsampling Number of subsampling to perform if the data are unbalanced
#'
#' @import data.table
#' @import magrittr
#' @importFrom stats median
#'
#' @return \code{s} List containing all the information about the AMOPLS model, organized in 2 groups :
#' \itemize{
#'     \item \code{general} General information about the parameters
#'     \item \code{output} Main outcomes from AMOPLS}
#' @export
#'
#' @examples
#'result <- run_AMOPLS(datamatrix = data_Ruiz2017$datamatrix,
#'                     samplemetadata = data_Ruiz2017$samplemetadata,
#'                     factor_names = c("Exposure time", "Dose"))
run_AMOPLS <- function(datamatrix,
                       samplemetadata,
                       factor_names,
                       interaction_level = 1,
                       scaling = T,
                       nb_perm = 100,
                       subsampling = NULL,
                       nb_compo_orthos = 1:3,
                       parallel = F,
                       debug = F) {
  # datamatrix = data_Ruiz2017$datamatrix
  # samplemetadata = data_Ruiz2017$samplemetadata
  # factor_names = c("Exposure time", "Dose")
  DEBUG <- Effect <- Iter <- Iteration <- Ortho <- Ortho_nb <- PermNb <- R2Y <- R2Y <- `p-value` <- R2Y_pval <- V_scores <- V_sign <- availableWorkers <- combn <- `cor.test` <- density <- id <- layout <- plot <- rn <- str <- tp_calc <- value <- variable <- variableid <- x <- y <- NULL
  if (debug) {
    DEBUG <- NULL
    requireNamespace("logger")
    logger::log_appender(logger::appender_console)
    logger::log_threshold(DEBUG)
  }

  if (debug) {logger::log_info("Starting function")}

  ## Format data
  ### Check all column are numeric
  Data <- fun_load_data(datamatrix, samplemetadata, factor_names)

  factor_index <- which(colnames(Data$factors) %in% factor_names)

  ## Generate formula with interaction levels
  equation.elements <- sapply(0:interaction_level, function(int) {
    temp <- utils::combn(factor_index, int + 1, simplify = F)
    sapply(temp, paste, collapse = "")
  }) %>% unlist() %>% paste(., collapse = ",")

  factors_element <- strsplit(equation.elements, ",") %>% unlist()
  nb_studied_factors <- length(factors_element)

  if (debug) {logger::log_info("Factors to study: {paste(factors_element, collapse = ', ')}")}

  ## Check if data are balanced
  if (!fun_is_balanced(Data, factor_names = factor_names, interaction_level = interaction_level)) {
    if (is.null(subsampling)) {
      stop("Data are unbalanced, set the subsampling argument to run subsampling stratification.")
    } else {
      message("Data are unbalanced, running stratified subsampling.")
      subsampling <- as.numeric(subsampling)
    }
  } else {
    subsampling <- 1
  }

  ## Run original model
  if (debug) {logger::log_info("Calculate full model for {length(nb_compo_orthos)} orthogonal components: ")}

  res_subsampling <- lapply(1:subsampling, function(zrf) {
    # zrf <- 1
    ## Balance the data only if subsampling is > 1
    if (subsampling > 1) {
      message("Run sub-sampling: ", zrf)
      temp_data <- fun_balance_data(Data)
    } else {
      temp_data <- Data
    }

    result_original <- lapply(1:length(nb_compo_orthos), function(x) {
      # x <- 1
      output <- fun_AMOPLS(Data = temp_data,
                           equation.elements = equation.elements,
                           scaling = scaling,
                           only.means.matrix = FALSE,
                           use.previous.model = NULL,
                           Nb_compo_ortho = nb_compo_orthos[[x]],
                           perm_t = NULL)

      if (is.null(output)) {stop("Resolve the collinearity problems in the data")} else {
        output <- fun_outputs(output)
      }

      if (debug) {logger::log_info("Ortho {nb_compo_orthos[x]}: R2Y={output$outputs$R2Y %>% formatC(., digits = 2)} Cp={ncol(output$outputs$block_saliences_norm)-1}")}
      output$outputs$summary <- data.table("Effect" = output$general$ee.names,
                                           Iter = 0,
                                           RSR = t(output$outputs$RSR),
                                           RSS = t(output$outputs$RSS),
                                           R2Y = t(output$kOPLS$R2Yhat %>% {.[length(.)]}))
      setnames(output$outputs$summary, c("Effect", "Iter", "RSR", "RSS", "R2Y"))
      return(output)
    })


    ## Application of AMOPLS for permuted data - for each factor + interaction - and calculation of the scores :
    ## Create iteration arguments
    iter_template <- CJ("Effect" = factors_element, "PermI" = 1:nb_perm, "Ortho" = nb_compo_orthos)
    iter_template[, Iter := 1:.N]
    apply_it <- nrow(iter_template)
    if (debug) {logger::log_info("Running {nb_perm} permutations for each factor and ortho cp: {nb_perm} x {nb_studied_factors} x {length(nb_compo_orthos)} = {nb_perm*nb_studied_factors*length(nb_compo_orthos)}")}

    if (!is.null(parallel) & !isFALSE(parallel)) {
      if(!requireNamespace("future")) {
        stop("You need to install future and furr packages to use parallelisation.")
      } else {requireNamespace("future")}
      if(!requireNamespace("furrr")) {
        stop("You need to install future and furr packages to use parallelisation.")
      } else {requireNamespace("furrr")}

      if (is.numeric(parallel)) {
        future::plan(future::multiprocess, workers = parallel)
      } else {
        future::plan(future::multiprocess, workers = (length(future::availableWorkers()) - 1))
      }

      temp <- furrr::future_map(1:apply_it, function(x) {
        # x <- 1
        temp_effect <- iter_template[Iter == x, Effect]
        temp_ortho <- iter_template[Iter == x, Ortho]
        P_Results <- fun_temp_perm(Data = temp_data,
                                   equation.elements = equation.elements,
                                   scaling = scaling,
                                   Nb_compo_ortho = temp_ortho,
                                   perm_t = temp_effect,
                                   iter = x)
        return(P_Results)
      }, .progress = TRUE)
    } else {
      pb <- progress::progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = apply_it)
      pb$tick(0)
      temp <- lapply(1:apply_it, function(x) {
        # x <- 1
        # message(x)
        pb$tick()
        if (debug) {logger::log_trace("Perm: {x}, Effect: {iter_template[x, Effect]}, Ortho: {iter_template[x, Ortho]}")}
        temp_effect <- iter_template[Iter == x, Effect]
        temp_ortho <- iter_template[Iter == x, Ortho]
        P_Results <- fun_temp_perm(Data = temp_data,
                                   equation.elements = equation.elements,
                                   scaling = scaling,
                                   Nb_compo_ortho = temp_ortho,
                                   perm_t = temp_effect,
                                   iter = x)
        return(P_Results)
      })
    }

    temp_dt <- temp %>%
      rbindlist() %>%
      merge(., iter_template, by = "Iter")

    ## P-value calculation
    output_pval <- lapply(1:length(nb_compo_orthos), function(cp) {
      # cp <- 1
      lapply(c("RSS", "RSR", "R2Y"), function(x) {
        # x <- "R2Y"
        lapply(iter_template[, unique(Effect)], function(effect) {
          # effect <- iter_template[, unique(Effect)][[2]]
          subset_perm <- temp_dt[Ortho == nb_compo_orthos[[cp]] & Effect == effect, x, with = F][[1]]
          subset_ori <- unlist(result_original[[cp]]$outputs$summary[Effect == effect, x, with = F])
          temp_subset <- length(which(subset_perm >= subset_ori))/length(subset_perm)
          if (temp_subset == 0) {temp_subset <- 1/length(subset_perm)}
          effect_name <- data.table(Effect_name = result_original[[cp]]$general$factor_names[which(factors_element %in% effect)])
          output <- data.table(effect, effect_name, temp_subset)
          setnames(output, c("Effect", "Effect_name", paste0(x, "_pvalue")))
          return(output)
        }) %>% rbindlist(use.names = T)
      }) %>% {Reduce(function(z, w) {merge(z, w, by = c("Effect", "Effect_name"))}, .)}
    })

    output_pval <- lapply(output_pval, function(x) {
      x[, PermNb := nb_perm]
      x[, Effect := factor(Effect, levels = factors_element)]
      x[order(Effect)]
    })

    output <- mapply(function(x, y, z) {
      # x <- result_original[[1]]
      # y <- output_pval[[1]]
      x$outputs$Permutation_result <- list("summary" = y,
                                           "details" = z)
      return(x)
    }, result_original, output_pval, split(temp_dt, temp_dt$Ortho), SIMPLIFY = F)
    names(output) <- paste0("orthoNb_", nb_compo_orthos)
    return(output)
  })

  if (!isFALSE(parallel)) {
    if (inherits(future::plan(), "multiprocess")) {future::plan(future::sequential)}
  }

  ## Aggregate subsampled models using median as in Boccard et al., 2019
  ### Extract data to combine in each subsampled results
  output <- lapply(res_subsampling, function(w) {
    # w <- res_subsampling[[1]]
    lapply(w, function(z) {
      # z <- w[[1]]
      temp_dt <- z
      output <- list(
        "general" = temp_dt$general[c("data",
                                      "factors",
                                      "ssq",
                                      "Nb_compo_pred",
                                      "Nb_compo_ortho",
                                      "ee.names",
                                      "studied_factors",
                                      "equation.elements",
                                      "order.to.evaluate.ee",
                                      "factor_names",
                                      "factors_data")],
        "decompo_ANOVA" = temp_dt$decompo_ANOVA,
        "kOPLS" = list("R2Yhat" = temp_dt$kOPLS$R2Yhat %>% {.[length(.)]}),
        "output" = list(
          "x_loadings" = fun_xloadings(temp_dt),
          "x_scores" = fun_xscores(temp_dt),
          "block_saliences_norm" = fun_block_saliences(temp_dt)[[2]],
          "RSS" = fun_get_RSS(temp_dt),
          "RSR" = fun_get_RSR(temp_dt),
          "SSQ" = fun_Sum_of_Squares(temp_dt)[[2]],
          "Permutation_result" = temp_dt$outputs$Permutation_result["summary"],
          "VIP" = fun_VIP(temp_dt),
          "Summary" = fun_AMOPLS_summary(temp_dt)
        )
      )
      return(output)
    })
  })

  ## Extract for each orthogonal component
  results_combined <- lapply(1:unique(sapply(output, length)), function(z) {
    # z <- 1
    temp_orthon <- lapply(output, function(w) {w[[z]]})

    ## Calculate median
    ### Need to check scores and loadings orientation of each component (may be arbitrarly reversed between models)

    #### SCORES
    if (temp_orthon %>% length() <= 1) {
      x <- temp_orthon[[1]]
      output_scores <- data.table("id" = rownames(x$general$factors), x$output$x_scores) %>% as.matrix(., rownames = "id")
    } else {
      output_scores <- lapply(temp_orthon, function(x) {
        # x <- temp_orthon[[1]]
        data.table("id" = rownames(x$general$factors), x$output$x_scores)
      }) %>%
        rbindlist(use.names = T, fill = TRUE, idcol = "rn")
      ## Check component correlation between each models
      scores_sign <- sapply(1:(ncol(output_scores)-2), function(x) {
        # x <- 1
        var_col <- names(output_scores)[x+2]
        col_sel <- c("id", "rn", var_col)
        output_scores[, col_sel, with = F] %>% dcast(., id ~ rn, value.var = var_col) %>% {.[, -1][, lapply(.SD, function(z) {if(stats::cor.test(z, .[, 2][[1]])$estimate < 0) {return(-1)} else {return(1)}} %>% round(., 1))]}
      }) %>% as.data.table(keep.rownames = 'rn')

      scores_sign[, rn := as.numeric(rn)]
      output_scores <- lapply(1:(ncol(output_scores)-2), function(x) {
        # x <- 1
        var_col <- names(output_scores)[x+2]
        ## Reverse axes with negative correlation by component
        temp_merge <- merge(output_scores[, .(id, rn, "V_scores" = get(var_col))],
                            scores_sign[, .(rn, "V_sign" = get(var_col))], by = "rn")
        row_nb <- temp_merge[, .N]
        temp_merge[, tp_calc := ifelse(any(is.na(V_scores), is.na(V_sign)), NA, as.numeric(V_scores) * as.numeric(V_sign)), by = 1:row_nb]
        ## Calculate median
        output <- temp_merge[, median(tp_calc), by = "id"]
        setnames(output, c("id", var_col))
        return(output)
      }) %>% {Reduce(function(x, y) {merge(x, y, by = "id", all = T)}, .)} %>%
        {.[rownames(Data$dataset)]} %>%
        as.matrix(rownames = "id")
    }

    ## LOADINGS
    if (temp_orthon %>% length() <= 1) {
      x <- temp_orthon[[1]]
      output_loadings <- data.table("id" = colnames(x$general$data), x$output$x_loadings) %>% as.matrix(., rownames = "id")
    } else {
      output_loadings <- lapply(temp_orthon, function(x) {
        # x <- temp_orthon[[1]]
        data.table("id" = colnames(x$general$data), x$output$x_loadings)
      }) %>%
        rbindlist(use.names = T, fill = TRUE, idcol = "rn")
      ## Check component correlation between each models
      loadings_sign <- sapply(1:(ncol(output_loadings)-2), function(x) {
        # x <- 1
        var_col <- names(output_loadings)[x+2]
        col_sel <- c("id", "rn", var_col)
        output_loadings[, col_sel, with = F] %>% dcast(., id ~ rn, value.var = var_col) %>% {.[, -1][, lapply(.SD, function(z) {if(stats::cor.test(z, .[, 2][[1]])$estimate < 0) {return(-1)} else {return(1)}} %>% round(., 1))]}
      }) %>%
        as.data.table(keep.rownames = 'rn')
      loadings_sign[, rn := as.numeric(rn)]
      output_loadings <- lapply(1:(ncol(output_loadings)-2), function(x) {
        # x <- 1
        var_col <- names(output_loadings)[x+2]
        ## Reverse axes with negative correlation by component
        temp_merge <- merge(output_loadings[, .(id, rn, "V_scores" = get(var_col))],
                            loadings_sign[, .(rn, "V_sign" = get(var_col))], by = "rn")
        row_nb <- temp_merge[, .N]
        temp_merge[, tp_calc := ifelse(any(is.na(V_scores), is.na(V_sign)), NA, as.numeric(V_scores) * as.numeric(V_sign)), by = 1:row_nb]
        ## Calculate median
        output <- temp_merge[, median(tp_calc), by = "id"]
        setnames(output, c("id", var_col))
        return(output)
      }) %>% {Reduce(function(x, y) {merge(x, y, by = "id", all = T)}, .)} %>%
        {.[colnames(Data$dataset)]} %>%
        as.matrix(rownames = "id")
    }
    output_saliences <- lapply(temp_orthon, function(x) {
      # x <- temp_orthon[[3]]
      data.table("id" = c(x$general$factor_names[x$general$order.to.evaluate.ee], "residual"), x$output$block_saliences_norm)
    }) %>%
      rbindlist(use.names = T, fill = TRUE) %>%
      {.[, lapply(.SD, median), keyby = "id"]} %>%
      as.matrix(rownames = "id")

    output_RSS <- lapply(temp_orthon, function(x) {
      # x <- temp_orthon[[3]]
      x$output$RSS
    }) %>%
      rbindlist(use.names = T, fill = TRUE) %>%
      {.[, lapply(.SD, median), keyby = c("Effect", "Effect Name")]}

    output_RSR <- lapply(temp_orthon, function(x) {
      # x <- temp_orthon[[3]]
      x$output$RSR
    }) %>%
      rbindlist(use.names = T, fill = TRUE) %>%
      {.[, lapply(.SD, median), keyby = c("Effect", "Effect Name")]}

    output_SSQ <- lapply(temp_orthon, function(x) {
      # x <- temp_orthon[[3]]
      x$output$SSQ
    }) %>%
      rbindlist(use.names = T, fill = TRUE) %>%
      {.[, lapply(.SD, median)]}

    output_Perm <- lapply(temp_orthon, function(x) {
      # x <- temp_orthon[[3]]
      x$output$Permutation_result$summary
    }) %>%
      rbindlist(use.names = T, fill = TRUE) %>%
      {.[, lapply(.SD, median), keyby = c("Effect", "Effect Name")]}

    output_R2Y <- sapply(temp_orthon, function(x) {
      # x <- temp_orthon[[1]]
      x$kOPLS$R2Yhat
    }) %>%
      median()

    output_VIP <- lapply(temp_orthon, function(x) {
      # x <- temp_orthon[[1]]
      x$output$VIP
    }) %>%
      rbindlist(use.names = T, fill = TRUE) %>%
      {.[, lapply(.SD, median), keyby = c("id")]}

    output_Summary <- lapply(temp_orthon, function(x) {
      # x <- temp_orthon[[1]]
      x$output$Summary
    }) %>%
      rbindlist(use.names = T, fill = TRUE) %>%
      {.[, lapply(.SD, function(x) {median(x, na.rm = T)}), keyby = c("Effect", "Effect Name")]}

    #### DEV
    output_residuals <- lapply(temp_orthon, function(x) {
      # x <- temp_orthon[[1]]
      x$decompo_ANOVA$residuals %>% as.data.table(keep.rownames = "sampleid")
    }) %>%
      rbindlist(use.names = T, fill = TRUE) %>%
      {.[, lapply(.SD, function(x) {median(x, na.rm = T)}), keyby = c("sampleid")]} %>%
      {as.data.frame(., row.names = "sampleid")}

    output_decompoANOVA <- lapply(1:(length(temp_orthon[[1]]$decompo_ANOVA) - 1), function(y) {
      # y <- 1
      lapply(temp_orthon, function(x) {
        # x <- temp_orthon[[1]]
        x$decompo_ANOVA %>%
          {.[which(!names(.) %in% "residuals")]} %>% .[[y]] %>% .[["means.matrix_res"]] %>%
          as.data.table(keep.rownames = "sampleid")
      }) %>%
        rbindlist(use.names = T, fill = TRUE) %>%
        {.[, lapply(.SD, function(x) {median(x, na.rm = T)}), keyby = c("sampleid")]} %>%
        as.data.frame(., row.names = "sampleid")
    })
    names(output_decompoANOVA) <- names(temp_orthon[[1]]$decompo_ANOVA) %>% {.[!. == "residuals"]}

    factors_data <- fun_factors_data(Data$factors, temp_orthon[[1]]$general$studied_factors)

    return(
      list(
        "general" = c(Data, list("factors_data" = factors_data), temp_orthon[[1]]$general %>% {.[!names(.) %in% c("dataset", "factors")]}),
        "decompo_ANOVA" = c(list("residuals" = output_residuals), output_decompoANOVA),
        "kOPLS" = list("R2Yhat" = output_R2Y),
        "output" = list(
          "x_loadings" = output_loadings,
          "x_scores" = output_scores,
          "block_saliences_norm" = output_saliences,
          "RSS" = output_RSS,
          "RSR" = output_RSR,
          "SSQ" = output_SSQ,
          "R2Y" = output_R2Y,
          "Permutation_result" = output_Perm,
          "Summary" = output_Summary,
          "VIP" = output_VIP
        )
      )
    )
  })

  names(results_combined) <- names(output[[1]])
  return(results_combined)
}



#' @title Score plot
#' @description Score plot of the x-scores from AMOPLS results, according to the selected factor and components.
#'
#' @param fact Studied factor
#' @param t_1 First component to project the x-scores
#' @param t_2 Second component to project the x-scores
#' @inheritParams fun_outputs
#'
#' @return 2D score plot of the x-scores from AMOPLS results according to the 2 selected components. Every datapoint is colored according to its factor-level. Every color group is surrounded by a convex hull.
#'
#' @import ggplot2
#' @import magrittr
#' @import data.table
#' @import ggpubr
#'
#' @references From grDevices chull
#' @export
fun_score_plot <- function(s, fact, t_1 = NULL, t_2 = NULL) {
  # s <- result_optimal
  # fact <- "Dose"
  x <- y <- NULL
  nb <- which(rownames(s$output$block_saliences_norm) == fact)
  nb_compo <- (s$general$Nb_compo_pred + 1)

  if (all(is.null(t_1), is.null(t_2))) {
    t_1 <- which(s$output$block_saliences_norm[nb, ] ==  max(s$output$block_saliences_norm[nb,-c(nb_compo)]))
    t_2 <- which(s$output$block_saliences_norm[nb, ] ==  max(s$output$block_saliences_norm[nb,-c(t_1, nb_compo)]))
  }

  temp_scores <- as.data.table(s$output$x_scores, keep.rownames = "sampleid")[, c(1, t_1+1, t_2+1), with = F]
  temp_plot <- data.table(temp_scores, s$general$factors_data)
  setnames(temp_plot, 2:3, c("x", "y"))

  ## Convex hulls :
  find_hull <- function(df) {df[grDevices::chull(df$x, df$y),]}
  hulls <- temp_plot[!is.na(x) | !is.na(y)][, find_hull(.SD), by = fact]

  sp <- temp_plot %>%
    ggplot2::ggplot(aes(x, y, color = factor(get(fact)))) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point() +
    geom_polygon(
      data = hulls,
      alpha = 0.2,
      aes(fill = factor(get(fact))),
      show.legend = F
    ) +
    labs(
      title = "",
      subtitle = s$general$factor_names[nb],
      x = paste('tp', as.character(t_1)),
      y = paste('tp', as.character(t_2)),
      color = fact
    ) +
    theme_bw() +
    labs(
      title = "AMOPLS score plot",
      subtitle = paste0("Colored by factor: ", fact),
      x = paste("tp", t_1),
      y = paste("tp", t_2)
    )

  return(sp)
}

#' Generate optimal score plots
#'
#' This function creates a plot for each factor considered
#' with the 2 best components for each factor.
#'
#' @inheritParams fun_outputs
#'
#' @import ggplot2
#' @import magrittr
#' @import data.table
#' @import ggpubr
#'
#' @export
#'
fun_plot_optimal_scores <- function(s) {
  # s <- result_optimal
  nb_factors <- length(s$general$factor_names)
  ncol <- ceiling(sqrt(nb_factors))
  nrow <- ceiling(sqrt(nb_factors))

  s$general$factor_names %>%
    lapply(., function (x) {
      fun_score_plot(s, x, NULL, NULL) + labs(title = NULL)
    }) %>%
    ggpubr::ggarrange(plotlist = ., ncol = ncol, nrow = nrow, align = "hv")
}

#' Loading plots
#'
#' Loading plot according to the selected components.
#'
#' @param fact Studied factor
#' @param t_1 First component to project the scores
#' @param VIP_nb Number of VIP to return (top n)
#' @param t_2 Second component to project the scores
#' @inheritParams fun_outputs
#'
#' @import ggplot2
#' @import magrittr
#' @import data.table
#'
#' @return 2D loading plot from the AMOPLS results (to complete).
#' @export
fun_loading_plot <- function(s,
                             fact,
                             t_1 = NULL,
                             t_2 = NULL,
                             VIP_nb = NULL) {
  variableid <- x <- y <- id <- NULL
  nb <- which(rownames(s$output$block_saliences_norm) == fact)
  nb_compo <- (s$general$Nb_compo_pred + 1)

  if (all(is.null(t_1), is.null(t_2))) {
    t_1 <- which(s$output$block_saliences_norm[nb, ] ==  max(s$output$block_saliences_norm[nb,-c(nb_compo)]))
    t_2 <- which(s$output$block_saliences_norm[nb, ] ==  max(s$output$block_saliences_norm[nb,-c(t_1, nb_compo)]))
  }

  temp_plot <- as.data.table(s$output$x_loadings, keep.rownames = "variableid")[, c(1, t_1+1, t_2+1), with = F]
  setnames(temp_plot, 2:3, c("x", "y"))

  if (is.null(VIP_nb)) {
    ## If null, show 10 variable names if available
    if (temp_plot[, length(unique(variableid))] > 10) {
      VIP_nb <- 10
    } else {
      VIP_nb <- temp_plot[, length(unique(variableid))]
    }
  }

  ## Convex hulls :
  sp <- temp_plot %>% {
    ggplot2::ggplot(., aes(x, y)) +
      geom_vline(xintercept = 0, linetype = 2) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_point(alpha = 0.6) +
      geom_point(data = .[variableid %in% s$output$VIP[, c("id", fact), with = F][order(-get(fact))][1:VIP_nb, id]], color = "red") +
      geom_text(data = .[variableid %in% s$output$VIP[, c("id", fact), with = F][order(-get(fact))][1:VIP_nb, id]], aes(label = variableid), color = "red", vjust = -0.5) +
      labs(
        title = "AMOPLS score plot",
        subtitle = paste0("Factor: ", fact),
        caption = paste0("In red: the top ", VIP_nb, " VIPs"),
        x = paste("tp", t_1),
        y = paste("tp", t_2)
      ) +
      theme_bw()
  }

  return(sp)
}

#' Generate optimal loading plot
#'
#' This function creates a plot for each factor considered
#' with the 2 best components for each factor.
#'
#' @inheritParams fun_loading_plot
#' @inheritParams fun_outputs
#'
#' @import ggplot2
#' @import magrittr
#' @import data.table
#' @import ggpubr
#' @export
fun_plot_optimal_loadings <- function(s, VIP_nb = NULL) {
  # s <- result_optimal
  nb_factors <- length(s$general$factor_names)
  ncol <- ceiling(sqrt(nb_factors))
  nrow <- ceiling(sqrt(nb_factors))

  s$general$factor_names %>%
    lapply(., function (x) {
      fun_loading_plot(s, x, NULL, NULL, VIP_nb = VIP_nb) + labs(title = NULL)
    }) %>%
    ggpubr::ggarrange(plotlist = ., ncol = ncol, nrow = nrow, align = "hv")
}

#' Plot VIP
#'
#' This function plot the VIP2 of all variables for each factors
#'
#' @param main_factor String to set the main_factor to order the plot
#' @param debugL Boolean to activate the debug mode
#' @inheritParams fun_outputs
#' @inheritParams fun_loading_plot
#'
#' @import ggplot2
#' @import magrittr
#' @import data.table
#'
#' @export
fun_plot_VIPs <- function(s, main_factor = NULL, VIP_nb = NULL, debugL = F) {
  # s <- result_optimal
  # main_factor <- "Dose"
  # VIP_nb <- NULL
  # debugL <- F
  str <- id <- variable <- value <- NULL
  temp_factors <- colnames(s$output$VIP[, -1])
  if (is.null(main_factor)) {
    main_factor <- temp_factors[[1]]
  } else if (!main_factor %in% temp_factors) {stop("The main_factor wasn't found in the dataset")}

  data_vips <- as.data.table(s$output$VIP)

  if (is.null(VIP_nb)) {
    VIP_nb <- data_vips[, .N]
  }

  if (all(VIP_nb != "Force", VIP_nb >= 200)) {
    message("There are a high number of variables (", VIP_nb, ") filter the 200 most significant. To force all variable, set VIP_nb argument to 'Force'.")
    VIP_nb <- 200
  }

  if (debugL) {message("Main fac: ", main_factor)}
  if (debugL) {message("Main fac (str): ", str(main_factor))}
  if (debugL) {message("VIP nb: ", VIP_nb)}
  if (debugL) {message("data_vips colnames: ", paste(names(data_vips), collapse = ", "))}
  ## Filter the most significant variables for the considered factor
  if (debugL) {message("data_vips (class): ", class(data_vips))}

  VIP <- data_vips[order(-get(main_factor))][1:VIP_nb]
  # Reorder variables order by decreasing order
  VIP[, id := factor(id, levels = unique(VIP$id))]
  # Melt data
  plot_data <- VIP %>% melt(id.vars = "id")
  ## Set factor order (first is factor of interest)
  plot_data[, variable := factor(variable, levels = rev(union(main_factor, temp_factors)))]

  ggplot2::ggplot(data = plot_data, aes(x = id, y = value, fill = variable)) +
    geom_bar(stat = "identity", col = 'black') +
    labs(title = "Variable Important in the Projection (VIP2)",
         subtitle = paste0("By decreasing order of importance for factor: ", main_factor),
         x = '',
         y = bquote(~VIP^2)) +
    theme_bw() +
    theme(legend.position = c("right"),
          axis.text.x = element_text(angle = 60, hjust = 1),
          legend.title = element_blank(),
          legend.key.size = unit(0.8, "cm"),
          legend.text = element_text(size = 11, hjust = 0.3, face = 'bold'))
}


#' Function to get RSS
#'
#' @inheritParams fun_outputs
#'
#' @export
fun_get_RSS <- function(s) {
  Effect <- NULL
  temp <- fun_Sum_of_Squares(s)[[2]] %>% t() %>% {as.data.table(., keep.rownames = "Effect")}
  temp[, "Effect Name" := c(s$general$factor_names, "residuals")[s$general$ee.names == Effect]]
  setnames(temp, c("Effect", "RSS", "Effect Name"))
  setcolorder(temp, c("Effect", "Effect Name", "RSS"))
  return(temp)
}

#' Function to get RSR
#'
#' @inheritParams fun_outputs
#'
#' @export
fun_get_RSR <- function(s) {
  Effect <- NULL
  temp <- fun_RSR(s) %>% t() %>% {as.data.table(., keep.rownames = "Effect")}
  temp[, "Effect Name" := c(s$general$factor_names, "residuals")[s$general$ee.names == Effect]]
  setnames(temp, c("Effect", "RSR", "Effect Name"))
  setcolorder(temp, c("Effect", "Effect Name", "RSR"))
  return(temp)
}

#' Function to get normalized block contribution
#'
#' @inheritParams fun_outputs
#'
#' @export
fun_get_blockcontrib <- function(s) {
  Effect <- NULL
  temp <- s$outputs$block_saliences_norm %>% as.data.table() %>% {data.table("Effect" = s$general$ee.names, .)}
  temp[, "Effect Name" := c(s$general$factor_names, "residuals")[s$general$ee.names == Effect]]
  setcolorder(temp, c("Effect", "Effect Name"))
  setnames(temp, c("Effect", "Effect Name", paste0("Tp", 1:s$general$Nb_compo_pred), paste0("To", 1:s$general$Nb_compo_ortho)))
  return(temp)
}

#' Function to get permutation results
#'
#' @inheritParams fun_outputs
#'
#' @export
fun_get_perm <- function(s) {
  temp <- s$outputs$Permutation_result$summary
  setnames(temp, c("Effect", "Effect Name", "RSS p-value", "RSR p-value", "R2Y p-value", "PermNb"))
  return(temp)
}

#' Summary of AMOPLS results
#'
#' This function retrieve different levels of summary from the output of AMOPLS.
#'
#' @param type String to select the summary to return (All, RSS, RSR, Permutation or Block contrib)
#' @inheritParams fun_outputs
#'
#' @export
fun_AMOPLS_summary <- function(s, type = c("All", "RSS", "RSR", "Permutation", "Block contrib")) {
  # s <- temp_dt
  # type <- 'All'
  switch(type[[1]],
         'RSS' = fun_get_RSS(s),
         'RSR' = fun_get_RSR(s),
         'Permutation' = fun_get_perm(s),
         'Block contrib' = fun_get_blockcontrib(s),
         {
           s %>% {
             list(fun_get_RSS(.),
                  fun_get_RSR(.),
                  fun_get_perm(.),
                  fun_get_blockcontrib(.)
             )
           } %>% {
             Reduce(function(x, y) {merge(x, y, by = c("Effect", "Effect Name"), all = T)}, .)
           }
         }
  )
}

### Remove MetStat dependency

#' Function from MetStat package
#'
#' Copy of the MetStat function to remove partial dependency
#'
#' @param x.input The data matrix that needs to be scaled.
#' @param center Boolean. If TRUE the data will also be centered per column (the mean of each column will become zero).
#' @param scale This Argument defines which type of scaling is to be applied. With the default value of TRUE, the data is autoscaled. When set to "pareto", pareto scaling is applied.
#' @param quietly Boolan. If TRUE, no intermediate text output concerning the centering and scaling methods is returned.
#'
#' @export
MetStaT.ScalePip <- function (x.input, center = TRUE, scale = TRUE, quietly = FALSE) {
  options(warn = -1)
  no.col.x.input <- ncol(x.input)
  if (is.null(no.col.x.input)) {
    no.col.x.input <- 1
  }
  tryCatch({
    x <- matrix(as.numeric(x.input), ncol = no.col.x.input)
  }, error = function(ex) {
    bad.matrix <- x.input
    stop(ex)
  })
  colnames(x) <- colnames(x.input)
  rownames(x) <- rownames(x.input)
  options(warn = 0)
  x.scaled <- list()
  nc <- ncol(x)
  if (is.null(center))
    center <- FALSE
  if (is.character(center) && center == "true")
    center <- TRUE
  if (is.character(center) && center == "false")
    center <- FALSE
  if (is.character(scale) && scale == "true")
    scale <- TRUE
  if (is.character(scale) && scale == "false")
    scale <- FALSE
  center.description <- center
  if (is.logical(center)) {
    if (center) {
      center.description <- "Around mean. "
      center <- colMeans(x, na.rm = TRUE)
      x <- sweep(x, 2L, center, check.margin = FALSE)
    }
    else {
      x.scaled$description <- paste(x.scaled$description,
                                    "Not centered. ", sep = "")
      not.centered <- matrix(rep(0, nc), nrow = 1)
      colnames(not.centered) <- colnames(x)
      x.scaled$center.vector <- not.centered
    }
  }
  else if (is.numeric(center) && (length(center) == nc)) {
    center.description <- "Manual input by user used. "
    x <- sweep(x, 2L, center, check.margin = FALSE)
  }
  else {
    stop("length of 'center' must equal the number of columns of 'x'")
  }
  if (is.numeric(center)) {
    x.scaled$description <- paste(x.scaled$description, "Centered: ",
                                  center.description, sep = "")
    center <- matrix(center, nrow = 1)
    colnames(center) <- colnames(x)
    x.scaled$center.vector <- center
  }
  if (is.null(scale))
    scale <- FALSE
  if (is.logical(scale)) {
    if (scale) {
      scale = "stdev"
    }
  }
  scale.description <- scale
  if (is.logical(scale)) {
    x.scaled$description <- paste(x.scaled$description, "Not scaled. ",
                                  sep = "")
    not.scaled <- matrix(rep(1, nc), nrow = 1)
    colnames(not.scaled) <- colnames(x)
    x.scaled$scale.vector <- not.scaled
  }
  else if (is.character(scale)) {
    scale <- tolower(scale)
    if (scale == "stdev" || scale == "auto") {
      f <- function(v) {
        v <- v[!is.na(v)]
        sqrt(sum(v^2)/max(1, length(v) - 1L))
      }
    }
    else if (scale == "pareto") {
      f <- function(v) {
        v <- v[!is.na(v)]
        sqrt(sqrt(sum(v^2)/max(1, length(v) - 1L)))
      }
    }
    scale <- apply(x, 2L, f)
    x <- sweep(x, 2L, scale, "/", check.margin = FALSE)
  }
  else if (is.numeric(scale) && length(scale) == nc) {
    scale.description <- "Manual input by user used."
    x <- sweep(x, 2L, scale, "/", check.margin = FALSE)
  }
  else {
    stop("length of 'scale' must equal the number of columns of 'x'")
  }
  if (is.numeric(scale)) {
    x.scaled$description <- paste(x.scaled$description, "Scaled: ",
                                  scale.description, ".", sep = "")
    scale <- matrix(scale, nrow = 1)
    colnames(scale) <- colnames(x)
    x.scaled$scale.vector <- scale
  }
  x.scaled$data <- x
  if (!quietly) {
    print(x.scaled$description)
  }
  x.scaled
}

#' Title
#'
#' @param data A datamatrix (sample x variables)
#'
#' @export
MetStat.PCA.Calculate <- function (data) {
  svd.result <- svd(data)
  svd.result$var.explained <- svd.result$d^2
  svd.result$var.explained <- svd.result$var.explained/(sum(svd.result$var.explained))
  svd.result$t <- svd.result$u %*% diag(svd.result$d)
  svd.result$u <- NULL
  svd.result
}


#' Cutome koplsModel function with tol param
#'
#' Correct the error returned by solve:  system is computationally singular
#'
#' @inheritParams base::solve
#' @inheritParams kopls::koplsModel
#'
koplsModel_custom <- function(K, Y, A, nox, preProcK = "mc", preProcY = "mc", tol = 1e-20) {
  if (!requireNamespace("kopls", quietly = TRUE)) {
    stop("Package \"kopls\" needed for this function to work. Please install it using install_kopls()",
         call. = FALSE)
  } else {requireNamespace("kopls")}
  n = ncol(K)
  I <- diag(rep(1, n))
  if (preProcK == "mc") {
    Kmc <- kopls::koplsCenterKTrTr(K)
  } else {
    Kmc <- K
  }
  K <- matrix(list(), ncol = nox + 1, nrow = nox + 1)
  K[1, 1] <- list(Kmc)
  Y.old <- Y
  scale.params <- list()
  if (preProcY == "mc" | preProcY == "uv" | preProcY == "pareto") {
    scale.params <- kopls::koplsScale(Y, center = "mc", scale = ifelse(preProcY == "mc", "none", preProcY))
    Y <- scale.params$x
  }
  to <- list()
  co <- list()
  so <- list()
  toNorm <- list()
  Tp <- list()
  Cp <- list()
  Bt <- list()
  tmp <- svd(t(Y) %*% K[1, 1][[1]] %*% Y, nu = A, nv = A)
  Cp <- tmp$u
  if (A > 1) {
    Sp <- diag(tmp$d[1:A])
    Sps <- diag(tmp$d[1:A]^(-1/2))
  } else {
    Sp <- tmp$d[1]
    Sps <- tmp$d[1]^(-1/2)
  }
  Up <- Y %*% Cp
  if (nox > 0) {
    for (i in 1:nox) {
      Tp[[i]] <- t(K[1, i][[1]]) %*% Up %*% Sps
      solve_res <- solve(t(Tp[[i]]) %*% Tp[[i]], tol = tol)
      Bt[[i]] <- solve_res %*% t(Tp[[i]]) %*% Up
      tmp <- svd(t(Tp[[i]]) %*% (K[i, i][[1]] - Tp[[i]] %*% t(Tp[[i]])) %*% Tp[[i]], nu = 1, nv = 1)
      co[[i]] <- tmp$u
      so[[i]] <- tmp$d[1]
      to[[i]] <- (K[i, i][[1]] - Tp[[i]] %*% t(Tp[[i]])) %*% Tp[[i]] %*% co[[i]] %*% so[[i]]^(-1/2)
      toNorm[[i]] <- c(sqrt(t(to[[i]]) %*% to[[i]]))
      to[[i]] <- to[[i]]/toNorm[[i]]
      K[1, i + 1][[1]] <- K[1, i][[1]] %*% (I - to[[i]] %*% t(to[[i]]))
      K[i + 1, i + 1][[1]] <- (I - to[[i]] %*% t(to[[i]])) %*% K[i, i][[1]] %*% (I - to[[i]] %*% t(to[[i]]))
    }
  }
  Tp[[nox + 1]] = t(K[1, nox + 1][[1]]) %*% Up %*% Sps
  Bt[[nox + 1]] = solve(t(Tp[[nox + 1]]) %*% Tp[[nox + 1]]) %*% t(Tp[[nox + 1]]) %*% Up
  sstotY <- sum(sum(Y * Y))
  F <- Y - Up %*% t(Cp)
  R2Y <- 1 - sum(sum(F * F))/sstotY
  EEprime <- K[nox + 1, nox + 1][[1]] - Tp[[nox + 1]] %*% t(Tp[[nox + 1]])
  sstotK <- sum(diag(K[1, 1][[1]]))
  R2X <- NULL
  R2XO <- NULL
  R2XC <- NULL
  R2Yhat <- NULL
  for (i in 1:(nox + 1)) {
    rss <- sum(diag(K[i, i][[1]] - Tp[[i]] %*% t(Tp[[i]])))
    R2X <- c(R2X, 1 - rss/sstotK)
    rssc <- sum(diag(K[1, 1][[1]] - Tp[[i]] %*% t(Tp[[i]])))
    R2XC <- c(R2XC, 1 - rssc/sstotK)
    rsso <- sum(diag(K[i, i][[1]]))
    R2XO <- c(R2XO, 1 - rsso/sstotK)
    Yhat <- Tp[[i]] %*% Bt[[i]] %*% t(Cp)
    R2Yhat <- c(R2Yhat, 1 - sum(sum((Yhat - Y)^2))/sstotY)
  }
  model <- list()
  model$Cp <- Cp
  model$Sp <- Sp
  model$Sps <- Sps
  model$Up <- Up
  model$Tp <- Tp
  model$T <- as.matrix(Tp[[nox + 1]])
  model$co <- co
  model$so <- so
  model$to <- to
  if (nox > 0) {
    model$To <- matrix(nrow = nrow(model$T), ncol = nox,
                       data = unlist(to), byrow = FALSE)
  }
  else {
    model$To <- NULL
  }
  model$toNorm <- toNorm
  model$Bt <- Bt
  model$A <- A
  model$nox <- nox
  model$K <- K
  model$EEprime <- EEprime
  model$sstot_K <- sstotK
  model$R2X <- R2X
  model$R2XO <- R2XO
  model$R2XC <- R2XC
  model$sstot_Y <- sstotY
  model$R2Y <- R2Y
  model$R2Yhat <- R2Yhat
  model$preProc <- list()
  model$preProc$K <- preProcK
  model$preProc$Y <- preProcY
  model$preProc$paramsY <- scale.params
  class(model) <- "kopls"
  return(model)
}

#' Get summary results of run_AMOPLS
#'
#' @inheritParams fun_outputs
#'
#' @import magrittr
#'
#' @export
fun_get_summary <- function(s) {
  s$output$Summary %>% {.[s$general$ee.names]}
}

