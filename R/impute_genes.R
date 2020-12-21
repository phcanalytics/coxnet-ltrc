#' Select genes from dataset based on criteria
#' 
#' To be included a gene must satisfy the `max_missing` and
#' `min_prop` thresholds, as well as have a veriance that is non-zero.
#' @param data a `data.frame` containing the mutation information.
#' @param max_missing Genes with a proportion missing greater than 
#' `max_missing` are excluded.
#' @param min_prop Minimum proportion of patients with gene alteration. 
#' @return A character vector containing the genes meeting the criteria.
#' @export
select_genes <- function(data, max_missing = 0.3, min_prop = 0.001){
  
  # Missing threshold
  missingness <- apply(data, 2, function(x) sum(is.na(x))/length(x))
  data <- data[, missingness < max_missing]
  
  # Variance is non-zero
  variance <- apply(data, 2, function (x) var(x, na.rm = TRUE))
  data <- data[, variance > 0]
  
  # Must occur in certain proportion of population
  prop <- apply(data, 2, function (x) mean(x, na.rm = TRUE))
  data <- data[, prop > min_prop]
  
  # Return
  return(colnames(data))
}


#' Impute missing data via KNN
#' 
#' @param data a `data.frame` containing the variables to inpute
#' @param k Passed to `impute::impute.knn()`.
#' @param seed Passed to `rng.seed` in `impute::impute.knn()`.
#'  
#' @return An object of class `tibble` having missing values replaced with 
#' imputed values via k nearest neighbors.The random seed is also returned
#' as an attribute of the data frame for reproducibility.
#' 
#' @export
impute_knn <- function(data, k = 20, seed = NULL){
  if (is.null(seed)){
    if(!exists(".Random.seed")) set.seed(NULL)
    seed <- .Random.seed
  }
  
  result <- t(impute::impute.knn(t(data), k = k, colmax = 1, 
                                 rng.seed = seed)$data) %>%
    as.data.frame()
  attr(result, "seed") <- seed
  return(result)
}

impute_genes <- function(train, test) {
  is_gene <- grepl("^((SV)|(RE)|(CN))_", colnames(train))
  selected_genes <- train[, is_gene] %>% 
    select_genes()
  
  # First impute genes for training set  
  imputed_genes_train <- train[, selected_genes] %>%
    impute_knn()
  
  # Then impute both training and test data
  imputed_genes_test <- impute_knn(
    bind_rows(
      train[, selected_genes],
      test[, selected_genes]
    )
  ) %>%
    tail(., n = -nrow(train))
  
  # Update training and test data
  add_imputed_genes <- function(data, imputed_genes, is_gene){
    bind_cols(
      data[, !is_gene],
      imputed_genes
    )
  }
  
  train <- add_imputed_genes(train, imputed_genes_train, is_gene)
  test <- add_imputed_genes(test, imputed_genes_test, is_gene)
  
  return(list(train = train, test = test, selected_genes = selected_genes))
}