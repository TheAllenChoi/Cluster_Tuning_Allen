resample_function <- function(data = data,
                              formula = ~ .,
                              k = 3,
                              number_of_resamples = 15,
                              proportion_resample = 0.9,
                              starting_seed = 599,
                              algorithm = "kmeans") {
  # tictoc::tic()
  data <- data |>
    drop_na()
  data$index <- 1:nrow(data)
  results <- list()

  # For loop
  results <- foreach(i = 1:number_of_resamples,
                     .packages = c("tidyclust", "dplyr", "tidymodels")) %dopar% {
                       # Reproducibility over parallel
                       set.seed(starting_seed + i)
                       result_matrix <- matrix(0, nrow = nrow(data), ncol = nrow(data))

                       random_sample <- data |>
                         filter(index %in% sample(index, proportion_resample * nrow(data)))

                       if (algorithm == "kmeans") {
                         cluster_assigned <- k_means(num_clusters = k) |>
                           fit({{formula}},
                               data = random_sample |>
                                 select(-index))
                       } else if (algorithm == "hier_clust") {
                         cluster_assigned <- hier_clust(num_clusters = k) |>
                           fit({{formula}},
                               data = random_sample |>
                                 select(-index))
                       } else {
                         stop("Algorithm is not supported. Please select one of (kmeans, hier_clust)")
                       }


                       intermediate <- data.frame(random_sample$index,
                                                  extract_cluster_assignment(cluster_assigned) |>
                                                    mutate(.cluster = as.character(.cluster)),
                                                  stringsAsFactors = FALSE)
                       colnames(intermediate) <- c("index", "cluster")

                       for (c in unique(intermediate$cluster)) {
                         idx <- intermediate[intermediate$cluster == c, ]$index

                         # Check that the index list is not length 1 for purposes of n > m for combn
                         # This also fixes a bug in the previous code where idx was wrongly interpreted as single numerical
                         # value as param inside combn
                         if (length(idx) > 1) {
                           idx <- sort(idx)
                           ones <- t(combn(idx, 2))
                           result_matrix[cbind(ones[, 1], ones[, 2])] <- 1
                         }

                         # Note for future self: This calculates the indices, even across y=x line
                         # which means we are doing double the work for this part.
                         neg_one_idx <- expand.grid(idx, setdiff(random_sample$index, idx))
                         result_matrix[cbind(neg_one_idx[, 1], neg_one_idx[, 2])] <- -1
                       }
                       result_matrix[lower.tri(result_matrix, diag = TRUE)] <- NA
                       result_matrix
                     }
  # tictoc::toc()
  return(results)
}

one_k_mean_matrix <- function(data = data,
                              formula = ~.,
                              k = k,
                              starting_seed = 599,
                              number_of_resamples = 15,
                              proportion_resample = 0.9) {
  result <- resample_function(data = data,
                              formula = formula,
                              k = k,
                              starting_seed = starting_seed,
                              number_of_resamples = number_of_resamples,
                              proportion_resample = proportion_resample)

  # Absolute value the entire matrix, for all matrices
  absolute_sum <- lapply(result, abs)

  # Get the sum of absolute values in the final matrix
  absolute_final <- absolute_sum |>
    reduce(`+`)

  # Get the sum of values in the final matrix
  numerator <- result |>
    reduce(`+`)

  # Returns the matrix of zero-removed means.
  # The zero-removed mean is the sum over the (length - the number of zeroes).
  # The length - the number of zeroes is equivalent to the sum over absolute values
  return(numerator / absolute_final)
}

mean_matrix <- function(list_of_matrices = list_of_matrices) {
  # Absolute value the entire matrix, for all matrices
  absolute_sum <- lapply(list_of_matrices, abs)

  # Get the sum of absolute values in the final matrix
  absolute_final <- absolute_sum |>
    reduce(`+`)

  # Get the sum of values in the final matrix
  numerator <- list_of_matrices |>
    reduce(`+`)

  # Returns the matrix of zero-removed means.
  # The zero-removed mean is the sum over the (length - the number of zeroes).
  # The length - the number of zeroes is equivalent to the sum over absolute values
  return(numerator / absolute_final)
}

squared_distance_from_one <- function(mean_matrix = mean_matrix) {
  res_vec <- as.vector(mean_matrix)
  res_vec <- abs(res_vec[!is.na(res_vec)])
  return(sum((1 - res_vec)^2))
}

image_helper <- function(list_res = list_res) {
  lapply(list_res, function(x) image(x))
}
