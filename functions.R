resample_function <- function(data = data,
                              formula = ~ .,
                              k = 3,
                              number_of_resamples = 15,
                              proportion_resample = 0.9,
                              starting_seed = 599) {
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

                       kmeans <- k_means(num_clusters = k) |>
                         fit({{formula}},
                             data = random_sample)

                       intermediate <- data.frame(random_sample$index,
                                                  extract_cluster_assignment(kmeans) |>
                                                    mutate(.cluster = as.character(.cluster)),
                                                  stringsAsFactors = FALSE)
                       colnames(intermediate) <- c("index", "cluster")

                       for (c in unique(intermediate$cluster)) {
                         idx <- intermediate[intermediate$cluster == c, ]$index

                         # with combinations (x, y) != (y, x), possibly breaking the matrix calc...
                         # need to set x = y and y = x for guarantee

                         idx <- sort(idx)
                         ones <- t(combn(idx, 2))
                         result_matrix[ones[, 1], ones[, 2]] <- 1

                         neg_one_idx <- expand.grid(idx, setdiff(random_sample$index, idx))
                         result_matrix[neg_one_idx[, 1], neg_one_idx[, 2]] <- -1
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

squared_distance_from_one <- function(mean_matrix = mean_matrix) {
  res_vec <- as.vector(mean_matrix)
  res_vec <- abs(res_vec[!is.na(res_vec)])
  return(sum((1 - res_vec)^2))
}

image_helper <- function(list_res = list_res) {
  lapply(list_res, function(x) image(x))
}
