resample_function <- function(data = data,
                              k = 3,
                              number_of_resamples = 15,
                              proportion_resample = 0.9,
                              starting_seed = 599) {
  tictoc::tic()
  data <- data |>
    drop_na()
  # data <- data[sample(nrow(data)), ]
  data$index <- 1:nrow(data)
  results <- list()

  # For loop
  results <- foreach(i = 1:number_of_resamples,
                     .packages = c("tidyclust", "dplyr", "tidymodels")) %dopar% {
                       # Reproducibility over parallel
                       set.seed(starting_seed + i)
                       result_matrix <- matrix(0, nrow = nrow(data), ncol = nrow(data))
                       # 80% resample
                       random_sample <- data |>
                         filter(index %in% sample(index, proportion_resample * nrow(data)))

                       # Run k-means on resample
                       # kmeans <- k_means(num_clusters = k) |>
                       #   fit(~ bill_length_mm + flipper_length_mm,
                       #       data = random_sample)

                       kmeans <- k_means(num_clusters = k) |>
                         fit(~ X1 + X2,
                             data = random_sample)

                       intermediate <- data.frame(random_sample$index,
                                                  extract_cluster_assignment(kmeans) |>
                                                    mutate(.cluster = as.character(.cluster)),
                                                  stringsAsFactors = FALSE)
                       colnames(intermediate) <- c("index", "cluster")

                       # I don't know how to vectorize this :(
                       for (n in 1:(nrow(intermediate) - 1)) {
                         for (m in (n + 1):nrow(intermediate)) {
                           # if the points are the same then set to NA
                           if (intermediate[n, ]$index != intermediate[m, ]$index) {
                             # get first pointer cluster
                             first <- intermediate[n, ]$cluster
                             # get second pointer cluster
                             second <- intermediate[m, ]$cluster

                             # if clusters are the same
                             if (first == second) {
                               result_matrix[intermediate[n, ]$index, intermediate[m, ]$index] <- 1
                             } else {
                               # if clusters are different
                               result_matrix[intermediate[n, ]$index, intermediate[m, ]$index] <- -1
                             }
                           }
                         }
                       }
                       result_matrix[lower.tri(result_matrix, diag = TRUE)] <- NA
                       result_matrix
                     }
  tictoc::toc()
  return(results)
}

one_k_mean_matrix <- function(data = data,
                              k = k,
                              starting_seed = 599,
                              number_of_resamples = 15,
                              proportion_resample = 0.9) {
  result <- resample_function(data = data,
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
