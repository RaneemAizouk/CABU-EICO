generate_sequential_ids <- function(start_id, num_individuals) {
  individual_ids <- numeric(num_individuals)
  for (i in 1:num_individuals) {
    individual_ids[i] <- start_id + i - 1
  }
  return(individual_ids)
}

