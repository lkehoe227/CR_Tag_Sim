# functions 

clipped_tag <- function(df) {
  # Select all L* and D* columns
  target_cols <- grep("^(L|D)\\d+$", names(df), value = TRUE)
  
  # Subset the matrix to L/D columns
  mat <- as.matrix(df[, target_cols])
  
  # Apply function row-wise
  mat_clipped <- t(apply(mat, 1, function(row) {
    first2 <- which(row == 2)[1]
    if (!is.na(first2) && first2 < length(row)) {
      row[(first2 + 1):length(row)][row[(first2 + 1):length(row)] == 2] <- 0
    }
    return(row)
  }))
  
  # Assign modified matrix back to the original dataframe
  df[, target_cols] <- mat_clipped
  
  return(df)
}
