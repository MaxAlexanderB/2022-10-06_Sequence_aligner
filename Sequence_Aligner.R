#####################################################
#######        COMPUTATIONAL BIOLOGY         ########
#######             HOMEWORK 1               ########
#####################################################
#                                                   #
# Implement the pairwise alignment algorithms       #
# Needleman-Wunsch and Smith-Waterman.              #
#                                                   #
#####################################################
#####################################################



init_score_matrix = function(nrow, ncol, local, score_gap) {
    # Initialize the score matrix with zeros.
  
    score_matrix <- matrix(0, nrow, ncol)
  
    # If the alignment is global, the leftmost column and the top row will have incremental gap scores,
    
    
    if (local == F){
      for (i in 1:nrow) {
        score_matrix[i,1] <- (i-1)*score_gap
      }
      
      for (i in 1:ncol) {
        score_matrix[1,i] <- (i-1)*score_gap
      }
    }

    # Return the initialized empty score matrix
    return(score_matrix)
}

init_path_matrix = function(nrow, ncol, local) {
    # Initialize the path matrix with empty values ("").
    
    path_matrix <- matrix("", nrow, ncol)
    
    # Additionally, for GLOBAL alignment (i.e. local==FALSE), make the first row
    # have "left" on all positions except 1st, and make the first column
    # have "up" on all positions except 1st.
    if (local == FALSE){
      for (i in 1:(nrow-1)) {
        path_matrix[i+1,1] <- "up"
      }
      
      for (i in 1:(ncol-1)) {
        path_matrix[1,i+1] <- "left"
      }
    }
    # Return the initialized empty path matrix
    return(path_matrix)
}

get_best_score_and_path = function(row, col, nucA, nucB, score_matrix, score_gap, score_match, score_mismatch, local) {
    # Compute the score and the best path for a particular position in the score matrix

    diag <- score_matrix[row-1, col-1] + if (nucA==nucB){score_match} else {score_mismatch}
    left <- score_matrix[row, col-1] + score_gap
    up <- score_matrix[row-1, col] + score_gap
    options <- c(diag, left, up)

    best_score <- max(diag, left, up)
    
    if (local==T & best_score < 0) {best_path <- "-"} else {
      best_path <- if (which(options==best_score)[1] == 1){"diag"} else if (which(options==best_score)[1] == 2) {"left"} else if (which(options==best_score)[1] == 3) {"up"} #does not pick randomly if 2 best scores (prefers diagonal)
    }
    # Return the best score for the particular position in the score matrix

        return(list("score"=best_score, "path"=best_path))
}

fill_matrices = function(seqA, seqB, score_gap, score_match, score_mismatch, local, score_matrix, path_matrix) {
    # Compute the full score and path matrices
    for (i in 2:(nchar(seqA)+1)){
      for (j in 2:(nchar(seqB)+1)){
        score_matrix[i,j] <- get_best_score_and_path(i,j,substr(seqA,i-1,i-1), substr(seqB,j-1,j-1), score_matrix, score_gap, score_match, score_mismatch, local)$score 
      }
    }
    score_matrix
    for (i in 2:(nchar(seqA)+1)){
      for (j in 2:(nchar(seqB)+1)){
        path_matrix[i,j] <- get_best_score_and_path(i,j,substr(seqA,i-1,i-1), substr(seqB,j-1,j-1), score_matrix, score_gap, score_match, score_mismatch, local)$path
      }
    }
    path_matrix
    
    # Return the full score and path matrices
    return(list("score_matrix"=score_matrix, "path_matrix"=path_matrix))
}

get_best_move = function(nucA, nucB, path, row, col) {
    # Compute the aligned characters at the given position in the score matrix and return the new position,
    if (path=="diag"){
      char1 <- nucA
      char2 <- nucB
      newrow <- row-1
      newcol <- col-1
    } else if (path=="up"){
      char1 <- nucA
      char2 <- "-"
      newrow <- row-1
      newcol <- col
    } else if (path=="left") {
      char1 <- "-"
      char2 <- nucB
      newrow <- row
      newcol <- col-1
    }
    # Return the new row and column and the aligned characters
    return(list("newrow"=newrow, "newcol"=newcol, "char1"=char1, "char2"=char2))
}

get_best_alignment = function(seqA, seqB, score_matrix, path_matrix, local) {
    # Return the best alignment from the pre-computed score matrix
    if (local==F){
      row <- nchar(seqA)+1
      col <- nchar(seqB)+1
      score <- score_matrix[row,col]
    } else {
      row <- which(score_matrix==max(score_matrix), arr.ind = T)[1,][1]
      col <- which(score_matrix==max(score_matrix), arr.ind = T)[1,][2]
      score <- max(score_matrix)
    row
    col  
    score
    }
  
    alignment <- c("","")
    
    path_matrix[3,3]
    
    while (col > 1 || row > 1){
      
      if (path_matrix[row,col]=="-"){break}
      #get_best_move (substr(seqA,row-1,row-1), substr(seqB,col-1,col-1), path_matrix[row,col], row, col)
      
      alignment[1] <- paste0(get_best_move (substr(seqA,row-1,row-1), substr(seqB,col-1,col-1), path_matrix[row,col], row, col)$char1, alignment[1])
      alignment[2] <- paste0(get_best_move (substr(seqA,row-1,row-1), substr(seqB,col-1,col-1), path_matrix[row,col], row, col)$char2, alignment[2]) 
      
      currentrow <- get_best_move (substr(seqA,row-1,row-1), substr(seqB,col-1,col-1), path_matrix[row,col], row, col)$newrow
      currentcol <- get_best_move (substr(seqA,row-1,row-1), substr(seqB,col-1,col-1), path_matrix[row,col], row, col)$newcol
      
      score <- score+score_matrix[row,col]
      
      row <- currentrow
      col <- currentcol 
      
    }

    # Return the best score and alignment (or one thereof if there are multiple with equal score)
    return(list("score"=score, "alignment"=alignment))
}

align = function(seqA, seqB, score_gap, score_match, score_mismatch, local) {
    # Align the two sequences given the scoring scheme
    # Initialize score and path matrices
    score_matrix <- init_score_matrix(nchar(seqA)+1, nchar(seqB)+1,local, score_gap)
    path_matrix <- init_path_matrix(nchar(seqA)+1, nchar(seqB)+1,local)
    
    score_matrix
    path_matrix
    
    # Fill in the matrices with scores and paths using dynamic programming
    score_matrix <- fill_matrices(seqA,seqB,score_gap, score_match, score_mismatch, local, score_matrix, path_matrix)$score_matrix
    path_matrix <- fill_matrices(seqA,seqB,score_gap, score_match, score_mismatch, local, score_matrix, path_matrix)$path_matrix
    
    score_matrix
    path_matrix
    
    # Get the best score and alignment (or one thereof if there are multiple with equal score)
    
    get_best_alignment(seqA, seqB, score_matrix, path_matrix, local)
    
    result <- get_best_alignment(seqA, seqB, score_matrix, path_matrix, local)
    
    # Return the best score and alignment (or one thereof if there are multiple with equal score)

    return(result)
}

test_align = function() {
    seqA = "TCACACTAC"
    seqB = "AGCACAC"
    score_gap = -2
    score_match = +3
    score_mismatch = -1
    local = T
    result = align(seqA, seqB, score_gap, score_match, score_mismatch, local)
    print(result$alignment)
    print(result$score)
}

test_align()








