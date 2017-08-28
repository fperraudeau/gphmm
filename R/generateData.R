#' Initial parameters for gphmm.
#' 
#' \code{initializeGphmm} returns initial set of parameters for gphmm.
#' 
initializeGphmm <- function(){
  qR = rep(.25, 4)
  parameters = list(qR = qR,
                    qX = rep(.25, 4),
                    qY = rep(.25, 4),
                    deltaX = c(-.2, -.2),
                    deltaY = c(-.2, -.2),
                    epsX = 0.1,
                    epsY = 0.1,
                    tau = 1/1500,
                    eta = 1/1500,
                    pp = matrix(c(0.91, 0.03, 0.03, 0.03,
                                  0.03, 0.91, 0.03, 0.03,
                                  0.03, 0.03, 0.91, 0.03,
                                  0.03, 0.03, 0.03, 0.91),
                                nrow = 4, ncol = 4))
  parameters
}

#' Generate a noisy read from a true sequence.
#' 
#' \code{generateRead} returns a list with the noisy read, the true sequence, the path, and the phred quality score.
#' 
#' @param seq        - character vector of true sequence.
#' @param qv         - integer, wanted phred quality score for the read.
#' @param seed       - integer, seed for reproducibility.
#' @param paramgphmm - list of parameters.

generateRead <- function(seq = 'ATGCGGATCG', qv = NULL, seed = NULL, 
                          paramgphmm = initializeGphmm()){
  #seed
  if (!is.null(seed)) set.seed(seed)
  
  # param gphmm
  paramgphmm[['deltaX']] = 1/(1+exp(-sum(paramgphmm$deltaX * c(1, qv))))
  paramgphmm[['deltaY']] = 1/(1+exp(-sum(paramgphmm$deltaY * c(1, qv))))
  stopifnot(1 - paramgphmm$deltaX - paramgphmm$deltaY - paramgphmm$tau > 0)
  stopifnot(1 - paramgphmm$epsX - paramgphmm$tau > 0)
  stopifnot(paramgphmm$deltaX > 0)
  stopifnot(paramgphmm$deltaY > 0)
  
  nucleotides = c("A", "C", "G", "T")
  colnames(paramgphmm$pp) = rownames(paramgphmm$pp) = names(paramgphmm$qX) = names(paramgphmm$qY) = nucleotides
  transition_mat = matrix(c(1 - paramgphmm$deltaX - paramgphmm$deltaY - paramgphmm$tau,
                            1 - paramgphmm$deltaX - paramgphmm$deltaY - paramgphmm$tau,
                            1 - paramgphmm$epsX - paramgphmm$tau,
                            1 - paramgphmm$epsY - paramgphmm$tau,
                            paramgphmm$deltaX, paramgphmm$deltaX,
                            paramgphmm$epsX, 0, 
                            paramgphmm$deltaY, paramgphmm$deltaY,
                            0, paramgphmm$epsY), ncol = 3)
  colnames(transition_mat) = c('M', 'I', 'D')
  rownames(transition_mat) = c('Ini', 'M', 'I', 'D')
  
  # seq
  if (is.null(qv)) qv = 20
  seq = as.character(seq)
  seq = strsplit(seq, '')[[1]]
  
  # initialization
  state = 'Ini'
  state = sample(colnames(transition_mat), 1, prob = transition_mat[state, ])
  states = c()
  read = c()
  i = j = 0
  
  # simulate
  while(j < length(seq)){
    if (state == 'M'){
      i = i + 1
      j = j + 1
      read[i] = sample(rownames(paramgphmm$pp), 1, prob = paramgphmm$pp[, seq[j]])
    } else if (state == 'I'){
      i = i + 1
      read[i] = sample(names(paramgphmm$qX), 1, prob = paramgphmm$qX)
    } else{
      j = j + 1
    }
    states = c(states, state)
    state = sample(colnames(transition_mat), 1, prob = transition_mat[state, ])
  }
  
  list(read = paste(read, collapse = ''), ref = paste(seq, collapse = ''),
       path = paste(states, collapse = ''), tabstates = table(states), qv = qv)
}

#' Split randomly any R object for which method length() has been defined into a train set and a test set.
#' 
#' \code{generateRandomSequences} returns a list with indices for train and test sets.
#' 
#' @param n          - int, number of sequences to generate.
#' @param meanLen    - float, mean of the length distribution (gaussian) of the generated sequences.
#' @param sdLen      - float, sd of the length distribution (gaussian) of the generated sequences.
#' @param seed       - int, when the same seed and parameters are used, exact same sequences are generated.
#' @param ncores     - int, number of cores to use.
#' @param prob       - vector of length 4 with the probability for the 4 nucleotides (A, C, G, T).
generateRandomSequences <- function(n = 2,  meanLen = 10, sdLen = 0,
                                    seed = NULL, ncores = NULL,
                                    prob = rep(.25, 4)){
  if (!is.null(seed)){
    stopifnot(class(seed) == 'numeric')
    set.seed(seed)
  }
  
  length = round(rnorm(n = n, mean = meanLen, sd = sdLen))
  stopifnot(sum(prob) != 0)
  if (sum(prob) != 1) prob = prob/sum(prob)
  
  if (!is.null(ncores)){
    stopifnot(class(ncores) == 'numeric')
    seqs = mclapply(1:n, function(i){
      paste(sample(c('A', 'C', 'G', 'T'), length[i], replace = T, prob = prob), collapse = '')
    }, mc.cores = ncores)
    seqs = unlist(seqs)
  } else{
    seqs = sapply(1:n, function(i){
      paste(sample(c('A', 'C', 'G', 'T'), length[i], replace = T), collapse = '')
    })
  }
  
  seqs = DNAStringSet(seqs)
  names(seqs) = paste0('s', 1:n)
  seqs
}
