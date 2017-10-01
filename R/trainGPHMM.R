#' Create a function to compute gphmm probabilities during the training.
#' 
#' \code{makeGphmmPerRead} returns a function to compute gphmm probabilities for each row of the csv file.
#' 
#' @param seqs  - DNAStringSet with DNA sequences used for the training.
#' @param csv   - data.frame with first column = queries, second column = reference sequences, third column = qv
#' @examples 
#' library(Biostrings)
#' seqs <- DNAStringSet(c(a='ATGC', b = 'ATGG', c = 'ATGT'))
#' csv <- data.frame(queries = c('a', 'b'), refs = c('c', 'c'))
#' makeGphmmPerRead(seqs, csv)
makeGphmmPerRead <- function(seqs, csv){
  if (ncol(csv) == 3){
    stopifnot(class(csv[, 3]) %in% c('numeric', 'character', 'integer'))
    if (class(csv[, 3]) == 'character') csv[, 3] = as.integer(csv[, 3])
  }
  gphmmPerRead <- function(i, parameters){
    read = as.character(seqs[csv[i, 1]])
    ref = as.character(seqs[csv[i, 2]])
    qv = ifelse(ncol(csv) == 3, csv[i, 3], 20)
    path = computegphmm(read, ref, qv = qv, parameters = parameters, output = 'long') 
    et = computeCounts(path)
    et$vit = path$V
    et$qv = qv
    et 
  }
  gphmmPerRead
}

#' Compute counts for emissions and transitions.
#' 
#' \code{computeCounts} returns a list with emission and transition counts.
#' 
#' @param gphmm - output of computegphmm with arg output='long'.
#' @examples 
#' gphmmOut <- computegphmm('ATCG', 'ATG', output = 'long')
#' computeCounts(gphmmOut)
computeCounts <- function(gphmm){
  read = strsplit(gphmm$read, "")[[1]]
  ref = strsplit(gphmm$ref, "")[[1]]
  path = gphmm$path
  path = strsplit(path, "")[[1]]
  states = c("A", "C", "G", "T", "N")
  emissionM = matrix(0, ncol = 5, nrow = 5)
  emissionI = emissionD = rep(0, 5)
  names(emissionI) = names(emissionD) = colnames(emissionM) =
    rownames(emissionM) = states
  
  # emission proba
  i = j = 0
  for (k in 1:length(path)){
    if (path[k] == "M"){
      i = i + 1
      j = j + 1
      emissionM[read[i], ref[j]] = emissionM[read[i], ref[j]] + 1
    } else if (path[k] == "I"){
      i = i + 1
      emissionI[read[i]] = emissionI[read[i]] + 1
    } else if (path[k] == "D"){
      j = j + 1
      emissionD[ref[j]] = emissionD[ref[j]] + 1
    }
    if (i>0 & j>0){
      if (read[i] == "N" | ref[j] == "N"){
        path[k] = ""
      }
    }
  }
  
  # transition proba
  patterns = c("MM", "II", "IM", "DD", "DM", "MI", "MD")
  transition = stri_count_fixed(paste0(path, collapse = ""), patterns, overlap = TRUE)
  names(transition) = patterns
  
  return(list(emissionM = emissionM, 
              emissionD = emissionD,
              emissionI = emissionI,
              transition = transition))
}

#' Compute insertion and deletion rates, i.e., probability to transition from
#' match/mismatch state to insertion and deletion state.
#' 
#' \code{.computeDelta} returns list of coefficients from logit regression for insertion and deletion rates.
#' 
#' @param mat  - matrix with counts for insertions and deletions for each read.
#' @param qv   - vector with quality values for the reads.
.computeDelta <- function(mat, qv){
  mat = as.data.frame(mat)
  mat$failI = mat$MM + mat$MD
  mat$failD = mat$MM + mat$MI
  lrI = glm(as.matrix(mat[,c('MI', 'failI')]) ~ qv, family = binomial())
  lrD = glm(as.matrix(mat[,c('MD', 'failD')]) ~ qv, family = binomial())
  list(I = coefficients(lrI), D = coefficients(lrD))
}

#' Compute normalized emission probabilities.
#' 
#' \code{.computeQ} returns normalized emission probabilities.
#' 
#' @param emission - vector of length at least 4.
.computeQ <- function(emission){
  emission = emission[1:4]
  if (sum(emission) == 0) emission else emission / sum(emission)
}

#' Compute normalized emission probabilities.
#' 
#' \code{.computeP} returns normalized emission probabilities.
#' 
#' @param emission - data frame of dim at least 4x4.
.computeP <- function(emission){
  emission = emission[1:4,1:4]
  emission / rowSums(emission)
}

#' Compute normalized transition probabilities.
#' 
#' \code{.computeEps} returns normalized transition probabilities.
#' 
#' @param transition - named vector with at least elements 'DD','DM','II','IM'.
#' @param state      - str, "X" or "Y".
.computeEps <- function(transition, state = "X"){
  if (state == "Y"){
    normD = sum(transition[c("DD", "DM")])
    if (normD != 0) transition["DD"] / normD else 0
  } else if (state == "X"){
    normI = sum(transition[c("II", "IM")])
    if (normI != 0) transition["II"] / normI else 0
  } 
}

#' Compute gphmm parameters from counts.
#' 
#' \code{computeGphmmParam} returns a list with gphmm parameters.
#' 
#' @param emiTrans - list with emission and transition counts.
#' @examples 
#' library(Biostrings)
#' seqs <- DNAStringSet(c(a='ATGC', b = 'ATGG', c = 'ATGT'))
#' csv <- data.frame(queries = c('a', 'b'), refs = c('c', 'c'))
#' gphmmPerRead <- makeGphmmPerRead(seqs, csv)
#' parameters <- initializeGphmm()
#' counts <- lapply(1:nrow(csv), function(i) gphmmPerRead(i, parameters))
#' computeGphmmParam(counts)
computeGphmmParam <- function(emiTrans){
  emi = lapply(lapply(c(1:4), function(i) lapply(emiTrans, '[[', i)), function(x) Reduce('+', x))
  names(emi) = c("emissionM", "emissionD", "emissionI","transition")
  qR = colSums(emi$emissionM)[1:4]
  qR = qR/sum(qR)
  qX = .computeQ(emi$emissionI)
  qY = .computeQ(emi$emissionD)
  pp = .computeP(emi$emissionM)
  tau = 1/1500
  epsX = .computeEps(emi$transition, state = 'X')
  epsY = .computeEps(emi$transition, state = 'Y') 
  transList = lapply(emiTrans, '[[', 4)
  transMat = do.call(rbind, transList)
  delta = .computeDelta(transMat, qv = sapply(emiTrans, '[[', 6))
  deltaX = delta[['I']]
  deltaY = delta[['D']]
  list(qR = qR, qX = qX, qY = qY, deltaX = deltaX,
       deltaY = deltaY, epsX = epsX, epsY = epsY, tau = tau, eta = tau, pp = pp)
}

