#' Create a function to compute gphmm probabilities during the training.
#' 
#' \code{makeGphmmPerRead} returns a function to compute gphmm probabilities for each row of the csv file.
#' 
#' @param seqs  - DNAStringSet with DNA sequences used for the training.
#' @param csv   - data.frame with first column = queries, second column = reference sequences, third column = qv
makeGphmmPerRead <- function(seqs, csv){
  gphmmPerRead <- function(i, parameters){
    read = as.character(seqs[csv[i, 1]])
    ref = as.character(seqs[csv[i, 2]])
    path = computegphmm(read, ref, parameters = parameters, output = 'long')
    et = computeCounts(path)
    et$vit = path$V
    et$qv = csv[i, "qv"]
    et 
  }
  gphmmPerRead
}


#' Initial parameters for gphmm.
#' 
#' \code{initializeGphmm} returns initial set of parameters for gphmm.
#' 
initializeGphmm <- function(){
  qR = rep(.25, 4)
  parameters = list(qR = qR,
                    qX = rep(.25, 4),
                    qY = rep(.25, 4),
                    deltaX = c(-10, -1),
                    deltaY = c(-10, -1),
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

#' Compute GPHMM log probability.
#' 
#' \code{computegphmm} returns GPHMM log probability for a read and a reference sequence.
#' 
#' @param read         - chr str.
#' @param ref          - chr str.
#' @param parameters   - list of GPHMM parameters, can be created by initializeGphmm().
#' @param qv           - float, quality value.
#' @param output       - if 'long', output is a list with path, read, ref and log GPHMM proba, else output is just the log GPHMM proba.
computegphmm <- function(read, ref, parameters, qv = 20, output = "short"){
  p = parameters
  p$pp = log(p$pp)
  p$qX = log(p$qX)
  p$qY = log(p$qY)
  p$pp = cbind(p$pp, matrix(0, ncol=1, nrow=4))
  p$pp = rbind(p$pp, matrix(0, ncol=5, nrow=1))
  p$qX = c(p$qX, 0)
  p$qY = c(p$qY, 0)
  p[['deltaX']] = 1/(1+exp(-sum(p$deltaX * c(1, qv))))
  p[['deltaY']] = 1/(1+exp(-sum(p$deltaY * c(1, qv))))
  
  states = c("A", "C", "G", "T", "N")
  colnames(p$pp) = rownames(p$pp) = names(p$qX) = names(p$qY) = states
  
  res = calculategphmm(read, ref, p$tau, p$pp, p$qX, p$qY, p$deltaX, p$deltaY, p$epsX, p$epsY)
  
  if (output == "long"){
    res$path = paste0(res$path, collapse = '')
    res$read = read
    res$ref = ref
    return(res)
  } else {
    return(res$V)
  } 
}


#' Compute counts for emissions and transitions.
#' 
#' \code{computeCounts} returns a list with emission and transition counts.
#' 
#' @param gphmm - output of computegphmm with arg output='long'.
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
#' \code{computeDelta} returns list of coefficients from logit regression for insertion and deletion rates.
#' 
#' @param mat  - matrix with counts for insertions and deletions for each read.
#' @param qv   - vector with quality values for the reads.
computeDelta <- function(mat, qv){
  mat = as.data.frame(mat)
  mat$failI = mat$MM + mat$MD
  mat$failD = mat$MM + mat$MI
  lrI = glm(as.matrix(mat[,c('MI', 'failI')]) ~ qv, family = binomial())
  lrD = glm(as.matrix(mat[,c('MD', 'failD')]) ~ qv, family = binomial())
  list(I = coefficients(lrI), D = coefficients(lrD))
}

#' Compute normalized emission probabilities.
#' 
#' \code{computeQ} returns normalized emission probabilities.
#' 
#' @param emission - vector of length at least 4.
computeQ <- function(emission){
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
#' \code{computeEps} returns normalized transition probabilities.
#' 
#' @param transition - named vector with at least elements 'DD','DM','II','IM'.
#' @param state      - str, "X" or "Y".
computeEps <- function(transition, state = "X"){
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
computeGphmmParam <- function(emiTrans){
  emi = lapply(lapply(c(1:4), function(i) lapply(emiTrans, '[[', i)), function(x) Reduce('+', x))
  names(emi) = c("emissionM", "emissionD", "emissionI","transition")
  qR = colSums(emi$emissionM)[1:4]
  qR = qR/sum(qR)
  qX = computeQ(emi$emissionI)
  qY = computeQ(emi$emissionD)
  pp = .computeP(emi$emissionM)
  tau = 1/1500
  epsX = computeEps(emi$transition, state = 'X')
  epsY = computeEps(emi$transition, state = 'Y') 
  transList = lapply(emiTrans, '[[', 4)
  transMat = do.call(rbind, transList)
  delta = computeDelta(transMat, qv = sapply(emiTrans, '[[', 6))
  deltaX = delta[['I']]
  deltaY = delta[['D']]
  list(qR = qR, qX = qX, qY = qY, pp = pp, deltaX = deltaX,
       deltaY = deltaY, epsX = epsX, epsY = epsY, tau = tau)
}

