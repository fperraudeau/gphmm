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


