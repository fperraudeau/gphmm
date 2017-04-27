#' Find last json files created.
#' 
#' \code{findLastJson} returns name of last json created.
#' 
#' @param sfx     - chr str, suffix of the json file, for example, 'parambt'.
#' @param path    - chr str, path of the directory to look for.
findLastJson <- function(sfx = 'gphmm',
                         path = paste0(system.file(package = 'gphmm'),
                                       '/training/')){
  dir = list.files(path, pattern = "^v.*", full.names = F)
  lastDir = paste0(path, 'v', max(as.numeric(gsub('v', '', dir))))
  json = list.files(lastDir)
  paste(lastDir, json[grepl(sprintf('.*%s.json', sfx), json)], sep = '/')
}

#' Find format of a file.
#' 
#' \code{findFormat} returns format of file.
#' 
#' @param file - chr str, name of file to extract the format from.
findFormat <- function(file){
  format = strsplit(file, '\\.')[[1]]
  format[length(format)]
}

#' Transform queries in format ShortRead to format ShortReadQ.
#' 
#' \code{ShortReadToShortReadQ} returns queries in ShortReadQ format.
#' 
#' @param queries - ShortRead object.
ShortReadToShortReadQ <- function(queries){
  if (class(queries) == 'ShortRead'){
    qv = sapply(width(queries), function(x) paste(rep('I', x), collapse=''))
    ShortReadQ(sread = sread(queries), quality = BStringSet(qv), id = id(queries))
  } else{
    queries
  }
}

#' Transform queries in format ShortReadQ to format ShortRead.
#' 
#' \code{ShortReadQToShortRead} returns queries in ShortRead format.
#' 
#' @param queries - ShortReadQ object.
ShortReadQToShortRead  <- function(queries){
  if (class(queries) == 'ShortReadQ'){
    ShortRead(sread(queries), id(queries))
  } else{
    queries
  } 
}

#' Read fasta or fastq file and returns ShortReadQ format.
#' 
#' \code{toFastq} returns queries in ShortReadQ format.
#' 
#' @param fPath - chr str, path to file with queries, can be a fasta or fastq file.
toFastq <- function(fPath){
  format = findFormat(fPath)
  if (format %in% c('fastq', 'fq')){
    readFastq(fPath)
  } else{
    queries = readFasta(fPath)
    ShortReadToShortReadQ(queries)
  }
}

#' Read fasta or fastq file and returns ShortRead format.
#' 
#' \code{toFasta} returns queries in ShortRead format.
#' 
#' @param fPath - chr str, path to file with queries, can be a fasta or fastq file.
toFasta <- function(fPath){
  format = findFormat(fPath)
  if (format %in% c('fastq', 'fq')){
    queries = readFastq(fPath)
    ShortReadQToShortRead(queries)
  } else{
    readFasta(fPath)
  }
}

#' Compute consensus sequence for a strain.
#' 
#' \code{computeConsensus} returns consensus sequence for a strain.
#' 
#' @param strain  - chr str, name of the wb strain contained in externalId.
#' @param ref     - DNAStringSet from dbFasta().
#' @param anno    - data frame from dbAnno().
computeConsensus <- function(strain, ref, anno, thr = .01){
  print(strain)
  id = anno[grepl(strain, anno$externalId), 'refGid']
  id = unique(id)
  if (length(id) == 0){
    id = anno[grepl(strain, paste(anno$genus, anno$species, sep = '_')), 'refGid']
    regex = paste(paste0('(?=.*', strsplit(strain, '_')[[1]], '.*)'), collapse='')
    idd = anno[grepl(regex, anno$desc, perl = T), 'refGid']
    id = unique(c(idd, id))
  }
  alleles = ref[id]
  if (length(alleles) == 1){
    as.character(alleles)
  }else{
    if (length(alleles) > 20){
      alleles = alleles[sample(1:length(alleles), 20)]
    }
    msa = msa(alleles)
    conMat = consensusMatrix(msa, as.prob = T)
    consensusString(conMat, ambiguityMap = "N", threshold = thr)
  } 
}


#' Randomly generate sequences.
#' 
#' \code{splitTrainTest} returns a DNAStringSet with generated sequences.
#' 
#' @param y          - any R object for which method length() has been defined.
#' @param proptrain  - float, percentage to be included in training set.
#' @param seed       - integer.
splitTrainTest <- function(y, proptrain = 0.8, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  toSample = 1:length(y)
  trainIdx = sample(toSample, round(proptrain*length(y)))
  testIdx = toSample[! toSample %in% trainIdx]
  list(train = trainIdx, test = testIdx)
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
generateRandomSequences <- function(n = 2,  meanLen = 10, sdLen = 0,
                                    seed = NULL, ncores = NULL){
  if (!is.null(seed)){
    stopifnot(class(seed) == 'numeric')
    set.seed(seed)
  }
  
  length = round(rnorm(n = n, mean = meanLen, sd = sdLen))
  
  if (!is.null(ncores)){
    stopifnot(class(ncores) == 'numeric')
    seqs = mclapply(1:n, function(i){
      paste(sample(c('A', 'C', 'G', 'T'), length[i], replace = T), collapse = '')
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
