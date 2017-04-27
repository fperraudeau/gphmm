#!/usr/bin/env Rscript

suppressMessages({
  require(docopt)
  require(gphmm)
  require(Biostrings)
  require(jsonlite)
})

'Estimate.

Usage: 
gphmm compute <fasta> <csv> [--param=<param>] [out=<out>] [options]
gphmm train <fasta> <csv> [options] [--llthreshold=<l>] [--maxit=<i>] 

Options:
<fasta>                  path to fasta file with DNA sequences.
<csv>                    path to csv file.
--param=<param>          path to json file with GPHMM parameters [Default: ].
--out=<out>              path to the output file [Default: ].
--llthreshold=<l>        algo stops when diff. of log ll between 2 iterations is below the thr [Default: 0.00001].
--maxit=<i>              maximum number of iterations [Default: 10].
--seed=<seed>            integer, seed.
--ncores=<n>             integer, number of cores [Default: 0].
--verbose                if TRUE, print things.' -> doc

## Read in args.
opt              = docopt(doc)
compute          = opt$compute
train            = opt$train
fasta            = opt$fasta
csvPath          = opt$csv

param            = opt$param
out              = opt$out 

llthreshold      = as.numeric(opt$llthreshold)
maxit            = as.numeric(opt$maxit)

seed             = opt$seed
ncores           = as.numeric(opt$ncores)
verbose          = opt$verbose


# todo
stopifnot(train | compute)
todo = ifelse(compute, 'compute', 'train')
if (verbose) print(sprintf('We are going to : %s', todo))

# number of cores
if (ncores == 0) ncores = max(1, detectCores() - 2)
if (verbose) print(sprintf('Number of cores used : %s', ncores))

# start timing computation
start = Sys.time()

if (compute){
  
  if (param == '') paramgphmm = fromJSON(findLastJson())
  seqs = readDNAStringSet(fasta)
  csv = read.table(csvPath, stringsAsFactors = F)
  if (ncol(csv) == 2) csv[, 'qv'] = 20
  gphmmProb = mclapply(1:nrow(csv), function(i){
    read = as.character(seqs[csv[i, 1]])
    ref = as.character(seqs[csv[i, 2]])
    qv = csv[i, 3]
    computegphmm(read, ref, qv = qv, parameters = paramgphmm)
  }, mc.cores = ncores)
  csv[,'gphmm'] = unlist(gphmmProb)
  if (out == '') out = gsub('.csv', '_gphmm.csv', csvPath)
  if (verbose) print(sprintf('Writting output to : %s', out))
  write.table(csv, out)  
  
} else if (train){
  
  seqs = readDNAStringSet(fasta)
  csv = read.table(csvPath, stringsAsFactors = F)
  stopifnot(ncol(csv) == 3)
  
  # initialization
  lldiff = 1000
  ll = 0
  it = 0
  llvect = NULL
  parameters = initializeGphmm()
  gphmmPerRead = makeGphmmPerRead(seqs, csv)
  
  #train
  while( lldiff > llthreshold & it < maxit){
    N = nrow(csv)
    it = it + 1
    counts = mclapply(1:N, function(i) gphmmPerRead(i, parameters), mc.cores = ncores)
    parameters = computeGphmmParam(counts)
    llnew = Reduce('+', lapply(counts, '[[', 5))/N 
    lldiff = abs(llnew - ll)
    ll = llnew
    llvect = c(llvect, ll)
    if (verbose) print(sprintf("Training step %s",it))
    if (verbose) print(c(lldiff, ll))
  }
  
  #output
  if (out == '') out = gsub('.csv', '_paramgphmm.json', csvPath)
  if (verbose) print(sprintf('Writting gphmm parameters in : %s', out))
  parameters = toJSON(parameters)
  write(parameters, file = out)
  outll = gsub('.csv', '_llgphmm.json', csvPath)
  if (verbose) print(sprintf('Writting ll in : %s', outll))
  llvect = toJSON(llvect)
  write(llvect, file = outll)
}

# stop timing computation
stop = round(Sys.time() - start)
if (verbose) print(stop)



