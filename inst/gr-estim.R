#!/usr/bin/env Rscript

suppressMessages({
  require(docopt)
  require(wbretriever)
  require(Biostrings)
})

'Estimate.

Usage: 
gr-estim bowtie <fasta> [options]
gr-estim gphmm <fasta> [--paramgphmm=<pgphmm>] [options]

Options:
<fasta>                  path to fasta file with queries.
--verbose                if TRUE, printttt things.
--bam                    if TRUE, only optim, does not run bowtie2.
--dirIdx=<o>             name of the directory to save index of db [Default: getwd()].
--ncores=<n>             integer, number of cores [Default: 0].
--dbName=<dbn>           name of reference database, should be identical to one of db names in s3 [Default: 16S-reannotated].
--parambt=<pbt>          path to json file with Bowtie parameters, if not specified, last added in inst/training [Default: ].
--paramgphmm=<pgphmm>    path to json file with GPHMM parameters, if not specified, last added in inst/training [Default: ].
--B=<B>                  integer, number of iterations in Monte Carlo cross-validation [Default: 10].
--validset=<vs>          float, proportion of reads in validation set for Monte Carlo cross-validation [Default: 0.8].
--seed=<seed>            integer, seed.' -> doc

## Read in args.
opt              = docopt(doc)
bowtie           = opt$bowtie
gphmm            = opt$gphmm
fasta            = opt$fasta
verbose          = opt$verbose
bam              = opt$bam
dirIdx           = opt$dirIdx
ncores           = as.numeric(opt$ncores)
dbName           = opt$dbName
parambt          = opt$parambt
paramgphmm       = opt$paramgphmm
B                = as.numeric(opt$B)
validset         = as.numeric(opt$validset)
seed             = opt$seed

if (ncores == 0) ncores = max(1, detectCores() - 2)
if (dirIdx == 'getwd()') dirIdx = eval(parse(text=dirIdx))
out = gsub('.fasta', '', fasta)
rName = strsplit(gsub('.fasta', '', dbName), '/')[[1]]
rName = rName[length(rName)]
out = paste(out, rName, sep = '_')

stopifnot(bowtie | gphmm)
model = ifelse(bowtie, 'bowtie', 'gphmm')
paramgphmm = ifelse(bowtie, '', paramgphmm)
if (model == 'gphmm' & paramgphmm == '') paramgphmm = findLastJson('paramgphmm')
if (parambt == '') parambt = findLastJson('parambt')

####################
# Bowtie + GPHMM
####################
if (!bam){
  print('Start alignment with Bowtie2...')
  startGPHMM = Sys.time()
  if (verbose) print(paste0('Model to compute normalized probabilities is: ', model))
  assignFx = makeAssign(fasta, parambt, dbName, dirIdx, ncores)
  bamdf = assignFx(model, paramgphmm)
  write.table(bamdf, sprintf('%s_%s_bam.csv', out, model), row.names = T, col.names = T)
  stopGPHMM = Sys.time() - startGPHMM
} else{
  bamdf = read.table(sprintf('%s_%s_bam.csv', out, model), stringsAsFactors = F)
}

####################
# Optimization
####################
if (dbName %in% c('16S-wballeles', '16S-species', '16S-reannotated')){
  xdb = fetchDb(dbName)
  annoAll = dbAnno(xdb)
  annoAll$longname = paste(annoAll$genus, annoAll$species, sep='_')
  anno = annoAll[,c('refGid', 'longname', 'genus', 'family')]
  rm(annoAll)
}else{
  dna = readDNAStringSet(dbName)
  anno = data.frame(refGid = names(dna), longname = names(dna),
                    genus = sapply(strsplit(names(dna), '_'), '[[', 1),
                    family = '')
}


#download_rrnDB()

Pout = sprintf('%s_P.RData', out)
mccvRidgePath = sprintf('%s_mccvRidge.RData', out)
mccvLassoPath = sprintf('%s_mccvLasso.RData', out)
pihatout = sprintf('%s_pihat.csv', out)
P = computeP(bamdf)
rm(bamdf)

# naive
if (verbose) print('Naive estimator')
naive = Matrix::colSums(P) / sum(Matrix::colSums(P))

# mle
if (verbose) print('MLE')
startMLE = Sys.time()
mleOptim = optimPi(P, seed = seed)
stopifnot(mleOptim$convergence == 0)
mle = mleOptim$pars
stopMLE = Sys.time() - startMLE
rm(mleOptim)

# ridge
if (verbose) print('Ridge estimator')
startRidge = Sys.time()
mccvRidge = processMCCV(P, B = B, validset = validset, piIni = mle,
                        ncores = ncores, seed = seed)
save(mccvRidge, file =  mccvRidgePath)
load(mccvRidgePath)
ridgeOptim = optimPi(P, lambda =  mccvRidge$lambdahat, piIni = mle)
ridge = ridgeOptim$pars
stopRidge = Sys.time() - startRidge
rm(mccvRidge)

# grouplasso, species
print('Group Lasso estimator, species')
startLasso = Sys.time()
mccvLasso = processMCCV(P, B = B, validset = validset, piIni = mle, seed =seed,
                        ncores = ncores, ridge = F, anno = anno)
save(mccvLasso, file = mccvLassoPath)
lassoOptim = optimPi(P, lambda = mccvLasso$lambdahat, anno=anno,ridge=F,
                     piIni= mle)
lasso = lassoOptim$pars
stopLasso = Sys.time() - startLasso
rm(mccvLasso)

pihat = data.frame(Naive = naive, MLE = mle, Ridge = ridge, LASSO = lasso)
rownames(pihat) = colnames(P)
write.csv(pihat, pihatout)
rm(pihat)

if (verbose){
  if (!bam) print(stopGPHMM)
  print(stopMLE)
  print(stopRidge)
  print(stopLasso)
}


stopifnot(lassoOptim$convergence == 0)
stopifnot(ridgeOptim$convergence == 0)


