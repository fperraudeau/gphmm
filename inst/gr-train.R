#!/usr/bin/env Rscript

suppressMessages({
  require(docopt)
  require(wbretriever)
})

'train

Usage: 
gr-train.R bowtie <fasta> <csv> [options] [--objFx=<obj>] [--cpdir=<cpd>] [--domain=<dmn>] [--k=<k>] [--smin=<sm>] [--np=<np>] [--ma=<ma>] [--mp=<mp>] [--ai=<ai>] [--si=<si>] [--sd=<sd>] [--ad=<ad>] [--L=<L>] [--i=<i>] [--D=<D>] [--R=<R>] [--alpha=<a>] [--tempini=<ti>] [--tempmin=<tm>] [--it=<it>]
gr-train.R gphmm <fasta> <csv> [options] [--llthreshold=<l>] [--maxit=<i>] 

Options:
<fasta>                  path to fasta file with queries.
<csv>                    path to csv file with information on training queries.
--verbose                print things.
--outdir=<o>             name of the output directory file [Default: getwd()].
--ncores=<n>             number of cores [Default: 0].
--dbName=<dbn>           name of reference database, should be identical to one of db names in s3 [Default: 16S-reannotated].

--objFx=<obj>            objective function either refInTruth or accuracy [Default: refInTruth].
--cpdir=<cpd>            name of directory to copy paste results [Default: ].
--domain=<dmn>           json file with domain [Default: ].
--k=<k>                  k [Default: c(1,1,0,1)].
--smin=<sm>              sMin [Default: c(0.9,0.9,0,0.9)].
--ma=<ma>                ma [Default: c(2,2,0,2)].
--mp=<mp>                mp [Default: c(6,6,0,6)].
--si=<si>                si [Default: c(5,5,0,5)].
--ai=<ai>                ai [Default: c(3,3,0,3)].
--sd=<sd>                sd [Default: c(5,5,0,5)].
--ad=<ad>                ad [Default: c(3,3,0,3)].
--i=<i>                  i [Default: c(0.75,0.75,0,0.75)].
--L=<L>                  L [Default: c(20,20,0,20)].
--D=<D>                  D [Default: c(15,15,0,15)].
--R=<R>                  R [Default: c(2,2,0,2)].
--np=<np>                np [Default: c(60,60,0,60)].
--alpha=<a>              alpha [Default: 0.9].
--tempini=<ti>           initial temperature [Default: 1000].
--tempmin=<tm>           minimum temperature [Default: 10].
--it=<it>                it per step [Default: 10]

--llthreshold=<l>        algo stops when difference of log likelihood is below the threshold [Default: 0.00001].
--maxit=<i>              maximum number of iterations [Default: 10].' -> doc

## Read in args.
opt              = docopt(doc)
bowtie           = opt$bowtie
gphmm            = opt$gphmm
fasta            = opt$fasta
csv              = opt$csv
verbose          = opt$verbose
outdir           = opt$outdir
ncores           = as.numeric(opt$ncores)
dbName           = opt$dbName

llthreshold      = as.numeric(opt$llthreshold)
maxit            = as.numeric(opt$maxit)

bowtie           = opt$bowtie
obj              = opt$objFx
cpdir            = opt$cpdir
domain           = opt$domain
k                = as.numeric(eval(parse(text = opt$k)))
smin             = as.numeric(eval(parse(text = opt$smin)))
ma               = as.numeric(eval(parse(text = opt$ma)))
mp               = as.numeric(eval(parse(text = opt$mp)))
si               = as.numeric(eval(parse(text = opt$si)))
ai               = as.numeric(eval(parse(text = opt$ai)))
sd               = as.numeric(eval(parse(text = opt$sd)))
ad               = as.numeric(eval(parse(text = opt$ad)))
L                = as.numeric(eval(parse(text = opt$L)))
i                = as.numeric(eval(parse(text = opt$i)))
R                = as.numeric(eval(parse(text = opt$R)))
D                = as.numeric(eval(parse(text = opt$D)))
np               = as.numeric(eval(parse(text = opt$np)))
alpha            = as.numeric(opt$alpha)
tempini          = as.numeric(opt$tempini)
tempmin          = as.numeric(opt$tempmin)
it               = as.numeric(opt$it)

if (ncores == 0) ncores = max(1, detectCores() - 2)
if (outdir == 'getwd()') outdir = eval(parse(text=outdir))

stopifnot(bowtie | gphmm)
model = ifelse(bowtie, 'bowtie', 'gphmm')
trainFx = makeTrain(model, fasta, csv, verbose, outdir, ncores)

if (verbose) print(paste0('Model to train is: ', model))
if (model == 'bowtie'){
  domain = getOrCreateDomain(domain, cpdir, makeDomain(ma, mp, si, ai, sd, ad, L, i, D, R, k, np, smin))
  stopifnot(sum(domain['min', ] != domain['max',]) != 0)
  trainFx(dbName, domain, obj, tempini, tempmin, alpha, it, cpdir)
} else if (model == 'gphmm'){
  trainFx(llthreshold, maxit)
}


