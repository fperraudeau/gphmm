require(gphmm)
require(Biostrings)


generateReads <- function(seq = 'ATGCGGATCG', qv = NULL, seed = NULL, 
                          paramgphmm = initializeGphmm()){
  #seed
  if (!is.null(seed)) set.seed(seed)
  
  # param gphmm
  paramgphmm[['deltaX']] = 1/(1+exp(-sum(paramgphmm$deltaX * c(1, qv))))
  paramgphmm[['deltaY']] = 1/(1+exp(-sum(paramgphmm$deltaY * c(1, qv))))
  
  states = c("A", "C", "G", "T")
  colnames(paramgphmm$pp) = rownames(paramgphmm$pp) = names(paramgphmm$qX) = names(paramgphmm$qY) = states
  transition_mat = matrix(c(1 - paramgphmm$deltaX - paramgphmm$deltaY - paramgphmm$tau,
                            1 - paramgphmm$deltaX - paramgphmm$deltaY - paramgphmm$tau,
                            1 - paramgphmm$epsX - paramgphmm$tau,
                            1 - paramgphmm$epsY - paramgphmm$tau,
                            paramgphmm$deltaX, paramgphmm$deltaY,
                            paramgphmm$epsX, 0, 
                            paramgphmm$deltaX, paramgphmm$deltaY,
                            0, paramgphmm$epsY, 
                            rep(paramgphmm$tau, 4)), ncol = 4)
  colnames(transition_mat) = c('M', 'X', 'Y', 'T')
  rownames(transition_mat) = c('I', 'M', 'X', 'Y')
  
  # seq
  if (is.null(qv)) qv = 20
  seq = as.character(seq)
  seq = strsplit(seq, '')[[1]]
  
  # initialization
  state = 'I'
  state = sample(colnames(transition_mat), 1, prob = transition_mat[state, ])
  states = c()
  read = c()
  i = 1
  
  # simulate
  while(state != 'T' & i < (length(seq) + 1)){
    if (state == 'M'){
      read[i] = sample(colnames(paramgphmm$pp), 1, prob = paramgphmm$pp[, seq[i]])
    } else if (state == 'X'){
      read[i] = sample(names(paramgphmm$qX), 1, prob = paramgphmm$qX)
    } else{
      read[i] = ''
    }
    states = c(states, state)
    i = i + 1
    state = sample(colnames(transition_mat), 1, prob = transition_mat[state, ])
  }
  
  list(read = paste(read, collapse = ''), ref = paste(seq, collapse = ''),
       path = paste(states, collapse = ''), tabstates = table(states), qv = qv)
}


paramgphmm = initializeGphmm()
paramgphmm
paramgphmm$deltaX = c(-.1, -.1)
paramgphmm$deltaY = c(-1, -1)
qv = seq(1, 40, length.out = 100)
plot(qv, sapply(qv, function(x) 1/(1+exp(-sum(c(-.01, -.1) * c(1, x))))))

n = 10
seqs = generateRandomSequences(n = n, meanLen = 100, sdLen = 10, prob = paramgphmm$qR)
plot(density(width(seqs)))

reads = lapply(as.character(seqs), function(x){
  generateReads(seq = x, paramgphmm = paramgphmm, qv = 20, seed = 8263)
})


train = c(seqs, DNAStringSet(sapply(reads, '[[', 1)))
names(train) = c(names(train)[1:n], gsub('s', 't', names(train)[1:n]))
writeXStringSet(train, 'train.fasta')
write.table(data.frame(reads = paste0('t', 1:n),
                       ref = paste0('s', 1:n),
                       qv = rnorm(n, 20, 1)), 'train.csv')

system('inst/gphmm.R train train.fasta train.csv --verbose')
estimator = fromJSON('train_paramgphmm.json')
estimator









