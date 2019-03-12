args = commandArgs(trailingOnly=TRUE)
library(qvalue)


eFDR = function(pobs,plist,nmax=100,qval=T)
{
  if(is.matrix(plist)) pnullmat = plist else
  if(is.data.frame(plist)) pnullmat = as.matrix(plist) else
  if(is.list(plist)) pnullmat = plist2pmat(pobs,plist)
  if(is.null(names(pobs)))
  {
    warning("Observed p vector, pobs, is not labelled - labelling SNP1..n")
    names(pobs) = paste0("SNP",1:length(pobs))
  }
  fdrdata = eFDR.mat(pobs,pnullmat,nmax,qval)
  return(fdrdata)
}

## turns list of pnull to matrix pnullmat

plist2pmat = function(pobs,pnulllist)
{
  maxnrow = max(unlist(lapply(pnulllist,length)))
  maxnrow = max(maxnrow,pobs)
  nsim = length(pnulllist)
  pnullmat = matrix(NA,maxnrow,nsim)
  for(cc in 1:nsim)
  {
    pvec = pveclist[[cc]]
    pnullmat[1:length(pvec),cc] = pvec
  }
  return(pnullmat)
}

## 

eFDR.mat = function(pobs,pnullmat,nmax,qval)
{
  n_pobs = sum(!is.na(pobs))
  n_null = sum(!is.na(pnullmat))
  print(n_pobs)
  print(n_null)
  lambda = median(pobs,na.rm=T) ## lambda is used to estimate pi0
  pi0 = ( sum(pobs>lambda,na.rm=T) / n_pobs ) /
    ( sum(pnullmat>lambda,na.rm=T) / n_null )
  pi0 = min(pi0,1)
  print(pi0)
  fdrvec = rep(NA,nmax)
  for(rr in 1:nmax)
    {
      gamma = pobs[rr]
      npless.null = sum(pnullmat<=gamma,na.rm=T)
      npless.obs = sum(pobs<=gamma,na.rm=T)
      fdr.cons = (npless.null/n_null) / (npless.obs/n_pobs)
      fdrvec[rr] = fdr.cons
    }

  fdrvec[fdrvec>1]=1
  for(rr in nmax:2)
  {
    if(fdrvec[rr]<fdrvec[rr-1]) fdrvec[rr-1] = fdrvec[rr]
  }

  fdrdata = data.frame(pobs=pobs[1:nmax])
  if(qval) fdrdata$FDR = qvalue(pobs)$qval[1:nmax]
  fdrdata = data.frame(fdrdata,eFDR = fdrvec, eFDR.pi0 = fdrvec * pi0)
  rownames(fdrdata) = names(pobs)[1:nmax]
  return(fdrdata)
}


real_file = args[1]
null_file1 = args[2]
output_file = args[3]


# Load in data
real_data <- read.table(real_file,header=FALSE)
null_data1 <- read.table(null_file1, header=FALSE)


real_pvalues <- real_data$V4
null_pvalues <- null_data1$V4

final <- eFDR(sort(real_pvalues), as.matrix(null_pvalues),nmax=length(real_pvalues))

write.table(final, file = output_file, quote = FALSE, sep = "\t")