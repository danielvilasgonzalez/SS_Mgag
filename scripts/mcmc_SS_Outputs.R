### Example script to run MCMC using R package adnuts for a stock
### syntehsis assessment model. Started Feb 2021 by Cole Monnahan
### (cole.monnahan@noaa.gov | AFSC)

## Note I am using SS version 3.30.16, which is freely available
## from: https://vlab.ncep.noaa.gov/web/stock-synthesis/home
## Non-windows users will need to download and put that
## executable in the hake folder
library(adnuts)
library(rstan)
packageVersion('adnuts')        # 1.1.2
library(shinystan)

# read the output files
fit <- readRDS(file=paste0('fits\\','run1.RDS'))
chain <- length(fit$time.total)
print(fit$cmd)

## Key information from run. Including the two recommended
## convergence diagnostics:
summary(fit)
## or from the rstan::monitor output
#str(fit$monitor)



## Interactive tools (must close out browser to regain console)
launch_shinyadmb(fit)
#plot_sampler_params(fit)                # NUTS adaptation

## Extract posterior samples as a data.frame
post <- extract_samples(fit)
#str(post)


## The key parameters
KeyPars <- which(fit$par_names=="lp__")
for (Ipar in 1:length(fit$par_names)) if (substr(fit$par_names[Ipar],1,6)=="MGparm") KeyPars <- c(KeyPars,Ipar)
for (Ipar in 1:length(fit$par_names)) if (substr(fit$par_names[Ipar],1,7)=="SR_parm") KeyPars <- c(KeyPars,Ipar)
print(fit$par_names[KeyPars])
pairs_admb(fit, pars=KeyPars)

## The key parameters
KeyPars <- which(fit$par_names=="lp__")
for (Ipar in 1:length(fit$par_names)) if (substr(fit$par_names[Ipar],1,6)=="init_F") KeyPars <- c(KeyPars,Ipar)
for (Ipar in 1:length(fit$par_names)) if (substr(fit$par_names[Ipar],1,6)=="Q_parm") KeyPars <- c(KeyPars,Ipar)
for (Ipar in 1:length(fit$par_names)) if (substr(fit$par_names[Ipar],1,7)=="selparm") KeyPars <- c(KeyPars,Ipar)
#print(fit$par_names[KeyPars])
pairs_admb(fit, pars=KeyPars)


## The most correlated parameters
Npars <- length(fit$par_names)-1
Parmat <- matrix(0,nrow=length(post[[1]]),ncol=Npars)
for (Ipar in 1:Npars) Parmat[,Ipar] <- post[[Ipar]]
covar <- cor(Parmat)
KeyPars <- NULL
for (II in 2:Npars)
 for (JJ in 1:(II-1))  
  if (abs(covar[II,JJ]) > 0.9) { print(covar[II,JJ]); KeyPars <- c(KeyPars,II,JJ); }
KeyPars <- sort(unique(KeyPars))
print(fit$par_names[KeyPars])
#if (length(KeyPars)>0) pairs_admb(fit, pars=KeyPars)

## parameters with low CVs based on MLE
KeyPars <- which(abs(fit$mle$se[1:Npars]/fit$mle$est[1:Npars])<0.1)
print(fit$par_names[KeyPars])
pairs_admb(fit, pars=KeyPars)

## The 6 slowest/fastest mixing parameters
pairs_admb(fit, pars=1:12)
pairs_admb(fit, pars=1:8, order='slow',diag='trace')
pairs_admb(fit, pars=1:8, order='slow', diag='hist')
pairs_admb(fit, pars=1:6, order='fast')
## Can also specify names or use grep
#pairs_admb(fit, pars=c('recdev_early[21]','recdev_early[22]',
#                       'recdev_early[23]'))
pairs_admb(fit, pars=grep('_parm', fit$par_names))

## Marginal MLE vs posterior
plot_marginals(fit, pars=1:16)
x <- plot_uncertainties(fit)
plot_marginals(fit, pars=which.max(x$sd.post))

## Marginal comparisons as multipage PDF for easy scrolling
pdf('marginals.pdf', onefile=TRUE, width=7,height=5)
plot_marginals(fit)
dev.off()


## Get sampler characteristics
for (Ichain in 1:chain)
  print(summary(fit$sampler_params[[Ichain]]))

### ------------------------------------------------------------
### Extracting posterior samples
## Get post-warmup, merged chains into data frame (these two are
## identical)
str(as.data.frame(fit))
str(extract_samples(fit))
## If you want it in list form, e.g., to put into coda package
str(extract_samples(fit, as.list=TRUE))
## If you want to see warmup samples, and the log-posterior (lp__ column)
str(extract_samples(fit, inc_warmup=TRUE, inc_lp=TRUE))
