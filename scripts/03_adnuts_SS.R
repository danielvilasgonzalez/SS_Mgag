#load libraries  
library(adnuts)
library(snowfall)
library(rstan)
library(shinystan)
library(r4ss)

#cores to be used
reps <- parallel::detectCores()-1 # chains to run in parallel
  
#Reproducible seeds are passed to ADMB
set.seed(352)
seeds <- sample(1:1e4, size=reps)
  
#folder and model example (SS model)
setwd('D:/WFS_DV/manuscripts/stocksynthesis/models/drive/example/')
m <- 'example'
  
#list to store results
Mycomparisonlist<-list()

## First optimize the model to make sure the Hessian is good.
#setwd(m); system('example -nox -iprint 200 -mcmc 15'); setwd('..')
  
#arguments
iter <- 1000 #number of samples
warmup <- iter/4 #warmup iterations
inits <- NULL #start chains from MLE

## First optimize the model to make sure the Hessian is good.
#system('example -nox -iprint 200 -mcmc 15'); setwd('..')

#####################
# initial
#####################
#run the model at its initial values to check how things look like, and plot results
system('example -maxfn 0 -nohess')
file.copy(from="Report.sso",to="Report_initial.sso",overwrite=T)

#save results into a list
mymodel<-SS_output(dir = getwd(),verbose = T,printstats = T,forecast = F,covar = FALSE)
Mycomparisonlist[['initial']]<-mymodel

# plots the results
SS_plots(mymodel,png = TRUE) #plot = c(1:2,4:7,9:26),

#####################
# optimized
#####################
#run model optimizing
system('example -nohess')
#file.copy(from="Report.sso",to="Report_optimized.sso",overwrite=T)
std_table<-read.table(paste0('./ss.std'),header = TRUE)
mymodel<-SS_output(dir = getwd(),verbose = T,printstats = T,forecast = F,covar = TRUE)
SS_plots(mymodel,uncertainty = TRUE,png = TRUE) #plot = c(1:2,4:7,9:26),


#save results into a list
Mycomparisonlist[['optimized']]<-mymodel



#####################
# optimized slope (modified manually)
#####################
#run model optimizing
system('example -nohess')
file.copy(from="Report.sso",to="Report_optimized_slope.sso",overwrite=T)
#save results into a list
mymodel<-SS_output(dir = getwd(),verbose = T,printstats = T,forecast = F,covar = TRUE)
Mycomparisonlist[['optimized fixed slopes']]<-mymodel

#create df
# df<-data.frame(matrix(NA,nrow=0,ncol=5))
# names(df)<-c('index','name','value','std.dev','run')

#run SS model
r4ss::run(dir = getwd(),
          exe = 'ss_win',skipfinished = FALSE,
          show_in_console = TRUE)

SS_output("D:/WFS_DV/manuscripts/stocksynthesis/models/drive/example/")
#r4ss::getADMBHessian()

h<-r4ss::getADMBHessian(hesfile = paste0('./admodel.hes'))$hes#<0

#save std table
std_table<-read.table(paste0('./ss.std'),header = TRUE)
#recruits
df<-std_table[which(std_table$name %in% c('recdev1')),] #,'recr_std','F_std','SSB_std'
df$year<-1985:2020
ggplot()+
  geom_point(data=df,aes(y=value,x=year))+
  geom_line(data=df,aes(y=value,x=year))+
  geom_ribbon(data = df,aes(x=year,ymax=value+std.dev,ymin=value-std.dev),alpha=0.2)
#recruits
df<-std_table[which(grepl('F_rate',std_table$name)),] #,'recr_std','F_std','SSB_std'
df$year<-1985:2020
ggplot()+
  geom_point(data=df,aes(y=value,x=year))+
  geom_line(data=df,aes(y=value,x=year))+
  geom_ribbon(data = df,aes(x=year,ymax=value+std.dev,ymin=value-std.dev),alpha=0.2)
#recruits
df<-std_table[which(grepl('SPR',std_table$name)),] #,'recr_std','F_std','SSB_std'
df$year<-1985:2020
ggplot()+
  geom_point(data=df,aes(y=value,x=year))+
  geom_line(data=df,aes(y=value,x=year))+
  geom_ribbon(data = df,aes(x=year,ymax=value+std.dev,ymin=value-std.dev),alpha=0.2)



Report <- read.table("REPORT.SSO",fill=T,col.names=c(1:200),comment.char="?")
Index <- which (Report[,1]=="Number_of_parameters:"); Npars <- as.numeric(Report[Index,2])
Index <- which (Report[,1]=="Active_count:"); Nests <- as.numeric(Report[Index,2])
Index <- which (Report[,1]=="Number_of_active_parameters_on_or_near_bounds:"); Nbnd <- as.numeric(Report[Index,2])
cat("Npar =",Npars,"; Nest=",Nests,";Nbnd =",Nbnd,"\n")
if (Nbnd >0) { print("Numbers of bounds; stopping")}







SS_plots(mymodel,uncertainty = TRUE,png = TRUE,plot = c(4)) #plot = c(1:2,4:7,9:26),




Mycomparisonlist$example$derived_quants


Mycomparisonsummary<-SSsummarize(biglist = Mycomparisonlist)
SSplotComparisons(summaryoutput=Mycomparisonsummary,uncertainty = TRUE,
                  png=T,plotdir = getwd(),legendlabels = c('initial','optimized'),
                  col = c('#009E73'),indexUncertainty = TRUE)


# read the model output and print diagnostic messages 
replist <- SS_output(dir = getwd(), 
                     verbose = TRUE,
                     printstats = TRUE,covar = TRUE)


SS_plots(replist,plot = c(2),uncertainty = TRUE,png = TRUE) #plot = c(1:2,4:7,9:26),


system('example -nohess')


# plots the results
SS_plots(replist,uncertainty = TRUE)




## First optimize the model to make sure the Hessian is good.
system('example -nox -iprint 200 -mcmc 15'); setwd('..')

## Then run parallel RWM chains as a first test
thin <- 100
iter <- 1000*thin; warmup <- iter/4
inits <- NULL ## start chains from MLE
pilot <- sample_rwm(m, iter=iter, thin=thin, seeds=seeds, init=inits,
                     chains=reps, warmup=warmup, #parallel=TRUE, 
                     path=m, cores=reps,control=list(max_treedepth=5))

#In sample_rwm(m, iter = iter, thin = thin, seeds = seeds, init = inits,  :
#Default init of MLE used for each chain. Consider using dispersed inits3

## Check convergence
mon <- monitor(pilot$samples, warmup=pilot$warmup, print=FALSE)
max(mon[,'Rhat']) #maximum to be <1.01 #potential scale reduction factor
min(mon[,'n_eff']) #effective samples sizes to be >200 
#n_eff indicates the correlation
#%going to zero indicates high autocorrelation
#it should be 0 divergences
## Examine the slowest mixing parameters
slow <- names(sort(mon[,'n_eff']))[1:8]
pairs_admb(fit=pilot, pars=slow)
## Or can specify them by name
pairs_admb(fit=pilot, pars=c('MGparm[1]', 'SR_parm[1]', 'SR_parm[2]'))




#####################################################################
#####################################################################
inits <- sample_inits(pilot, 1)


nuts.mle <-
  sample_nuts(model=m, iter=800, init=inits,  seeds=seeds[1],
              chains=1, warmup=100, path=m, cores=reps)#,
              #control=list(metric="mle", adapt_delta=0.8))

## After regularizing we can run NUTS chains. First reoptimize to get the
## correct mass matrix for NUTS. Note the -hbf 1 argument. This is a
## technical requirement b/c NUTS uses a different set of bounding
## functions and thus the mass matrix will be different.
setwd(m); system(paste(m, '-hbf 1 -nox -iprint 200 -mcmc 15')); setwd('..')
## Use default MLE covariance (mass matrix) and short parallel NUTS chains
## started from the MLE.
nuts.mle <-
  sample_nuts(model=m, iter=500, init=NULL,  seeds=seeds[1],
              parallel=TRUE, chains=1, warmup=100, path=m, cores=reps,
              ) #control=list(metric="mle", adapt_delta=0.8)

## Check for issues like slow mixing, divergences, max treedepths with
## ShinyStan and pairs_admb as above. Fix and rerun this part as needed.
launch_shinyadmb(nuts.mle)

## If good, run again for inference using updated mass matrix. Increase
## adapt_delta toward 1 if you have divergences (runs will take longer).
mass <- nuts.mle$covar.est # note this is in unbounded parameter space
inits <- sample_inits(nuts.mle, reps) ## use inits from pilot run
nuts.updated <-
  sample_nuts(model=m, iter=1000, init=inits,  seeds=seeds,
              chains=reps, warmup=100, path=m, cores=reps,
              mceval=TRUE, control=list(metric=mass, adapt_delta=0.9))
## Again check for issues of nonconvergence and other standard checks. Then
## use for inference.
mon <- monitor(nuts.updated$samples, warmup=nuts.updated$warmup, print=FALSE)
max(mon[,'Rhat'])
min(mon[,'n_eff'])

## We can calculate efficiecy as ess/time. Since there's multiple chains
## add the time together because the ESS is summed across chains too.
(eff <- min(mon[,'n_eff'])/sum(nuts.updated$time.total))
## Or how long to get 1000 effective samples
1000/eff                                # in seconds
1000/eff/60                             # in minutes

## NOTE: the mceval=TRUE argument tells ADMB to run -mceval on ALL chains
## combined AFTER discarding warmup period and thinning. Thus whatever your
## model outputs during mceval is ready for use in
## management. Alternatively you can run -mceval from the command
## line. sample_admb will merge samples into the .psv file in the main
## folder so either way works.

#if low n_eff and >Rhat, check the pairs, and
#if the problem is related to high correlation, then try this tune
cov0<-adnuts:::.get.admb.cov(d)$cov.unbounded
cov0<-0.5*(cov0+t(cov0)) #force symmetric
## the last 3 parameters are the new sigmas in this case
n0<-nrow(cov0)
cov.new<-diag(x=.2^2,nrow=n0) #guess for unbounded variances!!
cov.bew[(1:(n0-3)),(1:(n0-3))]<-cov0[(1:(n0-3)),(1:(n0-3))]
#fit w/ dense adaptation
fit_est_sigmas<-sample_nuts(m,d,iter=2000,chains=5,cores=5,seeds=1:5,warmup=1000,
                            control=list(metric=cov.new,adapt_delta=.9,adapt_mass_dense=TRUE,adapt_window=500))

#diagnostics
print(fit)
summary(fit$monitor$n_eff)
summary(fit$monitor$Rhat)
post <- extract_samples(fit)
str(post[,1:5])

fit$cmd[1] # this is the command line arguments used
## Merged chains after discarding warmup phase
post <- extract_samples(fit)
str(post)
## A list with MLE fit
str(fit$mle)

#start from MLE (RWM gave error-NA)
pilot <- sample_nuts(m, iter=800, seeds=seeds, init=inits,
                     chains=reps, warmup=200, #parallel=TRUE, 
                     path=m, cores=reps)

inits <- sample_inits(pilot, reps) 