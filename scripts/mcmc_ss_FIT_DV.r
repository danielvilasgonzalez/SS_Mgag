### Example script to run MCMC using R package adnuts for a stock
### synthesis assessment model. Started Feb 2021 by Cole Monnahan
### (cole.monnahan@noaa.gov | AFSC)

## Note I am using SS version 3.30.16, which is freely available
## from: https://vlab.ncep.noaa.gov/web/stock-synthesis/home
## Non-windows users will need to download and put that
## executable in the hake folder
library(adnuts)
packageVersion('adnuts')        # 1.1.2
library(shinystan)

#library(adnuts)
library(snowfall)
library(rstan)
#library(shinystan)
### ------------------------------------------------------------
### Task 0: Set up and test model for running. This requires
### pointing to a folder and executable. The folder needs to
### contain all sufficient input files and assumes optimization
### has occurred and produced all necessary outputs. Temporary
### copies will be made in the working directory during execution

## Define the path and model name (without .exe extension)
# m <- 'ss'                               # model name
# p <- 'Run2'                             # path to folder
m <- 'example'                               # model name
p <- 'adnuts'                             # path to folder
## Assumes current working directory is where this R script is
#(wd <- getwd())
setwd('D:/WFS_DV/manuscripts/stocksynthesis/models/example/')
### Optimize model from R. Want to optimize w/ -mcmc flag b/c of
### bias adjustment.
#setwd(p);

#file.copy(from="STARTER1.SS", to="STARTER.SS",overwrite=T)
file.copy(from="starter.SS", to="starter1.ss",overwrite=T)
system("ss -nox");
file.copy(from="report.SSO",to="report0.SSO",overwrite=T)

Report <- read.table("report.SSO",fill=T,col.names=c(1:200),comment.char="?")
Index <- which (Report[,1]=="Number_of_parameters:"); Npars <- as.numeric(Report[Index,2])
Index <- which (Report[,1]=="Active_count:"); Nests <- as.numeric(Report[Index,2])
Index <- which (Report[,1]=="Number_of_active_parameters_on_or_near_bounds:"); Nbnd <- as.numeric(Report[Index,2])
cat("Npar =",Npars,"; Nest=",Nests,";Nbnd =",Nbnd,"\n")
if (Nbnd >0) { print("Numbers of bounds; stopping"); AAA}

#AAA

file.copy(from="STARTER2.SS", to="STARTER.SS",overwrite=T)
system("ss -hbf -nox -mcmc 10"); 
file.copy(from="REPORT.SSO",to="REPORT.SSO",overwrite=T)
Report <- read.table("REPORT.SSO",fill=T,col.names=c(1:200),comment.char="?")
Index <- which (Report[,1]=="Number_of_parameters:"); Npars <- as.numeric(Report[Index,2])
Index <- which (Report[,1]=="Active_count:"); Nests <- as.numeric(Report[Index,2])
Index <- which (Report[,1]=="Number_of_active_parameters_on_or_near_bounds:"); Nbnd <- as.numeric(Report[Index,2])
cat("Npar =",Npars,"; Nest=",Nests,";Nbnd =",Nbnd,"\n")
if (Nbnd >0) { print("Numbers of bounds; stopping")}

#setwd('..')

#AAAA

# Run specs
print(Sys.time())
chains <- 5
thin <- 1
iter <- 1000*thin
# fit <- sample_nuts(model=m, path=p, iter=iter, warmup=iter*0.5,
#                    chains=chains, thin=thin,
#                    control=list(max_treedepth=15,metric='mle'))
fit <- sample_nuts(model=m, path=, iter=iter, warmup=iter*0.5,
                   chains=chains, thin=thin,
                   control=list(max_treedepth=5)) #,metric='mle'
#adnuts::
print(fit$cmd)
print(mean(fit$time.total)/60)
print(Sys.time())

## Check convergence
mon <- monitor(fit$samples, warmup=fit$warmup, print=FALSE)
max(mon[,'Rhat']) #maximum to be <1.01
min(mon[,'n_eff']) #effective samples sizes to be >200
#%going to zero indicates high autocorrelation
#it should be 0 divergences
## Examine the slowest mixing parameters
slow <- names(sort(mon[,'n_eff']))[1:8]
## We can calculate efficiecy as ess/time. Since there's multiple chains
## add the time together because the ESS is summed across chains too.
(eff <- min(mon[,'n_eff'])/sum(fit$time.total))
## Or how long to get 1000 effective samples
1000/eff                                # in seconds
1000/eff/60                             # in minutes

## Good idea to save the output, I recommend RDS format.
cat("saving",paste0('fits/',p,'.RDS'),"\n")
saveRDS(fit, file=paste0('fits/',p,'.RDS'))
#fit <- readRDS(file=paste0('fits/',p,'.RDS'))

## Key information from run. Including the two recommended
## convergence diagnostics:
summary(fit)

## Interactive tools (must close out browser to regain console)
launch_shinyadmb(fit)

