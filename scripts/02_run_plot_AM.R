####################################################################
####################################################################
##
##    RUN ECOSPACE CONSOLE AND INVESTIGATE RESULTS FROM IBM
##    Daniel Vilas (danielvilasgonzalez@gmail.com)
##
####################################################################
####################################################################

#set library
#.libPaths( c( .libPaths(), "C:/R/win-library/") )

#load packages
# library(ggplot2)
# library(ggthemes)
# library(scales)
#library(colorRamps)
#install.packages("r4ss")
library(r4ss)
#library(reshape2)

#set working directory
wd<-'/Users/daniel/Work/SS time-varying M/'
setwd(wd)

##########################################################
# RUN MODEL
##########################################################

#vector of models
models<-list.dirs(path = './models/',full.names = FALSE,recursive = FALSE)
models<-models[!grepl('old',models)]

#create list models
Mycomparisonlist<-list()

#not working for 4, 8 
for (m in models) {
  
  #m<-models[1]
  
  #print scenario to check progress
  cat(paste(" #############   Model", m, match(m,models), 'out of',length(models),  "  #############\n"))
  
  #get ss exe file in the model folder
  r4ss::get_ss3_exe(paste0(getwd(),'/models/',m))
  
  #run SS model
  r4ss::run(dir = paste0(getwd(),'/models/',m),skipfinished = FALSE,show_in_console=TRUE)
  
  #check Hessian matrix
  r4ss::getADMBHessian(hesfile = paste0('./models/',m,'/admodel.hes'))$hes<0 #negative------ non-convergence
  
  #save results into a list
  mymodel<-SS_output(dir = paste0('./models/',m),verbose = T,printstats = T,forecast = F,covar = TRUE)
  Mycomparisonlist[[m]]<-mymodel
  
  # read the model output and print diagnostic messages 
  replist <- SS_output(dir = paste0('./models/',m), 
                       verbose = TRUE,
                       printstats = TRUE,covar = TRUE)
  
  # plots the results
  SS_plots(replist,plot = c(1:2,4:7,9:26),uncertainty = TRUE)
}

#check models in the list
names(Mycomparisonlist)

#save list
saveRDS(Mycomparisonlist, file="./models/summary_list.RData")
#readRDS(file = "./models/summary_list.RData")

#summary list
Mycomparisonsummary<-SSsummarize(biglist = Mycomparisonlist)

#####################
# maximum final parameter gradient
#####################




#####################
# RMSE - relative error for survey indices and age comp
#####################




#####################
# Likelihood components
#####################



#####################
# Retrospective (bias) analysis
#####################



#####################
# Hindcast cross-validation
#####################




#####################
# Hindcast cross-validation
#####################



#####################
# correlation MRT and recruitment and SSB
#####################



#####################
# timeseries comparison of recruitment, M, F, SSB and SPR relative to OBS values
#####################

#color palette
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
names<-names(palettes)
pal <- ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value

summaryoutput<-replist
sum<-r4ss::SSsummarize(Mycomparisonlist)
sum$likelihoods_by_fleet
sum$likelihoods
sum(sum$likelihoods[1,])



# run the retrospective analyses
retro(
  dir = new_mod_path, # wherever the model files are
  oldsubdir = "", # subfolder within dir
  newsubdir = "retrospectives", # new place to store retro runs within dir
  years = 0:-5, # years relative to ending year of model
  exe = "ss3"
)

# load the 6 models
retroModels <- SSgetoutput(dirvec = file.path(
  new_mod_path, "retrospectives",
  paste("retro", 0:-5, sep = "")
))
# summarize the model results
retroSummary <- SSsummarize(retroModels)
# create a vector of the ending year of the retrospectives
endyrvec <- retroSummary[["endyrs"]] + 0:-5
# make plots comparing the 6 models
# showing 2 out of the 17 plots done by SSplotComparisons
SSplotComparisons(retroSummary,
                  endyrvec = endyrvec,
                  legendlabels = paste("Data", 0:-5, "years"),
                  subplot = 2, # only show one plot in vignette
                  print = TRUE, # send plots to PNG file
                  plot = FALSE, # don't plot to default graphics device
                  plotdir = new_mod_path
)

# calculate Mohn's rho, a diagnostic value
rho_output <- SSmohnsrho(
  summaryoutput = retroSummary,
  endyrvec = endyrvec,
  startyr = retroSummary[["endyrs"]] - 5,
  verbose = FALSE
)

read.table(paste0('./models/',m,'/Report.sso'),header = TRUE)









save(Mycomparisonlist,file=paste0(dr,':/WFS_DV/manuscripts/stocksynthesis/outputs/comparison_list.RData'))

##########################################################
# COMPARISON PLOT
##########################################################

    Mycomparisonsummary<-SSsummarize(biglist = Mycomparisonlist)

#color palette
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
names<-names(palettes)
pal <- ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value

summaryoutput<-replist
sum<-r4ss::SSsummarize(Mycomparisonlist)
sum$likelihoods_by_fleet
sum$likelihoods
sum(sum$likelihoods[1,])


# get stuff from summary output
n <- summaryoutput[["n"]]
likelihoods <- summaryoutput[["likelihoods"]]
if (is.null(likelihoods)) {
  stop(
    "Input 'summaryoutput' needs to be a list output from SSsummarize\n",
    "and have an element named 'likelihoods'."
  )
}
pars <- summaryoutput[["pars"]]
par_prior_likes <- summaryoutput[["par_prior_likes"]]

r4ss::SS_output()


SSplotComparisons(summaryoutput=Mycomparisonsummary,
                  png=T,uncertainty = TRUE,show_equilibrium = FALSE,
                  plotdir = paste0(dr,':/WFS_DV/manuscripts/stocksynthesis/outputs/model_comparison/'),
                  legendlabels = c(models),indexUncertainty = TRUE,
                  col = c(pal[1:length(Mycomparisonlist)])) #,subplots = c(1:2,4:7,9:26) #1,5,7,9,11,13,14 #'#009E73','#F0E442','#0072B2'



##########################################################
# INDEX
##########################################################

SSplotComparisons(summaryoutput=Mycomparisonsummary, png=T,plotdir=paste0(dr,':/WFS_DV/manuscripts/stocksynthesis/outputs/model_comparison/'))

for (m in models) {
  m<-models[6]
  
  #save std table
  std_table<-read.table(paste0('./models/',m,'/ss.std'),header = TRUE)

#std_table[,which('index'==)]
unique(std_table$name)
df<-std_table[which(std_table$name %in% c('SSB_std')),] #,'recr_std','F_std','SSB_std'
df$idx<-1:nrow(df)
df$std.dev<-df$value
df$value<-Mycomparisonlist[[m]]$sprseries$SSB
df$value+df$std.dev

#recruits
df<-std_table[which(std_table$name %in% c('recr_std')),] #,'recr_std','F_std','SSB_std'
df$year<-1985:2020
Mycomparisonlist[[m]]$recruit

#F_rate
df<-std_table[which(grepl('F_rate',std_table$name)),] #,'recr_std','F_std','SSB_std'
df$year<-1985:2020

ggplot()+
  geom_line(data=df,aes(x=year,y=value))+
  geom_ribbon(data=df,aes(x=year,ymax=value+std.dev,ymin=value-std.dev),alpha=0.5)


grepl('spr',names(Mycomparisonlist[[m]]))

read.table(paste0('./models/',m,'/Report.sso'),header = TRUE)
  
SStableComparisons(Mycomparisonsummary)

Mycomparisonsummary

r4ss::SSplotIndices(replist)

##############################################
##############################################
##############################################


Mycomparisonlist$AM_blk.a$derived_quants
replist$derived_quants$StdDev






covarfile<-'E:/WFS_DV/manuscripts/stocksynthesis/models/AM_blk.a/covar.sso'
covarskip <- grep("active-i", covarhead) - 1
CoVar <- read.table(covarfile, header = TRUE, colClasses = c(rep("numeric", 4), rep("character", 4), "numeric"), skip = covarskip)
covarhead <- readLines(con = covarfile, n = 10)

covartime <- findtime(covarhead)

stdtable <- CoVar[CoVar[["Par..j"]] == "Std", c(7, 9, 5)]
names(stdtable) <- c("name", "std", "type")
N_estimated_parameters2 <- sum(stdtable[["type"]] == "Par")

###########################################################
###########################################################
###########################################################

replist <- SS_output(dir = 'E:/software courses/StockSynthesis/2020 SS Class/In Class Example', 
                     verbose = TRUE,
                     printstats = TRUE,covar = TRUE)

'E:/software courses/StockSynthesis/2020 SS Class/In Class Example'





















SSplotPars(replist,showpost = F)

SSplotIndices(replist)


mymodel$Kobe$F.Fmsy
mymodel$sprseries
mymodel$timeseries

mymodel$sprseries$SPR_report #B/Bmsy
mymodel$Kobe$F.Fmsy
mymodel$maximum_gradient_component

Mycomparisonsummary$SPRratio

Mycomparisonsummary$Bratio
#rm(.SavedPlots,pos=1)

SSplotComparisons(summaryoutput=Mycomparisonsummary,
                  png=T,
                  plotdir = mydir,
                  #legendlabels = c('fix_scn','env_scn','blk_scn'),
                  col = c('#009E73','#F0E442','#0072B2'),subplots = c(1)) #1,5,7,9,11,13,14


#######################################################
#
# RESIDUAL INDEX
#
########################################################

head(Mycomparisonsummary$indices)

ind<-Mycomparisonsummary$indices
ind$Yr<-as.numeric(ind$Yr)

df<-data.frame(matrix(NA,nrow = 6,ncol = 3))
colnames(df)<-c('scenario','index','RMSE')
df$scenario<-c(unique(ind$name),unique(ind$name))
df$index<-rep( unique(ind$Fleet_name),each=3)

for (i in unique(ind$name)) {
  for (j in unique(ind$Fleet_name)) {
    
    #j='SURVEYADULT'
    x<-subset(ind,name==i & Fleet_name==j)
    y<-sqrt(mean((x[,'Exp'] - x[,'Obs'])^2))
    y1<-sqrt(sum((x[,'Obs'] - x[,'Exp'])^2)/length(x[,'Exp']))
    print(paste0(i," ",j," - RMSE ",y))
    print(paste0(i," ",j," - RMSE ",y1))
    #print(rmse<-sqrt(mean((ind$Exp - ind$Obs)^2)))
    df[which(df$scenario==i & df$index==j),'RMSE']<-y
    
    
    
  }
}

df$index<-rep( c('juvenile survey','adult survey'),each=3)
colnames(df)[2]<-'survey'
df$RMSE<-round(df$RMSE,digits = 2)
ind$label<-paste0(ind$name,' - ',ind$Fleet_name)


ind1<-subset(ind,Fleet_name=='SURVEYJUVENILE')
ind2<-subset(ind,Fleet_name=='SURVEYADULT')

col_grid <- rgb(235, 235, 235, 100, maxColorValue = 500)

p1<-ggplot()+
  geom_linerange(data=ind1,aes(x=Yr,y=Dev,color=label,ymin=0,ymax=Dev),position = position_dodge(.5),show.legend = F) +
  #annotate("text", x=2015, y=0.7, label= paste0("RMSE = ",round(rmse,digits = 2))) + 
  geom_hline(yintercept=0,linetype='dashed')+
  #stat_smooth(data=ind,aes(x=Yr,y=Dev),method = 'loess',se=FALSE, colour="black",fill=NA)+
  geom_boxplot(data=ind1,aes(x=Yr,y=Dev,group=Yr))+
  geom_point(data=ind1,aes(x=Yr,y=Dev,fill=label),position = position_dodge(.5),color='black',size=1.5,alpha=0.8,shape=21)+
  
  #scale_fill_tableau()+
  #scale_color_manual(values=c('fix_scn'='#009E73','env_scn'='#F0E442','blk_scn'='#0072B2'))+
  scale_fill_manual(values=c('fix_scn - SURVEYJUVENILE'='#009E73',
                             'env_scn - SURVEYJUVENILE'='#F0E442',
                             'blk_scn - SURVEYJUVENILE'='#0072B2'),
                    label=c('fix_SS',
                            'env_SS',
                            'blk_SS'))+
  scale_color_manual(values=c('fix_scn - SURVEYJUVENILE'='#009E73',
                              'env_scn - SURVEYJUVENILE'='#F0E442',
                              'blk_scn - SURVEYJUVENILE'='#0072B2'),
                     label=c('fix_SS',
                             'env_SS',
                             'blk_SS'))+
  ggtitle('juvenile index')+
  theme_bw()+
  ylab('')+
  xlab('')+
  guides(fill = guide_legend(override.aes = list(size = 3)))+
  theme(aspect.ratio = 0.33,legend.position = 'none',strip.text.x = element_text(size = 9,color='black'), 
        strip.background = element_rect(color="black", fill="white", size=1.2, linetype="solid"),plot.title = element_text(hjust = 0.5),
        #plot.title = element_text(hjust = 0.5,size = 14),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.minor = element_line(color = col_grid,
                                        size = 0.3,
                                        linetype = 'dashed'),
        panel.grid.major = element_line(color = col_grid,
                                        size = 0.5))+
  scale_y_continuous(limits=c(-1.75,+1.75),breaks = c(-1,0,1))+
  
  #coord_cartesian(ylim=c(-1,+1)) +
  #geom_smooth(method="lm")+
  scale_x_continuous(breaks = c(1985,1990,1995,2000,2005,2010,2015,2020),labels = c(1985,1990,1995,2000,2005,2010,2015,2020),
                     minor_breaks = c(1986:1989,1991:1994,1996:1999,2001:2004,2006:2009,2011:2014,2015:2019))+
  
  expand_limits(y=0)


p2<-ggplot()+
  geom_linerange(data=ind2,aes(x=Yr,y=Dev,color=label,ymin=0,ymax=Dev),position = position_dodge(.5),show.legend = F) +
  #annotate("text", x=2015, y=0.7, label= paste0("RMSE = ",round(rmse,digits = 2))) + 
  geom_hline(yintercept=0,linetype='dashed')+
  #stat_smooth(data=ind,aes(x=Yr,y=Dev),method = 'loess',se=FALSE, colour="black",fill=NA)+
  geom_boxplot(data=ind2,aes(x=Yr,y=Dev,group=Yr))+
  geom_point(data=ind2,aes(x=Yr,y=Dev,fill=label),position = position_dodge(.5),color='black',size=1.5,alpha=0.8,shape=21)+
  
  #scale_fill_tableau()+
  #scale_color_manual(values=c('fix_scn'='#009E73','env_scn'='#F0E442','blk_scn'='#0072B2'))+
  scale_fill_manual(values=c('fix_scn - SURVEYADULT'='#009E73',
                             'env_scn - SURVEYADULT'='#F0E442',
                             'blk_scn - SURVEYADULT'='#0072B2'),
                    label=c('fix_SS',
                            'env_SS',
                            'blk_SS'))+
  scale_color_manual(values=c('fix_scn - SURVEYADULT'='#009E73',
                              'env_scn - SURVEYADULT'='#F0E442',
                              'blk_scn - SURVEYADULT'='#0072B2'),
                     label=c('fix_SS',
                             'env_SS',
                             'blk_SS'))+
  ggtitle('adult index')+
  theme_bw()+
  ylab('deviance residuals')+
  xlab('')+
  guides(fill = guide_legend(override.aes = list(size = 3)))+
  theme(aspect.ratio = 0.33,legend.position = 'none',strip.text.x = element_text(size = 9,color='black'), 
        strip.background = element_rect(color="black", fill="white", size=1.2, linetype="solid"),plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        #plot.title = element_text(hjust = 0.5,size = 14),
        panel.grid.minor = element_line(color = col_grid,
                                        size = 0.3,
                                        linetype = 'dashed'),
        panel.grid.major = element_line(color = col_grid,
                                        size = 0.5))+
  scale_y_continuous(limits=c(-1.75,+1.75),breaks = c(-1,0,1))+
  
  #coord_cartesian(ylim=c(-1,+1)) +
  #geom_smooth(method="lm")+
  scale_x_continuous(breaks = c(1985,1990,1995,2000,2005,2010,2015,2020),labels = c(1985,1990,1995,2000,2005,2010,2015,2020),
                     minor_breaks = c(1986:1989,1991:1994,1996:1999,2001:2004,2006:2009,2011:2014,2015:2019))+
  
  expand_limits(y=0)
