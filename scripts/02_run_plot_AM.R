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
library(ggplot2)
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
  
  #m<-models[3]
  
  #print scenario to check progress
  cat(paste(" #############   Model", m, match(m,models), 'out of',length(models),  "  #############\n"))
  
  # #get ss exe file in the model folder
  r4ss::get_ss3_exe(paste0(getwd(),'/models/',m))
  # 
  # #run SS model
  r4ss::run(dir = paste0(getwd(),'/models/',m),skipfinished = FALSE,show_in_console=TRUE)
  # 
  # #check Hessian matrix
  # r4ss::getADMBHessian(hesfile =paste0('./models/'models,m,'/admodel.hes'))$hes<0 #negative------ non-convergence
  
  #save results into a list
  mymodel<-SS_output(dir = paste0('./models/',m),verbose = T,printstats = T,forecast = F,covar = TRUE)
  Mycomparisonlist[[m]]<-mymodel
  
  # read the model output and print diagnostic messages
  replist <- SS_output(dir = paste0('./models/',m),
                       verbose = TRUE,
                       printstats = TRUE,covar = TRUE)

  # plots the results
  SS_plots(replist,uncertainty = TRUE,png = TRUE) #plot = c(1:2,4:7,9:22),
  #mymodel$maximum_gradient_component
  }

#check models in the list
names(Mycomparisonlist)

#save list
save(Mycomparisonlist, file="./models/summary_list.RData")
load("./models/summary_list.RData") #Mycomparisonlist

#summary list
Mycomparisonsummary<-SSsummarize(biglist = Mycomparisonlist)

#####################
# maximum final parameter gradient
#####################

#create a df with data
df<-
data.frame('model'=Mycomparisonsummary$modelnames,
           'max gradient'=Mycomparisonsummary$maxgrad)

# Set the levels first
df$model <- factor(df$model, levels = c('fixed', 
                                        'env.add', 'env.mul', 
                                        'blk.add.all','blk.add.sev',
                                        'blk.add3.all', 'blk.add3.sev', 
                                        'pred.all', 'pred.sev', 
                                        'bycatch.all', 'bycatch.sev', 
                                        'bycatchF.all', 'bycatchF.sev'))

# Assign new names to the factor levels
levels(df$model) <- c('fixed', 
                      'env_add', 'env_mul', 
                      'blk_all',  'blk_sev',
                      'blk3_all', 'blk3_sev', 
                      'pred_all', 'pred_sev', 
                      'bycatch_all', 'bycatch_sev', 
                      'bycatchF_all', 'bycatchF_sev')


#plot
p<-
ggplot(data=df,aes(y=max.gradient,x=model,fill=model))+
  geom_hline(yintercept = 0.0001,linetype='dashed',color='grey70')+
  geom_point(size=3,shape=21)+
  theme_bw()+
  labs(y='maximum gradient')+
  scale_y_continuous(limits = c(0,0.0003),expand=c(0,0))+
  scale_fill_manual(values = c('fixed'='black',
                                'env_add'='#01AD2C','env_mul'='#6FBD84',
                                'blk_all'='#B30101','blk3_all'='#E66000','blk_sev'='#C25460','blk3_sev'='#E09353',
                                'pred_all'='#BC9912','pred_sev'='#DECE8B',
                                'bycatch_all'='#154ABD','bycatch_sev'='#6DA4C2',
                                'bycatchF_all'='#7A29E5','bycatchF_sev'='#8A5FC1'),name='Assessment\nmodel')+
  guides(fill = guide_legend(override.aes = list(size=5,shape = rep(21,times=13)))) +  # Custom shapes in legend
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1),
        axis.title.x = element_blank())

#plot
ragg::agg_png('./figures/max_gradient.png',  width = 6, height = 5, units = "in", res = 300)
print(
  p
)
dev.off()

#####################
# deviance residuals and calculating the RMSE of abundance indices and age composition time series
#####################

#RMSE, p<value, mean residuals

#loop getting the 
Mycomparisonsummary$indices$Dev-(log(Mycomparisonsummary$indices$Obs)-log(Mycomparisonsummary$indices$Exp))
Mycomparisonsummary$indices$Dev-(log(Mycomparisonsummary$indices$Obs)-log(Mycomparisonsummary$indices$Exp))^2
Mycomparisonsummary$indices$Dev-(Mycomparisonsummary$indices$Obs-Mycomparisonsummary$indices$Exp)

df_diag<-
data.frame(models=names(Mycomparisonlist),
           RMSE_nearshore=NA,
           RMSE_offshore=NA,
           p_value_nearshore=NA,
           p_value_offshore=NA)


for (imodel in names(Mycomparisonlist)) {
  
  #imodel<-names(Mycomparisonlist)[1]
  
  # Extract model data
  m <- Mycomparisonlist[[imodel]]
  
  # Clean and fix m$cpue dataframe
  cpue_clean <- m$cpue
  colnames(cpue_clean) <- as.character(unlist(cpue_clean[1, ]))  # Use first row as header
  cpue_clean <- cpue_clean[-1, ]                                # Remove first row
  m$cpue <- cpue_clean
  m$cpue$Obs <- as.numeric(m$cpue$Obs)
  m$cpue$Exp <- as.numeric(m$cpue$Exp)
  m$cpue$Time <- as.numeric(m$cpue$Yr)
  
  # Calculate RMSE with SSplotJABBAres
  rmse <- ss3diags::SSplotJABBAres(ss3rep = m)
  
  # Store RMSE values in diagnostics dataframe
  df_diag[df_diag$models == imodel, c('RMSE_nearshore', 'RMSE_offshore')] <- rmse$RMSE.perc[1:2]
  
  # Extract and clean ind_data from summary for this model
  ind_data <- Mycomparisonsummary$indices[Mycomparisonsummary$indices$name == imodel, ]
  colnames(ind_data) <- as.character(unlist(ind_data[1, ]))  # Set header from first row
  ind_data <- ind_data[-1, ]                                 # Remove header row
  ind_data$Obs <- as.numeric(ind_data$Obs)
  ind_data$Exp <- as.numeric(ind_data$Exp)
  ind_data$Time <- as.numeric(ind_data$Yr)
  
  # Extract residuals and convert to numeric
  residuals_juv_num <- as.numeric(as.character(ind_data$Dev[ind_data$Fleet_name == 'SURVEYJUVENILE']))
  residuals_adult_num <- as.numeric(as.character(ind_data$Dev[ind_data$Fleet_name == 'SURVEYADULT']))
  
  # Calculate and store p-values from runs.test
  df_diag[df_diag$models == imodel, 'p_value_nearshore'] <- randtests::runs.test(residuals_juv_num)$p.value
  df_diag[df_diag$models == imodel, 'p_value_offshore'] <- randtests::runs.test(residuals_adult_num)$p.value
  
}


# Set the levels first
df_diag$models <- factor(df_diag$models, levels = c('fixed', 
                                            'env.add', 'env.mul', 
                                            'blk.add.all','blk.add.sev',
                                            'blk.add3.all', 'blk.add3.sev', 
                                            'pred.all', 'pred.sev', 
                                            'bycatch.all', 'bycatch.sev', 
                                            'bycatchF.all', 'bycatchF.sev'))

# Assign new names to the factor levels
levels(df_diag$models) <- c('fixed', 
                        'env_add', 'env_mul', 
                        'blk_all',  'blk_sev',
                        'blk3_all', 'blk3_sev', 
                        'pred_all', 'pred_sev', 
                        'bycatch_all', 'bycatch_sev', 
                        'bycatchF_all', 'bycatchF_sev')
ggplot()+
  geom_point(data=reshape2::melt(df_diag[,c('models','RMSE_nearshore','RMSE_offshore')]),aes(y=value,x=models,fill=models,shape=variable))+
  theme_bw()+
  scale_shape_manual(values = c('RMSE_offshore'=21,'RMSE_nearshore'=24),labels=c('offshore','nearshore'))+
  scale_fill_manual(values = c('fixed'='black',
                               'env_add'='#01AD2C','env_mul'='#6FBD84',
                               'blk_all'='#B30101','blk3_all'='#E66000','blk_sev'='#C25460','blk3_sev'='#E09353',
                               'pred_all'='#BC9912','pred_sev'='#DECE8B',
                               'bycatch_all'='#154ABD','bycatch_sev'='#6DA4C2',
                               'bycatchF_all'='#7A29E5','bycatchF_sev'='#8A5FC1'),name='Assessment\nmodel')+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1),axis.title.x = element_blank())+
  guides(shape = guide_legend(override.aes = list(shape = c(1,2))),
         fill = guide_legend(override.aes = list(size=5,shape = rep(22,times=13)))) +  # Custom shapes in legend
  scale_y_continuous(limits = c(0,NA))




#The runs tests were implemented using the function runs.test in the R package tseries (Trapletti, 2011). 
#This function calculates the 2-sided p-value of the Wald-Wolfowitz runs test, which is a nonparametric statistical test 
#that checks a randomness hypothesis for a data sequence. 

#####################
# RMSE - relative error for survey indices and age comp
#####################

#df<-Mycomparisonsummary$quants
#df_subset <- subset(df, grepl("SURVEY", Label))

df<-Mycomparisonsummary$indices

# Extract first row as header
new_header <- as.character(unlist(df[1, ]))

# Fix last two column names manually (adjust names as needed)
new_header[length(new_header)-1] <- "name"
new_header[length(new_header)] <- "imodel"

# Assign fixed header as colnames
colnames(df) <- new_header

# Remove the first row now that it’s used as header
df <- df[-1, ]

#numeric
df$Obs <- as.numeric(df$Obs)
df$Exp <- as.numeric(df$Exp)

df$obs_est2<-(df$Obs-df$Exp)^2

out<-data.frame(matrix(NA,nrow=0,ncol = 4))
names(out)<-c('imodel','Fleet','RMSE','RRMSE')

for (m in unique(df$name)) {
  for (sur in c('SURVEYJUVENILE','SURVEYADULT')) {
   
#sur<-unique(df$Fleet_name)[1];m<-unique(df$name)[1]
    
    df1<-subset(df,Fleet_name==sur & name==m)
    
    iout<-data.frame('imodel'=m,
               'Fleet'=sur,
               'RMSE'=(sqrt(sum(df1$obs_est2)/nrow(df1))),
               'RRMSE'=sqrt(sum(df1$obs_est2)/nrow(df1))/mean(df1$Exp))
    
    out<-rbind(out,iout)
     
  }
}

# Set the levels first
out$imodel <- factor(out$imodel, levels = c('fixed', 
                                        'env.add', 'env.mul', 
                                        'blk.add.all','blk.add.sev',
                                        'blk.add3.all', 'blk.add3.sev', 
                                        'pred.all', 'pred.sev', 
                                        'bycatch.all', 'bycatch.sev', 
                                        'bycatchF.all', 'bycatchF.sev'))

# Assign new names to the factor levels
levels(out$imodel) <- c('fixed', 
                      'env_add', 'env_mul', 
                      'blk_all',  'blk_sev',
                      'blk3_all', 'blk3_sev', 
                      'pred_all', 'pred_sev', 
                      'bycatch_all', 'bycatch_sev', 
                      'bycatchF_all', 'bycatchF_sev')


diag_pvalue<-reshape2::melt(df_diag[,c('models','p_value_nearshore','p_value_offshore')])


levels(diag_pvalue$variable)<-c('SURVEYJUVENILE','SURVEYADULT')

out<-merge(out,diag_pvalue,by.x = c('imodel','Fleet'),by.y=c('models','variable'))

#plot
# p<-
# ggplot(data=out,aes(y=RRMSE,x=imodel,fill=imodel,shape=Fleet))+
#   geom_point(size=3)+
#   theme_bw()+
#   scale_shape_manual(values = c('SURVEYADULT'=21,'SURVEYJUVENILE'=24),labels=c('offshore','nearshore'))+
#   scale_fill_manual(values = c('fixed'='black',
#                                 'env_add'='#01AD2C','env_mul'='#6FBD84',
#                                 'blk_all'='#B30101','blk3_all'='#E66000','blk_sev'='#C25460','blk3_sev'='#E09353',
#                                 'pred_all'='#BC9912','pred_sev'='#DECE8B',
#                                 'bycatch_all'='#154ABD','bycatch_sev'='#6DA4C2',
#                                 'bycatchF_all'='#7A29E5','bycatchF_sev'='#8A5FC1'),name='Assessment\nmodel')+
#   theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1),axis.title.x = element_blank())+
#   guides(shape = guide_legend(override.aes = list(shape = c(1,2))),
#          fill = guide_legend(override.aes = list(size=5,shape = rep(22,times=13)))) +  # Custom shapes in legend
#   scale_y_continuous(limits = c(0,NA))

p1<-
ggplot(data=out, aes(y=RRMSE, x=imodel, fill=imodel, shape=Fleet)) +
  geom_point(size=3) +
  geom_text(data = subset(out, value < 0.05), # Adjust condition for significance
            aes(label = "*"), 
            nudge_x = 0.3, # Adjust position of asterisks
            size = 6, color = "black", fontface = "bold") +
  theme_bw() +
  scale_shape_manual(values = c('SURVEYADULT' = 21, 'SURVEYJUVENILE' = 24),
                     labels = c('offshore', 'nearshore')) +
  scale_fill_manual(values = c('fixed' = 'black',
                               'env_add' = '#01AD2C', 'env_mul' = '#6FBD84',
                               'blk_all' = '#B30101', 'blk3_all' = '#E66000',
                               'blk_sev' = '#C25460', 'blk3_sev' = '#E09353',
                               'pred_all' = '#BC9912', 'pred_sev' = '#DECE8B',
                               'bycatch_all' = '#154ABD', 'bycatch_sev' = '#6DA4C2',
                               'bycatchF_all' = '#7A29E5', 'bycatchF_sev' = '#8A5FC1'),
                    name = 'Assessment\nmodel') +
  theme(#axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),axis.text.x = element_blank()) +
  guides(shape = guide_legend(override.aes = list(shape = c(1, 2))),
         fill = guide_legend(override.aes = list(size = 4, shape = rep(22, times = 13)))) +  
  #scale_y_continuous(limits = c(0.05,0.35))+
  annotate("text", x = 1, y = 0.35, label = "Survey indices", hjust = 0, size = 5,vjust=0.9)



# ggplot(data=out,aes(y=RMSE,x=imodel,fill=imodel,shape=Fleet))+
#   geom_point(size=3)+
#   theme_bw()+
#   scale_shape_manual(values = c('SURVEYADULT'=21,'SURVEYJUVENILE'=24),labels=c('offshore','nearshore'))+
#   scale_fill_manual(values = c('fixed'='black',
#                                'env_add'='#01AD2C','env_mul'='#6FBD84',
#                                'blk_all'='#B30101','blk3_all'='#E66000','blk_sev'='#C25460','blk3_sev'='#E09353',
#                                'pred_all'='#BC9912','pred_sev'='#DECE8B',
#                                'bycatch_all'='#154ABD','bycatch_sev'='#6DA4C2',
#                                'bycatchF_all'='#7A29E5','bycatchF_sev'='#8A5FC1'),name='Assessment\nmodel')+
#   theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1),axis.title.x = element_blank())+
#   guides(shape = guide_legend(override.aes = list(shape = c(1,2))),
#          fill = guide_legend(override.aes = list(size=5,shape = rep(22,times=13)))) +  # Custom shapes in legend
#   scale_y_continuous(limits = c(0,NA))

  
#plot
ragg::agg_png('./figures/RRMSE_surveys.png',  width = 6, height = 5.5, units = "in", res = 300)
print(
  p1
)
dev.off()

##AGE COMP

age_df<-
data.frame(imodel=names(Mycomparisonlist),
           rrmse_fish=NA,
           rrmse_juv=NA,
           rrmse_adu=NA,
           p_value_fish=NA,
           p_value_juv=NA,
           p_value_adu=NA)


for (imodel in names(Mycomparisonlist)) {
  
  #imodel<-names(Mycomparisonlist)[1]
  
  m<-Mycomparisonlist[[imodel]]
  
  rmse<-ss3diags::SSplotJABBAres(ss3rep = m,subplots = 'age')
  
  agedata<-ss3diags::SScompsTA1.8(m,type='age',plotit = FALSE)$runs_dat
  agedata$dev<-agedata$Obs-agedata$Exp
  agedata$obs_est2<-(agedata$Obs-agedata$Exp)^2
    
    agedata1<-subset(agedata,Fleet_name=='FISHERY')
    age_df[which(age_df$imodel==imodel),'rrmse_fish']<-sqrt(mean(agedata1$obs_est2))/mean(agedata1$Exp)

    agedata1<-subset(agedata,Fleet_name=='SURVEYADULT')
    age_df[which(age_df$imodel==imodel),'rrmse_adu']<-sqrt(mean(agedata1$obs_est2))/mean(agedata1$Exp)
    
    agedata1<-subset(agedata,Fleet_name=='SURVEYJUVENILE')
    age_df[which(age_df$imodel==imodel),'rrmse_juv']<-sqrt(mean(agedata1$obs_est2))/mean(agedata1$Exp)
  
  
  residuals_adult<-agedata[which(agedata$Fleet_name=='SURVEYADULT'),'dev']
  residuals_juv<-agedata[which(agedata$Fleet_name=='SURVEYJUVENILE'),'dev']
  residuals_fish<-agedata[which(agedata$Fleet_name=='FISHERY'),'dev']
  age_df[which(age_df$imodel==imodel),'p_value_juv']<-randtests::runs.test(residuals_juv)$p.value
  age_df[which(age_df$imodel==imodel),'p_value_adu']<-randtests::runs.test(residuals_adult)$p.value
  age_df[which(age_df$imodel==imodel),'p_value_fish']<-randtests::runs.test(residuals_fish)$p.value
  
  
}



# Set the levels first
age_df$imodel <-factor(age_df$imodel, levels = c('fixed', 
                                            'env.add', 'env.mul', 
                                            'blk.add.all','blk.add.sev',
                                            'blk.add3.all', 'blk.add3.sev', 
                                            'pred.all', 'pred.sev', 
                                            'bycatch.all', 'bycatch.sev', 
                                            'bycatchF.all', 'bycatchF.sev'))

# Assign new names to the factor levels
levels(age_df$imodel) <- c('fixed', 
                        'env_add', 'env_mul', 
                        'blk_all',  'blk_sev',
                        'blk3_all', 'blk3_sev', 
                        'pred_all', 'pred_sev', 
                        'bycatch_all', 'bycatch_sev', 
                        'bycatchF_all', 'bycatchF_sev')

diag_pvalue<-reshape2::melt(age_df[,c('imodel','p_value_juv','p_value_adu','p_value_fish')])
diag_rrmse<-reshape2::melt(age_df[,c('imodel','rrmse_juv','rrmse_adu','rrmse_fish')])

levels(diag_rrmse$variable)<-c('SURVEYJUVENILE','SURVEYADULT','FISHERY')
levels(diag_pvalue$variable)<-c('SURVEYJUVENILE','SURVEYADULT','FISHERY')

# Merge diag_rrmse and diag_pvalue based on 'imodel' (or another appropriate key)
#merged_data <- merge(diag_rrmse, diag_pvalue[, c("imodel", "value")], by = "imodel", suffixes = c("_rrmse", "_pvalue"))


diag_age<-cbind(diag_rrmse,pvalue=diag_pvalue$value)

p2<-
ggplot(data=diag_age, aes(y=value, x=imodel, fill=imodel, shape=variable)) +
  geom_point(size=3) +
  geom_text(data = subset(diag_age, pvalue < 0.01), # Adjust condition for significance
            aes(label = "*"), 
            nudge_x = 0.3, # Adjust position of asterisks
            size = 6, color = "black", fontface = "bold") +
  theme_bw() +
  labs(y='RRMSE')+
  scale_shape_manual(values = c('SURVEYADULT' = 21, 'SURVEYJUVENILE' = 24, 'FISHERY'=23),
                     labels = c('offshore', 'nearshore','fishery'),name='Fleet') +
  scale_fill_manual(values = c('fixed' = 'black',
                               'env_add' = '#01AD2C', 'env_mul' = '#6FBD84',
                               'blk_all' = '#B30101', 'blk3_all' = '#E66000',
                               'blk_sev' = '#C25460', 'blk3_sev' = '#E09353',
                               'pred_all' = '#BC9912', 'pred_sev' = '#DECE8B',
                               'bycatch_all' = '#154ABD', 'bycatch_sev' = '#6DA4C2',
                               'bycatchF_all' = '#7A29E5', 'bycatchF_sev' = '#8A5FC1'),
                    name = 'Assessment\nmodel') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank()) +
  guides(shape = guide_legend(override.aes = list(shape = c(1, 2,5)),order=2),
         fill = guide_legend(override.aes = list(size = 4, shape = rep(22, times = 13)),order=1)) +  
  #scale_y_continuous(limits = c(0.05,0.35))+
  annotate("text", x = 1, y = 0.23, label = "Age composition", hjust = 0, size = 5,vjust=0.9)

#plot
ragg::agg_png('./figures/RRMSE_agecomp.png',  width = 6, height = 5.5, units = "in", res = 300)
print(
  p2
)
dev.off()


# Extract the legend from one of the plots
shared_legend <- cowplot::get_legend(
  p2 + theme(legend.position = "right")
)

# Combine plots and the shared legend
final_plot <- cowplot::plot_grid(
  cowplot::plot_grid(p1 + theme(legend.position = "none"),
            p2 + theme(legend.position = "none"),
            ncol = 1,
            align = 'v',rel_heights = c(0.5, 0.6)),
  shared_legend,
  ncol = 2,
  rel_widths = c(0.8, 0.2) # Adjust the width ratio for the plots and the legend
)

ragg::agg_png(paste0('./figures/age_ind_rrmse.png'), width = 8, height = 6, units = "in", res = 300)
final_plot
dev.off()

#####################
# Likelihood components
#####################

#reshape
df<-reshape2::melt(Mycomparisonsummary$likelihoods)

#subset
df<-subset(df,Label %in% c('Age_comp','Catch','Recruitment','Survey','TOTAL'))

# Set the levels first
df$variable <- factor(df$variable, levels = c('fixed', 
                                            'env.add', 'env.mul', 
                                            'blk.add.all','blk.add.sev',
                                            'blk.add3.all', 'blk.add3.sev', 
                                            'pred.all', 'pred.sev', 
                                            'bycatch.all', 'bycatch.sev', 
                                            'bycatchF.all', 'bycatchF.sev'))

# Assign new names to the factor levels
levels(df$variable) <- c('fixed', 
                        'env_add', 'env_mul', 
                        'blk_all',  'blk_sev',
                        'blk3_all', 'blk3_sev', 
                        'pred_all', 'pred_sev', 
                        'bycatch_all', 'bycatch_sev', 
                        'bycatchF_all', 'bycatchF_sev')

#rename facets
df$Label<-as.factor(df$Label)
levels(df$Label)
levels(df$Label)<-c('Age composition','Catch','Recruitment',
                    'Survey','Total')
#plot
p<-
ggplot(data=df,aes(y=value,x=variable,fill=variable))+
  geom_point(size=3,shape=21)+
  theme_bw()+
  labs(y='negative log(likelihood)')+
  scale_fill_manual(values = c('fixed'='black',
                               'env_add'='#01AD2C','env_mul'='#6FBD84',
                               'blk_all'='#B30101','blk3_all'='#E66000','blk_sev'='#C25460','blk3_sev'='#E09353',
                               'pred_all'='#BC9912','pred_sev'='#DECE8B',
                               'bycatch_all'='#154ABD','bycatch_sev'='#6DA4C2',
                               'bycatchF_all'='#7A29E5','bycatchF_sev'='#8A5FC1'),name='Assessment\nmodel')+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1),axis.title.x = element_blank(),strip.background = element_rect(fill='white'))+
  facet_wrap(~Label,scales = 'free_y',ncol = 1)
#scale_y_continuous(limits = c(0,0.0001))

#Large negative values indicate components that are well-fit by the model.

#plot
ragg::agg_png('./figures/NLLs.tiff',  width = 6, height = 8, units = "in", res = 300)
print(
  p
)
dev.off()

#####################
# timeseries comparison of recruitment, M, F, SSB and SPR relative to OBS values
#####################

df<-Mycomparisonsummary$quants
unique(df$Label)

df<-subset(df,Label %in% c(paste0('SSB_',1985:2020),
                           paste0('SPRratio_',1985:2020),
                           paste0('F_',1985:2020),
                           paste0('Recr_',1985:2020),
                           paste0('SURVEYJUVENILE_',1985:2020),
                           paste0('SURVEYADULT_',1985:2020)))


df$Label<-gsub('_\\d{4}','',df$Label)

df<-reshape2::melt(df,id.vars=c('Label','Yr'))

# Set the levels first
df$variable <- factor(df$variable, levels = c('fixed', 
                                              'env.add', 'env.mul', 
                                              'blk.add.all','blk.add.sev',
                                              'blk.add3.all', 'blk.add3.sev', 
                                              'pred.all', 'pred.sev', 
                                              'bycatch.all', 'bycatch.sev', 
                                              'bycatchF.all', 'bycatchF.sev'))

# Assign new names to the factor levels
levels(df$variable) <- c('fixed', 
                         'env_add', 'env_mul', 
                         'blk_all',  'blk_sev',
                         'blk3_all', 'blk3_sev', 
                         'pred_all', 'pred_sev', 
                         'bycatch_all', 'bycatch_sev', 
                         'bycatchF_all', 'bycatchF_sev')

#rename facets
df$Label<-as.factor(df$Label)
levels(df$Label)
# levels(df$Label)<-c('Age composition',
#                     'Catch',
#                     'Recruitment',
#                     'Survey',
#                     'Total')
#plot
#p<-
  ggplot(data=df,aes(y=value,x=Yr,color=variable))+
  geom_line()+
  theme_bw()+
  labs(y='')+
  scale_color_manual(values = c('fixed'='black',
                               'env_add'='#01AD2C','env_mul'='#6FBD84',
                               'blk_all'='#B30101','blk3_all'='#E66000','blk_sev'='#C25460','blk3_sev'='#E09353',
                               'pred_all'='#BC9912','pred_sev'='#DECE8B',
                               'bycatch_all'='#154ABD','bycatch_sev'='#6DA4C2',
                               'bycatchF_all'='#7A29E5','bycatchF_sev'='#8A5FC1'),name='Assessment\nmodel')+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1),axis.title.x = element_blank(),strip.background = element_rect(fill='white'))+
  facet_wrap(~Label,scales = 'free_y',ncol = 2)
#scale_y_continuous(limits = c(0,0.0001))

#RECDEVS
recdevs<-Mycomparisonsummary$recdevs[-nrow(Mycomparisonsummary$recdevs),]
# recdevslow<-Mycomparisonsummary$recdevsLower[-nrow(Mycomparisonsummary$recdevsLower),]
# recdevsup<-Mycomparisonsummary$recdevsUpper[-nrow(Mycomparisonsummary$recdevsUpper),]
recdevsSD<-Mycomparisonsummary$recdevsSD[-nrow(Mycomparisonsummary$recdevsSD),]
recdevs_df<-cbind(reshape2::melt(recdevs,id.vars=c('Label','Yr')),
                  'SD'=reshape2::melt(recdevsSD,id.vars=c('Label','Yr'))[,'value'])

# #Fs
# fs<-Mycomparisonsummary$Fvalue[-nrow(Mycomparisonsummary$Fvalue),]
# fslow<-Mycomparisonsummary$FvalueLower[-nrow(Mycomparisonsummary$FvalueLower),]
# fsup<-Mycomparisonsummary$FvalueUpper[-nrow(Mycomparisonsummary$FvalueUpper),]
# fs_df<-cbind(reshape2::melt(fs,id.vars=c('Label','Yr')),
#                   'lower'=reshape2::melt(fslow,id.vars=c('Label','Yr'))[,'value'],
#                   'upper'=reshape2::melt(fsup,id.vars=c('Label','Yr'))[,'value'])
# 
# #RECRUITS
# recruits<-Mycomparisonsummary$recruits[-nrow(Mycomparisonsummary$recruits),]
# recruitslow<-Mycomparisonsummary$recruitsLower[-nrow(Mycomparisonsummary$recruitsLower),]
# recruitsup<-Mycomparisonsummary$recruitsUpper[-nrow(Mycomparisonsummary$recruitsUpper),]
# recruits_df<-cbind(reshape2::melt(recruits,id.vars=c('Label','Yr')),
#              'lower'=reshape2::melt(recruitslow,id.vars=c('Label','Yr'))[,'value'],
#              'upper'=reshape2::melt(recruitsup,id.vars=c('Label','Yr'))[,'value'])
# 
# #SPBs
# SPRratios<-Mycomparisonsummary$SPRratio[-nrow(Mycomparisonsummary$SPRratio),]
# SPRratioslow<-Mycomparisonsummary$SPRratioLower[-nrow(Mycomparisonsummary$SPRratioLower),]
# SPRratiosup<-Mycomparisonsummary$SPRratioUpper[-nrow(Mycomparisonsummary$SPRratioUpper),]
# SPRratios_df<-cbind(reshape2::melt(SPRratios,id.vars=c('Label','Yr')),
#                    'lower'=reshape2::melt(SPRratioslow,id.vars=c('Label','Yr'))[,'value'],
#                    'upper'=reshape2::melt(SPRratiosup,id.vars=c('Label','Yr'))[,'value'])
# 
# #SPRs
# SpawnBios<-Mycomparisonsummary$SpawnBio[-nrow(Mycomparisonsummary$SpawnBio),]
# SpawnBioslow<-Mycomparisonsummary$SpawnBioLower[-nrow(Mycomparisonsummary$SpawnBioLower),]
# SpawnBiosup<-Mycomparisonsummary$SpawnBioUpper[-nrow(Mycomparisonsummary$SpawnBioUpper),]
# SpawnBios_df<-cbind(reshape2::melt(SpawnBios,id.vars=c('Label','Yr')),
#                    'lower'=reshape2::melt(SpawnBioslow,id.vars=c('Label','Yr'))[,'value'],
#                    'upper'=reshape2::melt(SpawnBiosup,id.vars=c('Label','Yr'))[,'value'])
# SpawnBios_df<-subset(SpawnBios_df,Yr %in% 1985:2020)
# 
# #Indices
# index<-Mycomparisonsummary$indices[,c('Fleet_name')]

quants_df<-reshape2::melt(Mycomparisonsummary$quants,id.vars=c('Label','Yr'))
quantsSD_df<-reshape2::melt(Mycomparisonsummary$quantsSD,id.vars=c('Label','Yr'))

quants_df<-subset(quants_df,Label %in% c(paste0('SSB_',1985:2020),
                           paste0('SPRratio_',1985:2020),
                           paste0('F_',1985:2020),
                           paste0('Recr_',1985:2020),
                           paste0('SURVEYJUVENILE_',1985:2020),
                           paste0('SURVEYADULT_',1985:2020)))

quantsSD_df<-subset(quantsSD_df,Label %in% c(paste0('SSB_',1985:2020),
                                         paste0('SPRratio_',1985:2020),
                                         paste0('F_',1985:2020),
                                         paste0('Recr_',1985:2020),
                                         paste0('SURVEYJUVENILE_',1985:2020),
                                         paste0('SURVEYADULT_',1985:2020)))

quants<-cbind(quants_df,'SD'=quantsSD_df$value)

library(dplyr)
library(reshape2)
library(ggplot2)

# Assuming recdevs_df and quants are ready and have columns: Label, Yr, variable, value, SD, Observed (if joined)

# Combine data frames (use bind_rows for safety)
quants1 <- bind_rows(recdevs_df, quants)

# Set factor levels for 'variable' (once)
quants1$variable <- factor(quants1$variable, levels = c(
  'fixed', 'env.add', 'env.mul',
  'blk.add.all','blk.add.sev',
  'blk.add3.all','blk.add3.sev',
  'pred.all','pred.sev',
  'bycatch.all','bycatch.sev',
  'bycatchF.all','bycatchF.sev'
))

# Rename factor levels for 'variable'
levels(quants1$variable) <- c(
  'fixed', 'env_add', 'env_mul',
  'blk_all', 'blk_sev',
  'blk3_all', 'blk3_sev',
  'pred_all', 'pred_sev',
  'bycatch_all', 'bycatch_sev',
  'bycatchF_all', 'bycatchF_sev'
)

# Clean Label by removing year suffix "_YYYY"
quants1$Label <- gsub('_\\d{4}', '', quants1$Label)

# Set factor levels for Label to control plotting order
quants1$Label <- factor(quants1$Label, levels = c(
  "Main_RecrDev",
  "Recr",
  "SPRratio",
  "SSB",
  "SURVEYJUVENILE",
  "SURVEYADULT",
  "F"
))

# Subset survey data for ratio calculation
quants2 <- subset(quants1, Label %in% c('SURVEYJUVENILE','SURVEYADULT'))
quants2$ratio <- quants2$value / quants2$SD

# Facet labels for nicer plot labels
facet_labels <- c(
  "Main_RecrDev" = "Recruitment deviations",
  "SSB" = "SSB",
  "Recr" = "Recruits",
  "SPRratio" = "SPR Ratio",
  "F" = "Fishing Mortality (F)",
  "SURVEYJUVENILE" = "Nearshore Survey",
  "SURVEYADULT" = "Offshore Survey"
)


#sort levels
quants1$variable<-factor(quants1$variable,levels = c('fixed',
                                                     'env.add','env.mul',
                                                     'blk.add.all','blk.add.sev','blk.add3.all','blk.add3.sev',
                                                     'pred.all','pred.sev',
                                                     'bycatch.all','bycatch.sev',
                                                     'bycatchF.all','bycatchF.sev'))

# Set the levels first
quants1$variable <- factor(quants1$variable, levels = c('fixed', 
                                                        'env.add', 'env.mul', 
                                                        'blk.add.all','blk.add.sev',
                                                        'blk.add3.all', 'blk.add3.sev', 
                                                        'pred.all', 'pred.sev', 
                                                        'bycatch.all', 'bycatch.sev', 
                                                        'bycatchF.all', 'bycatchF.sev'))

# Assign new names to the factor levels
levels(quants1$variable) <- c('fixed', 
                              'env_add', 'env_mul', 
                              'blk_all',  'blk_sev',
                              'blk3_all', 'blk3_sev', 
                              'pred_all', 'pred_sev', 
                              'bycatch_all', 'bycatch_sev', 
                              'bycatchF_all', 'bycatchF_sev')

quants1$Label<-gsub('_\\d{4}','',quants1$Label)
#quants1<-subset(quants1,Label %in% c('SSB','SPRratio','Recr','Main_RecrDev','F'))

quants2<-subset(quants1,Label %in% c('SURVEYJUVENILE','SURVEYADULT'))
quants2$ratio<-quants2$value/quants2$SD

# Create a named vector for relabeling
facet_labels <- c(
  "Main_RecrDev" = "Recruitment deviations",
  "SSB" = "SSB",
  "Recr" = "Recruits",
  "SPRratio" = "SPR Ratio",
  "F" = "F",
  "SURVEYJUVENILE" = "nearshore survey",
  "SURVEYADULT" = "offshore survey"
)

quants1$Label<-factor(quants1$Label,levels = c(  "Main_RecrDev" ,
                                                 "Recr" ,
                                                 "SPRratio",
                                                 "SSB" ,
                                                 "SURVEYJUVENILE" ,
                                                 "SURVEYADULT","F"  ))

quants1 <- quants1 %>%
  mutate(
    predicted_backtrans = ifelse(Label %in% c("SURVEYJUVENILE", "SURVEYADULT"),
                                 exp(value),   # back-transform only surveys
                                 value)        # keep others unchanged
  )


# Create a new data frame for observed points with aesthetics to appear in legend
observed_points <- quants1 %>% 
  filter(!is.na(Observed)) %>%
  mutate(variable = "Observed")  # Assign a new 'variable' factor for observed points

# Add Observed to factor levels of variable for proper legend handling
quants1$variable <- factor(quants1$variable, levels = c(levels(quants1$variable), "Observed"))
observed_points$variable <- factor(observed_points$variable, levels = levels(quants1$variable))

quants1 <- quants1 %>%
  mutate(
    SD_real = ifelse(Label %in% c("SURVEYJUVENILE", "SURVEYADULT"),
                     predicted_backtrans * sqrt(exp(SD^2) - 1),
                     SD)  # Keep original SD for non-survey data
  )

# Example plot: Value and Observed over Years faceted by Label, colored by variable
ggplot() +
  geom_line(data=quants1, aes(y=predicted_backtrans, x=Yr, color=variable, group=interaction(variable,Label))) +
  geom_ribbon(
    data = quants1,
    aes(
      ymin = predicted_backtrans - SD_real,
      ymax = predicted_backtrans + SD_real,
      x = Yr,
      fill = variable,
      group = interaction(variable, Label)
    ),
    alpha = 0.2
  )+
  geom_point(data = quants1 %>% filter(!is.na(Observed)),
             aes(x = Yr, y = Observed, color = "Observed", shape = "Observed"),
             size = 1.5, alpha = 0.7) +
  
  scale_color_manual(
    values = c(
      'fixed' = 'black',
      'env_add' = '#01AD2C', 'env_mul' = '#6FBD84',
      'blk_all' = '#B30101', 'blk3_all' = '#E66000', 'blk_sev' = '#C25460', 'blk3_sev' = '#E09353',
      'pred_all' = '#BC9912', 'pred_sev' = '#DECE8B',
      'bycatch_all' = '#154ABD', 'bycatch_sev' = '#6DA4C2',
      'bycatchF_all' = '#7A29E5', 'bycatchF_sev' = '#8A5FC1',
      'Observed' = 'black'  # observed points
    ),
    name = 'Assessment model'
  ) +
  
  scale_shape_manual(
    values = c('Observed' = 16),
    name = ''
  )+
  theme_bw() +
  labs(y = "") +
  scale_color_manual(
    values = c(
      'fixed' = 'black',
      'env_add' = '#01AD2C', 'env_mul' = '#6FBD84',
      'blk_all' = '#B30101', 'blk3_all' = '#E66000', 'blk_sev' = '#C25460', 'blk3_sev' = '#E09353',
      'pred_all' = '#BC9912', 'pred_sev' = '#DECE8B',
      'bycatch_all' = '#154ABD', 'bycatch_sev' = '#6DA4C2',
      'bycatchF_all' = '#7A29E5', 'bycatchF_sev' = '#8A5FC1'
    ),
    name = 'Assessment model'
  ) +
  scale_fill_manual(
    values = c(
      'fixed' = 'black',
      'env_add' = '#01AD2C', 'env_mul' = '#6FBD84',
      'blk_all' = '#B30101', 'blk3_all' = '#E66000', 'blk_sev' = '#C25460', 'blk3_sev' = '#E09353',
      'pred_all' = '#BC9912', 'pred_sev' = '#DECE8B',
      'bycatch_all' = '#154ABD', 'bycatch_sev' = '#6DA4C2',
      'bycatchF_all' = '#7A29E5', 'bycatchF_sev' = '#8A5FC1'
    ),
    name = 'Assessment model'
  ) +
  theme(
    axis.title.x = element_blank(),
    legend.position = c(0.7, 0.15),
    strip.background = element_rect(fill = 'white'),
    text = element_text(size = 12)
  ) +
  facet_wrap(~Label, scales = 'free_y', nrow = 3, labeller = labeller(Label = facet_labels)) +
  guides(color = guide_legend(ncol = 3), fill = guide_legend(ncol = 3))





#plot
#p<-
ggplot()+
  #geom_point(size=1)+
  #geom_point(data=quants1,aes(y=value,x=Yr,color=variable,group=interaction(variable,Label)))+
  #geom_errorbar(data=quants1,aes(ymin=value-SD,ymax=value+SD,x=Yr,color=variable,group=interaction(variable,Label)))+
  geom_line(data=quants1,aes(y=value,x=Yr,color=variable,group=interaction(variable,Label)))+
  #geom_line(data=subset(quants1,variable=='bycatch.all'),aes(y=value,x=Yr,color=variable,group=interaction(variable,Label)))+
  
  geom_ribbon(data=quants1,aes(ymin=value-SD,ymax=value+SD,x=Yr,fill=variable,group=interaction(variable,Label)),alpha=0.2)+
  #geom_point(data=obs,aes(x=Yr,y=Obs))+
  theme_bw()+
  labs(y='')+
  scale_color_manual(values = c('fixed'='black',
                                'env_add'='#01AD2C','env_mul'='#6FBD84',
                                'blk_all'='#B30101','blk3_all'='#E66000','blk_sev'='#C25460','blk3_sev'='#E09353',
                                'pred_all'='#BC9912','pred_sev'='#DECE8B',
                                'bycatch_all'='#154ABD','bycatch_sev'='#6DA4C2',
                                'bycatchF_all'='#7A29E5','bycatchF_sev'='#8A5FC1'),name='Assessment model')+
  scale_fill_manual(values = c('fixed'='black',
                                'env_add'='#01AD2C','env_mul'='#6FBD84',
                                'blk_all'='#B30101','blk3_all'='#E66000','blk_sev'='#C25460','blk3_sev'='#E09353',
                                'pred_all'='#BC9912','pred_sev'='#DECE8B',
                                'bycatch_all'='#154ABD','bycatch_sev'='#6DA4C2',
                                'bycatchF_all'='#7A29E5','bycatchF_sev'='#8A5FC1'),name='Assessment model')+
  theme(axis.title.x = element_blank(),
        legend.position = c(0.7,0.15),
        strip.background = element_rect(fill='white'),text = element_text(size=12))+
  facet_wrap(~Label, scales = 'free_y', nrow = 3, labeller = labeller(Label = facet_labels))+
  guides(color = guide_legend(ncol = 3), fill = guide_legend(ncol = 3))
#scale_y_continuous(limits = c(0,0.0001))

#plot
ragg::agg_png('./figures/derived_quants.png',  width =7, height = 6, units = "in", res = 300)
print(
  p
)
dev.off()



m<-Mycomparisonlist$fixed

# Clean and fix m$cpue dataframe
cpue_clean <- m$cpue
colnames(cpue_clean) <- as.character(unlist(cpue_clean[1, ]))  # Use first row as header
cpue_clean <- cpue_clean[-1, ]                                # Remove first row
m$cpue <- cpue_clean
m$cpue$Obs <- as.numeric(m$cpue$Obs)
m$cpue$Exp <- as.numeric(m$cpue$Exp)
m$cpue$Time <- as.numeric(m$cpue$Yr)

# Example observed values for recruitment and biomass
obs_bio <- data.frame(
  Label = paste0("SSB_", m$recruit$Yr),
  Yr = m$recruit$Yr,
  Observed = m$recruit$SpawnBio  # or replace with actual observed biomass
)

obs_rec <- data.frame(
  Label = paste0("Recr_", m$recruit$Yr),
  Yr = m$recruit$Yr,
  Observed = m$recruit$pred_recr  # or actual observed recruitment if available
)

obs_f <- data.frame(
  Label = paste0("F_", m$timeseries$Yr),
  Yr = m$timeseries$Yr,
  Observed = m$timeseries$F  # fishing mortality
)

# Survey observed data
survey_obs <- data.frame(
  Label = paste0("SURVEY", m$cpue$Fleet, "_", m$cpue$Yr),
  Yr = m$cpue$Yr,
  Observed = m$cpue$Obs
)

survey_obs$Label <- ifelse(m$cpue$Fleet_name == 'SURVEYJUVENILE',
                           paste0("SURVEYJUVENILE_", m$cpue$Yr),
                           paste0("SURVEYADULT_", m$cpue$Yr))

# Combine observed data into one dataframe
obs_quants <- rbind(obs_bio, obs_rec, obs_f, survey_obs)

obs_quants$Yr <- as.numeric(obs_quants$Yr)


# Join the Observed values into your quants dataframe by Label and Yr
library(dplyr)
quants <- quants %>%
  left_join(obs_quants %>% select(Label, Yr, Observed), by = c("Label", "Yr"))

# Remove Observed and SD columns if they exist in quants_df
quants_df$Observed <- NULL
quants_df$SD <- NULL

# Extract unique observed values by Label and Yr
obs_unique <- unique(obs_quants[, c("Label", "Yr", "Observed")])

# Merge observed values by Label and Yr into quants_df
quants_df_obs <- merge(quants_df, obs_unique, by = c("Label", "Yr"), all.x = TRUE)

# Add SD column from quantsSD_df (make sure it aligns row-wise)
quants <- cbind(quants_df_obs, SD = quantsSD_df$value)

# Prepare recdevs_df: add Observed = NA and reorder columns to match quants
recdevs_df$Observed <- NA
recdevs_df <- recdevs_df[, c("Label", "Yr", "variable", "value", "Observed", "SD")]

# Combine recdevs and quants for plotting
quants1 <- rbind(recdevs_df, quants)

# Plot
library(ggplot2)

ggplot() +
  geom_line(data = quants1, aes(y = value, x = Yr, color = variable, group = interaction(variable, Label))) +
  geom_ribbon(data = quants1, aes(ymin = value - SD, ymax = value + SD, x = Yr, fill = variable, group = interaction(variable, Label)), alpha = 0.2) +
  #geom_point(data = quants1, aes(x = Yr, y = Observed, shape = variable, group = interaction(variable, Label)), color = "black", size = 1.5) +
  theme_bw() +
  labs(y = '') +
  scale_color_manual(values = c(
    'fixed' = 'black',
    'env_add' = '#01AD2C', 'env_mul' = '#6FBD84',
    'blk_all' = '#B30101', 'blk3_all' = '#E66000',
    'blk_sev' = '#C25460', 'blk3_sev' = '#E09353',
    'pred_all' = '#BC9912', 'pred_sev' = '#DECE8B',
    'bycatch_all' = '#154ABD', 'bycatch_sev' = '#6DA4C2',
    'bycatchF_all' = '#7A29E5', 'bycatchF_sev' = '#8A5FC1'
  ), name = 'Assessment model') +
  scale_fill_manual(values = c(
    'fixed' = 'black',
    'env_add' = '#01AD2C', 'env_mul' = '#6FBD84',
    'blk_all' = '#B30101', 'blk3_all' = '#E66000',
    'blk_sev' = '#C25460', 'blk3_sev' = '#E09353',
    'pred_all' = '#BC9912', 'pred_sev' = '#DECE8B',
    'bycatch_all' = '#154ABD', 'bycatch_sev' = '#6DA4C2',
    'bycatchF_all' = '#7A29E5', 'bycatchF_sev' = '#8A5FC1'
  ), name = 'Assessment model') +
  theme(
    axis.title.x = element_blank(),
    legend.position = c(0.7, 0.15),
    strip.background = element_rect(fill = 'white'),
    text = element_text(size = 12)
  ) +
  facet_wrap(~Label, scales = 'free_y', nrow = 3, labeller = labeller(Label = facet_labels)) +
  guides(color = guide_legend(ncol = 3), fill = guide_legend(ncol = 3))


#####################
# Index
#####################

#obs
obs<-Mycomparisonsummary$indices[c('Fleet_name','Yr','Obs',"Exp",'name',"SE")]
obs<-subset(obs,Fleet_name %in% c('SURVEYADULT','SURVEYJUVENILE'))
names(obs)[c(1,5)]<-c('Label','variable')

# Set the levels first
obs$variable <- factor(obs$variable, levels = c('fixed', 
                                                        'env.add', 'env.mul', 
                                                        'blk.add.all','blk.add.sev',
                                                        'blk.add3.all', 'blk.add3.sev', 
                                                        'pred.all', 'pred.sev', 
                                                        'bycatch.all', 'bycatch.sev', 
                                                        'bycatchF.all', 'bycatchF.sev'))

# Assign new names to the factor levels
levels(obs$variable) <- c('fixed', 
                              'env_add', 'env_mul', 
                              'blk_all',  'blk_sev',
                              'blk3_all', 'blk3_sev', 
                              'pred_all', 'pred_sev', 
                              'bycatch_all', 'bycatch_sev', 
                              'bycatchF_all', 'bycatchF_sev')





obs1<-merge(obs,quants2,by=c("Label",'Yr','variable'))

meanvalues<-aggregate(value ~ Label,obs1,FUN=mean)
meanvalues<-reshape2::melt(meanvalues)[c(1,3)]
names(meanvalues)[2]<-'meanvalues'

obs1<-merge(obs1,meanvalues,by=c('Label'))

meanvalues<-aggregate(Exp ~ Label,obs1,FUN=mean)
meanvalues<-reshape2::melt(meanvalues)[c(1,3)]
names(meanvalues)[2]<-'meanExp'

obs1<-merge(obs1,meanvalues,by=c('Label'))

obs1$valueRatioSD<-obs1$SD/obs1$meanvalues
obs1$ExpSD<-obs1$valueRatioSD*obs1$meanExp

obs1$Label<-as.factor(obs1$Label)
levels(obs1$Label)<-c('offshore','nearshore')

#sort levels
obs1$variable<-factor(obs1$variable,levels = c('fixed', 
                                               'env_add', 'env_mul', 
                                               'blk_all',  'blk_sev',
                                               'blk3_all', 'blk3_sev', 
                                               'pred_all', 'pred_sev', 
                                               'bycatch_all', 'bycatch_sev', 
                                               'bycatchF_all', 'bycatchF_sev'))


# Assign new names to the factor levels
levels(obs1$variable) <- c('fixed', 
                          'env_add', 'env_mul', 
                          'blk_all',  'blk_sev',
                          'blk3_all', 'blk3_sev', 
                          'pred_all', 'pred_sev', 
                          'bycatch_all', 'bycatch_sev', 
                          'bycatchF_all', 'bycatchF_sev')


#plot
ggplot()+
  #geom_point(size=1)+
  geom_point(data=obs1,aes(y=Obs,x=Yr,group=interaction(variable,Label),shape='observations'),color='grey30',alpha=0.5)+
  geom_errorbar(data=obs1,aes(ymin=Obs-(SE*mean(Obs)),ymax=Obs+(SE*mean(Obs)),x=Yr,group=interaction(variable,Label)),color='grey30',alpha=0.5)+
  geom_line(data=obs1,aes(y=Exp,x=Yr,color=variable,group=interaction(variable,Label)),linewidth=0.7)+
  #geom_ribbon(data=obs1,aes(ymin=Explower,ymax=Expupper,x=Yr,fill=variable,group=interaction(variable,Label)),alpha=0.2)+
  #geom_point(data=obs,aes(x=Yr,y=Obs))+
  theme_minimal()+
  labs(y='biomass')+
  scale_color_manual(values = c('fixed'='black',
                                'env_add'='#01AD2C','env_mul'='#6FBD84',
                                'blk_all'='#B30101','blk3_all'='#E66000','blk_sev'='#C25460','blk3_sev'='#E09353',
                                'pred_all'='#BC9912','pred_sev'='#DECE8B',
                                'bycatch_all'='#154ABD','bycatch_sev'='#6DA4C2',
                                'bycatchF_all'='#7A29E5','bycatchF_sev'='#8A5FC1'),name='Assessment\nmodel')+  
  #theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1),axis.title.x = element_blank())+
  scale_shape_discrete(name='')+
  facet_wrap(~Label,nrow = 2)

#####################
# AgeComp
#####################

unique(Mycomparisonsummary$quants$Label)
#View(data.frame(unique(Mycomparisonsummary$pars$Label)))
Mycomparisonsummary$pars

#####################
# Natural mortality
#####################

unique(Mycomparisonsummary$quants$Label)
#View(data.frame(unique(Mycomparisonsummary$pars$Label)))

natM<-c('NatM','M2',paste0('F_fleet_',4:36))

natM1<-paste(natM, collapse = "|")

natM_df<-reshape2::melt(Mycomparisonsummary$pars[grep(natM1,Mycomparisonsummary$pars$Label),],id.vars=c('Yr','Label'))

#######
#fixed
#######
natM_fix<-subset(natM_df,variable=='fixed')[!is.na(subset(natM_df,variable=='fixed')["value"]),]

natM_fix1<-data.frame('Yr'=rep(1985:2020,each=32),
                      'age'=rep(1:32,times=length(1985:2020)),
                      'RTM'=0,
                      'NatM'=rep(natM_fix$value,times=length(1985:2020)),
                      'model'='fixed')

######
#env
######
envdata <- read.delim("./outputs/mrt_as_env.txt",header = FALSE)
names(envdata)<-c("Yr",'age','zmrt')
df<-subset(natM_df,variable %in% c('env.add','env.mul'))[!is.na(subset(natM_df,variable %in% c('env.add','env.mul'))["value"]),]

natM_env<-data.frame(matrix(NA,nrow = 0,ncol = 5))
names(natM_env)<-c('Yr','age','RTM','NatM','model')

#ADDITIVE
for (m in unique(df$variable)[1]) {
  
  #m<-unique(df$variable)[1]
  
  df1<-subset(df,variable==m)
  df1<-df1[grep('ENV_add',df1$Label),]
  
  for (age in 1:32) {
    
    #age<-1
    
    sc<-df1[age,'value']
    
    for (y in 1985:2020) {
      
      #y<-1985
      
      age1<-ifelse(age %in% 1:6 ,age,6)
      
      RTM<-sc*envdata[which(envdata$Yr==y & envdata$age==age1),'zmrt']
      
      totalM<-RTM+natM_fix[age,'value']
      
      
      inatM<-data.frame('Yr'=y,
                        'age'=age,
                        'RTM'=RTM,
                        'NatM'=totalM,
                        'model'=m)
      
      natM_env<-rbind(natM_env,inatM)
      
    }
  }
}


ggplot()+
  geom_point(data=natM_env,aes(x=Yr,y=NatM,color=NatM))+
  scale_y_continuous(limits = c(0,NA))+
  theme_bw()+
  facet_wrap(~age)

#MULTIPLICATIVE
for (m in unique(df$variable)[2]) {
  
  #m<-unique(df$variable)[2]
  
  df1<-subset(df,variable==m)
  df1<-df1[grep('ENV_mult',df1$Label),]
  
  for (age in 1:32) {
    
    #age<-1
    
    sc<-df1[age,'value']
    
    for (y in 1985:2020) {
      
      #y<-1985
      
      age1<-ifelse(age %in% 1:6 ,age,6)
      
      RTM<-exp(sc*envdata[which(envdata$Yr==y & envdata$age==age1),'zmrt'])
      
      totalM<-RTM*natM_fix[age,'value']
      
      RTM<-totalM-natM_fix[age,'value']
      inatM<-data.frame('Yr'=y,
                        'age'=age,
                        'RTM'=RTM,
                        'NatM'=totalM,
                        'model'=m)
      
      natM_env<-rbind(natM_env,inatM)
      
    }
  }
}


ggplot()+
  geom_point(data=natM_env,aes(x=Yr,y=RTM,color=NatM))+
  scale_y_continuous(limits = c(0,NA))+
  scale_color_viridis_b()+
  theme_bw()+
  facet_wrap(~age + model)

 ######
#blk
######

blk_mods<-c('blk.add.all','blk.add3.all','blk.add.sev','blk.add3.sev')

df<-subset(natM_df,variable %in% blk_mods)[!is.na(subset(natM_df,variable %in% blk_mods)["value"]),]

natM_blk<-data.frame(matrix(NA,nrow = 0,ncol = 5))
names(natM_blk)<-c('Yr','age','RTM','NatM','model')

#ADDITIVE
for (m in unique(df$variable)) {
  
  #m<-unique(df$variable)[1]
  
  df1<-subset(df,variable==m)
  df1<-df1[grep('BLK',df1$Label),]
  
  if (m %in% c('blk.add.all','blk.add3.all')) {
    yrs<-c(1985:2020)
  } else {
    yrs<-c(2005,2006,2018,2019)
  }
  
  
  df1$age<-rep(1:32,each=length(yrs))
  
  for (age in 1:32) {
    
    #age<-1
  
    
    for (y in 1985:2020) {
      
      #y<-1985
      sc<-df1[which(df1$age==age & df1$Yr==y),'value']
      if (length(sc)==0) { sc<-0}
      totalM<-sc+natM_fix[age,'value']
      
      
      inatM<-data.frame('Yr'=y,
                        'age'=age,
                        'RTM'=sc,
                        'NatM'=totalM,
                        'model'=m)
      
      natM_blk<-rbind(natM_blk,inatM)
      
    }
  }
}


ggplot()+
  geom_point(data=natM_blk,aes(x=Yr,y=NatM,color=model))+
  scale_y_continuous(limits = c(0,NA))+
  theme_bw()+
  facet_wrap(~age)

######
#bycatch
######

byc_mods<-c('bycatch.all','bycatch.sev',
            'bycatchF.all','bycatchF.sev')

df<-subset(natM_df,variable %in% byc_mods)[!is.na(subset(natM_df,variable %in% byc_mods)["value"]),]

natM_byc<-data.frame(matrix(NA,nrow = 0,ncol = 5))
names(natM_byc)<-c('Yr','age','RTM','NatM','model')

#ADDITIVE
for (m in unique(df$variable)) {
  
  #m<-unique(df$variable)[3]
  
  df1<-subset(df,variable==m)
  df1<-df1[grep('F_fleet',df1$Label),]
  df1$age<-as.numeric(substr(df1$Label,9,9))-3
  
  if (m %in% c('bycatch.all','bycatchF.all')) {
    yrs<-c(1985:2020)
  } else {
    yrs<-c(2005,2006,2018,2019)
  }
  
  for (age in 1:32) {
    
    #age<-1
    
    if (m %in% c('bycatchF.all','bycatchF.sev')) {
      
      age1<-ifelse(age %in% c(1:6),age,6)
      df2<-df1[which(df1$age==age1),]
    } else {
      df2<-df1
     } 
    
    
    for (y in 1985:2020) {
      
      
      #y<-1985
      if (y %in% unique(df2$Yr)) {
        sc<-df2[which(df2$Yr==y),'value']
      } else {sc<-0}
      
      totalM<-sc+natM_fix[age,'value']
      
      
      inatM<-data.frame('Yr'=y,
                        'age'=age,
                        'RTM'=sc,
                        'NatM'=totalM,
                        'model'=m)
      
      natM_byc<-rbind(natM_byc,inatM)
      
    }
  }
}


ggplot()+
  geom_point(data=natM_byc,aes(x=Yr,y=RTM,color=model))+
  scale_y_continuous(limits = c(0,NA))+
  theme_bw()+
  facet_wrap(~age)


######
#predator
######

pred_mods<-c('pred.all','pred.sev')

df<-subset(natM_df,variable %in% pred_mods)[!is.na(subset(natM_df,variable %in% pred_mods)["value"]),]

natM_pred<-data.frame(matrix(NA,nrow = 0,ncol = 5))
names(natM_pred)<-c('Yr','age','RTM','NatM','model')

#ADDITIVE
for (m in unique(df$variable)) {
  
  #m<-unique(df$variable)[1]
  
  df1<-subset(df,variable==m)
  df1<-df1[grep('M2_pred1',df1$Label),]
  df1<-subset(df1, Yr %in% 1985:2020)
  
  if (m %in% c('pred.all')) {
    yrs<-c(1985:2020)
  } else {
    yrs<-c(2005,2006,2018,2019)
  }
  
  for (age in 1:32) {
    
    #age<-1
    
    for (y in 1985:2020) {
      
      
      #y<-1985
      if (y %in% unique(df1$Yr)) {
        sc<-df1[which(df1$Yr==y),'value']
      } else {sc<-0}
      
      totalM<-sc+natM_fix[age,'value']
      
      
      inatM<-data.frame('Yr'=y,
                        'age'=age,
                        'RTM'=sc,
                        'NatM'=totalM,
                        'model'=m)
      
      natM_pred<-rbind(natM_pred,inatM)
      
    }
  }
}


ggplot()+
  geom_point(data=natM_pred,aes(x=Yr,y=RTM,color=model))+
  scale_y_continuous(limits = c(0,NA))+
  theme_bw()+
  facet_wrap(~age)

natM_all<-
rbind(natM_fix1,
      natM_env,
      natM_blk,
      natM_byc,
      natM_pred)


#sort levels
natM_all$model<-factor(natM_all$model,levels = c('fixed',
                                                     'env.add','env.mul',
                                                     'blk.add.all','blk.add.sev','blk.add3.all','blk.add3.sev',
                                                     'pred.all','pred.sev',
                                                     'bycatch.all','bycatch.sev',
                                                     'bycatchF.all','bycatchF.sev'))

#plot
ggplot()+
  geom_line(data=natM_all,aes(x=Yr,y=NatM,color=model),alpha=0.6)+
  scale_y_continuous(limits = c(0,NA))+
  theme_bw()+
  labs(x='',y='M (M0+RTM)')+
  scale_color_manual(values = c('fixed'='black',
                                'env.add'='#01AD2C','env.mul'='#6FBD84',
                                'blk.add.all'='#B30101','blk.add3.all'='#E66000','blk.add.sev'='#C25460','blk.add3.sev'='#E09353',
                                'pred.all'='#BC9912','pred.sev'='#DECE8B',
                                'bycatch.all'='#154ABD','bycatch.sev'='#6DA4C2',
                                'bycatchF.all'='#7A29E5','bycatchF.sev'='#8A5FC1'),name='Assessment\nmodel')+
  facet_wrap(~age,ncol = 10)

#plot
ggplot()+
  geom_line(data=natM_all,aes(x=Yr,y=RTM,color=model),alpha=0.6)+
  scale_y_continuous(limits = c(0,NA))+
  theme_bw()+
  labs(x='',y='M (M0+RTM)')+
  scale_color_manual(values = c('fixed'='black',
                                'env.add'='#01AD2C','env.mul'='#6FBD84',
                                'blk.add.all'='#B30101','blk.add3.all'='#E66000','blk.add.sev'='#C25460','blk.add3.sev'='#E09353',
                                'pred.all'='#BC9912','pred.sev'='#DECE8B',
                                'bycatch.all'='#154ABD','bycatch.sev'='#6DA4C2',
                                'bycatchF.all'='#7A29E5','bycatchF.sev'='#8A5FC1'),name='Assessment\nmodel')+
  facet_wrap(~age,ncol = 10)


#add obs
natM_obs<-Mycomparisonsummary$par_prior_likes[,c("Label","blk.add.all")][46:nrow(Mycomparisonsummary$par_prior_likes[,c("Label","blk.add.all")]),]
natM_obs<-natM_obs[grepl("Fem_GP_1_BLK", natM_obs$Label), ]
natM_obs$age<-gsub('NatM_break_',"",natM_obs$Label)
natM_obs$age<-gsub('Fem_GP_1_BLK1add_',"",natM_obs$age)
natM_obs$Yr<-sub(".*?_(.*)", "\\1", natM_obs$age)
natM_obs$age<-sub("_.*", "", natM_obs$age)
natM_obs$Yr<-as.numeric(natM_obs$Yr)
summary(natM_obs)
natM_obs[is.na(natM_obs)] <- 0
natM_obs$Yr<-as.numeric(natM_obs$Yr)

natM_fix$age<-1:32

natM_obs1<-
merge(natM_obs,
      natM_fix,
      by='age')

natM_obs1$NatM<-natM_obs1$blk.add.all+natM_obs1$value

natM_all$age<-factor(natM_all$age,1:32)
natM_obs1$age<-factor(natM_obs1$age,1:32)
natM_obs1<-natM_obs1[,c("age","blk.add.all","Yr.x","NatM")]
names(natM_obs1)<-c('age','obs','Yr','M')

#write csv mrt obs
obs<-read.csv2('./tables/mrt_annual.csv')
yr<-obs$yr
rmt<-obs$value
obs$age<-as.numeric(substr(obs$age,nchar(obs$age)-1,nchar(obs$age)))+1
obs[is.na(obs$age),]<-'6'
obs$yr<-yr
obs$value<-rmt
obs$age<-obs$age-1

# Define a named vector for facet labels
age_labels <- c('1' = 'gag 0', '2' = 'gag 1', '3' = 'gag 2',
                '4' = 'gag 3', '5' = 'gag 4', '6' = 'gag 5 plus')


p<-
ggplot() +
   geom_line(data = subset(natM_all, age %in% 1:6),
            aes(x = Yr, y = RTM, color = model), linewidth = 0.5) +
  geom_point(data = subset(obs, age %in% 1:6),
             aes(x = yr, y = value, shape = "observation"), stroke = 0.7,size=1) +
   scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = 'white')) +
  labs(x = '', y = expression(MRT~(year^{-1})))+
  scale_color_manual(values = c(
    'fixed' = 'black',
    'env.add' = '#01AD2C', 'env.mul' = '#6FBD84',
    'blk.add.all' = '#B30101', 'blk.add3.all' = '#E66000', 'blk.add.sev' = '#C25460', 'blk.add3.sev' = '#E09353',
    'pred.all' = '#BC9912', 'pred.sev' = '#DECE8B',
    'bycatch.all' = '#154ABD', 'bycatch.sev' = '#6DA4C2',
    'bycatchF.all' = '#7A29E5', 'bycatchF.sev' = '#8A5FC1'
  ), name = 'Assessment\nmodel',
  labels = c('fixed', 
             'env_add', 'env_mul', 
             'blk_all',  'blk_sev',
             'blk3_all', 'blk3_sev', 
             'pred_all', 'pred_sev', 
             'bycatch_all', 'bycatch_sev', 
             'bycatchF_all', 'bycatchF_sev'),
  guide = guide_legend(override.aes = list(linewidth = 1))) +  # Increase legend line width
  guides(
    shape = guide_legend(override.aes = list(size = 3))  # Increase size of shape legend symbol
  ) +
  scale_shape_manual(values = c("observation" = 1), name = " ") +
  facet_wrap(~age, ncol = 2, labeller = labeller(age = age_labels))


ragg::agg_png('./figures/MRT.png',  width = 7, height = 5, units = "in", res = 300)
print(
  p
)
dev.off()


df<-Mycomparisonsummary$indices
df$obs_est2<-(df$Obs-df$Exp)^2

out<-data.frame(matrix(NA,nrow=0,ncol = 4))
names(out)<-c('imodel','Fleet','RMSE','RRMSE')

for (m in unique(df$name)) {
  for (sur in c('SURVEYJUVENILE','SURVEYADULT')) {
    
    #sur<-unique(df$Fleet_name)[1];m<-unique(df$name)[1]
    
    df1<-subset(df,Fleet_name==sur & name==m)
    
    iout<-data.frame('imodel'=m,
                     'Fleet'=sur,
                     'RMSE'=sqrt(sum(df1$obs_est2)/nrow(df1)),
                     'RRMSE'=sqrt(sum(df1$obs_est2)/nrow(df1))/mean(df1$Exp))
    
    out<-rbind(out,iout)
    
  }
}








natM_all1<-merge(natM_all,natM_obs1,by=c('age','Yr'))
natM_all1$sqdiff<-(natM_all1$obs-natM_all1$RTM)^2

out<-data.frame(matrix(NA,nrow = 0,ncol = 4))
names(out)<-c('imodel','age','RMSE','RRMSE')

for (m in unique(natM_all1$model)) {
  for (a in unique(natM_all1$age)) {
    
  df2<-subset(natM_all1,model==m & age==a)
  
  #for (y in 1985:2020) {
      
      #df1<-subset(natM_all1,model==m & Yr==y)
      
      iout<-data.frame('imodel'=m,
                       'age'=a,
                       'RMSE'=sqrt(sum(df2$sqdiff))/nrow(df2),
                       'RRMSE'=sqrt(sum(df2$sqdiff))/nrow(df2)/mean(df2$RTM))
      
      out<-rbind(out,iout)
  }
}


# Set the levels first
out$imodel <- factor(out$imodel, levels = c('fixed', 
                                                'env.add', 'env.mul', 
                                                'blk.add.all','blk.add.sev',
                                                'blk.add3.all', 'blk.add3.sev', 
                                                'pred.all', 'pred.sev', 
                                                'bycatch.all', 'bycatch.sev', 
                                                'bycatchF.all', 'bycatchF.sev'))

# Assign new names to the factor levels
levels(out$imodel) <- c('fixed', 
                          'env_add', 'env_mul', 
                          'blk_all',  'blk_sev',
                          'blk3_all', 'blk3_sev', 
                          'pred_all', 'pred_sev', 
                          'bycatch_all', 'bycatch_sev', 
                          'bycatchF_all', 'bycatchF_sev')


#inf to NA
out[sapply(out, is.infinite)] <- NA#plot

p<-
ggplot() +
  geom_point(data = subset(out, age %in% 1:6 & imodel != 'fixed'), aes(x = imodel, y = RRMSE, fill = imodel),color='black', size=4,shape=21) +
  geom_vline(xintercept = c(2.5, 6.5, 8.5,10.5), linetype = "dashed", color = "gray50") + # Add vertical lines
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(axis.text.x = element_blank(),#element_text(angle=45,hjust=1,vjust=1),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill='white'))+
  labs(x = '', y = 'RRMSE of the MRT') +
  scale_fill_manual(
    values = c(
      #'fixed' = 'black',
      'env_add' = '#01AD2C', 'env_mul' = '#6FBD84',
      'blk_all' = '#B30101', 'blk3_all' = '#E66000', 'blk_sev' = '#C25460', 'blk3_sev' = '#E09353',
      'pred_all' = '#BC9912', 'pred_sev' = '#DECE8B',
      'bycatch_all' = '#154ABD', 'bycatch_sev' = '#6DA4C2',
      'bycatchF_all' = '#7A29E5', 'bycatchF_sev' = '#8A5FC1'
    ),
    name = 'Assessment\nmodel'
  )+
  facet_wrap(~age, ncol = 2, labeller = labeller(age = age_labels))


ragg::agg_png('./figures/MRT RRMSE.png',  width = 7, height = 5, units = "in", res = 300)
print(
  p
)
dev.off()

p<-
  ggplot() +
  geom_boxplot(data = subset(out, imodel != 'fixed'), aes(x = imodel, y = RRMSE, fill = imodel),outliers = FALSE,color='black',) +
  geom_vline(xintercept = c(2.5, 4.5, 8.5,10.5), linetype = "dashed", color = "gray50") + # Add vertical lines
  scale_y_continuous(limits = c(0, NA),expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1),axis.title.x = element_blank(),strip.background = element_rect(fill='white'))+
  labs(x = '', y = 'RRMSE of the MRT') +
  scale_fill_manual(
    values = c(
      'fixed' = 'black',
      'env_add' = '#01AD2C', 'env_mul' = '#6FBD84',
      'blk_all' = '#B30101', 'blk3_all' = '#E66000', 'blk_sev' = '#C25460', 'blk3_sev' = '#E09353',
      'pred_all' = '#BC9912', 'pred_sev' = '#DECE8B',
      'bycatch_all' = '#154ABD', 'bycatch_sev' = '#6DA4C2',
      'bycatchF_all' = '#7A29E5', 'bycatchF_sev' = '#8A5FC1'
    ),
    name = 'Assessment\nmodel'
  )


ragg::agg_png('./figures/MRT RRMSE distribution.png',  width = 6, height = 5, units = "in", res = 300)
print(
  p
)
dev.off()

#####################
# Retrospective (bias) analysis
#####################
install_github("jabbamodel/ss3diags")

library(ss3diags)



mase_df<-data.frame(matrix(NA,ncol=4,nrow = 0))
names(mase_df)<-c('Index','MASE','model','bias_rhoSSB')

for (m in models) {
  
  #m<-models[1]
  
  dir.create(paste0('./models/',m,'/retrosp/'))
  
  # run the retrospective analyses
  retro(
    dir = paste0('./models/',m), # wherever the model files are
    #oldsubdir = "", # subfolder within dir
    newsubdir = "retrosp", # new place to store retro runs within dir
    years = 0:-5, # years relative to ending year of model
    exe = "ss3",overwrite = TRUE
  )
  
  #r4ss::retro(dir = dir_retro, exe = "ss3", years = 0:-5, verbose = FALSE)
  
  retro_mods <- r4ss::SSgetoutput(dirvec = file.path(paste0('./models/',m,'/retrosp/'), 
                                                     paste0("retro", seq(0, -5, by = -1))), 
                                  verbose = F)
  
  retroSummary <- r4ss::SSsummarize(retro_mods, verbose = F)
  
  #evaluate.bias
  rho_output <- SSmohnsrho(
    summaryoutput = retroSummary,
    endyrvec = retroSummary$endyrs,
    startyr = retroSummary[["endyrs"]] - 5,
    verbose = FALSE
  )
  
  
  r4ss::sspar(mfrow = c(1, 2))
  mase<-ss3diags::SSplotHCxval(retroSummary, subplots = "cpue", add = TRUE)
  
  mase_df<-
  rbind(mase_df,
        data.frame('Index'=mase$Index,
                   'MASE'=mase$MASE, #cross-validation
                   'bias_rhoSSB'=rho_output$SSB, 
                   'model'=m))
  
}


saveRDS(object = mase_df,file = './outputs/mase_rho.RData')
mase_df<-readRDS(file = './outputs/mase_rho.RData')

#MASE 
#Mean Absolute Scaled Error (MASE) in a retrospective analysis is a metric used to evaluate the accuracy of a model by comparing its error to a baseline model, often a naive forecast. Here's what it typically tells you:
#Relative Forecast Accuracy: MASE measures forecast error relative to a naive or baseline forecast (often the "last observed value" forecast). A MASE value below 1 suggests that the model performs better than the naive approach, while a value above 1 indicates worse performance.
#Model Reliability: In retrospective analyses, particularly for time series, MASE can help determine how consistently a model would have performed if it had been used in past scenarios. Low MASE values suggest that the model reliably captures patterns in historical data.
#Robust Error Assessment: MASE is scale-independent, so it allows for comparisons across different time series without the errors being influenced by the scale of the data. This makes it a useful metric in retrospective studies with multiple series or different measurement units.
#A MASE score > 1 indicates that the average model forecasts are worse than a random walk. Conversely, a MASE score of 0.5 indicates that the model forecasts twice as accurately as a naïve baseline prediction; thus, the model has prediction skill.



#suggests values of 
#that fall outside (-0.15 to 0.20) for SSB for longer-lived species, or outside (-0.22 to 0.30) for shorter-lived species indicates an undesirable retrospective pattern. In addition, the direction of the retrospective bias has implications for characterizing risk associated with management advice. A positive 
#for SSB is of particular concern because it implies a systemic overestimation of biomass, which would lead to over-optimistic quota advice if not taken into consideration

mase_df1<-subset(mase_df,Index %in% c('SURVEYJUVENILE', 'SURVEYADULT'))

# Set the levels first
mase_df1$model <- factor(mase_df1$model, levels = c('fixed', 
                                            'env.add', 'env.mul', 
                                            'blk.add.all','blk.add.sev',
                                            'blk.add3.all', 'blk.add3.sev', 
                                            'pred.all', 'pred.sev', 
                                            'bycatch.all', 'bycatch.sev', 
                                            'bycatchF.all', 'bycatchF.sev'))

# Assign new names to the factor levels
levels(mase_df1$model) <- c('fixed', 
                        'env_add', 'env_mul', 
                        'blk_all',  'blk_sev',
                        'blk3_all', 'blk3_sev', 
                        'pred_all', 'pred_sev', 
                        'bycatch_all', 'bycatch_sev', 
                        'bycatchF_all', 'bycatchF_sev')

index_labels <- c(
  "SURVEYADULT" = "offshore survey",
  "SURVEYJUVENILE" = "nearshore survey")

p<-
ggplot()+
  geom_point(data=mase_df1,aes(x=model,y=MASE,fill=model),size=3,shape=21)+
  theme_bw()+
  facet_wrap(~Index,ncol=1, labeller = labeller(Index = index_labels))+
  geom_hline(yintercept = 1,linetype='dashed')+
theme_bw() +
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1),axis.title.x = element_blank(),strip.background = element_rect(fill='white'))+
  #labs(x = '', y = 'RRMSE (MRT)') +
  scale_fill_manual(
    values = c(
      'fixed' = 'black',
      'env_add' = '#01AD2C', 'env_mul' = '#6FBD84',
      'blk_all' = '#B30101', 'blk3_all' = '#E66000', 'blk_sev' = '#C25460', 'blk3_sev' = '#E09353',
      'pred_all' = '#BC9912', 'pred_sev' = '#DECE8B',
      'bycatch_all' = '#154ABD', 'bycatch_sev' = '#6DA4C2',
      'bycatchF_all' = '#7A29E5', 'bycatchF_sev' = '#8A5FC1'
    ),
    name = 'Assessment\nmodel'
  )

ragg::agg_png('./figures/MASE.png',  width = 6, height = 5, units = "in", res = 300)
print(
  p
)
dev.off()

p<-
ggplot()+
  geom_point(data=mase_df1,aes(x=model,y=bias_rhoSSB,fill=model),size=3,shape=21)+
  theme_bw()+
  #facet_wrap(~Index,ncol=1, labeller = labeller(Index = index_labels))+
  geom_hline(yintercept = -0.15,linetype='dashed')+
  geom_hline(yintercept = 0.30,linetype='dashed')+
  theme_bw() +
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1),axis.title.x = element_blank(),strip.background = element_rect(fill='white'))+
  #labs(x = '', y = 'RRMSE (MRT)') +
  scale_fill_manual(
    values = c(
      'fixed' = 'black',
      'env_add' = '#01AD2C', 'env_mul' = '#6FBD84',
      'blk_all' = '#B30101', 'blk3_all' = '#E66000', 'blk_sev' = '#C25460', 'blk3_sev' = '#E09353',
      'pred_all' = '#BC9912', 'pred_sev' = '#DECE8B',
      'bycatch_all' = '#154ABD', 'bycatch_sev' = '#6DA4C2',
      'bycatchF_all' = '#7A29E5', 'bycatchF_sev' = '#8A5FC1'
    ),
    name = 'Assessment\nmodel'
  )

ragg::agg_png('./figures/bias_rhoSSB.png',  width = 6, height = 5, units = "in", res = 300)
print(
  p
)
dev.off()

#####################
# correlation MRT and recruitment and SSB
#####################


#scatterplot X MRT vs SSB vs Recruitment (color intensity by MRT)

natM_all1<-natM_all[,-4]


Mycomparisonsummary$quants[grepl('SSB',Mycomparisonsummary$quants$Label),]
ssb<-melt(Mycomparisonsummary$quants[which(Mycomparisonsummary$quants$Label %in% paste0("SSB_",1985:2020)),])
recr<-melt(Mycomparisonsummary$quants[which(Mycomparisonsummary$quants$Label %in% paste0("Recr_",1985:2020)),])
names(recr)[2]<-'model'
names(ssb)[2]<-'model'
ssb$Yr<-sub(".*_", "", ssb$Label)
recr$Yr<-sub(".*_", "", recr$Label)

recr$ref<-sub("_.*", "", recr$Label)
ssb$ref<-sub("_.*", "", ssb$Label)


natM_all2<-merge(natM_all1,recr,by=c('Yr','model'))
natM_all3<-merge(natM_all2,ssb,by=c('Yr','model'))
natM_all3<-natM_all3[,c(1:4,6,9)]
names(natM_all3)[c(5,6)]<-c("recr",'ssb')


ggplot()+
  geom_point(data=natM_all3,aes(x=ssb,y=recr,color=RTM))+
  #scale_y_continuous(limits = c(0,NA))+
  theme_bw()+
  facet_wrap(~model,scales='free')+
  scale_color_continuous(limits = c(0, 0.5),              # Set color scale limits
                        oob = scales::oob_squish)
  #labs(x='',y='M (M0+RTM)')+
  # scale_color_manual(values = c('fixed'='black',
  #                               'env.add'='#01AD2C','env.mul'='#6FBD84',
  #                               'blk.add.all'='#B30101','blk.add3.all'='#E66000','blk.add.sev'='#C25460','blk.add3.sev'='#E09353',
  #                               'pred.all'='#BC9912','pred.sev'='#DECE8B',
  #                               'bycatch.all'='#154ABD','bycatch.sev'='#6DA4C2',
  #                               'bycatchF.all'='#7A29E5','bycatchF.sev'='#8A5FC1'),name='Assessment\nmodel')






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
