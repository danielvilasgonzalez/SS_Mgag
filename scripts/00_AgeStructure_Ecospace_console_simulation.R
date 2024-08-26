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
library(ggthemes)
library(scales)
library(colorRamps)
#library(r4ss)
library(reshape2)
library(raster)

#set working directory
#wd<-"C:/WFS EwE/age structure ecospace/"
wd<-"E:/WFS_DV/manuscripts/stocksynthesis/"
#setwd(wd)

##########################################################
#
# RUN ECOSPACE CONSOLE
# 
##########################################################

#number of years to simulate
nyrs<-36

#read commandfile
commandfile<-read.delim(paste0(wd,'CommandFile.txt'))
  
# modify commandfile txt file
commandfile[43,]<-paste0("<EWE_MODEL_FILE>, ",wd,"WFS EwE 2.09 chapter3/WFS EwE 2.09 chapter3.eweaccdb, System.String")
commandfile[44,]<-paste0("<SPATIAL_CONFIG_FILE>, ",wd,"WFS EwE 2.09 chapter3/chapter3_PhD.xml, System.String")
commandfile[46,]<-paste0("<N_ECOSPACE_YEARS>, ",nyrs,", System.Int32, Run length in years")
commandfile[58,]<-paste0('<ECOSPACE_USE_ANNUAL_OUTPUT>, False, System.Boolean, Output monthly maps and spatially averaged .csv files')
commandfile[59,]<-paste0("<ECOSPACE_OUTPUT_DIR>, ",wd,"outputs/, System.String")
commandfile[61,]<-paste0('<ECOSPACE_SAVE_MAP_BIOMASS>, False, System.Boolean')
commandfile[62,]<-paste0('<MO_SAVE_MAP_BIOMASS>, False, System.Boolean')
commandfile[65,]<-paste0("<ECOSPACE_ENVIRONMENTAL_RESPONSE_STRING>, , System.String")
commandfile[67,]<-paste0("<ECOSPACE_SPINUP_LENGTH>, 2, System.Int32")
commandfile[69,]<-paste0('<ECOINDICATORS_SAVE_ECOSPACE>, False, System.Boolean, True will save all the indicators to file. There is no fine grained manipulation of variable')
commandfile[74,]<-paste0("//<ECOSPACE_DISPERSAL_RATE>,  300 300 300 300 300 300 300 50 50 500 100 10 10 300 300 1, System.Single[]")
commandfile[80,]<-paste0("//<ECOSIM_VULNERABILITIES_BY_PREY_PRED>, , System.Single[]")
commandfile[86,]<-paste0('<IBM_SAVE_AGE_STRUCTURE>, True, System.Boolean, True to save the IBM age structure file')

#write new commandfile
write.table(commandfile,file = paste0(wd,"/outputs/CommandFile.txt"),quote = FALSE,row.names = FALSE)

#line including commandfile and ecospace console to run through command prompt
cmd<-paste(paste0('"','C:/Program Files/EcoSpace Console 1.3.8.1 64-bit/EwEClientConsole.exe','"'),
             paste0('"',wd,'outputs/CommandFile.txt','"'))

#starting time to check how much time took to run
start_time<-Sys.time()

#run command prompt
sys::exec_internal(cmd)

#ending time to check how much time took to run
end_time<-Sys.time()

##########################################################
#
# INVESTIGATE RESULTS (EXTRA TESTING)
# 
##########################################################

#filter excel file for gag
lf<-list.files(paste0(wd,'/outputs'))
lf.gag<-lf[grep("gag",lf)]

##########################################################
# Annual average
##########################################################

lf.annual<-lf[grep("Annual_Average_Biomass",lf)]
filen<-read.csv(paste(wd,'outputs',lf.annual,sep = '/'),skip = 31)
file.gag<-cbind(year=filen$Year,filen[,grep('gag',colnames(filen))])
file.gag$year<-c(1985:(1985+nyrs-1))
colnames(file.gag)<-c('year','0','1','2','3','4','5')
average_gag<-reshape2::melt(file.gag,id.vars='year')
dim(average_gag)
colnames(average_gag)<-c('year','age','x')
average_gag$file<-'annual_avg'

##########################################################
# All region (as file)
##########################################################

#open number and weight gag files
filen<-read.csv(paste(wd,'outputs',lf.gag[grep("Region_All_Number",lf.gag)],sep = '/'),skip = 27)
filew<-read.csv(paste(wd,'outputs',lf.gag[grep("Region_All_Weight",lf.gag)],sep = '/'),skip = 27)

#Joe Bukwoski (age structure developer) comments
#The IBM Age Structure files are weight and number at age by month averaged by the Ecospace regions. 
#Each file is for one region and multi-stanza group. The file naming is “AgeStructure _ [multistanza group name] _Region_[x]_Weight or Number.csv” i.e. “AgeStructure_gag_Region_1_Number.csv”. 
#The total area is in the “Region_All” file.   The first row of data is the baseline Ecosim weight and number. The units are the Ecopath t/km2 or number/km2. The Ecopath biomass is the sum of weight * number over the ages for a group.
#For this round I didn’t sum the monthly ages into an annual because we can’t be sure how someone else would use the data. I think it would be easy to post-process the data for plotting into any kind of monthly grouping you like.

#multiply weight by number
file<-filen[,3:ncol(filen)]*filew[,3:ncol(filew)]
#file<-filen[,3:ncol(filen)]
file$Timestep<-filen$Timestep

#reshape and prepare dataframe
df<-reshape2::melt(file,id='Timestep')
df<-df[order(df$Timestep),]

#add month
df$month<-1:(ncol(file)-1)

#add age
age<-sort(rep(0:(round((ncol(file)-1)/12)),times=12))
age<-age[1:(ncol(file)-1)]
df$age<-age

#remove timestep0 because it is the ecopath baseline
df<-df[which(df$Timestep!=0),]

#add year
df$year<-sort(rep(1985:2020,times=(12*377)))
#length(sort(rep(1:nyrs,times=(12*375))))

#aggregate by year timestep by averaging monthly steps
df1<-aggregate(df$value,by=list(year=df$year,Timestep=df$Timestep,age=df$age),FUN=sum) #FUN=mean - it looks mean because if not those estimates look pretty high

#packets are by months (age) and monthly timesteps 
#lets aggregate monthly numbers or weights by year (age) and monthly timesteps
sm0<-aggregate(df1$x, by=list(year=df1$year,age=df1$age), FUN=sum)
#mn1<-aggregate(df1$x, by=list(year=df1$year,age=df1$age), FUN=mean)

#aggregate by year timestep by averaging monthly steps
sm0<-aggregate(df$value,by=list(year=df$year,age=df$age),FUN=mean)

sm1<-subset(sm0,age<=5)

#plot
ggplot()+
  geom_bar(data=sm1,aes(x=age,y=x),stat='identity',position = position_dodge())+
  facet_wrap(~year,scales="free")

ggplot()+
  geom_bar(data=sm1,aes(x=year,y=x),stat='identity',position = position_dodge())+
  facet_wrap(~age,scales="free")

#comparison annual average and region All
sm1$file<-'age_structure'
x<-rbind(average_gag,sm1[c(1:216),])

ggplot()+
  geom_point(data=x,aes(x=year,y=x,color=file))+
  facet_wrap(~age,scales="free")

ggplot()+
  geom_point(data= average_gag,aes(x=year,y=x,color=file))+
  facet_wrap(~age,scales="free")

ggplot()+
  geom_point(data= sm1,aes(x=year,y=x,color=file))+
  facet_wrap(~age,scales="free")

##########################################################
# region1 - MPA
##########################################################

#open number and weight gag files
filen<-read.csv(paste(wd,'outputs',lf.gag[grep("Region_1_Number",lf.gag)],sep = '/'),skip = 27)
filew<-read.csv(paste(wd,'outputs',lf.gag[grep("Region_1_Weight",lf.gag)],sep = '/'),skip = 27)

#Joe Bukwoski (age structure developer) comments
#The IBM Age Structure files are weight and number at age by month averaged by the Ecospace regions. 
#Each file is for one region and multi-stanza group. The file naming is “AgeStructure _ [multistanza group name] _Region_[x]_Weight or Number.csv” i.e. “AgeStructure_gag_Region_1_Number.csv”. 
#The total area is in the “Region_All” file.   The first row of data is the baseline Ecosim weight and number. The units are the Ecopath t/km2 or number/km2. The Ecopath biomass is the sum of weight * number over the ages for a group.
#For this round I didn’t sum the monthly ages into an annual because we can’t be sure how someone else would use the data. I think it would be easy to post-process the data for plotting into any kind of monthly grouping you like.

#multiply weight by number
file<-filen[,3:ncol(filen)]*filew[,3:ncol(filew)]
#file<-filew[,3:ncol(filew)]
file$Timestep<-filen$Timestep
  
#reshape and prepare dataframe
df<-reshape2::melt(file,id='Timestep')
df<-df[order(df$Timestep),]

#add month
df$month<-1:(ncol(file)-1)

#add age
age<-sort(rep(0:(round((ncol(file)-1)/12)),times=12))
age<-age[1:(ncol(file)-1)]
df$age<-age

#remove timestep0 because it is the ecopath baseline
df<-df[which(df$Timestep!=0),]
  
#add year
df$year<-sort(rep(1985:2020,times=(12*377)))
#length(sort(rep(1:nyrs,times=(12*375))))

#aggregate by year timestep by averaging monthly steps
df1<-aggregate(df$value,by=list(year=df$year,Timestep=df$Timestep,age=df$age),FUN=sum) #FUN=mean - it looks mean because if not those estimates look pretty high

#packets are by months (age) and monthly timesteps 
#lets aggregate monthly numbers or weights by year (age) and monthly timesteps
sm1<-aggregate(df1$x, by=list(year=df1$year,age=df1$age), FUN=sum)
#mn1<-aggregate(df1$x, by=list(year=df1$year,age=df1$age), FUN=mean)
#aggregate by year timestep by averaging monthly steps
sm0<-aggregate(df$value,by=list(year=df$year,age=df$age),FUN=mean)

r1<-subset(sm0,age<=5)

#plot
ggplot()+
  geom_bar(data=r1,aes(x=age,y=x),stat='identity',position = position_dodge())+
  facet_wrap(~year,scales="free")

ggplot()+
  geom_bar(data=r1,aes(x=year,y=x),stat='identity',position = position_dodge())+
  facet_wrap(~age,scales="free")

#maybe there is a need to standardize values using timestep=0

##########################################################
# region2 - NON_MPA
##########################################################
  
#open number and weight gag files
filen<-read.csv(paste(wd,'outputs',lf.gag[grep("Region_2_Number",lf.gag)],sep = '/'),skip = 27)
filew<-read.csv(paste(wd,'outputs',lf.gag[grep("Region_2_Weight",lf.gag)],sep = '/'),skip = 27)

#multiply weight by number
file<-filen[,3:ncol(filen)]*filew[,3:ncol(filew)]
#file<-filew[,3:ncol(filew)]
file$Timestep<-filen$Timestep
  
#reshape and prepare dataframe
df<-reshape2::melt(file,id='Timestep')
df<-df[order(df$Timestep),]

#add month
df$month<-1:(ncol(file)-1)
age<-sort(rep(0:(round((ncol(file)-1)/12)),times=12))

#add age
age<-age[1:(ncol(file)-1)]
df$age<-age  

#remove timestep0 because it is the ecopath baseline
df<-df[which(df$Timestep!=0),]
  
#add year
df$year<-sort(rep(1985:2020,times=(12*377)))
#length(sort(rep(1:nyrs,times=(12*375))))

#aggregate by year timestep by averaging monthly steps
df1<-aggregate(df$value,by=list(year=df$year,Timestep=df$Timestep,age=df$age),FUN=sum) #FUN=mean - it looks mean because if not those estimates look pretty high
#packets are by months (age) and monthly timesteps 

#aggregate monthly numbers or weights by year (age) and monthly timesteps
sm2<-aggregate(df1$x, by=list(year=df1$year,age=df1$age), FUN=sum)
#mn2<-aggregate(df1$x, by=list(year=df1$year,age=df1$age), FUN=mean)
#aggregate by year timestep by averaging monthly steps
sm0<-aggregate(df$value,by=list(year=df$year,age=df$age),FUN=mean)

r2<-subset(sm0,age<=5)
#plot
ggplot()+
  geom_bar(data=r2,aes(x=age,y=x),stat='identity',position = position_dodge())+
  facet_wrap(~year,scales="free")
  
ggplot()+
  geom_bar(data=r2,aes(x=year,y=x),stat='identity',position = position_dodge())+
  facet_wrap(~age,scales="free")

#add region column and join two df
sm1$region<-'mpa'
sm2$region<-'other'
sm<-rbind(sm1,sm2)


ggplot()+
  geom_bar(data=sm,aes(x=age,y=x,fill=region),stat='identity',position = position_dodge())+
  facet_wrap(~year,scales="free")

##########################################################
# regionsum - sum region1 and region2
##########################################################
rsum<-r2
rsum$file<-'regionsum' 
rsum$x<-sum(rsum$x,r1$x) #rowMeans(cbind(rsum$x,r1$x))

r1$file<-'region1'
r2$file<-'region2'


xx<-rbind(x,r1,r2,rsum)

ggplot()+
  geom_point(data=xx,aes(x=age,y=x,color=file))+
  facet_wrap(~year,scales="free")

ggplot()+
  geom_point(data=xx,aes(x=year,y=x,color=file),alpha=0.5)+
  facet_wrap(~age,scales="free")
