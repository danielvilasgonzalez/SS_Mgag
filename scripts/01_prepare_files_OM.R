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
  #library(ggthemes)
  #library(scales)
  library(wesanderson)
  #library(colorRamps)
  #library(r4ss)
  #library(dplyr)
  library(reshape2)
  library(raster)
  library(foreach)
  library(doParallel)

  #set working directory
  wd<-"/Users/daniel/Work/SS time-varying M/"
  setwd(wd)
  
  #calculate df form ascii files?
  df.asci<-'no'

##########################################################
#
# INVESTIGATE AND PREPARE DATA FROM OM
# 
##########################################################

  #########################################################
  # BIOMASS with age structure number files
  #########################################################

  #filter file for gag
  lf<-list.files(paste0(wd,'/outputs/OM'))
  lf.gag<-lf[grep("gag",lf)]
  
  if (df.asci=='yes') {

  #list to store
  list.gag<-list()
  
  #check starting time
  start.time<-Sys.time()
  
  #loop if ascii to speed up un parallel
  list.gag<- 
  foreach(gag = c('gag 0','gag 1','gag 2','gag 3','gag 4','gag 5+')) %dopar% {
  #for (gag in c('gag 0','gag 1','gag 2','gag 3','gag 4','gag 5+')) {
    #gag<-'gag 5+'
    
    lf.gag.age<-lf.gag[grep(gag,lf.gag)]
    lf.gag.age.asc<-lf.gag.age[grep('.asc',lf.gag.age)]
    
    mean.gag<-c()

    for (r in 1:length(lf.gag.age.asc)) {
       #r<-1
       r.gag<-raster::raster(paste0('./outputs/OM/',lf.gag.age.asc[r]))
       print(names(r.gag))
       r.gag.mean<-raster::cellStats(r.gag,mean)
       mean.gag<-c(mean.gag,r.gag.mean)  
    }
    
    list.gag[[gag]]<-mean.gag
    
  }
  
  #check end time
  end.time<-Sys.time()
  
  #how much time to run this
  t2<-end.time-start.time
  
  #list to dataframe
  df <- do.call(cbind, list.gag)
  df<-as.data.frame(df)
  df<-as.data.frame(scale(df))
  colnames(df)<-c(0,1,2,3,4,5)
  #add month time step
  df$Timestep<-c(1:nrow(df))

  #reshape and prepare dataframe
  df<-reshape2::melt(df,id='Timestep')
  }
  
  #df to store results
  df1.IBM.mean<-data.frame(matrix(data=NA,nrow = 0,ncol=4))
  colnames(df1.IBM.mean)<-c('year','age','value','region')
  df1.IBM.sd<-df1.IBM.mean
  df2.IBM<-data.frame(matrix(data=NA,nrow = 0,ncol=7))
  colnames(df2.IBM)<-c('Timestep','variable','value','month.age','age','year','region')
  
  #each region a time series of numbers from age structure file
  for (region in c('1','2','All')) {
    #region='2'
  
    lf.gag.csv<-lf.gag[grep(paste0(region,'_Number.csv'),lf.gag)]
    csv.gag<-read.csv(paste0('./outputs/OM/',lf.gag.csv),skip = 27) 
    csv.gag<-csv.gag[,-2]
    
    #reshape and prepare dataframe
    df.IBM<-reshape2::melt(csv.gag,id='Timestep')
    df.IBM<-df.IBM[order(df.IBM$Timestep),]
    
    #add month
    df.IBM$month.age<-1:(ncol(csv.gag)-1)
    
    #add age
    age<-sort(rep(0:(round((ncol(csv.gag)-1)/12)),times=12))
    age<-age[1:(ncol(csv.gag)-1)]
    df.IBM$age<-age
    
    #remove timestep0 because it is the ecopath baseline
    df.IBM<-df.IBM[which(df.IBM$Timestep!=0),]
    
    #add year
    df.IBM$year<-sort(rep(1985:2020,times=(12*377)))
    #length(sort(rep(1:nyrs,times=(12*375))))
    
    #append
    df1.IBM<-df.IBM
    df1.IBM$region<-region
    df2.IBM<-rbind(df2.IBM,df1.IBM)
    
    #aggregate by year timestep by averaging monthly steps
    IBM.mean<-aggregate(list(value=df.IBM$value),
                            by=list(year=df.IBM$year,
                                    age=df.IBM$age),
                            FUN=mean) 
  
    #standardize index value
    #IBM.mean$value<-scale(IBM.mean$value)
    
    #region
    IBM.mean$region<-ifelse(is.na(as.numeric(region)),
                            3,
                            as.numeric(region))
    
    #aggregate by year timestep by averaging monthly steps
    IBM.sd<-aggregate(list(value=df.IBM$value),
                       by=list(year=df.IBM$year,
                               age=df.IBM$age),
                       FUN=function(x) sd(x) / sqrt(length(x))) 
    
    #region
    IBM.sd$region<-ifelse(is.na(as.numeric(region)),
                          3,
                          as.numeric(region))
    
    #append mean and sd data
    df1.IBM.mean<-rbind(df1.IBM.mean,IBM.mean)
    df1.IBM.sd<-rbind(df1.IBM.sd,IBM.sd)
    
  }
  
  
  #aggregate by age to create index
  df2.IBM<-aggregate(list(value=df1.IBM.mean$value),
                     by=list(year=df1.IBM.mean$year,
                             region=df1.IBM.mean$region),
                     FUN=sum) 
  
  #get sd
  df2.IBM.sd<-aggregate(list(value=df1.IBM.sd$value),
                        by=list(year=df1.IBM.sd$year,
                                region=df1.IBM.sd$region),
                        FUN=mean) 
  
  #biomass plot by region
  ggplot(df2.IBM)+
    geom_line(aes(x=year,y=value,color=region))+
    facet_wrap(~region)  
  
  #remove All region
  df3.IBM<-df2.IBM[which(df2.IBM$region!=3),]
  
  #add obs error
  oerr<-rnorm(length(unique(df3.IBM$value)),0,0.1) #SD=0.2
  
  #_yr month fleet obs stderr
  cpue<-data.frame('yr'=df3.IBM$year,
                   'month'=7,
                   'fleet'=df3.IBM$region+1, #or 1; 2 because 2 fisheries (bycatch and fishery)
                   'obs'=df3.IBM$value*188688.717455821,  #modeled area (km2) 188688.717455821,
                   'srderr'=0.1) #IBM.sd$value
  
  cpue<-rbind(cpue,c(-9999,1,1,1,1))
  
  #save data
  write.table(cpue, 
              file = "./outputs/abundance_indeces_agestructure.txt", sep = "\t",
              row.names = FALSE,
              col.names = FALSE)
  
  #_yr month fleet obs stderr
  cpue<-data.frame('yr'=df3.IBM$value*exp(oerr),
                   'month'=7,
                   'fleet'=df3.IBM$region+1, #or 1; 2 because 2 fisheries (bycatch and fishery)
                   'obs'=df3.IBM$value*188688.717455821,  #modeled area (km2) 188688.717455821,
                   'srderr'=0.1) #IBM.sd$value
  
  cpue<-rbind(cpue,c(-9999,1,1,1,1))
  
  #save data
  write.table(cpue, 
              file = "./outputs/abundance_indeces_agestructure_error.txt", sep = "\t",
              row.names = FALSE,
              col.names = FALSE)
  
  #############################################
  # BIOMASS data from biomass files
  #############################################
  
  #create df to store
  csv.gag3<-data.frame(matrix(data=NA,nrow = 0,ncol=4))
  colnames(csv.gag3)<-c('Timestep','Biomass','year','region')
  
  #each region a time series of biomass
  for (region in c('1','2')) {
    #region='1'
    
    lf.gag.csv<-lf[grep(paste0(region,'_Biomass.csv'),lf)]
    csv.gag<-read.csv(paste0('./outputs/OM/',lf.gag.csv[2]),skip = 31) 
    csv.gag1<-csv.gag[,grep('gag',colnames(csv.gag))]
    csv.gag2<-data.frame(cbind(TimeStep=csv.gag$TimeStep,
                         Biomass=rowSums(csv.gag1)))
    csv.gag2$year<-rep(1985:2020,each=12)
    csv.gag2$region<-region
    csv.gag3<-rbind(csv.gag3,csv.gag2)
  
  }
  
  #add region 
  csv.gag3$region<-as.factor(csv.gag3$region)
  bio.gag<-aggregate(list(Biomass=csv.gag3$Biomass), by = list(region=csv.gag3$region,year=csv.gag3$year), mean)
  bio.gag<-bio.gag[order(bio.gag$region),]
  
  #add obs error
  oerr<-rnorm(nrow(bio.gag),0,0.2) #SD=0.2
  
  #_yr month fleet obs stderr
  cpue<-data.frame('yr'=bio.gag$year,
                   'month'=7,
                   'fleet'=as.numeric(bio.gag$region)+1, #or 1; 2 because 2 fisheries (bycatch and fishery)
                   'obs'=round(bio.gag$Biomass*188688.717455821,digits = 3),  #modeled area (km2) 188688.717455821,
                   'srderr'=0.1) #IBM.sd$value
  
  cpue<-rbind(cpue,c(-9999,1,1,1,1))
  
  #save data
  write.table(cpue, 
              file = "./outputs/abundance_indeces.txt", sep = "\t",
              row.names = FALSE,
              col.names = FALSE)
  
  #_yr month fleet obs stderr
  cpue<-data.frame('yr'=bio.gag$year,
                   'month'=7,
                   'fleet'=as.numeric(bio.gag$region)+1, #or 1; 2 because 2 fisheries (bycatch and fishery)
                   'obs'=round(bio.gag$Biomass*exp(oerr)*188688.717455821,digits=3),  #modeled area (km2) 188688.717455821,
                   'srderr'=0.1) #IBM.sd$value
  
  cpue<-rbind(cpue,c(-9999,1,1,1,1))
  
  #save data
  write.table(cpue, 
              file = "./outputs/abundance_indeces_error.txt", sep = "\t",
              row.names = FALSE,
              col.names = FALSE)
  
  #change levels
  levels(csv.gag3$region)<-c('b','c')
  
  #annotate
  ann_text <- data.frame(year = c(2020,2020),Biomass = c(0.13803667,0.05121769),lab = c('b','c'),
                         region = as.factor(c('b','c')))
  
    #plots total biomass
p1c<-
  ggplot(subset(csv.gag3,region=='c'))+
    geom_boxplot(aes(x=year,y=Biomass*188688.717455821,group=year,fill=region))+
    facet_wrap(~region,scales='free')+
    scale_fill_manual(values = c('b'='#CE906C','c'='#448373'),name='region')+
    theme(aspect.ratio = 1,strip.text.x = element_text(size = 10,color='black'), 
          strip.background = element_rect(color="black", fill="white", linewidth = 1.2, linetype="solid"),
          legend.position = 'right')+ #legend.position=c(1,1)plot.title = element_text(margin = margin(t = 10, b = -20))
    scale_x_continuous(breaks=c(1985,1990,1995,2000,2005,2010,2015,2020),expand = c(0.01,0.01))+ #,position = "top"
    scale_y_continuous(expand=c(0,0.001),limits = c(0,10500),breaks = c(0,5000,10000))+
    theme_bw()+
    xlab('year')+
    ylab('biomass (tones)')+
    theme(aspect.ratio = 1,text = element_text(size=10),axis.title = element_blank(),
          plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines"),legend.background =  element_rect(fill = "transparent", colour = "transparent"),
          legend.key = element_rect(fill = NA, colour = NA),legend.position = 'none',
          panel.border = element_rect(fill = NA, colour = 'black'),panel.grid = element_line(linetype = 'dashed'))+
    theme(strip.background = element_blank(),strip.text = element_blank())#+

p1b<-
  ggplot(subset(csv.gag3,region=='b'))+
  geom_boxplot(aes(x=year,y=Biomass*188688.717455821,group=year,fill=region))+
  scale_fill_manual(values = c('b'='#CE906C','c'='#448373'),name='region')+
  theme(aspect.ratio = 1,strip.text.x = element_text(size = 10,color='black'), 
        strip.background = element_rect(color="black", fill="white", linewidth=1.2, linetype="solid"),
        legend.position = 'right',
        plot.title = element_text(margin = margin(t = 10, b = -20),hjust = 0.9))+ #legend.position=c(1,1)plot.title = element_text(margin = margin(t = 10, b = -20))
  scale_x_continuous(breaks=c(1985,1990,1995,2000,2005,2010,2015,2020),expand = c(0.01,0.01))+ #,position = "top"
  scale_y_continuous(expand=c(0,0.001),limits = c(0,27500),breaks = c(0,10000,20000))+
  theme_bw()+
  xlab('year')+
  ylab('biomass (tones)')+
  theme(aspect.ratio = 1,text = element_text(size=10),axis.title.x = element_blank(), 
        plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines"),legend.background =  element_rect(fill = "transparent", colour = "transparent"),
        legend.key = element_rect(fill = NA, colour = NA),legend.position = 'none',
        panel.border = element_rect(fill = NA, colour = 'black'),panel.grid = element_line(linetype = 'dashed'))+
  theme(strip.background = element_blank(),strip.text = element_blank())#+

p1bc<-
  ggplot(csv.gag3)+
    geom_boxplot(aes(x=year,y=Biomass*188688.717455821,group=interaction(region,year),fill=region,color=region),outlier.colour = 'transparent',size=0.3)+
    #facet_wrap(~region,scales='free')+
    scale_fill_manual(values = c('b'='#CE906C80','c'='#44837380'),name='region')+
    scale_color_manual(values = c('b'='#CE906C','c'='#448373'),name='region')+
    theme(aspect.ratio = 1,strip.text.x = element_text(size = 10,color='black'), 
          strip.background = element_rect(color="black", fill="white", linewidth=1.2, linetype="solid"),
          legend.position = 'right',
          plot.title = element_text(margin = margin(t = 10, b = -20),hjust = 0.9))+ #legend.position=c(1,1)plot.title = element_text(margin = margin(t = 10, b = -20))
    scale_x_continuous(breaks=c(1985,1990,1995,2000,2005,2010,2015,2020),minor_breaks = 1985:2020,
                       expand = c(0.01,0.01))+ #,position = "top"
    scale_y_continuous(expand=c(0,0.001),limits = c(0,27500),breaks = c(0,10000,20000))+
    theme_bw()+
    xlab('year')+
    ylab('biomass (tones)')+
    theme(text = element_text(size=12),axis.title.x = element_blank(), 
          plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines"),legend.background =  element_rect(fill = "transparent", colour = "transparent"),
          legend.key = element_rect(fill = NA, colour = NA),legend.position = 'none',
          panel.border = element_rect(fill = NA, colour = 'black'),panel.grid.minor =  element_line(linetype = 'dashed'))+
    theme(strip.background = element_blank(),strip.text = element_blank())#+

  #########################################################
  # MAPS regions
  #########################################################

#map plot
region<-raster('./processed data/2regions.asc')
region[region==0]<-NA
plot(region)

#color transparent grid
col_grid <- rgb(235, 235, 235, 100, maxColorValue = 255)

#blank polygon
x_coord <- c(-82,  -82, -80, -80)
y_coord <- c(30.5823, 28, 28, 30.5823)
xym <- cbind(x_coord, y_coord)

#spatial arrangements
p = Polygon(xym)
ps = Polygons(list(p),1)
sps = SpatialPolygons(list(ps))
plot(sps)

#florida polygon for crop raster
#us <- readOGR(dsn="E:/GAM surveys/General data/.", layer="cb_2018_us_nation_5m")
land<-rnaturalearth::ne_countries('large')
land_sp <- as(land, "Spatial")  # Convert to Spatial object
bbox = c(latN = 30.5823, latS = 24.91563, lonW = -87.5833, lonE = -80.91663)
FL <- extent(bbox[3],bbox[4], bbox[2],bbox[1])
fl <- crop(land_sp, FL)
plot(fl)

#plot map regions
p2 <- rasterVis::gplot(region) +
  geom_tile(aes(fill = factor(value), color = factor(value))) +
  coord_sf(crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", 
           xlim = c(bbox[3], bbox[4]), ylim = c(bbox[2], bbox[1])) +
  geom_polygon(data = sps, aes(x = long, y = lat), fill = 'white') +
  geom_polygon(data = fl, aes(x = long, y = lat, group = group),
               fill = 'grey80', size = 0.5, color = 'black') +
  theme(aspect.ratio = 1,
        strip.background = element_rect(color = "black", fill = "white", size = 1.2, linetype = "solid")) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-86, -84, -82), position = 'top') +
  scale_y_continuous(expand = c(0, 0), breaks = c(26, 28, 30), position = "right") +
  scale_fill_manual(values = c('1' = '#CE906C80', '2' = '#44837380'),
                    labels = c('nearshore', 'offshore'), 
                    name = 'region', na.translate = FALSE) +
  scale_color_manual(values = c('1' = '#CE906C80', '2' = '#44837380'),
                     labels = c('nearshore', 'offshore'), 
                     name = 'region', na.translate = FALSE) +
  theme(aspect.ratio = 1,
        panel.grid.major = element_line(color = col_grid, linetype = 'dashed', size = 0.5),
        panel.background = element_rect(fill = NA),
        panel.ontop = TRUE,
        text = element_text(size = 10),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines"),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.title = element_blank(),
        legend.position = c(0.25, 0.38),
        panel.border = element_rect(fill = NA, colour = 'black')) +
  theme(axis.text.y = element_text(margin = margin(0, -25, 0, 25, unit = 'points'), color = 'black'),
        axis.text.x = element_text(vjust = 5, margin = margin(-7, 0, 7, 0, unit = 'points'), color = 'black'),
        axis.ticks.length = unit(-5, "points")) +
  guides(fill = guide_legend(override.aes = list(color = "black")))

  #save multiplots
  ragg::agg_png(filename = './figures/regions_surveys2.tif',res = 300,width = 11,height = 3,units = 'in')
  cowplot::plot_grid(p2,p1bc,nrow = 1,rel_widths = c(1,2),rel_heights = c(1,0.8),labels = c('a','b'),
                     label_x = c(0.92,0.96),label_y = c(0.98,0.97),byrow = TRUE,align = 'v')
  dev.off()
  
  #save multiple plot
  ragg::agg_png(filename = './figures/regions_surveys.tif',res = 300,width = 11,height = 4,units = 'in')
  cowplot::plot_grid(p2,p1b,p1c,nrow = 1,rel_widths = c(1,0.93,0.93),rel_heights = c(1,0.93,0.93),labels = c('a','b','c'),
                     label_x = 0.92,label_y = 0.89,byrow = TRUE,align = 'v')
  dev.off()
  
 
#########################################################
# AGE COMPOSITION
#########################################################
 
  #to store
  df5.IBM<-data.frame(matrix(nrow = 0,ncol = 4))
  colnames(df5.IBM)<-c('year','variable','value','region')
  sel_data3<-df5.IBM
  
  #list to store
  agecomp_list<-list()
  
  df5.IBM_error<-data.frame(matrix(nrow = 0,ncol = 4))
  colnames(df5.IBM_error)<-c('year','variable','value','region')
  agecomp_list_error<-list()

  #loop over regions to get age comp based on age structure numbers
for (region in c('1','2','All')) {
  
  #region<-'All'
  
  #read files
  lf.gag.csv<-lf.gag[grep(paste0(region,'_Number.csv'),lf.gag)]
  csv.gag<-read.csv(paste0('./outputs/OM/',lf.gag.csv),skip = 27) 
  csv.gag<-csv.gag[,-2]
  
  #reshape and prepare dataframe
  df.IBM<-reshape2::melt(csv.gag,id='Timestep')
  df.IBM<-df.IBM[order(df.IBM$Timestep),]
  
  #add month
  df.IBM$month.age<-1:(ncol(csv.gag)-1)
  
  #add age
  age<-sort(rep(0:(round((ncol(csv.gag)-1)/12)),times=12))
  age<-age[1:(ncol(csv.gag)-1)]
  df.IBM$age<-age
  
  #remove timestep0 because it is the ecopath baseline
  df.IBM<-df.IBM[which(df.IBM$Timestep!=0),]
  
  #add year
  df.IBM$year<-sort(rep(1985:2020,times=(12*377)))
  #length(sort(rep(1:nyrs,times=(12*375))))
  
  #aggregate by year timestep by averaging monthly steps
  IBM.mean<-aggregate(list(value=df.IBM$value),
                      by=list(year=df.IBM$year,
                              age=df.IBM$age),
                      FUN=mean) #FUN=mean - it looks mean because if not those estimates look pretty high
  
  #select year age and value
  df1.IBM.mean<-IBM.mean[,c('year','age','value')]
  
  #add obs error
  oerr<-rnorm(nrow(df1.IBM.mean),0,0.2) #SD=0.2
  
  #select year age and value
  df1.IBM.mean_error<-IBM.mean[,c('year','age','value')]
  df1.IBM.mean_error$value<-IBM.mean$value*exp(oerr)
  
  #from long to wide
  df2.IBM.mean<-reshape(df1.IBM.mean, idvar = "year", timevar = "age", direction = "wide")
  df2.IBM.mean_error<-reshape(df1.IBM.mean_error, idvar = "year", timevar = "age", direction = "wide")
  yrs<-df2.IBM.mean$year
  
  #total numbers
  total<-rowSums(df2.IBM.mean[,-1])
  df3.IBM.mean<-df2.IBM.mean[,-1]/total

  #change col names
  colnames(df3.IBM.mean)<-c(0:(ncol(df3.IBM.mean)-1))

  #append years
  df4.IBM.mean<-cbind(year=yrs,df3.IBM.mean)
  
  #selectivity data
  sel_data<-df2.IBM.mean[,-1]
  colnames(sel_data)<-c(0:(ncol(df3.IBM.mean)-1))
  sel_data1<-cbind(year=yrs,sel_data)
  sel_data2<-reshape2::melt(sel_data1,id='year')
  sel_data2$region<-region
  sel_data3<-rbind(sel_data3,sel_data2)
  
  #append into list
  agecomp_list[[region]]<-df4.IBM.mean
  
  #reshape and region
  df5.IBM.mean<-reshape2::melt(df4.IBM.mean,id='year')
  df5.IBM.mean$region<-region
  
  #append
  df5.IBM<-rbind(df5.IBM,df5.IBM.mean)
  
  #total numbers with errors
  total_error<-rowSums(df2.IBM.mean_error[,-1])
  df3.IBM.mean_error<-df2.IBM.mean_error[,-1]/total_error
  colnames(df3.IBM.mean_error)<-c(0:(ncol(df3.IBM.mean_error)-1))
  df4.IBM.mean_error<-cbind(year=yrs,df3.IBM.mean_error)
  agecomp_list_error[[region]]<-df4.IBM.mean_error
  
  #reshape
  df5.IBM.mean_error<-reshape2::melt(df4.IBM.mean_error,id='year')
  df5.IBM.mean_error$region<-region
  
  #append
  df5.IBM_error<-rbind(df5.IBM_error,df5.IBM.mean_error)
  }  
    
  
  #example age comp data from 1985
  ggplot()+
    geom_bar(data=subset(df5.IBM,year==1985),aes(x=variable,y=value),stat='identity')+
    facet_wrap(~region,scales='free_y')
  
  #selectivity data vt region
  region1<-subset(sel_data3,region=='1')
  region2<-subset(sel_data3,region=='2')
  regionall<-subset(sel_data3,region=='All')
  
  #get proportion relative to total
  region1$value<-region1$value/regionall$value
  region2$value<-region2$value/regionall$value
  #y1985<-subset(region1,year==1985)
  
  #plot annual selectivity region 1
  ragg::agg_png(filename = './figures/sel_survey1.png',res = 300,width = 10,height = 10,units='in')
  ggplot()+
    geom_point(data=region1,aes(x=variable,y=value,group=year),stat='identity',alpha=0.7,color='#CE906C')+
    geom_line(data=region1,aes(x=variable,y=value,group=year),stat='identity',alpha=0.7,color='#CE906C')+
    facet_wrap(~year)+
    theme_bw()+
    theme(strip.background = element_blank())+
    scale_y_continuous(name='selectivity')+
    scale_x_discrete(breaks = c(0,5,10,15,20,25,30),name='age')
  dev.off()
  
  #plot annual selectivity region 2
  ragg::agg_png(filename = './figures/sel_survey2.png',res = 300,width = 10,height = 10,units='in')
  ggplot()+
    geom_point(data=region2,aes(x=variable,y=value,group=year),stat='identity',alpha=0.7,color='#448373')+
    geom_line(data=region2,aes(x=variable,y=value,group=year),stat='identity',alpha=0.7,color='#448373')+
    facet_wrap(~year)+
    theme_bw()+
    theme(strip.background = element_blank())+
    scale_y_continuous(name='selectivity')+
    scale_x_discrete(breaks = c(0,5,10,15,20,25,30),name='age')
  dev.off()
  
  ggplot()+
    geom_point(data=region1,aes(x=variable,y=value,group=year),stat='identity',alpha=0.7,color='#CE906C')+
    geom_line(data=region1,aes(x=variable,y=value,group=year),stat='identity',alpha=0.7,color='#CE906C')+
    geom_point(data=region2,aes(x=variable,y=value,group=year),stat='identity',alpha=0.7,color='#448373')+
    geom_line(data=region2,aes(x=variable,y=value,group=year),stat='identity',alpha=0.7,color='#448373')+
    facet_wrap(~year)+
    theme_bw()+
    theme(strip.background = element_blank())+
    scale_y_continuous(name='selectivity')+
    scale_x_discrete(breaks = c(0,5,10,15,20,25,30),name='age')
  
  #mean across years
  region11<-aggregate(value ~ variable, region1,FUN=mean)
  region11$region<-'nearshore'
  region22<-aggregate(value ~ variable, region2,FUN=mean)
  region22$region<-'offshore' 
  region<-rbind(region11,region22)
  region$variable<-as.numeric(region$variable)
  
  #plot selectivity surveys
  ragg::agg_png(filename = './figures/sel_survey.png',res = 300,width = 6,height = 3,units='in')
  ggplot()+
    geom_point(data=region,aes(x=variable,y=value,group=region,color=region),stat='identity',alpha=0.7)+
    geom_line(data=region,aes(x=variable,y=value,group=region,color=region),stat='identity',alpha=0.7)+
    theme_bw()+
    scale_color_manual(values=c('nearshore'='#CE906C','offshore'='#448373'))+
    theme(strip.background = element_blank(),legend.position = c(0.8,0.5),legend.key.size  = unit(0.3,'in'))+
    scale_y_continuous(name='selectivity')+
    scale_x_continuous(breaks = c(1,6,11,16,21,26,31),labels=c(0,5,10,15,20,25,30),name='age',limits=c(1,20))
  dev.off()
  
  #age-comp offshore
  p2<-ggplot()+
    geom_point(data = subset(df5.IBM,variable %in% 0:15 & region=='2'),aes(y=variable,x=year,size=value),color='#448373',alpha=0.5)+
    scale_size(range=c(0,10),
               breaks=c(0,0.001,0.01,0.1,0.2,0.5),
               name = "proportion",
               guide="legend")+
    scale_x_continuous(breaks=c(1985:2020),labels=c(1985:2020))+ #,expand = c(0.01,0.01)
    theme_bw()+
    xlab('year')+
    ylab('age')+
    theme(aspect.ratio = 0.5,panel.grid.minor = element_blank(),
          panel.grid = element_line(linetype = 'dashed'),panel.grid.major.y = element_line(linetype = 'dashed'),
          axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))#+

  #age-comp nearshore
  p1<-ggplot()+
    geom_point(data = subset(df5.IBM,variable %in% 0:15 & region=='1'),aes(y=variable,x=year,size=value),color='#CE906C',alpha=0.5)+
    scale_size(range=c(0,10),
               breaks=c(0.001,0.01,0.1,0.2,0.5),
               #labels=c("0","5","10","20","25+"),
               name = "proportion",
               guide="legend")+
    scale_x_continuous(breaks=c(1985:2020),labels=c(1985:2020))+ #,expand = c(0.01,0.01)
    theme_bw()+
    xlab('year')+
    ylab('age')+
    theme(aspect.ratio = 0.5,panel.grid.minor = element_blank(),
          panel.grid = element_line(linetype = 'dashed'),panel.grid.major.y = element_line(linetype = 'dashed'),
          axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1),axis.title.x.bottom = element_blank())#+

  #multiplot age comp
  tiff(filename = './outputs/age_comp1.tif',res = 300,width = 2200,height = 2000)
  cowplot::plot_grid(p1,p2,nrow = 2,rel_heights = c(0.945,1),labels = c('a','b'),
                     label_x = 0.80,label_y = 0.98,byrow = TRUE,align = 'v')
  dev.off()
  
  #save age comp 2
  tiff(filename = './outputs/age_comp2.tif',res = 300,width = 1500,height = 2000)
   ggplot(data = subset(df5.IBM.mean,variable %in% 0:9),aes(x=variable,y=year,group=year,height=value)) + 
    ggridges::geom_density_ridges2(stat = "identity", scale = 1,fill='#007bf3',color='#007bf3',alpha=0.5,panel_scaling = FALSE)+
    theme_bw()+
    scale_y_reverse(breaks=c(1985:2020),labels=c(1985:2020),expand = c(0.01,0.01))+
    scale_x_discrete(expand = c(0.01,0.01))+
    xlab('age')+
    ylab('year')+
    theme(aspect.ratio = 1.5,panel.grid.minor.y = element_blank(),
          panel.grid = element_line(linetype = 'dashed'),panel.grid.major.y = element_line(linetype = 'dashed'))
  dev.off()
  
  #_yr month fleet sex part ageerr Lbin_lo Lbin_hi Nsamp datavector(female-male)
  agecomp<-data.frame('yr'=c(df4.IBM.mean$year,df4.IBM.mean$year,df4.IBM.mean$year),
                      'month'=7,
                      'fleet'=rep(c(2,3,'All'),each=length(df4.IBM.mean$year)),
                      'sex'=0,
                      'part'=0,
                      'ageerr'=1,
                      'Lbin_lo'=-1,
                      'Lbin_hi'=-1,
                      'Nsamp'=100,
                      round(rbind(agecomp_list$`1`[,-1],agecomp_list$`2`[,-1],agecomp_list$All[,-1]),digits = 3))
  
  #add end line
  agecomp1<-rbind(agecomp,c(-9999,rep(0,ncol(agecomp)-1)))
  agecomp2<-agecomp1[agecomp1$fleet!='All',]
  
  #save data
  write.table(agecomp2, 
              file = "./outputs/age_composition.txt", sep = "\t",
              row.names = FALSE,
              col.names = FALSE)

  #_yr month fleet sex part ageerr Lbin_lo Lbin_hi Nsamp datavector(female-male)
  agecomp<-data.frame('yr'=c(df4.IBM.mean_error$year,df4.IBM.mean_error$year,df4.IBM.mean_error$year),
                      'month'=7,
                      'fleet'=rep(c(2,3,'All'),each=length(df4.IBM.mean_error$year)),
                      'sex'=0,
                      'part'=0,
                      'ageerr'=1,
                      'Lbin_lo'=-1,
                      'Lbin_hi'=-1,
                      'Nsamp'=100,
                      round(rbind(agecomp_list_error$`1`[,-1],agecomp_list_error$`2`[,-1],agecomp_list_error$All[,-1]),digits = 3))
  
  #add end line
  agecomp1<-rbind(agecomp,c(-9999,rep(0,ncol(agecomp)-1)))
  agecomp2<-agecomp1[agecomp1$fleet!='All',]

  #save data
  write.table(agecomp2, 
              file = "./outputs/age_composition_error.txt", sep = "\t",
              row.names = FALSE,
              col.names = FALSE)
  
  
#########################################################
# OBS ERROR
#########################################################

  #oerr<-rnorm(length(unique(samplesadu2$year)),0,sigp2)
  #b*exp(oerr)
  
#########################################################
# CATCH DATA - AGECOMP
#########################################################
#_Catch data in SS file: yr, seas, fleet, catch, catch_se
  
  #agecomp3<-agecomp1[agecomp1$fleet=='All',paste0('X',0:31)] #1985:2020
  agecomp3<-agecomp_list$All[,as.character(5:31)]
  agecomp4<-sweep(agecomp3,1,rowSums(agecomp3),FUN = '/')
  
  #agecomp3<-agecomp1[agecomp1$fleet=='All',paste0('X',0:31)] #1985:2020
  agecomp3_error<-agecomp_list_error$All[,as.character(5:31)]
  agecomp4_error<-sweep(agecomp3_error,1,rowSums(agecomp3_error),FUN = '/')
    
#########################################################
# CATCH DATA
#########################################################
#_Catch data in SS file: yr, seas, fleet, catch, catch_se
  
  #filter catch file
  lf<-list.files(paste0(wd,'/outputs/OM'))
  lf.catch<-lf[grep("Annual_Average_Catch",lf)]
  
  #read csv catch file
  catch<-read.csv(paste0(wd,'/outputs/OM/',lf.catch),skip = 31)
  
  #filter gag
  catch.gag<-catch[,grep('gag',colnames(catch))]
  
  #df to store values
  df.catch<-data.frame(matrix(nrow=36,ncol = 0))
  #df to store values
  df.catch_error<-data.frame(matrix(nrow=36,ncol = 0))
  
  #loops over stanzas
  for (g in paste0('gag.',1:5)) {
    
    #g<-'gag.1'
    
    df.catch<-cbind(df.catch,rowSums(catch.gag[,grep(g,colnames(catch.gag))]))
    colnames(df.catch)[ncol(df.catch)]<-g
    
    df.catch_error<-cbind(df.catch_error,rowSums(catch.gag[,grep(g,colnames(catch.gag))]))
    colnames(df.catch_error)[ncol(df.catch_error)]<-g
    
    #add obs error
    oerr<-rnorm(nrow(df.catch_error),0,0.2) #SD=0.2
    
    df.catch_error[,g]<-df.catch_error[,g]*exp(oerr)
  }
  
  #df to store values
  df.catch.all<-data.frame(matrix(nrow=36,ncol = 0))
  #df to store values
  df.catch.all_error<-data.frame(matrix(nrow=36,ncol = 0))
  
  #loop over age 5+
  for (g in as.character(5:31)) {
  
    #g<-6
    
    #without error
    df.catch.all<-cbind(df.catch.all,df.catch[,'gag.5']*agecomp4[,paste0(g)])
    colnames(df.catch.all)[ncol(df.catch.all)]<-paste0('X',g)
    
    #with error
    df.catch.all_error<-cbind(df.catch.all_error,df.catch_error[,'gag.5']*agecomp4_error[,paste0(g)])
    colnames(df.catch.all_error)[ncol(df.catch.all_error)]<-paste0('X',g)
    
  }
  
  #append and colnames
  catch.age<-cbind(0,df.catch[1:4],df.catch.all)
  colnames(catch.age)<-as.character(0:31)
  
  #catch apecomp
  catch.comp<-sweep(catch.age,1,rowSums(catch.age),FUN = '/')
  catch.age_error<-cbind(0,df.catch_error[1:4],df.catch.all_error)
  colnames(catch.age_error)<-as.character(0:31)
  catch.comp_error<-sweep(catch.age_error,1,rowSums(catch.age_error),FUN = '/')
  
  #_yr month fleet sex part ageerr Lbin_lo Lbin_hi Nsamp datavector(female-male)
  agecomp<-data.frame('yr'=c(df4.IBM.mean$year,df4.IBM.mean$year,df4.IBM.mean$year),
                      'month'=7,
                      'fleet'=rep(c(2,3,1),each=length(df4.IBM.mean$year)),
                      'sex'=0,
                      'part'=0,
                      'ageerr'=1,
                      'Lbin_lo'=-1,
                      'Lbin_hi'=-1,
                      'Nsamp'=100,
                      round(rbind(agecomp_list$`1`[,-1],agecomp_list$`2`[,-1],catch.comp),digits = 3))
  
  #add end line
  agecomp1<-rbind(agecomp,c(-9999,rep(0,ncol(agecomp)-1)))
  
  #save data
  write.table(agecomp1, 
              file = "./outputs/age_composition_wfishery.txt", sep = "\t",
              row.names = FALSE,
              col.names = FALSE)
  
  #reshape and colnames
  rownames(catch.comp)<-1985:2020
  catch_plot<-reshape2::melt(as.matrix(catch.comp))
  
  #fishery selectivity smoother
  ggplot()+
    geom_point(data=catch_plot,aes(x=Var2,y=value),stat='identity')+
    geom_smooth(data=catch_plot,aes(x=Var2,y=value))

  #_yr month fleet sex part ageerr Lbin_lo Lbin_hi Nsamp datavector(female-male)
  agecomp<-data.frame('yr'=c(df4.IBM.mean_error$year,df4.IBM.mean_error$year,df4.IBM.mean_error$year),
                      'month'=7,
                      'fleet'=rep(c(2,3,1),each=length(df4.IBM.mean_error$year)),
                      'sex'=0,
                      'part'=0,
                      'ageerr'=1,
                      'Lbin_lo'=-1,
                      'Lbin_hi'=-1,
                      'Nsamp'=100,
                      round(rbind(agecomp_list_error$`1`[,-1],agecomp_list_error$`2`[,-1],catch.comp_error),digits = 3))
  
  #add end line
  agecomp1<-rbind(agecomp,c(-9999,rep(0,ncol(agecomp)-1)))
  
  #save data
  write.table(agecomp1, 
              file = "./outputs/age_composition_wfishery_error.txt", sep = "\t",
              row.names = FALSE,
              col.names = FALSE)
  
  #fleets
  fleets<-c('private.rec.nearshore','private.rec.offshore','charter.inshore','charter.offshore',
            'headboat','shore','reef.fish.vertical.line','reef.fish.bottom.longline')
  
  #df to store values
  df.catch<-data.frame(matrix(nrow=36,ncol = 0))
  
  #loops over fleets
  for (f in fleets) {
    #f<-fleets[6]
    df.catch<-cbind(df.catch,rowSums(catch.gag[,grep(paste0('^',f),colnames(catch.gag))]))
    colnames(df.catch)[ncol(df.catch)]<-f
  }
  
  #absolute catch
  df.catch<-df.catch*188688.717455821
  df.catch$total<-rowSums(df.catch)
  df.catch$year<-c(1985:2020)
  
  #fishing mortality
  bio.csv<-read.csv('./outputs/OM/Ecospace_Average_Biomass.csv',skip = 31) 
  catch.csv<-read.csv('./outputs/OM/Ecospace_Average_Catch.csv',skip = 31) 
  head(bio.csv)
  bio.gag<-rowSums(bio.csv[,c('gag.0','gag.1','gag.2','gag.3','gag.4','gag.5.')])
  head(catch.csv)
  catch.gag<-rowSums(catch.csv[,grep('gag',colnames(catch.csv))])
  catch.gag/bio.gag
  
  catch<-data.frame(year=rep(1985:2020,each=12),
                    catch=catch.gag,
                    f=catch.gag/bio.gag)
  
  #absolute catch
  p1<-
    ggplot(data = catch,aes(x=year,y=catch*188688.717455821,group=year)) +
    geom_boxplot(fill='grey70',color='grey40',outlier.colour = 'transparent',size=0.3)+ #,size=1.5 #007bf3
    theme_bw()+
    scale_x_continuous(breaks=c(1985,1990,1995,2000,2005,2010,2015,2020),minor_breaks = 1985:2020,expand = c(0.01,0.01))+ #,position = "top"
    scale_y_continuous(limits = c(0,5100),expand = c(0,0))+
    theme_bw()+
    xlab('year')+
    ylab('catch (tones)')+
    theme(text = element_text(size=12),
          axis.title.x = element_blank(), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"),
          legend.background =  element_rect(fill = "transparent", colour = "transparent"),
          legend.key = element_rect(fill = NA, colour = NA),
          legend.position = 'none',
          panel.border = element_rect(fill = NA, colour = 'black'),
          panel.grid.minor =  element_line(linetype = 'dashed'))+
    theme(strip.background = element_blank(),strip.text = element_blank())#+
  
  #fishing mortality
  p2<-ggplot(data = catch,aes(x=year,y=f,group=year)) +
    geom_boxplot(fill='grey70',color='grey40',outlier.colour = 'transparent',size=0.3)+ #,size=1.5 #007bf3
    theme_bw()+
    scale_x_continuous(breaks=c(1985,1990,1995,2000,2005,2010,2015,2020),minor_breaks = 1985:2020,expand = c(0.01,0.01))+ #,position = "top"
    scale_y_continuous(limits = c(0,0.66),expand = c(0,0))+
    theme_bw()+
    xlab('year')+
    ylab(expression (paste("F (",year^{-1},')')))+
    theme(text = element_text(size=12),
          axis.title.x = element_blank(), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"),
          legend.background =  element_rect(fill = "transparent", colour = "transparent"),
          legend.key = element_rect(fill = NA, colour = NA),
          legend.position = 'none',
          panel.border = element_rect(fill = NA, colour = 'black'),
          panel.grid.minor =  element_line(linetype = 'dashed'))+
    theme(strip.background = element_blank(),strip.text = element_blank())#+
  
  #save multiple plot
  tiff(filename = './outputs/catch_f.png',res = 300,width =5,height = 5,units = 'in')
  cowplot::plot_grid(p1,p2,nrow = 2,labels = c('a','b'),rel_heights = c(1,1), 
                     label_x = 0.93,label_y = 0.97,byrow = TRUE,align = 'v')
  dev.off()
  
  #trend absolute catch
  ggplot(data = df.catch,aes(x=year,y=total)) +
    geom_line(color='#007bf3')+ #,size=1.5
    theme_bw()+
    #theme(panel.grid.minor = element_blank())
    scale_x_continuous(breaks=c(1985:2020),labels=c(1985:2020),expand = c(0.01,0.01))+ #,expand = c(0.01,0.01)
    scale_y_continuous(limits = c(0,5000),expand = c(0.01,0.01))+
    theme_bw()+
    xlab('year')+
    ylab('catch (tones)')+
    theme(aspect.ratio = 0.5,panel.grid.minor = element_blank(),
          panel.grid = element_line(linetype = 'dashed'),panel.grid.major.y = element_line(linetype = 'dashed'),
          axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))
  
  #_Catch data: yr, seas, fleet, catch, catch_se
  df.catch1<-data.frame('yr'=df.catch$year,
                        'seas'=1,
                        'fleet'=1,
                        'catch'=df.catch$total,
                        'catch_se'=0.1)
  
  #initial state
  df.catch1<-rbind(c(-999,1,1,mean(df.catch1$catch[1:5]),0.1),
                   df.catch1,
                   c(-9999,0,0,0,0))
  
  #save data
  write.table(df.catch1, 
              file = "./outputs/yield.txt", sep = "\t",
              row.names = FALSE,
              col.names = FALSE)
  
  #add obs error
  oerr<-rnorm(nrow(df.catch),0,0.2) #SD=0.2
  
  #_Catch data: yr, seas, fleet, catch, catch_se
  df.catch1<-data.frame('yr'=df.catch$year,
                        'seas'=1,
                        'fleet'=1,
                        'catch'=df.catch$total*exp(oerr),
                        'catch_se'=0.1)
  
  #initial state
  df.catch1<-rbind(c(-999,1,1,mean(df.catch1$catch[1:5]),0.1),
                   df.catch1,
                   c(-9999,0,0,0,0))
  
  #save data
  write.table(df.catch1, 
              file = "./outputs/yield_error.txt", sep = "\t",
              row.names = FALSE,
              col.names = FALSE)
  
############################################################
# GET RED TIDE MORTALITY DATA
###########################################################

  #read csv loss and bio files
  loss<-read.csv(paste0('./outputs/OM/Ecospace_Average_OtherMortalityLoss.csv'),skip = 31)
  loss.gag<-as.data.frame(loss[,grepl("gag", colnames(loss))])
  max(loss.gag)
  
  bio<-read.csv(paste0('./outputs/OM/Ecospace_Average_Biomass.csv'),skip = 31)
  bio.gag<-as.data.frame(bio[,grepl("gag", colnames(bio))])
  max(bio.gag)
  
  #divide biomass by mortality loss to get red tide mortality rate in each year
  df.mrt<-(loss.gag/bio.gag)
  loss.gag/bio.gag
  max(df.mrt$gag.0)
  
  #df.mrt<-loss.gag
  mf <- function(x){sd(df.mrt[,x])/sqrt(length(df.mrt[,x]))}
  sapply(1:length(df.mrt[1,]),mf)
  
  #add year
  df.mrt$yr<-rep(1985:2020,each=12)
  
  #reshape
  df.mrt1<-reshape2::melt(df.mrt,id.vars='yr')
  
  #mean, sd and se
  df.mrt2.m<-aggregate(df.mrt1$value,by=list(year=df.mrt1$yr,age=df.mrt1$variable),FUN=mean) #FUN=mean - it looks mean because if not those estimates look pretty high
  df.mrt2.sd<-aggregate(df.mrt1$value,by=list(year=df.mrt1$yr,age=df.mrt1$variable),FUN=sd)
  df.mrt2.l<-aggregate(df.mrt1$value,by=list(year=df.mrt1$yr,age=df.mrt1$variable),FUN=length)
  df.mrt2<-cbind(df.mrt2.m,df.mrt2.sd$x,df.mrt2.l$x)
  colnames(df.mrt2)<-c('year','age','mean','sd','length')
  df.mrt2$se <- df.mrt2$sd / sqrt(df.mrt2$length)
  
  #read csv loss and bio files
  loss<-read.csv(paste0('./outputs/OM/Ecospace_Annual_Average_OtherMortalityLoss.csv'),skip = 31)
  loss.gag<-as.data.frame(loss[,grepl("gag", colnames(loss))])
  max(loss.gag)
  bio<-read.csv(paste0('./outputs/OM/Ecospace_Annual_Average_Biomass.csv'),skip = 31)
  bio.gag<-as.data.frame(bio[,grepl("gag", colnames(bio))])
  max(bio.gag)
  df.mrt.annual<-(loss.gag/bio.gag)
  max(df.mrt.annual)
  
  #reshape
  df.mrt.annual$yr<-rep(1985:2020)
  df.mrt1.annual<-reshape2::melt(df.mrt.annual,id.vars='yr')
  names(df.mrt1.annual)[2]<-'age'
  levels(df.mrt1.annual$age)<-c('gag 0','gag 1','gag 2','gag 3','gag 4','gag 5 plus')
  df.mrt1.annual$se<-df.mrt2$se
  
  #palette
  pal <- wesanderson::wes_palette("Zissou1", 6, type = "continuous")
  
  #plot MRT
  tiff(filename = './outputs/mrt3.tif',res = 250,width = 2200,height = 1200)
  ggplot() +
    geom_line(data = df.mrt1.annual,aes(x=yr,y=value,color=age),linewidth=1,alpha=0.7)+ #,size=1.5
    theme_bw()+
    scale_x_continuous(breaks=c(1985,1990,1995,2000,2005,2010,2015,2020),minor_breaks = 1985:2020,expand = c(0,0))+ #,position = "top"
    scale_y_continuous(expand = c(0.01,0),limits=c(0,1.10))+
    theme_bw()+
    scale_color_manual(values = pal) +
    xlab('year')+
    ylab(expression (paste("MRT (",year^{-1},')')))+
    theme(text = element_text(size=12),
          axis.title.x = element_blank(), 
          plot.margin = unit(c(1,1,1,1), "lines"),
          legend.key = element_rect(fill = NA, colour = NA),
          panel.border = element_rect(fill = NA, colour = 'black'),
          panel.grid.minor =  element_line(linetype = 'dashed'),
          legend.position = c(0.1,0.77),legend.background = element_rect(color = "black", fill = NA),
          legend.title = element_blank(),legend.spacing.y = unit(0, "lines"))
  dev.off()
  
  
  ragg::agg_png('./figures/MRT_age_uncertainty.png',  width = 8, height = 5, units = "in", res = 300)
  ggplot(data = df.mrt1.annual, aes(x = yr, y = value, color = age)) +
    geom_line(linewidth = 1, alpha = 0.7) +  # Main line plot
    geom_ribbon(aes(ymin = value - se, ymax = value + se, fill = age), alpha = 0.3,color=NA) +  # Uncertainty ribbon
    theme_bw() +
    scale_x_continuous(
      breaks = c(1990, 1995, 2000, 2005, 2010, 2015,2020),
      minor_breaks = c(1985:2020),
      expand = c(0, 0)) +
    #scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +  # Use the same color palette for fill
    xlab('Year') +
    ylab(expression(paste("MRT (", year^{-1}, ')'))) +
    theme(
      text = element_text(size = 12),
      axis.title.x = element_blank(),
      axis.text.x =element_text(angle=30,hjust=1,vjust=1),
      #plot.margin = unit(c(1, 1, 1, 1), "lines"),
      #legend.key = element_rect(fill = NA, colour = NA),
      #panel.border = element_rect(fill = NA, colour = 'black'),
      panel.grid.minor = element_line(linetype = 'dashed'),
      legend.position = 'none',#c(0.1, 0.77),
      #legend.background = element_rect(color = "black", fill = NA),
      legend.title = element_blank(),strip.background = element_rect(fill='white'))+
    facet_wrap(~age)
  dev.off()
  
  
  #scaled function
  scale_value<- function(x){
    (x-mean(x,na.rm=TRUE))/sd(x)
  }
  
  #scale
  df.mrt1.annual$scaled<-scale_value(df.mrt1.annual$value)
  
  
  #save
  write.csv2(df.mrt1.annual,file='./tables/mrt_annual.csv')
  
  
  #selected severe years
  sev_yrs<-unique(df.mrt1.annual[which(df.mrt1.annual$scaled>=1),'yr'])
  # sev_yrs1<-sort(unique(df.mrt3[which(df.mrt3$scaled>=0),'year']))
  # sev_yrs2<-sort(unique(df.mrt3[which(df.mrt3$scaled>=0.5),'year']))
  
  #plot scaled MRT
  tiff(filename = './outputs/mrt5.tif',res = 250,width = 2200,height = 1200)
  ggplot() +
    #geom_hline(yintercept = c(-1,1),linetype='dotted',alpha=0.7)+
    geom_rect(aes(xmin=1985,xmax=2020,ymin=-1,ymax=+1),alpha=0.2)+
    geom_hline(yintercept = c(0),linetype='dotted',alpha=0.7)+
    geom_vline(xintercept = c(2005,2006,2018,2019),linetype='dashed',alpha=0.7)+
    geom_line(data = df.mrt1.annual,aes(x=yr,y=scaled,color=age),size=1,alpha=0.7)+ #,size=1.5
    theme_bw()+
    #theme(panel.grid.minor = element_blank())
    #scale_x_continuous(breaks=c(1985:2020),labels=c(1985:2020),expand = c(0.01,0.01))+ #,expand = c(0.01,0.01)
    scale_x_continuous(breaks=c(1985,1990,1995,2000,2005,2010,2015,2020),minor_breaks = 1985:2020,expand = c(0,0))+ #,position = "top"
    scale_y_continuous(expand = c(0.01,0),limits = c(-1,6))+ #,limits=c(-0.5,1.10)
    #scale_y_continuous(limits = c(0,0.15),expand=c(0,0.001),breaks = c(0,0.05,0.10,0.15))+
    theme_bw()+
    scale_color_manual(values = pal) +
    #facet_wrap(~variable,nrow = 3)+
    xlab('year')+
    ylab(expression (paste("scaled MRT (",year^{-1},')')))+
    theme(text = element_text(size=12),
          axis.title.x = element_blank(), 
          plot.margin = unit(c(1,1,1,1), "lines"),
          #legend.background =  element_rect(fill = "transparent", colour = "transparent"),
          legend.key = element_rect(fill = NA, colour = NA),
          #legend.position = 'none',
          panel.border = element_rect(fill = NA, colour = 'black'),
          panel.grid.minor =  element_line(linetype = 'dashed'),
          legend.position = c(0.1,0.77),legend.background = element_rect(color = "black", fill = NA),
          legend.title = element_blank(),legend.spacing.y = unit(0, "lines"))
  dev.off()
  
  
  
  
  #build MRT for all ages and years  
  mrt5<-df.mrt1.annual[which(df.mrt1.annual$age=='gag 5 plus'),]
  mrt5$age<-'gag 5'
  mrt51<-dplyr::bind_rows(replicate(length(6:32),mrt5, simplify=FALSE))
  mrt51$age<-rep(paste0('gag',c(6:32)),each=length(1985:2020))
  
  #append
  df.mrt22<-rbind(df.mrt1.annual,mrt51)
  
  #build txt
  df.mrt22$value<-ifelse(df.mrt22$value<0.00001,0,df.mrt22$value)
  df.mrt22$se<-ifelse(df.mrt22$se<0.00001,0,df.mrt22$se)
  df.mrt22$pr_type<-ifelse(df.mrt22$se<=0.00001,0,6)
  
  ################
  # SAVE BYCATCH DATA
  ################
  
  #save
  write.table(df.mrt22, 
              file = "./mrt_at_age.txt", sep = "\t",
              row.names = FALSE,
              col.names = FALSE)
  
  #create table
  df.mrt222<-data.frame('LO'=0,
                        'HI'=1.21,
                        'INIT'=df.mrt22$value,
                        'PRIOR'=df.mrt22$value  ,
                        'PR_SD'=df.mrt22$se,
                        'PR_type'=df.mrt22$pr_type,
                        'PHASE'=4,
                        "desc"=paste0("# NatM_break_1_Fem_GP_",as.numeric(df.mrt22$age),'_BLK1add_',df.mrt22$yr))
  
  #create txt
  age<-gsub('NatM_break_1_Fem_GP_','',df.mrt222$desc)
  age1<-gsub("\\_.*","",age)
  age1<-gsub('# ','',age1)
  
  #fxn
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  #modify year
  yr<-1985:2020
  
  #txt file
  df.mrt333<-data.frame('yr'=ifelse(round(df.mrt22$value,5)!=0,yr,paste0('#',yr)),
                        'month'=7,
                       'fleet'=as.double(age1)+3,
                       'catch'=round(df.mrt22$value,5),
                       'catch_se'=df.mrt22$se,
                       '#info'=paste0('# MRT_Age',age1))
  
  #file
  write.csv(df.mrt333,file = './tables/mrt_bycatch_F.csv',row.names = FALSE)
  
  #df
  df.mrt444<-data.frame('#_LO'=0,
                        'HI'=10,
                        'INIT'=aggregate(value ~ yr,df.mrt22,FUN=mean)[,'value'],
                        'PRIOR'=0.1,
                        'PR_SD'=1,
                        'PR_type'=0,
                        'PHASE'=3)
  
  #save
  write.csv(df.mrt444,file = './tables/mrt_pred.csv',row.names = FALSE)
  write.table(df.mrt222, 
              file = "./tables/mrt_at_age_blk.txt", sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  write.csv(df.mrt222,file = './tables/mrt_blk_add.csv',row.names = FALSE)
  
  #select for severe years
  df.mrt222sev<-df.mrt222[which(grepl(pattern = '2005|2006|2018|2019',df.mrt222$desc)),]
  
  #save
  write.table(df.mrt222sev, 
              file = "./mrt_at_age_blk_add_sev.txt", sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  #save
  write.csv(df.mrt222sev,file = './tables/mrt_blk_add_sev.csv',row.names = FALSE)
  
  #SD * 3
  df.mrt222$PR_SD<-(df.mrt222$PR_SD*3)
  
  #save
  write.table(df.mrt222, 
              file = "./mrt_at_age_blk3.txt", sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  write.csv(df.mrt222,file = './tables/mrt_blk_add3.csv',row.names = FALSE)
  
  #severe SD *3
  df.mrt222sev<-df.mrt222[which(grepl(pattern = '2005|2006|2018|2019',df.mrt222$desc)),]
  
  #save
  write.table(df.mrt222sev, 
              file = "./mrt_at_age_blk_add3_sev.txt", sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  write.csv(df.mrt222sev,file = './tables/mrt_blk_add3_sev.csv',row.names = FALSE)
  
  ################
  # SAVE ENV DATA
  ################
  #Yr Variable Value (gag red tide mortality rates) 
  #Z-transformation
  levels(df.mrt1.annual$age)<-c(1:6)
  df.mrt1.annual$age<-as.numeric(df.mrt1.annual$age)
  #z-scored
  df.mrt1.annual$value<-(df.mrt1.annual$value)/mean(df.mrt1.annual$value)
  df.mrt1.annual$value<-round(df.mrt1.annual$value,digits=3)
  df.mrt1<-rbind(df.mrt1.annual[,c("yr","age","value")],c(-9999,0,0))
  write.table(df.mrt1, 
              file = "./outputs/mrt_as_env.txt", sep = "\t",
              row.names = FALSE,
              col.names = FALSE) 
  
  ################
  # SAVE BYCATCH DATA
  ################
  #_Catch data: yr, seas, fleet, catch, catch_se
  
  #read txt yield file
  yield<-read.table('./outputs/yield.txt')
  yield<-yield[-nrow(yield),]
  colnames(yield)<-c('yr','seas','fleet','catch','catch_se')
  yield_error<-yield
  yield_error$catch<-yield_error$catch*exp(oerr)
  
  #read csv loss and bio files
  loss<-read.csv(paste0('./outputs/OM/Ecospace_Annual_Average_OtherMortalityLoss.csv'),skip = 31)
  loss.gag<-rowSums(loss[,grepl("gag", colnames(loss))])
  loss.gag1<-loss.gag*188688.717455821
  
  #_Catch data: yr, seas, fleet, catch, catch_se
  df.catch2<-data.frame('yr'=c(1985:2020),
                        'seas'=1,
                        'fleet'=4,
                        'catch'=loss.gag1*exp(oerr),
                        'catch_se'=0.1)
  
  
  df.catch2<-rbind(yield_error,df.catch2,c(-9999,0,0,0,0))
  
  #save
  write.table(df.catch2, 
              file = "./outputs/yield_bycatch_error.txt", sep = "\t",
              row.names = FALSE,
              col.names = FALSE) 

  ###################################################
  # FISHERY SELECTIVITY
  ###################################################
  
  #catch at age
  sel.catch<-catch.age
  sel.catch1<-data.frame('yr'=1985:2020,sel.catch)
  
  #reshape
  sel.catch2<-reshape2::melt(sel.catch1,id.vars='yr')
  sel.catch2$variable<-gsub('X','',sel.catch2$variable)
  
  dim(regionall);dim(sel.catch2)
  
  #proportion catch relative biomass
  sel.catch2$prop<-sel.catch2$value/regionall$value
  
  #relative to max
  sel.catch2$prop<-sel.catch2$prop/max(sel.catch2$prop)
  
  #factors
  sel.catch2$variable<-factor(sel.catch2$variable,levels=c(0:31))
  
  #ragg::agg_png(filename = './outputs/sel_survey2.png',res = 300,width = 10,height = 10,units='in')
  ggplot()+
    geom_bar(data=sel.catch2,aes(x=variable,y=value),stat='identity',alpha=0.7)+
    facet_wrap(~yr,scales='free_y')+
    theme_bw()+
    theme(strip.background = element_blank())+
    scale_y_continuous(name='proportion')+
    scale_x_discrete(breaks = c(0,5,10,15,20,25,30),name='age')
  #dev.off()
  
  ggplot()+
    geom_bar(data=sel.catch2,aes(x=variable,y=prop),stat='identity',alpha=0.7)+
    facet_wrap(~yr,scales='free_y')+
    theme_bw()+
    theme(strip.background = element_blank())+
    scale_y_continuous(name='proportion')+
    scale_x_discrete(breaks = c(0,5,10,15,20,25,30),name='age')
  
  #read csv catch file
  bio<-read.csv(paste0(wd,'/outputs/OM/Ecospace_Annual_Average_Biomass.csv'),skip = 31)
  
  #filter gag
  bio.gag<-bio[,grep('gag',colnames(bio))]
  
  #df to store values
  df.bio<-data.frame(matrix(nrow=36,ncol = 0))
  
  #df to store values
  df.bio.all<-data.frame(matrix(nrow=36,ncol = 0))
  
  #loop over 5 plus
  for (g in as.character(5:31)) {
    
    
    df.bio.all<-cbind(df.bio.all,bio.gag[,'gag.5.']*agecomp4[,paste0(g)])
    colnames(df.bio.all)[ncol(df.bio.all)]<-paste0('X',g)
  }
  
  #append
  bio.age<-cbind(bio.gag[1:5],df.bio.all)
  colnames(bio.age)<-as.character(0:31)
  
  #reshape
  bio.age2<-reshape2::melt(bio.age)
  
  #proportion
  sel.catch2$prop2<-sel.catch2$value/bio.age2$value
  
  #relative to max to get 1 as max selectivity
  sel.catch2$variable<-as.numeric(sel.catch2$variable)
  
  #check plot - 
  ggplot()+
    geom_bar(data=sel.catch2,aes(x=variable,y=prop2),stat='identity',alpha=0.7)+
    facet_wrap(~yr,scales='free_y')+
    theme_bw()+
    theme(strip.background = element_blank())+
    scale_y_continuous(name='proportion')+
    scale_x_continuous(name='age',limits=c(0,5))
  
  #read csv catch file
  bio<-read.csv(paste0(wd,'/outputs/OM/Ecospace_Annual_Average_Biomass.csv'),skip = 31)
  
  #filter gag
  bio.gag<-bio[,grep('gag',colnames(bio))]
  
  #read csv catch file
  catch<-read.csv(paste0(wd,'/outputs/OM/Ecospace_Annual_Average_Catch.csv'),skip = 31)
  
  #filter gag
  catch.gag<-catch[,grep('gag',colnames(catch))]
  
  #df to store values
  df.catch<-data.frame(matrix(nrow=36,ncol = 0))
  
  #loops over stanzas
  for (g in paste0('gag.',1:5)) {
    
    #g<-'gag.1'
    
    df.catch<-cbind(df.catch,rowSums(catch.gag[,grep(g,colnames(catch.gag))]))
    colnames(df.catch)[ncol(df.catch)]<-g
    
  }
  
  #for age 0 - 0 catches
  df.catch<-cbind(0,df.catch)
  
  #proportion and names
  sel.catch5<-df.catch/bio.gag
  rownames(sel.catch5)<-c(1985:2020)
  colnames(sel.catch5)<-as.character(c(0:5))
  sel.catch5<-reshape2::melt(as.matrix(sel.catch5))
  sel.catch5$value<-sel.catch5$value/max(sel.catch5$value)
  
  #aggregate to get the mean over years
  sel.catch6<-aggregate(value ~ Var2,sel.catch5,FUN=mean)
  sel.catch6$region<-'fishery'
  names(sel.catch6)<-names(region)
  sel.catch6$variable<-1:6
  
  #plot line plot
  ggplot()+
    geom_point(data=sel.catch5,aes(x=Var2,y=value),stat='identity',alpha=0.7)+
    geom_line(data=sel.catch5,aes(x=Var2,y=value),stat='identity',alpha=0.7)+
    facet_wrap(~Var1)+
    theme_bw()+
    theme(strip.background = element_blank())+
    scale_y_continuous(name='selectivity',breaks = c(0,0.5,1))+
    scale_x_continuous(breaks = c(0:5),name='age',labels=c('0','1','2','3','4','5+'))
  
  #join all selectivities (surveys and fishery)
  sel.catch7<-rbind(region,sel.catch6)
  
  #plot all selectivities
  ragg::agg_png(filename = './outputs/sel_fishery_survey.png',res = 300,width = 6,height = 3,units='in')
  ggplot()+
    geom_point(data=sel.catch7,aes(x=variable,y=value,group=region,color=region),stat='identity',alpha=0.7)+
    geom_line(data=sel.catch7,aes(x=variable,y=value,group=region,color=region),stat='identity',alpha=0.7)+
    theme_bw()+
    scale_color_manual(values=c('nearshore'='#CE906C','offshore'='#448373','fishery'='grey20'),labels=c('fishery','survey nearshore','survey offshore'))+
    theme(strip.background = element_blank(),legend.position = c(0.8,0.5),legend.key.size  = unit(0.3,'in'),legend.title = element_blank())+
    scale_y_continuous(name='selectivity')+
    scale_x_continuous(breaks = c(1,6,11,16,21,26,31),labels=c(0,5,10,15,20,25,30),name='age',limits=c(1,20))
  dev.off()
  
###################################################
# ADDITIONAL EXAMINATIONS
###################################################
 
#CHECK FG BIOMASS PREDICTIONS
   
  bio<-read.csv(paste0('./outputs/OM/Ecospace_Annual_Average_Biomass.csv'),skip = 31)
  bio<-bio[,-1]
  bio1<-scale(bio)
  bio1<-data.frame(cbind('year'=1985:2020,bio1))
  bio2<-reshape2::melt(bio1,id.vars='year')

  ggplot()+
    geom_line(data=bio2,aes(x=year,y=value))+ #,color=variable
    facet_wrap(~variable)
  
#CHECK WEIGHT - ERROR
  
  #all weight file
  w<-read.csv(paste0('./outputs/OM/AgeStructure_gag_Region_All_Weight.csv'),skip = 27)
  
  w.base<-w[1,-(1:2)]
  w.base1<-reshape2::melt(w.base)
  w.base1$months<-1:nrow(w.base1)
  
  ggplot()+
    geom_point(data=w.base1,aes(y=value,x=months))
  
  w1<-w[-1,3:ncol(w)]
  
  #file<-filew[,3:ncol(filew)]
  w1$Timestep<-w$Timestep[-1]
  
  #reshape and prepare dataframe
  df<-reshape2::melt(w1,id='Timestep')
  df<-df[order(df$Timestep),]
  
  #add month
  df$month<-1:(ncol(w1)-1)
  summary(df)
  age<-sort(rep(0:(round((ncol(w1)-1)/12)),times=12))
  
  #add age
  age<-age[1:(ncol(w1)-1)]
  df$age<-age  
  
  #add year
  df$year<-sort(rep(1985:2020,times=(12*377)))
  #length(sort(rep(1:nyrs,times=(12*375))))
  
  #aggregate by year timestep by averaging monthly steps
  df1<-aggregate(df$value,by=list(year=df$year,Timestep=df$Timestep,age=df$age),FUN=mean) #FUN=mean - it looks mean because if not those estimates look pretty high

  #plot
  ggplot()+
    geom_line(data=df1,aes(x=Timestep,y=x))+
    facet_wrap(~age)
  
  

  
  
  #####plot mrt with uncertainty
  
  
  x <- read.delim('./tables/mrt_at_age_blk.txt', header = FALSE)
  head(x)
  
  data.frame(mean=x$V3,
             se=x$V5) #sqrt(12)

