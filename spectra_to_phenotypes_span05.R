#!/usr/bin/env Rscript

#SCRIPT TO OUTPUT PHENOTYPES FROM SPECTRA
# Based on CLR v. 1.05, R. Montgomerie

library(argparser,quietly=T)
library(stringr)

p<-arg_parser('Calculate brightness, hue, and saturation from measured spectra with 1 nm bin width.')

p<- add_argument(p,'folder',help='folder to find .txt files to run')
p<- add_argument(p, 'lambda_min', help = 'min wavelength')
p<- add_argument(p, 'lambda_max', help = 'max wavelength')

argv<-parse_args(p)

lambda_min = as.numeric(argv$lambda_min) # min wavelength
lambda_max = as.numeric(argv$lambda_max) # max wavelength
strans = argv$folder # folder to search

today = format(Sys.Date(),"%Y%m%d")
binsize = 1 #bin size in nm
pattern = "*.txt" #which files it picks to process from 'fold'


out_file = file.path(strans,paste0('output_',lambda_min,'_to_',lambda_max,'_binsize_',binsize,'_',today,'_',basename(strans),'.csv'))
print('saving in file')
print(out_file)


#define the functions to read the spectrum files
read_txt_to_table<-function(f){
  lns<-readLines(f)
  brkpt<-grep('>',lns) #find the last line of the header
  head<-grepl('>',lns) #a vector the length of the file with 'true' if brkpt
  head[1:brkpt]=TRUE
  spc<-read.table(text=lns[!head])
  r<-list('header' = lns[head],'dat'=spc)
  return(r)
}

refine_and_smooth<-function(spc,xs){
  names(spc)=c('wavelength','refl')
  spc<-subset(spc,wavelength>=max(min(xs)-25,min(wavelength))&wavelength<=min(max(xs)+25,max(wavelength)))
  loessM<-loess(refl~wavelength,data=spc,span=0.05)
  smoothed<-predict(loessM,xs)
  smoothed[smoothed<0]=0
  dat<-data.frame('x'=xs,'y'=smoothed)
  return(dat)
}

calc_phenos <- function(df,xmin){
  df$yprop = df$y/100
  B1 = sum(df$yprop)
  if(xmin == 400){
    Rmax = 700
    RminYmax = 625
    YminGmax = 550
    GminBmax = 475
    Bmin = 400
  }
  if(xmin == 320){
    Rmax = 700
    RminYmax = 605
    YminGmax = 510
    GminBmax = 415
    Bmin = 320
  }
  if(xmin == 300){
    Rmax = 700
    RminYmax = 600
    YminGmax = 500
    GminBmax = 400
    Bmin = 300
  }
  R = sum(subset(df,x>=RminYmax & x<=Rmax)$yprop)
  Y = sum(subset(df,x>=YminGmax & x<=RminYmax)$yprop)
  G = sum(subset(df,x>=GminBmax & x<=YminGmax)$yprop)
  B = sum(subset(df,x>=Bmin & x<=GminBmax)$yprop)
  S5 = sqrt((R-G)**2 + (Y-B)**2)
  H4 = atan2((Y-B)/B1, (R-G)/B1)
  sat = S5/B1
  hue = H4 * 180/pi
  
  return(c(B1,hue,sat))
  
}


#Read in the list of files to process
fs<-list.files(strans,pattern='*.txt')
print('Including data from files:')
print(fs)
#define the wavelengths/bins:
xs = seq(lambda_min,lambda_max,binsize)

#Get rid of any preexisting data frame named 'out_df'
if(exists('out_df')){rm(out_df)}

#Set up the output file
out_df <- data.frame()

#For each spectrum file:
for(i in 1:length(fs)){
  
  #read in the data and separate the header and data
  fl<-read_txt_to_table(file.path(strans,fs[i]))
  spc<-fl$dat #the data only
  
  #Smooth the data with loess method, span = 0.05; very little smoothing.
  Sm<-refine_and_smooth(spc,xs)
  Y<-Sm$y
  
  # Calculate the phenotypes, based on Montgomerie CLR v 1.05 
  dout <- calc_phenos(Sm, min(Sm$x))
  
  #Use the header to get the file name and Date (Time = Date, as in CLR Vars)
  filename = str_sub(strsplit(fl$header[1],' ')[[1]][3],1,-5)
  d = strsplit(fl$header[3],'Date: ')[[1]][2]
  
  #Add the smoothed data to the output data frame
  out_df<-rbind(out_df,data.frame('file'=filename,'folder'=strans,'date'=d,
                                  'Brightness'=dout[1], 'Hue'=dout[2], 'Saturation' = dout[3],
                                  'min_lambda' = min(Sm$x)))
}

#Write the combined data frame to the output folder as a csv
write.csv(out_df,out_file,row.names=F)



