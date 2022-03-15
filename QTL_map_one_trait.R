library(qtl)
library(argparser)
library(ggplot2)
library(stringr)

# Parse command line arguments.
ap <- arg_parser('Run QTL mapping and save output for one trait.')
ap <- add_argument(ap, "--rqtl_object",help="path to rqtl object with genotypes and phenotypes to use for mapping")
ap <- add_argument(ap, "--var", help = "text name of the column corresponding to the trait you wish to map.")
ap <- add_argument(ap, "--save_folder", help = "path to the folder to save the resulting output",default='.')
ap <- add_argument(ap, "--model", help = "model to use for mapping. default is normal",default='normal')
ap <- add_argument(ap, "--method", help = "method to use for mapping, default is ehk",default='ehk')
ap <- add_argument(ap, "--nperm", help = "number of permutations to use",default='1000')
ap <- add_argument(ap, "--xperm", help = "whether to treat x separately; must be 'true','True','T','false','False' or 'F'.",default='T')
ap <- add_argument(ap, "--covs", help = "additive covariates to include (names of columns in phenotype data; note, these columns must be numeric and covariates are only supported for normal and binary models)", nargs = Inf, default = NA)
argv <- parse_args(ap)

# Load rqtl object.
if(str_to_upper(str_sub(argv$rqtl_object,start=-3))=='RDA'){
  dat<-get(load(argv$rqtl_object)) # this renames it as 'dat'
}
if(str_to_upper(str_sub(argv$rqtl_object,start=-5))=='RDATA'){ #
  dat<-get(load(argv$rqtl_object)) # this renames it as 'dat'
}
if(str_to_upper(str_sub(argv$rqtl_object,start=-3))=='RDS'){
  dat<-readRDS(argv$rqtl_object)
}

# Set varname
varname = argv$var

# Set savefolder
savefold <- argv$save_folder

# Set covariates
if(is.na(argv$covs)){
  covs = NULL
}
if(!is.na(argv$covs)){
  covs = cbind(dat$pheno[argv$covs])
}


# Set nperms
nperm = as.numeric(as.character(argv$nperm))

# Set Xperm
if(argv$xperm %in% c('T','true','True')){
  xperm=T
} else if(argv$xperm %in% c('F','false','False')){
  xperm=F
} else{
  print(paste0('xperm is ',argv$xperm, ' which is not allowed. Going with default: xperm is T'))
  xperm=T
}

# Report the plan.
writeLines(paste0("Performing QTL mapping using the file '",
             argv$rqtl_object,
             "', the variable '",
             argv$var,
             "', with covariates '",
            argv$covs,
            "', the model '",
             argv$model, 
             "', and the method '",
             argv$method,
             "'. Finding LOD threshold with nperms '",
             nperm,
             "'. Xsp is '",
             xperm,
             "'. Saving output in the folder '",
             argv$save_folder),
           file.path(savefold,paste0('configuration_',varname,'.txt')))

# Plot the trait distribution.
pl <- ggplot(dat$pheno)+theme_classic()+geom_histogram()+
  aes_string(x=varname)+ggtitle(varname)
ggsave(file.path(savefold,paste0('histogram_',varname,'.pdf')),plot=pl)

# Perform QTL map. Note, if model is np, method is ignored.
out <- scanone(dat,pheno.col=varname,model=argv$model,method=argv$method,addcovar = covs)
saveRDS(out,file=file.path(savefold,paste0('qtl_map_',varname,'.RDS')))

# Save the plot with no threshold.
pdf(file.path(savefold,paste0('qtl_plot_',varname,'.pdf')),width=10,height=6)
plot(out,main=varname,bandcol='aliceblue')
dev.off()


# Perform permutations.
perms <- scanone(dat,pheno.col=varname, model=argv$model, method=argv$method,
                 n.perm = nperm, perm.Xsp = xperm, addcovar = covs)
saveRDS(perms,file=file.path(savefold,paste0('perms_',varname,'.RDS')))

# Save the plot.
pdf(file.path(savefold,paste0('qtl_plot_',varname,'.pdf')),width=10,height=6)
plot(out,main=varname,bandcol='aliceblue')
add.threshold(out,perms=perms,alpha=0.05)
dev.off()

# Save the summary
peaks<-summary(out,perms=perms,alpha=0.05)
write.csv(peaks,file.path(savefold,paste0('peaks_',varname,'.csv')),row.names=T)

# Intervals
print('Creating interval information...')
qtlmap <- out
intervaldf <- data.frame(chr = factor(), lod = double(), pos = double())
intervaldf_bayes <- intervaldf
intervaldf_lod <- intervaldf
for(chr in peaks$chr){
  intervaldf_bayes<-rbind(intervaldf_bayes,bayesint(qtlmap,chr,0.95))
  intervaldf_lod<-rbind(intervaldf_lod,lodint(qtlmap,chr,drop=1.8))
}
if((dim(intervaldf_bayes))[1]>0){
  intervaldf_bayes$interval_type = 'bayes_95'
  intervaldf_lod$interval_type = 'lod_1.8'
  intervaldf = rbind(intervaldf_bayes,intervaldf_lod)
  
  print('saving interval data...')
  
  
  write.csv(intervaldf,file.path(savefold,paste0('intervals_',varname,'.csv')),row.names=T)
}


# plot and save intervals
print('Making interval plots...')
for(chr in peaks$chr){
  pdf(file=file.path(savefold,paste0('chr_',chr,'_with_bayes_interval.pdf')))
  plot(qtlmap,chr=chr)
  add.threshold(qtlmap,chr,perms,alpha=0.05)
  rect(xleft=intervaldf[intervaldf$interval_type=='bayes_95' & intervaldf$chr==chr,][1,'pos'],xright=intervaldf[intervaldf$interval_type=='bayes_95' & intervaldf$chr==chr,][3,'pos'],ybottom=0,ytop=100,col=rgb(0,0.5,1,0.1),border='transparent')
  dev.off()
  
  pdf(file=file.path(savefold,paste0('chr_',chr,'_with_1.8lod_interval.pdf')))
  plot(qtlmap,chr=chr)
  add.threshold(qtlmap,chr,perms,alpha=0.05)
  rect(xleft=intervaldf[intervaldf$interval_type=='lod_1.8' & intervaldf$chr==chr,][1,'pos'],xright=intervaldf[intervaldf$interval_type=='lod_1.8' & intervaldf$chr==chr,][3,'pos'],ybottom=0,ytop=100,col=rgb(0,0.5,1,0.1),border='transparent')
  dev.off()
}


# Save effect plots
print('Making effect plots...')
for(chr in peaks$chr){
  pdf(file=file.path(savefold,paste0('chr_',chr,'_pxg_',varname,'.pdf')))
  plotPXG(dat,pheno.col = varname, marker = row.names(peaks[peaks$chr==chr,]),infer=F)
  dev.off()
}

