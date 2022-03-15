#!/usr/bin/Rscript

# Goal: run HZAR cline fitting for each genotype, in the cascades-only and full dataset for the cline. 

library(hzar)
library(ggplot2)
library(argparser)

ap <- arg_parser('Run MCMC fit for clines with HZAR.')
ap <- add_argument(ap, "--var", help = "text name of the column corresponding to the trait to use for the cline.")
ap <- add_argument(ap, "--dtype", help = "data type for cline: must be soil or mouse",default='mouse')
ap <- add_argument(ap, "--region", help = "region to get data from; must be full (whole transect) or mid (cascades only)",default='full')
ap <- add_argument(ap, "--chain_length", help = "chain length for MCMC (burnin is chain_length/10, thin is chain_length/1000)",default=1e6)
argv <- parse_args(ap)

chain_length <- argv$chain_length
pheno <- argv$var
dtype <- argv$dtype
region <- argv$region


##### FUNCTION #####

run_hzar_pheno <- function(df, pheno.col, model_tails, chain_length, model_scales){ 
  
  print(paste0('Running HZAR cline fit for ',pheno.col,'. Fitting ',length(model_tails), ' models, with tails ',
               paste0(model_tails, collapse = ' and '),', and fitting ',length(model_scales),' models with scaling ',paste0(model_scales,collapse= ' and '), '. Chain length will be ',chain_length, ' with burnin ',chain_length/10,' and thin ',chain_length/1e3,'. Performing 3 runs of each of 3 independent chains for each model.'))
  
  # Sets thin to 1/1000 of chain length (ie save 1000 estimates)
  # Sets burnin to 1/10 of chain length
  
  # Set up place to hold output
  output <- list()
  output$obs <- list()
  output$models <- list()
  output$fitRs <- list()
  output$runs <- list()
  output$analysis <- list()
  
  # Set up observations
  print('Setting up observations')
#   output$obs <- hzar.doNormalData1DRaw(hzar.mapSiteDist(distdf$Site, distdf$dist_east)
#                                            df$Site, 
#                                            df[[pheno.col]])
  
  af_var_name0<-paste('SW',pheno.col,sep='_'); af_var_name<-paste(af_var_name0,'af',sep='_') #get name of column that holds allele frequency for locus of interest
  ntotal_name<-paste(pheno.col,'ntotal_alleles',sep='_') #get name of column that holds number of sampled alleles for locus of interest
  output$obs <-hzar.doMolecularData1DPops(df$dist_east,df[,af_var_name],df[,ntotal_name])
  
  # Set up the model(s): 
  print('Setting up models')
  for(tail_setting in model_tails){ 
    for(scale_setting in model_scales){
      output$models[[paste0('model.tails.',tail_setting,'.model.scales.',scale_setting)]] <- hzar.makeCline1DFreq(output$obs,scale_setting,tail_setting)
    }
  }

  
  # Focus the MCMC sampling in the experimental region. 
  # This sets upper and lower bounds on the cline center, and an upper bound on cline width.
  output$models <- lapply(output$models,
                              hzar.model.addBoxReq, 
                              low = min(df$dist_east)-10,
                              high = max(df$dist_east)+10)
  
  
  # Make first fit request; this sets up the parameters for the MCMC chain.
  output$fitRs$init <- lapply(output$models,
                              hzar.first.fitRequest.old.ML,
                                  obsData = output$obs)
  
  
  # Set seeds to be different and update MCMC parameters:
  print('Updating MCMC parameters')
  for(i in 1:length(output$fitRs$init)){
    output$fitRs$init[[i]]$mcmcParam$seed[[1]] <- sample(1e2:1e4,6) 
    output$fitRs$init[[i]]$mcmcParam$chainLength <- chain_length
    output$fitRs$init[[i]]$mcmcParam$thin <- chain_length/1e3 # obtain 1000 estimates from each chain
    output$fitRs$init[[i]]$mcmcParam$burnin <- chain_length/10 # burnin is 10% of chain length
  }
 
  # Use these initial fit requests as a base to generate
  # a total of three independent MCMC chains for each model.
  # The seeds are the same for each chain but the channel differentiates them.
  output$fitRs$init <- hzar.multiFitRequest(output$fitRs$init,
                                                each=3,
                                                baseSeed = NULL)
  
  # Do a first MCMC run for each of the chains.
  print('Running first chains')
  output$runs$init <- lapply(output$fitRs$init,
                                 hzar.doFit)
  
  # Set up next fit request for each chain.
  print('Setting up fit requests')
  
  output$fitRs$chains <- lapply(output$runs$init,
                                    hzar.next.fitRequest)
  
  
  # Repeat the cycle ( output of doFit --> next.fitRequest --> doFit etc)
  # 3x for each chain. This would be run for one chain with hzar.chain.doSeq.
  print('Running all chains')
  output$runs$chains <- hzar.doChain.multi(output$fitRs$chains,
                                               collapse = F,
                                               count = 3)
  
  # Perform analysis: group data from all runs for each model, pick best model from AICc
  print('Chains complete; running analysis')
  output$analysis$initDGs <- list()
  
  # Groups multiple fits of the same model into a single object
  for(tail_setting in model_tails){
    for(scale_setting in model_scales){
      model.name = paste0('model.tails.',tail_setting,'.model.scales.',scale_setting)
      output$analysis$initDGs[[model.name]] <- hzar.dataGroup.add(output$runs$init[[model.name]])
    }
    }
  
  # Collect data from multiple models fit to the same data into a single object with obsDataGroup.
  output$analysis$oDG <- hzar.make.obsDataGroup(output$analysis$initDGs)
  output$analysis$oDG <- hzar.copyModelLabels(output$analysis$initDGs, output$analysis$oDG)
  
  # Add all the runs to the obsDataGroup:
  output$analysis$oDG <- hzar.make.obsDataGroup(lapply(output$runs$chains, hzar.dataGroup.add), output$analysis$oDG)
  
  # Collect AICc into table
  output$analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(output$analysis$oDG)
  
  # Identify and save the best model
  output$analysis$model.name <- rownames(output$analysis$AICcTable)[[which.min(output$analysis$AICcTable$AICc)]]
  output$analysis$model.selected <- output$analysis$oDG$data.groups[[output$analysis$model.name]]
  
  output$analysis$model.LL.param.limits <- hzar.getLLCutParam(output$analysis$model.selected, names(output$analysis$model.selected$data.param))
  output$analysis$ML.cline <- hzar.get.ML.cline(output$analysis$model.selected)
  
  print('Analysis complete')
  return(output)
  
}


###### DATA #######
source_dir = ''

geno_data <- read.csv(paste0(source_dir,'/data/cline_allele_freqs_by_site_with_distance.csv'))

if(region=='mid'){
  dat <- geno_data[!geno_data$Site %in% c('Baker City A','Baker City B','Baker City C','Baker City D','Siuslaw A','Siuslaw B','Siuslaw C','Siuslaw D','Siuslaw E'),]
  model_tails = 'none'
  model_scales = c('fixed','free')
}
if(region=='full'){
  dat <- geno_data
  model_tails = c('none','both','right','left','mirror')
  model_scales = c('fixed','free')
}

# Run HZAR
mcmc_output <- run_hzar_pheno(dat, pheno, model_tails, chain_length, model_scales)

# Save output
saveRDS(mcmc_output,file=paste0(source_dir, '/analysis/hzar_fits/mcmc_results_',pheno,'_',region,'.RDS'), version = 2)

# Make some plots
# 1. MCMC output of final chain for each run of each model type
n_scales<-length(model_scales)
for(mod in 1:length(model_tails)){
  for(scl in 1:length(model_scales)){
  pdf(paste0(source_dir,'/analysis/hzar_fits/chains_',pheno,'_',region,'_model_',model_tails[mod],'_',model_scales[scl],'.pdf'),width = 8, height = 10.5)
  plot(do.call(mcmc.list,
	lapply(mcmc_output$runs$chains[(3*scl+3*n_scales*(mod-1)-2):(3*scl+3*n_scales*(mod-1))],
		function(x) hzar.mcmc.bindLL(x[[3]]))))
  dev.off()
}}

# 2. fZCline best fit model from hzar.
pdf(paste0(source_dir,'/analysis/hzar_fits/fZCline_',pheno,'_',region,'_',mcmc_output$analysis$model.name,'.pdf'), width = 5, height = 4)
hzar.plot.fzCline(mcmc_output$analysis$model.selected)
dev.off()



