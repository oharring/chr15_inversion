#!/usr/bin/Rscript

# Goal: run HZAR cline fitting for each phenotype (and soil trait), in the cascades-only and full dataset for the cline. 

library(hzar)
library(coda)
library(ggplot2)
library(argparser)

ap <- arg_parser('Run MCMC fit for clines with HZAR.')
ap <- add_argument(ap, "--var", help = "text name of the column corresponding to the trait to use for the cline.")
ap <- add_argument(ap, "--dtype", help = "data type for cline: must be soil or mouse",default='mouse')
ap <- add_argument(ap, "--region", help = "region to get data from; must be full (whole transect) or mid (cascades only)",default='full')
ap <- add_argument(ap, "--chain_length", help = "chain length for MCMC (burnin is chain_length/10, thin is chain_length/1000)",default=1e5)
ap <- add_argument(ap, '--tune', help = 'tune parameter for MCMCmetrop1R, applied to all parms. Default is 1.5 (in hzar)', default = NA)
ap <- add_argument(ap, '--subfolder','-f', help = 'sub folder to put results, inside analysis/hzar_fits; defaults to todays date', default = format(Sys.Date(),"%Y%m%d"))
argv <- parse_args(ap)

chain_length <- argv$chain_length
pheno <- argv$var
dtype <- argv$dtype
region <- argv$region
tune <- as.numeric(as.character(argv$tune))
mid_folder = argv$subfolder




##### FUNCTION #####

run_hzar_pheno <- function(df, dtype, distdf, pheno.col, model_tails, chain_length, tune = NA){
  
  print(paste0('Running HZAR cline fit for ',pheno.col,'. Fitting ',length(model_tails), ' models, with tails ',
               paste0(model_tails, collapse = ' and '), '. Chain length will be ',chain_length, ' with burnin ',chain_length/10,
               ' and thin ',chain_length/1e3,'. Performing 3 runs of each of 3 independent chains for each model.'))
  
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
  if(dtype=='mouse'){
    output$obs <- hzar.doNormalData1DRaw(hzar.mapSiteDist(distdf$Site, distdf$dist_east),
                                           df$Site, 
                                           df[[pheno.col]])
  }
  if(dtype=='soil'){
    print('Working with pre-summarized data; ignoring distdf.')
    output$obs <- hzar.doNormalData1DPops(distance = df$dist_east,
                                          siteID = df$Site, 
                                          muObs = df[[pheno.col]],
                                          nEff = df$nEff,
                                          varObs = df$var)
    
  }
  if(!(dtype %in% c('mouse','soil'))){
    print('Error: dtype must be "mouse" or "soil".')
  }
  
  # Set up the model(s): 
  print('Setting up models')
  for(tail_setting in model_tails){
    output$models[[paste0('model.tails.',tail_setting)]] <- hzar.makeCline1DNormal(output$obs, tail_setting)
  }
  
  # Focus the MCMC sampling in the experimental region. 
  # This sets upper and lower bounds on the cline center, and an upper bound on cline width.
  output$models <- lapply(output$models,
                              hzar.model.addBoxReq, 
                              low = min(distdf$dist_east)-10,
                              high = max(distdf$dist_east)+10)
  
  # Place max value on varH, varL and varR, using the total phenotypic variance
  output$models <- lapply(output$models,
                          hzar.model.addMaxVariance,
                          maxValue = 1.5 * var(df[[pheno.col]],na.rm=T)
                          )
  
  # Update tune parameter if requested
  if(!is.na(tune)){
    print(paste0('Setting tune: ',tune))
    for(m in names(output$models)){
      hzar.meta.tune(output$models[[m]]) <- tune
    }
  }
  
  # Make first fit request; this sets up the parameters for the MCMC chain.
  output$fitRs$init <- lapply(output$models,
                                  hzar.first.fitRequest.gC,
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
    model.name = paste0('model.tails.',tail_setting)
    output$analysis$initDGs[[model.name]] <- hzar.dataGroup.add(output$runs$init[[model.name]])
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
source_dir = '/n/holylfs03/LABS/hoekstra_lab/Users/hager/papers/cline_paper'

# Load the distance data:
distdf <- read.csv(paste0(source_dir,'/data/transect_distance.csv'))

# Load the mouse data:
mousedf <- read.csv(paste0(source_dir,'/data/adult_transect_data.csv'))
mousedf_mid <- droplevels(subset(mousedf, cascades == TRUE))

# Load the soil data (already subset to 1 km, see set_up_soil_1D_pops.R':
soildf <- read.csv(paste0(source_dir,'/data/soil_summary_1DPops.csv'))

if(region=='mid'){
  distdat <- subset(distdf, Site %in% levels(mousedf_mid$Site))
  model_tails = c('none','both','right','left','mirror')
  if(dtype == 'soil'){
    dat <- subset(soildf, Site %in% levels(mousedf_mid$Site))
  }
  if(dtype == 'mouse'){
    dat <- mousedf_mid
  }
}
if(region=='full'){
  distdat <- distdf
  model_tails = c('none','both','right','left','mirror')
  if(dtype=='soil'){
    dat <- soildf
  }
  if(dtype=='mouse'){
    dat <- mousedf
  }
}

# Run HZAR
mcmc_output <- run_hzar_pheno(dat, dtype, distdat, pheno, model_tails, chain_length, tune)

# Save output
dir.create(paste0(source_dir,'/analysis/hzar_fits/',mid_folder))
saveRDS(mcmc_output,file=paste0(source_dir, '/analysis/hzar_fits/',mid_folder,'/mcmc_results_',dtype,'_',pheno,'_',region,'.RDS'), version = 2)

# Make some plots
# 1. MCMC output of final chain for each run of each model type
for(mod in 1:length(model_tails)){
  png(paste0(source_dir,'/analysis/hzar_fits/',mid_folder,'/chains_',dtype,'_',pheno,'_',region,'_model_',model_tails[mod],'.png'),width = 720, height = 1280)
  plot(do.call(mcmc.list,
               lapply(mcmc_output$runs$chains[(3*mod-2):(3*mod)], 
                      function(x) hzar.mcmc.bindLL(x[[3]]))))
  dev.off()
}

# 2. fZCline best fit model from hzar.
png(paste0(source_dir,'/analysis/hzar_fits/',mid_folder,'/fZCline_',dtype,'_',pheno,'_',region,'.png'), width = 600, height = 480)
hzar.plot.fzCline(mcmc_output$analysis$model.selected)
dev.off()

# 3. Autocorrelation and effective sample size for the final chain for each run
# of each model type.
for(mod in 1:length(model_tails)){
    A = autocorr(do.call(mcmc.list,
        lapply(mcmc_output$runs$chains[(3*mod-2):(3*mod)], 
               function(x) hzar.mcmc.bindLL(x[[3]]))),lags = 1, relative = TRUE)
    for(i in seq(1,3,1)){
      #print(dim(A[[i]][1,,]))
      write.csv(A[[i]][1,,],file = paste0(source_dir,'/analysis/hzar_fits/',mid_folder,'/autocorr_lag1rel_',dtype,'_',pheno,'_',region,'_model_',model_tails[mod],'_chain_',i,'.csv'),row.names=T)
    }
    print(paste0('Max autocorrelation, model ',model_tails[mod],' : ',
                 max(max(A[[1]]),max(A[[2]]),max(A[[3]]))))
}

dn = list()
for(mod in 1:length(model_tails)){
  n = effectiveSize(do.call(mcmc.list,
                          lapply(mcmc_output$runs$chains[(3*mod-2):(3*mod)], 
                                 function(x) hzar.mcmc.bindLL(x[[3]]))))
  print(paste0('Min effective sample for model ',model_tails[mod],' : ',min(n)))
  dn[model_tails[mod]]=list(n)
}
saveRDS(dn,file = paste0(source_dir,'/analysis/hzar_fits/',mid_folder,'/effectiveSize_',dtype,'_',pheno,'_',region,'_all_models.RDS'),version=2)
    


#pheno_list = c('tail_mm.museum','foot_mm.museum','body_mm.museum',
#               'D.Hue.427','D.Brightness.427','D.Saturation',
#               'F.Hue.427','F.Brightness.427','F.Saturation',
#               'tail','body')

#soil_trait_list = c('mean_Munsell_hue','mean_Munsell_value','mean_Munsell_chroma')





