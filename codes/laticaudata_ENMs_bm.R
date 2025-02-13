##########  Ensemble ENMs for 3 Laticauda species
##########  Ensemble ENMs for L. laticaudata

# clear working env
rm(list = ls(all.names = T))
gc()

# set seed for repeatability
set.seed(111)

# load packages
library(biomod2)
library(ecospat)
library(terra)
library(dplyr)


#####  part 1 ::: input data prep ----------

### occurrence data
l.occs <- read.csv('data/occs/laticaudata.csv') %>% select(1,3,2)
colnames(l.occs) = c('species', 'long', 'lat')
l.occs$pa <- rep(1, nrow(l.occs)) # define coordinates as presences
head(l.occs)


### present environmental data
p.envs <- rast(list.files(path = 'data/envs/geotiff//present/', pattern = '.tif$', full.names = T))
p.envs <- p.envs[[c('bathy', 'dtp', 'dts', 'slope', 'sst_lt.max', 'sst_range', 'sss_mean')]]
print(p.envs)
plot(p.envs[[1]])


### future environmental data == ssp245
f.envs1 <- rast(list.files(path = 'data/envs/geotiff/2090s_SSP245/', pattern = '.tif$', full.names = T))
f.envs1 <- f.envs1[[names(p.envs)]]
print(f.envs1)
plot(f.envs1[[1]])


### future environmental data == ssp585
f.envs2 <- rast(list.files(path = 'data/envs/geotiff/2090s_SSP585/', pattern = '.tif$', full.names = T))
f.envs2 <- f.envs2[[names(p.envs)]]
print(f.envs2)
plot(f.envs2[[1]])


### format the dataset // sample pseudoabsences here // random selection of 10 sets
l.data.bm <- BIOMOD_FormatingData(resp.name = 'laticaudata',
                                  resp.var = vect(l.occs[, c(2,3)], geom = c('long', 'lat'), crs = 'EPSG:4326'),
                                  expl.var = p.envs,
                                  dir.name = 'outputs/models/laticaudata_bm/',
                                  PA.nb.rep = 2,                                # 2 replicates
                                  PA.nb.absences = 2000,                        # 2000 pseudo-absences
                                  PA.strategy = 'random',
                                  na.rm = T)

# check formatted data
summary(l.data.bm)
plot(l.data.bm)


### cross-validation dataset prep // use random crossvalidation
l.cv.bm <- bm_CrossValidation(bm.format = l.data.bm,
                              strategy = 'random',
                              nb.rep = 2,                                       # 2 replicates
                              perc = 0.7,
                              do.full.models = T)



#####  part 2 ::: single models ----------

### use the below code if running a large number of replications and model tuning is computationally intractable == which is probably the case for our purposes
# run single models == use the pre-defined model parameters (the "bigboss" strategy) for computational tractability
# the individual tuning of each modeling algorithm (the "tuned" strategy) take way too long for the models to run
# tuning may not be doable with higher number of replications
l.mods.single_bb <- BIOMOD_Modeling(bm.format = l.data.bm,
                                    modeling.id = 'laticaudata_singles_bigboss',
                                    models = c('GBM','RF','MAXENT'),
                                    CV.strategy = 'random',
                                    CV.nb.rep = 2,                              # 2 replicates
                                    CV.perc = 0.7,
                                    CV.do.full.models = T,
                                    OPT.data.type = 'binary',
                                    OPT.strategy = 'bigboss',
                                    metric.eval = c('ROC','TSS','BOYCE'),
                                    seed.val = 123,
                                    do.progress = T)

# check model results
print(l.mods.single_bb)

# get model evaluation results
get_evaluations(l.mods.single_bb)

# look at mean response curves
bm_PlotResponseCurves(bm.out = l.mods.single_bb, models.chosen = 'all', fixed.var = 'mean')

### retrieve various evaluation metrics
# get observed data, and corresponding calibration lines, and predictions
tab <- get_species_data(l.data.bm)
calib <- get_calib_lines(l.mods.single_bb, as.data.frame = T)
pred <- get_predictions(l.mods.single_bb, model.as.col = T)

# get evaluation data
eval <- get_evaluations(l.mods.single_bb)

mod1 = 'laticaudata_allData_allRun_GBM'
mod2 = 'laticaudata_allData_allRun_RF'
mod3 = 'laticaudata_allData_allRun_MAXENT'

ref1 = eval[which(eval$full.name == mod1 & eval$metric.eval == 'TSS'), ]
ref2 = eval[which(eval$full.name == mod2 & eval$metric.eval == 'TSS'), ]
ref3 = eval[which(eval$full.name == mod3 & eval$metric.eval == 'TSS'), ]

ind1 = which(calib[, 1] == TRUE)
bm_FindOptimStat(metric.eval = 'TSS', obs = tab[ind1, 1], fit = pred[ind1, mod1])
bm_FindOptimStat(metric.eval = 'TSS', obs = tab[ind1, 1], fit = pred[ind1, mod1])
bm_FindOptimStat(metric.eval = 'TSS', obs = tab[ind1, 1], fit = pred[ind1, mod1])

ind2 = which(calib[, 1] == FALSE)
bm_FindOptimStat(metric.eval = 'TSS', obs = tab[ind2, 1], fit = pred[ind2, mod1], threshold = ref1$cutoff)


#####  part 3 ::: ensemble models ----------
### look at eval metrics to dicede on the cutoff value
eval.metrics <- get_evaluations(l.mods.single_bb)[, c('metric.eval', 'validation')]
print(eval.metrics)

# TSS
eval.tss <- eval.metrics %>% filter(metric.eval == 'TSS') %>% na.omit()
mean(eval.tss$validation)

# ROC
eval.roc <- eval.metrics %>% filter(metric.eval == 'ROC') %>% na.omit()
mean(eval.roc$validation)

# BOYCE
eval.boyce <- eval.metrics %>% filter(metric.eval == 'BOYCE') %>% na.omit()
mean(eval.boyce$validation)


### ensemble
l.mods.em <- BIOMOD_EnsembleModeling(bm.mod = l.mods.single_bb,
                                     models.chosen = 'all',
                                     em.by = 'all',
                                     em.algo = c('EMmean','EMmedian'),
                                     metric.select = c('TSS','ROC','BOYCE'),
                                     metric.select.thresh = c(mean(eval.tss$validation), mean(eval.roc$validation), mean(eval.boyce$validation)),
                                     var.import = 5,
                                     seed.val = 123,
                                     do.progress = T)

# save the ensemble object
saveRDS(l.mods.em, 'outputs/models/laticaudata_bm/laticaudata/rds_out/laticaudata_ensemble_out.rds')

### get evaluation data
eval.em <- get_evaluations(l.mods.em)

mod.em = 'laticaudata_EMmeanByTSS_mergedData_mergedRun_mergedAlgo'
ref.em = eval.em[which(eval.em$full.name == mod.em & eval.em$metric.eval == 'TSS'), ]


#####  part 4 ::: project ensemble models to current & future envs ----------

### project ensemble models to the current envs
l.em.proj.cur <- BIOMOD_EnsembleForecasting(bm.em = l.mods.em,
                                            bm.proj = NULL,
                                            proj.name = 'laticaudata_em_proj',
                                            new.env = p.envs,
                                            models.chosen = 'all',
                                            metric.binary = c('TSS','ROC','BOYCE'),
                                            metric.filter = c('TSS','ROC','BOYCE'))

# save the current ensemble projection object
saveRDS(l.em.proj.cur, 'outputs/models/laticaudata_bm/laticaudata/rds_out/laticaudata_ensemble_proj_current_out.rds')


### project ensemble models to the future envs (ssp245)
l.em.proj.fut1 <- BIOMOD_EnsembleForecasting(bm.em = l.mods.em,
                                             bm.proj = NULL,
                                             proj.name = 'laticaudata_em_proj_fut_245',
                                             new.env = f.envs1,
                                             models.chosen = 'all',
                                             metric.binary = c('TSS','ROC','BOYCE'),
                                             metric.filter = c('TSS','ROC','BOYCE'))

# save the ssp245 ensemble projection object
saveRDS(l.em.proj.fut1, 'outputs/models/laticaudata_bm/laticaudata/rds_out/laticaudata_ensemble_proj_future245_out.rds')


### project ensemble models to the future envs (ssp585)
l.em.proj.fut2 <- BIOMOD_EnsembleForecasting(bm.em = l.mods.em,
                                             bm.proj = NULL,
                                             proj.name = 'laticaudata_em_proj_fut_585',
                                             new.env = f.envs2,
                                             models.chosen = 'all',
                                             metric.binary = c('TSS','ROC','BOYCE'),
                                             metric.filter = c('TSS','ROC','BOYCE'))

# save the ssp585 ensemble projection object
saveRDS(l.em.proj.fut2, 'outputs/models/laticaudata_bm/laticaudata/rds_out/laticaudata_ensemble_proj_future585_out.rds')


### see model predictions == note that the scale is (predicted habitat suitability) * 1000 == thus on the scale of 0-1000
### simply divide this raster by 1000 to make the scale 0-1
### current
l.preds.cur <- rast(list.files(path = 'outputs/models/laticaudata_bm/laticaudata/proj_laticaudata_em_proj', pattern = '.tif$', full.names = T))
print(l.preds.cur)
plot(l.preds.cur[[1]]/1000)

# export current esnsemble raster
writeRaster(l.preds.cur[[1]]/1000, 'outputs/predictions/laticaudata/cont/laticaudata_current_ensembleProj.tif', overwrite = T)


### future == ssp 245
l.preds.fut1 <- rast(list.files(path = 'outputs/models/laticaudata_bm/laticaudata/proj_laticaudata_em_proj_fut_245/', pattern = '.tif$', full.names = T))
print(l.preds.fut1)
plot(l.preds.fut1[[1]]/1000)

# export future (ssp585) ensemble raster
writeRaster(l.preds.fut1[[1]]/1000, 'outputs/predictions/laticaudata/cont/laticaudata_future245_ensembleProj.tif', overwrite = T)


### future == ssp 585
l.preds.fut2 <- rast(list.files(path = 'outputs/models/laticaudata_bm/laticaudata/proj_laticaudata_em_proj_fut_585/', pattern = '.tif$', full.names = T))
print(l.preds.fut2)
plot(l.preds.fut2[[1]]/1000)

# export future (ssp585) ensemble raster
writeRaster(l.preds.fut2[[1]]/1000, 'outputs/predictions/laticaudata/cont/laticaudata_future585_ensembleProj.tif', overwrite = T)


#####  part 5 ::: binary rasters ----------
# function to calculate p10 == from == https://babichmorrowc.github.io/post/2019-04-12-sdm-threshold/
sdm_threshold <- function(sdm, occs, type = "mtp", binary = FALSE){
  occPredVals <- raster::extract(sdm, occs)
  if(type == "mtp"){
    thresh <- min(na.omit(occPredVals))
  } else if(type == "p10"){
    if(length(occPredVals) < 10){
      p10 <- floor(length(occPredVals) * 0.9)
    } else {
      p10 <- ceiling(length(occPredVals) * 0.9)
    }
    thresh <- rev(sort(occPredVals))[p10]
  }
  sdm_thresh <- sdm
  sdm_thresh[sdm_thresh < thresh] <- NA
  if(binary){
    sdm_thresh[sdm_thresh >= thresh] <- 1
  }
  return(sdm_thresh)
}


# calculate threshold from the current prediction
l.cur_th <- sdm_threshold(sdm = raster::raster(l.preds.cur[[1]]/1000), occs = l.occs[, c('long', 'lat')], type = 'p10', binary = F)
l.th <- raster::minValue(l.cur_th)

# binary transition
l.cur.bin <- ecospat.binary.model(Pred = l.preds.cur[[1]]/1000, Threshold = l.th)
l.fut.bin1 <- ecospat.binary.model(Pred = l.preds.fut1[[1]]/1000, Threshold = l.th)
l.fut.bin2 <- ecospat.binary.model(Pred = l.preds.fut2[[1]]/1000, Threshold = l.th)

# export binary
writeRaster(l.cur.bin, 'outputs/predictions/laticaudata/bin/laticaudata_current_ensembleProj_bin.tif', overwrite = T)
writeRaster(l.fut.bin1, 'outputs/predictions/laticaudata/bin/laticaudata_future245_ensembleProj_bin.tif', overwrite = T)
writeRaster(l.fut.bin2, 'outputs/predictions/laticaudata/bin/laticaudata_future585_ensembleProj_bin.tif', overwrite = T)


#####  part 6 ::: calculate range size change ----------
# calculate range size change from current to 245 & current to 585
l.cur_to_245 <- BIOMOD_RangeSize(proj.current = l.cur.bin, proj.future = l.fut.bin1)
l.cur_to_585 <- BIOMOD_RangeSize(proj.current = l.cur.bin, proj.future = l.fut.bin2)


