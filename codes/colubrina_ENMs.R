##########  Ensemble ENMs for 3 Laticauda species
##########  Ensemble ENMs for L. colubrina

# clear working env
rm(list = ls(all.names = T))
gc()

# load packages
library(biomod2)
library(terra)


#####  part 1 ::: input data prep ----------

### occurrence data
c.occs <- read.csv('data/occs/colubrina.csv')
c.occs$pa <- rep(1, nrow(c.occs)) # define coordinates as presences
head(c.occs)

### present environmental data
p.envs <- rast(list.files(path = 'data/envs/present/', pattern = '.asc$', full.names = T))
print(p.envs)
plot(p.envs[[1]])

### future environmental data
f.envs <- rast(list.files(path = 'data/envs/SSP585-2090-2100/', pattern = '.asc$', full.names = T))
print(f.envs)
plot(f.envs[[1]])

### format the dataset // sample pseudoabsences here // random selection of 10 sets
c.data.bm <- BIOMOD_FormatingData(resp.name = 'colubrina',
                                  resp.var = vect(c.occs, geom = c('Long', 'Lat'), crs = 'EPSG:4326'),
                                  expl.var = p.envs,
                                  dir.name = 'outputs/models/colubrina/',
                                  PA.nb.rep = 10,
                                  PA.nb.absences = 10000,
                                  PA.strategy = 'random',
                                  na.rm = T)

# check formatted data
summary(c.data.bm)
