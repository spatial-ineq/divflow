# ws ---------------------------------------------------------------------------
rm(list = ls())
require(tidyverse)
require(sf)
devtools::load_all()

# option setting
sf_use_s2(T)
options(tigris_use_cache = TRUE)
Sys.setenv("VROOM_SHOW_PROGRESS"="false")



# set directories ---------------------------------------------------------


prx.dir <-
  '/scratch/gpfs/km31/adjacencies+proximities/'

save.dir <-
  paste0(prx.dir, 'flow-weights/')
  #'~/R/all sharkey geoseg work/divflow/R/creating SpWs/spWs/'


sfg.dir <-
  '/projects/SHARKEY/safegraph/processed/orig_dest_annual/'



# gen for tracts --------------------------------------------------------------

czs <- geox::rx$cz %>% unique()
tract.params <-
  tibble(
     cz = czs
    ,agg2tracts = T
    ,weight.floor = 0.001
    ,drop.loops = T
    ,sfg.dir = sfg.dir
    ,year = '2019'
    ,save.dir = save.dir
  )

# within-distance list

# I am still generating code for BGs so will use tracts for now
prx.dir %>% list.files()

spws <- prx.dir %>%
  list.files(pattern = 'tract-adjacencies.rds$'
             ,full.names = T) %>%
  read_rds()
spws <- spws %>% select(geoid, below.cutoff)

# tracts test run
Della.wrapper_flow.weights(
  cz = czs[2]
  ,agg2tracts = T
  ,weight.floor = 0.001
  ,drop.loops = T
  ,sfg.dir = sfg.dir
  ,year = '2019'
  ,save.dir = '/projects/SHARKEY/flww-tests/'
)
