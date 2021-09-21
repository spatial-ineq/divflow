# ws ---------------------------------------------------------------------------
rm(list = ls())
require(tidyverse)
require(sf)
# devtools::load_all()

# option setting
sf_use_s2(T)
options(tigris_use_cache = TRUE)
Sys.setenv("VROOM_SHOW_PROGRESS"="false")



# set directories ---------------------------------------------------------

prx.dir <-
  '/scratch/gpfs/km31/adjacencies+proximities/'

sfg.dir <-
  '/projects/SHARKEY/safegraph/processed/orig_dest_annual/'



# test run ----------------------------------------------------------------

devtools::load_all()
Della.wrapper_flow.weights_by.rt

tmp <- Della.wrapper_flow.weights_by.rt(cz_id = '19700'
                                        ,agg2tracts = T
                                        ,drop.loops = T
                                        ,sfg.dir = sfg.dir
                                        ,year = '2019'
                                        ,save.dir = '/scratch/gpfs/km31/tests/flww-by-region/')
tsts <- list.files('/scratch/gpfs/km31/tests/flww-by-region/'
                   ,full.names = T)

# gen for tracts --------------------------------------------------------------

# within-distance tract list

## params for job ##

# save dir
save.dir <-
  paste0(prx.dir, 'flow-weights/by-region/')
# '~/R/all sharkey geoseg work/divflow/R/creating SpWs/spWs/'
save.dir

czs <- geox::rx$cz %>% unique()

tract.params <-
  tibble(
    cz_id = czs
    ,cbsa_id = NULL
    ,agg2tracts = T
    ,weight.floor = 0.001
    ,drop.loops = T
    ,sfg.dir = sfg.dir
    ,year = '2019'
    ,save.dir = save.dir
  )


# send job
library(rslurm)
tract.flwws.dellajob <-
  slurm_apply(f =
                Della.wrapper_flow.weights_by.rt,
              params = tract.params,
              jobname = 'tract flow weights by cz',
              nodes = 19,
              cpus_per_node = 1,
              slurm_options = list(time = '10:00:00',
                                   'mem-per-cpu' = '30G',
                                   'mail-type' = list('begin', 'end', 'fail'),
                                   'mail-user' = 'km31@princeton.edu')
  )


# by cbsa
cbsas <- geox::rx$cbsa %>% unique()

tract.params.cbsas <-
  tibble(
    cz_id = NULL
    ,cbsa_id = cbsas
    ,agg2tracts = T
    ,weight.floor = 0.001
    ,drop.loops = T
    ,sfg.dir = sfg.dir
    ,year = '2019'
    ,save.dir = save.dir
  )


# check generation --------------------------------------------------------

# (also i should move to tract-specific directory)
save.dir
fns <- list.files(save.dir
                  ,full.names = T
                  ,pattern = 'rds$')

smpl <- fns[1:2] %>%
  map_dfr(read_rds)

ctfw <- fns %>% map_dfr(read_rds)

# looks great and complete!
ctfw$geoid %>% duplicated() %>% sum()

# change save directory ---------------------------------------------------
save.dir
ct.save.dir <- paste0(save.dir
                      ,'cts-20mi-limit/')
dir.create(ct.save.dir)
fns <- fns %>% gsub('//', '/', .)

fn.dests <-
  fns %>% gsub(save.dir,ct.save.dir, .
               ,fixed = T)

file.rename(fns, fn.dests)

# bg flow weights ---------------------------------------------------------

# save dir
save.dir <-
  paste0(prx.dir, 'flow-weights/bgs-by-region/')

# within-distance tract list
prx.dir %>% list.files()


# params for job
czs <- geox::rx$cz %>% unique()

bg.params <-
  tibble(
    cz_id = czs
    ,cbsa_id = NULL
    ,agg2tracts = F
    ,weight.floor = 0.001
    ,drop.loops = T
    ,sfg.dir = sfg.dir
    ,year = '2019'
    ,save.dir = save.dir
  )


# send job
library(rslurm)
cbg.flwws.dellajob <-
  slurm_apply(f =
                Della.wrapper_flow.weights_by.rt,
              params = bg.params,
              jobname = 'cbg flow weights by cz',
              nodes = 19,
              cpus_per_node = 1,
              slurm_options = list(time = '10:00:00',
                                   'mem-per-cpu' = '40G',
                                   'mail-type' = list('begin', 'end', 'fail'),
                                   'mail-user' = 'km31@princeton.edu')
  )


## by cbsa ##
cbsas <- geox::rx$cbsa %>% unique()

bg.params <-
  tibble(
    cz_id = NULL
    ,cbsa_id = cbsas
    ,agg2tracts = F
    ,weight.floor = 0.001
    ,drop.loops = T
    ,sfg.dir = sfg.dir
    ,year = '2019'
    ,save.dir = save.dir
  )


# send job
library(rslurm)
cbg.flwws.dellajob <-
  slurm_apply(f =
                Della.wrapper_flow.weights_by.rt,
              params = bg.params,
              jobname = 'cbg flow weights by cbsa',
              nodes = 19,
              cpus_per_node = 1,
              slurm_options = list(time = '10:00:00',
                                   'mem-per-cpu' = '40G',
                                   'mail-type' = list('begin', 'end', 'fail'),
                                   'mail-user' = 'km31@princeton.edu')
  )


