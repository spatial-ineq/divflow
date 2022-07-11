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

sfg.dir <-
  '/projects/SHARKEY/safegraph/processed/orig_dest_annual/'


# Della fcn ---------------------------------------------------------------

#' Della.wrapper_flow.weights_by.rt
#'
#' Wrapper function to generate all flow weights by either CZ or CBSA.
#'
#' @inheritParams get.flow.weights.within.distance
#' @param list.colms whether to return list columns (or tibble colms)
#'
Della.wrapper_flow.weights_by.rt <- function(
  cz_id = NULL
  ,cbsa_id = NULL
  ,agg2tracts = T
  ,drop.loops = T
  ,weight.floor = .001
  ,sfg.dir
  ,year = '2019'
  ,save.dir = NULL
  #,list.colms = T
) {

  if(is.null(cbsa_id))
    czs2load <- cz_id
  else
    czs2load <- geox::rx %>%
      filter(cbsa %in% cbsa_id) %>%
      pull(cz) %>% unique()

  # sfg load
  sfg <- sfg.seg::read.sfg.CZs(czs2load,
                               sfg.dir = sfg.dir,
                               year = year)

  # subset to trips w/in region (cz/cbsa)
  sfg <- geox::geosubset(sfg
                         ,subset.cols = c('origin', 'dest')
                         ,cz = cz_id
                         ,cbsa = cbsa_id)

  # agg2 tracts
  if(agg2tracts)
    sfg <- sfg.seg::cbg.flows2tracts(sfg)

  # drop loops if appropriate
  if(drop.loops)
    sfg <- sfg %>% filter(origin != dest)


  # get inc and visited weights
  inc.flwws <- sfg %>%
    group_by(dest) %>%
    mutate(flww = n / sum(n)) %>%
    ungroup() %>%
    filter(flww >= weight.floor) %>%
    select(geoid = dest, j = origin, flww) %>%
    nest(inc.flwws = c(j, flww))

  vis.flwws <- sfg %>%
    group_by(origin) %>%
    mutate(flww = n / sum(n)) %>%
    ungroup() %>%
    filter(flww >= weight.floor) %>%
    select(geoid = origin, j = dest, flww) %>%
    nest(vis.flwws = c(j, flww))

  flwws <- full_join(inc.flwws, vis.flwws)

  #browser()

  if(!is.null(save.dir)) {
    if(!dir.exists(save.dir))
      dir.create(save.dir, recursive = T)

    rids <- geox::get.region.identifiers(cz = cz_id
                                         ,cbsa = cbsa_id) %>%
      paste0(collapse = '-') %>%
      gsub('\\/', '-', .)

    write_rds(flwws
              ,file = paste0(save.dir
                             ,rids
                             ,'-flow-weights.rds'))
  }

  return(flwws)
}


# test run ----------------------------------------------------------------
"
devtools::load_all()
Della.wrapper_flow.weights_by.rt

phl <- sfg.seg::read.sfg.CZs('19700'
                             ,sfg.dir)
geox::rpops %>% filter(rt=='cbsa') %>% geox::add.rns()
tmp <- Della.wrapper_flow.weights_by.rt(cz_id = '10100'
                                        ,agg2tracts = F
                                        ,drop.loops = T
                                        ,sfg.dir = sfg.dir
                                        ,weight.floor = .01
                                        ,year = '2019'
                                        ,save.dir = '/scratch/gpfs/km31/tests/flww-by-region/')


tmpcbsa <- Della.wrapper_flow.weights_by.rt(cbsa_id = '31140'
                                        ,agg2tracts = F
                                        ,drop.loops = T
                                        ,sfg.dir = sfg.dir
                                        ,weight.floor = .01
                                        ,year = '2019'
                                        ,save.dir = '/scratch/gpfs/km31/tests/flww-by-region/')

tsts <- list.files('/scratch/gpfs/km31/tests/flww-by-region/'
                   ,full.names = T)
tmpcbsa
tmp
tmp
# saved 2disk version
tmp2 <- tsts[2] %>% read_rds()
tmp2$geoid %>% head() %>% nchar()

tmp2
"

# weight floor descriptives -----------------------------------------------

# for a sample BG in Philly, a weight floor of .01 filters out 727/732
# connections (5 remaining with a share of >1%). Those remaining connections
# comprise ~33% of visits. With no weight floor, flww quartiles were:
# .0000 0.0005 0.0007 0.0013 0.0443


# gen lists ---------------------------------------------------------------

genczs <- geox::rpops %>% filter(pop > 25e3) %>% filter(rt == 'cz') %>% pull(rid)
gencbsas <- geox::rpops %>% filter(pop > 25e3)  %>% filter(rt == 'cbsa') %>% pull(rid)

# gen for tracts --------------------------------------------------------------

# within-distance tract list

## params for job ##

# save dir
save.dir <-
  paste0(prx.dir, 'flow-weights/cts-by-region/')
# '~/R/all sharkey geoseg work/divflow/R/creating SpWs/spWs/'
save.dir

## tract x CZ
tract.params <-
  tibble(
    cz_id = genczs
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
                                   'mem-per-cpu' = '20G',
                                   'mail-type' = list('begin', 'end', 'fail'),
                                   'mail-user' = 'km31@princeton.edu')
  )


## tract x CBSA

tract.params.cbsas <-
  tibble(
    cz_id = NULL
    ,cbsa_id = gencbsas
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
              params = tract.params.cbsas,
              jobname = 'tract flow weights by cbsa',
              nodes = 19,
              cpus_per_node = 1,
              slurm_options = list(time = '10:00:00',
                                   'mem-per-cpu' = '20G',
                                   'mail-type' = list('begin', 'end', 'fail'),
                                   'mail-user' = 'km31@princeton.edu')
  )


# check generation --------------------------------------------------------
# (also i should move to tract-specific directory)
save.dir
fns <- list.files(save.dir
                  ,full.names = T
                  ,pattern = 'rds$')

smpl <- fns[1:3] %>%
  map_dfr(read_rds)

smpl$inc.flwws[[1]]
geox::rpops %>% filter(pop > 25e3) %>% count(rt)
fns[grepl('cbsa', fns)] %>% length()
fns[grepl('cz', fns)] %>% length()
# looks great and complete!
ctfw$geoid %>% duplicated() %>% sum()

# change save directory ---------------------------------------------------
"save.dir
ct.save.dir <- paste0(save.dir
                      ,'cts-20mi-limit/')
dir.create(ct.save.dir)
fns <- fns %>% gsub('//', '/', .)

fn.dests <-
  fns %>% gsub(save.dir,ct.save.dir, .
               ,fixed = T)

file.rename(fns, fn.dests)"

# bg flow weights ---------------------------------------------------------

# save dir
save.dir <-
  paste0(prx.dir, 'flow-weights/bgs-by-region/')

# within-distance tract list
prx.dir %>% list.files()


# params for job
bg.params <-
  tibble(
    cz_id = genczs
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
                                   'mem-per-cpu' = '20G',
                                   'mail-type' = list('begin', 'end', 'fail'),
                                   'mail-user' = 'km31@princeton.edu')
  )


## by cbsa ##

bg.params <-
  tibble(
    cz_id = NULL
    ,cbsa_id = gencbsas
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
                                   'mem-per-cpu' = '20G',
                                   'mail-type' = list('begin', 'end', 'fail'),
                                   'mail-user' = 'km31@princeton.edu')
  )



# checking genration ------------------------------------------------------



# save dir
save.dir <-
  paste0(prx.dir, 'flow-weights/bgs-by-region/')

fns <- list.files(save.dir) %>%
  str_extract('[a-z]+-[0-9]{5}')

bgczs <- fns[grepl('cz', fns)]
bgczs <- paste0(save.dir,bgczs) %>% map_dfr(read_rds)
bgczs %>% nrow()
bgczs$geoid[1:5] %>% nchar()
bgczs[1,]$inc.flwws[[1]] %>% filter( j %in% '470190701002')

