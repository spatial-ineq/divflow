# ws ---------------------------------------------------------------------------
rm(list = ls())
require(tidyverse)
require(sf)

# option setting
sf_use_s2(T)
options(tigris_use_cache = TRUE)

# dropbox dir or della base dir
ddir <- # Sys.getenv('drop_dir')
  '/scratch/gpfs/km31/'

# notes ------------------------------------------------------------------------

#' the prepped dataset "full-spatial-composites-by-rt" (generated through della
#' scripts) has everything in percents. I need counts and regional totals for flows
#' and res, so have to calcute from tract-level residence- & flow- & spatial-
#' weights.



# load spatial weights ------------------------------------------------------


spw.dir <- paste0(ddir, 'adjacencies+proximities/')
spwfn <- spw.dir %>% list.files(pattern = 'tract-adjacencies.rds')
spws <- paste0(spw.dir, spwfn) %>% readRDS()

spws

# make a tible( it should already be a tibble...)
spws <- spws %>%
  mutate(dists =
           map(dists,
               ~tibble(
                  j = names(.)
                 ,dist = .
               ))
  )

# get flow/res --------------------------------------------------------------

# which has flow totals

flrdir <- paste0(ddir
                 ,'seg-measures/by tract/broader ineq flows/')

flrdir %>% list.files()

flr <- flrdir %>%
  list.files(pattern = 'flow-res'
             ,full.names = T) %>%
  vroom::vroom()

flr <- flr %>%
  geox::abv.rcols()


tots <- flr %>%
  select(geoid, rt, rid,
         matches('^total'))

# tots <- tots %>% select(geoid, matches('^total'))


# get long res vars -----------------------------------------------------------

# also get all residential values, long by variable
resl <- flrdir %>%
  list.files(pattern = 'res-chars-long'
             ,full.names = T) %>%
  vroom::vroom()

resl


# proximity counts --------------------------------------------------------


devtools::load_all()

#' calc.proxim.counts
#'
#'
calc.proxim.counts <- function(cz = NULL, cbsa = NULL
                               #spatial.weights = spws,
                               #long.res.chars = trimmed.resl
                               ,save.csv = F) {

  require(tidyverse)

  if(! (exists('trimmed.resl') & exists('spws')) )
    warning('Wrapper expects trimmed.resl and spws in global env')

  # local decay fcns
  bisq.dist2weights <- function(d
                                , cutoff = dst.ceiling
                                , ... ) {

    if_else( d >= cutoff
             ,0
             ,(1-as.numeric(d/cutoff)^2)^2)

  }
  neg.exp <- function(d, ...) {
    exp(-d)
  }

  # trim to region
  lw <- spws %>%
    geox::geosubset(cz = cz, cbsa = cbsa)

  #  browser()

  # unnest, apply decay fcns
  lw <- lw %>%
    unnest(dists) %>%
    mutate( negexp.decay =
              neg.exp( as.numeric(dist) )
            ,bisq.8km.decay =
              bisq.dist2weights( as.numeric(dist), cutoff = 8 )
    )

  # browser()

  # process: get distance-from-i-weighted avg pop/hh of neighboring tracts. Then
  # multiply dist/pop weights by value in those tracts, and sum over i's.

  # merge to get % value and total pop/hh at each j
  js <- lw %>%
    left_join(trimmed.resl
              , by = c('j' = 'geoid'))

  # remove where i == j, and trim to js in region
  js <- js %>%
    select(i = geoid, j,
           matches('dist|decay'),
           var, value, weight
    ) %>%
    filter(i != j ) %>%
    geox::geosubset(subset.cols = 'j',
                    cz = cz, cbsa = cbsa)

  # row standardise decays
  js <- js %>%
    group_by(i, var) %>%
    mutate( across(matches('decay$')
                   , ~ .x / sum(.x)
    ))

  # "weight" colm is either population or hhs, depending on the variable
  js <-
    js %>%
    mutate( across(matches('decay$')
                   , ~ .x * weight * value
    ))


  n.dsts <- js %>%
    group_by(i, var) %>%
    summarise(across(matches('decay$')
                     , list(n = ~sum(.x, na.rm = T) )
                     , .names = '{.fn}.{.col}')) %>%
    ungroup() %>%
    rename(geoid = i)

  if(save.csv) {

    save.dir <- '/scratch/gpfs/km31/adjacencies+proximities/spatial-composites/composite-counts/'

    if(!exists(save.dir))
      dir.create(save.dir)

    rids <- geox::get.region.identifiers(cz = cz, cbsa = cbsa)
    rids <- paste0(rids, collapse = '-')

    save.pth <- paste0(save.dir, rids, '.csv')

    write.csv(n.dsts
              ,save.pth)
  }

  return(n.dsts)
}

#trimmed.resl <- resl %>% filter(var == 'perc_bl')
#calc.proxim.counts(cz = '19700')
calc.proxim.counts(cz = '00100')
#rm(js)


# slurm apply -------------------------------------------------------------

czs <- geox::rx %>% pull(cz) %>% unique()
# missing_czs <- c('19400', '17900', '17700', '16500', '18100', '17800')

resl$var %>% unique()

trimmed.resl <- resl %>%
  filter(var %in% c('perc_bl', 'perc_wh'
                    ,'perc_hsp', 'perc_asian'
                    ,'perc_below.pov', 'perc_above.pov'))
spws

cz.pms <- tibble(
  cz =  #missing_czs
     czs
  ,cbsa = NULL
  ,save.csv = T
  )

library(rslurm)
cz.dstcomposites.counts.dellajob <-
  slurm_apply(f =
                calc.proxim.counts,
              params = cz.pms,
              jobname = 'czcalc.proxim.counts',
              nodes = 6,
              cpus_per_node = 1,
              slurm_options = list(time = '15:00:00',
                                   'mem-per-cpu' = '30G',
                                   'mail-type' = list('begin', 'end', 'fail'),
                                   'mail-user' = 'km31@princeton.edu'),
              global_objects = c( 'trimmed.resl'
                                 ,'spws')
  )

# cbsas
cbsas <- geox::rx %>% pull(cbsa) %>% unique()

trimmed.resl
spws

cbsa.pms <- tibble(
  cz = NULL
  ,cbsa =
    cbsas
  ,save.csv = T)

library(rslurm)
cbsa.dstcomposites.counts.dellajob <-
  slurm_apply(f =
                calc.proxim.counts,
              params = cbsa.pms,
              jobname = 'cbsacalc.proxim.counts',
              nodes = 1, #22
              cpus_per_node = 1,
              slurm_options = list(time = '15:00:00',
                                   'mem-per-cpu' = '30G',
                                   'mail-type' = list('begin', 'end', 'fail'),
                                   'mail-user' = 'km31@princeton.edu'),
              global_objects = c( 'trimmed.resl'
                                  ,'spws')
  )

# check genr --------------------------------------------------------------
sdir <-
  '/scratch/gpfs/km31/adjacencies+proximities/spatial-composites/composite-counts/'

fns <- sdir %>% list.files()
fns[grepl('cz', fns)] %>% length()
fns[grepl('cbsa', fns)] %>% length()
# cool!


# do adjacencies ----------------------------------------------------------

# this is the same but all the weights are equal



devtools::load_all()

#' calc.adj.counts
#'
#'
calc.adj.counts <- function(cz = NULL, cbsa = NULL
                               #spatial.weights = spws,
                               #long.res.chars = trimmed.resl
                               ,save.csv = F) {

  require(tidyverse)

  if(! (exists('trimmed.resl') & exists('spws')) )
    warning('Wrapper expects trimmed.resl and spws in global env')

  # browser()

  # trim to region
  lw <- spws %>%
    geox::geosubset(cz = cz, cbsa = cbsa)


  # unnest. I'm gonna do just rook adj.
  js <- lw %>%
    unnest(rook.adj.nbs) %>%
    select(i = geoid, j = rook.adj.nbs)

  # process: get i-adjacent-weighted pop/hh of neighboring tracts. Then
  # multiply dist/pop weights by value in those tracts, and sum over i's.

  # merge to get % value and total pop/hh at each j
  js <- js %>%
    left_join(trimmed.resl
              , by = c('j' = 'geoid'))

  # remove where i == j, and trim to js in region
  js <- js %>%
    select(i, j,
           var, value, weight
    ) %>%
    filter(i != j ) %>%
    geox::geosubset(subset.cols = 'j',
                    cz = cz, cbsa = cbsa)

  # row standardised adjacency weights
  js <- js %>%
    group_by(i, var) %>%
    mutate( adj.weight =
              1 / n()
    )

  # "weight" colm is either population or hhs, depending on the variable
  js <- js %>%
    mutate( n.adj =
              adj.weight * weight * value
            )

  n.adjs <- js %>%
    group_by(i, var) %>%
    summarise( n.adj = sum(n.adj, na.rm = T) ) %>%
    ungroup() %>%
    rename(geoid = i)

  if(save.csv) {

    save.dir <- '/scratch/gpfs/km31/adjacencies+proximities/spatial-composites/adj-composite-counts/'

    if(!exists(save.dir))
      dir.create(save.dir)

    rids <- geox::get.region.identifiers(cz = cz, cbsa = cbsa)
    rids <- paste0(rids, collapse = '-')

    save.pth <- paste0(save.dir, rids, '.csv')

    write.csv(n.adjs
              ,save.pth)
  }

  return(n.adjs)
}

#rm(js)

calc.adj.counts(cz = '00402')
calc.adj.counts(cz = '00402', save.csv = T)


# slurm adj apply ---------------------------------------------------------

czs <- geox::rx %>% pull(cz) %>% unique()

resl$var %>% unique()

trimmed.resl <- resl %>%
  filter(var %in% c('perc_bl', 'perc_wh'
                    ,'perc_hsp', 'perc_asian'
                    ,'perc_below.pov', 'perc_above.pov'))
spws

cz.pms <- tibble(
  cz = czs
  ,cbsa = NULL
  ,save.csv = T
)

library(rslurm)
cz.adjcomposites.counts.dellajob <-
  slurm_apply(f =
                calc.adj.counts,
              params = cz.pms,
              jobname = 'czcalc.adj.counts',
              nodes = 22,
              cpus_per_node = 1,
              slurm_options = list(time = '15:00:00',
                                   'mem-per-cpu' = '20G',
                                   'mail-type' = list('begin', 'end', 'fail'),
                                   'mail-user' = 'km31@princeton.edu'),
              global_objects = c( 'trimmed.resl'
                                  ,'spws')
  )

# cbsas
cbsas <- geox::rx %>% pull(cbsa) %>% unique()

trimmed.resl
spws

cbsa.pms <- tibble(
  cz = NULL
  ,cbsa = cbsas
  ,save.csv = T)

library(rslurm)
cbsa.adjcomposites.counts.dellajob <-
  slurm_apply(f =
                calc.adj.counts,
              params = cbsa.pms,
              jobname = 'cbsacalc.adj.counts',
              nodes = 22,
              cpus_per_node = 1,
              slurm_options = list(time = '15:00:00',
                                   'mem-per-cpu' = '10G',
                                   'mail-type' = list('begin', 'end', 'fail'),
                                   'mail-user' = 'km31@princeton.edu'),
              global_objects = c( 'trimmed.resl'
                                  ,'spws')
  )

# check genr --------------------------------------------------------------
sdir <-
  '/scratch/gpfs/km31/adjacencies+proximities/spatial-composites/adj-composite-counts/'
fns <- sdir %>% list.files()
fns[grepl('cz', fns)] %>% length()
fns[grepl('cbsa', fns)] %>% length()
# cool!


