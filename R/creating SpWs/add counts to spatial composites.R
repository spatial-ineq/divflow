# ---------------------------------------------------------------------------
rm(list = ls())
require(tidyverse)
require(sf)

# option setting
sf_use_s2(T)
options(tigris_use_cache = TRUE)

# dropbox dir or della base dir
ddir <- Sys.getenv('drop_dir')

# notes ------------------------------------------------------------------------

#' the prepped dataset "full-spatial-composites-by-rt" (generated through della
#' scripts) has everything in percents. I need counts and regional totals for flows
#' and res, so have to calcute from tract-level residence- & flow- & spatial-
#' weights.


# load residence values --------------------------------------------------------



# load existing composite ------------------------------------------------------

# get flow composites
lag.dir <- paste0(ddir,
                  'adjacencies+proximities/spatial-composites/' )
lag.dir %>% list.files()

splags <- lag.dir %>%
  list.files(pattern = 'full-spatial-composites-by-rt.csv'
             ,full.names = T) %>%
  vroom::vroom()

splags
splags$geoid %>% head() %>% nchar()

# trim to scratch cz
"splags <- splags %>%
  filter(rt == 'cz' &
           rid %in% smplcz)

splags <- splags %>%
  filter(var == 'perc_bl')
"

# avail vars
splags %>% count(var) %>% pull(var)
# drop the fam vars / any others ?


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

# also get all residential values, long by variable
resl <- flrdir %>%
  list.files(pattern = 'res-chars-long'
             ,full.names = T) %>%
  vroom::vroom()

resl

# gen res/flow counts by var -------------------------------------------------------------

# these are easy, i just multiply %flow composite by total flows or total pop/hh
ccs <- splags %>%
  left_join(tots)

ccs <- ccs %>%
  mutate(
    n.value = value * weight
    ,n.inc.flw =
      inc.flww.composite * total.incoming
    ,n.vis.flw =
      vis.flww.composite * total.visits.from
  )

# a check
ccs %>%
  geox::geosubset(cz = '24701') %>%
  filter(var == 'perc_bl') %>%
  select(1,var,value, matches('^n|flw'))


# prx counts -------------------------------------------------------------------

# these are more of a pain, because the weights should be: row-standardized
# population/hh * decayed distance or adjacency. And i don't have those to just merge
# in.


# prx process scratch ----------------------------------------------------------

# process is:

# get spatial weights
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

# for scratch, trim to sample
smpl <- spws %>% geox::geosubset(cz = '24701')

smpl[1,]$dists
splags

smpl <- smpl %>%
  unnest(dists) %>%
  mutate( negexp.decay =
           neg.exp( as.numeric(dist) )
         ,bisq.8km.decay =
           bisq.dist2weights( as.numeric(dist), cutoff = 8 )
         )


# get distance-from-i-weighted avg pop/hh of neighboring tracts. Then multiply
# dist/pop weights by value in those tracts, and sum over i's.


# merge to get % value and total pop/hh at each j
js <- smpl %>%
  select(i = geoid, j, matches('dist|decay')) %>%
  left_join(resl
            , by = c('j' = 'geoid'))

smpl %>%
  filter(i != j ) %>%
  geox::geosubset(subset.cols = 'j', cz = '24701') %>%
  mutate( across(matches('decay$')
                 , ~ .x / sum(.x)
  ))

js

# remove where i == j
js <- js %>%
  filter(i != j ) %>%
  geox::geosubset(subset.cols = 'j', cz = '24701')

# row standardise decays
js <- js %>%
  group_by(i, var) %>%
  mutate( across(matches('decay$')
                 , ~ .x / sum(.x)
  ))

# check
js %>% group_by(i, var) %>% summarise(across(matches('decay$')
                                             , sum ))

n.dstsst <-
  js %>%
  mutate( across(matches('decay$')
                 , ~ .x * weight * value
  )) %>%
  group_by(i, var) %>%
  summarise(across(matches('decay$')
                   , list(n = sum)
                   , .names = '{.fn}.{.col}'))


# laugh test:
n.dstsst %>%
  filter(substr(var,1,4) == 'perc')

ccs %>%
  filter(geoid == '17133600101') %>%
  filter(substr(var,1,4) == 'perc') %>%
  select(1:4,var, n.value)



# functionalize ----------------------------------------------------------------

devtools::load_all()

#' calc.proxim.counts
#'
#'
calc.proxim.counts <- function(cz = NULL, cbsa = NULL,
                               spatial.weights = spws,
                               long.res.chars = resl) {

  lw <- spatial.weights %>%
    geox::geosubset(cz = cz, cbsa = cbsa)

 #  browser()

  lw <- lw %>%
    unnest(dists) %>%
    mutate( negexp.decay =
              neg.exp( as.numeric(dist) )
            ,bisq.8km.decay =
              bisq.dist2weights( as.numeric(dist), cutoff = 8 )
    )

  # get distance-from-i-weighted avg pop/hh of neighboring tracts. Then multiply
  # dist/pop weights by value in those tracts, and sum over i's.

  # merge to get % value and total pop/hh at each j
  js <- lw %>%
    left_join(long.res.chars
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

  n.dsts <-
    js %>%
    mutate( across(matches('decay$')
                   , ~ .x * weight * value
    )) %>%
    group_by(i, var) %>%
    summarise(across(matches('decay$')
                     , list(n = sum)
                     , .names = '{.fn}.{.col}')) %>%
    ungroup() %>%
    rename(geoid = i)

  return(n.dsts)
}

#rm(js)

tmp <- calc.proxim.counts(cz= '24701'
                          ,long.res.chars =
                            filter(resl
                                   ,var %in% c('perc_bl', 'perc_wh'
                                               ,'perc_hsp', 'perc_asian'))
                          )
tmp
n.dstsst %>% filter(var %in% c('perc_bl', 'perc_wh'
                   ,'perc_hsp', 'perc_asian'))

ccs %>% filter(geoid == '17133600101') %>% filter(var %in% c('perc_bl', 'perc_wh'
                                                         ,'perc_hsp', 'perc_asian'))


# map over regions czs ---------------------------------------------------------------

# merge into ccs.



calc.proxim.counts

# a more straightforward way: --------------------------------------------------


n.dsts
n.dsts2

across(matches('decay$')
                , list(decayed.popw =
                         stats::weighted.mean(weight, ))
         )



# scratch ----------------------------------------------------------------------

# merge to get % value and total pop/hh at each j
js <- smpl %>%
  select(i = geoid, j, matches('decay')) %>%
  left_join(
    select(resl
           , geoid, var, value)
    , by = c('j' = 'geoid')) %>%
  left_join(distinct(select(resl
                            , geoid, var, pop, hh, weight))
            , by = c('i' = 'geoid', 'var'))

# form and row-standardize wegiths
js <- js %>%
  mutate( negexp.w =
            negexp.decay * weight
          ,bisq.8km.w =
            bisq.8km.decay * weight)

js <- js %>%
  group_by(i, var) %>%
  mutate( snegexp.w =
            negexp.w / sum(negexp.w)
          ,sbisq.8km.w =
            bisq.8km.w / sum(bisq.8km.w)
  )

# checks
js %>% group_by(i, var) %>% summarise(across(matches('.w$')
                                             , sum ))
js %>% filter(var == 'perc_bl')
splags
# multiply standardized flow weights for each j by variable % at each j and the
# total incoming to or visited from at i
js %>%
  group_by(i, var) %>%
  mutate(bisq.8km.n =
           sbisq.8km.w * )

