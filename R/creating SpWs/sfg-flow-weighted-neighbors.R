# ws ---------------------------------------------------------------------------
rm(list = ls())
require(tidyverse)
require(sf)

# option setting
sf_use_s2(T)
options(tigris_use_cache = TRUE)


wdir <-
  '~/R/all sharkey geoseg work/divflow/R/creating SpWs/spWs/'


# load a sample sfg cz ---------------------------------------------------------
Sys.getenv('drop_dir') %>%
  paste0()
  list.files(pattern = 'sfg-processed/annual'
                                      ,recursive = T)
sfg.dir <-
  paste0(Sys.getenv('drop_dir'),
         'sfg-processed/orig_dest_annual/'
         )

sample.cz <- '19700'
sfg <- sfg.seg::read.sfg.CZs(czs2load = sample.cz
                             ,sfg.dir)


# descriptives -----------------------------------------------------------------

# i did live in CBGs
sfg$origin[1:1000] %>% nchar() %>% unique()

deciles <- seq(0,1,.1)

sfg$n %>% quantile(deciles)
# 70% are fewer than 1 avg trips

sfg

# what do trips outside of CZ look like?
coids <- geox::x2cos(cz = sample.cz)
outside.of.cz <- sfg %>%
  filter(!substr(dest, 1,5) %in%
           coids)
outside.of.cz$n %>% quantile(deciles)
where.ppl.going <- outside.of.cz %>%
  group_by(dest) %>%
  summarise(visits.to = sum(n)) %>%
  ungroup()
where.ppl.going %>% arrange(desc(visits.to))
# They're going to disney world, florida, places just barely outside CZ boundaries.

# I think a good approach is subsetting to tracts within 20 miles away...

# get tracts 20 mi away --------------------------------------------------------

# I am still generating code for BGs so will use tracts for now
prx.dir <- paste0(Sys.getenv('drop_dir'),
                    'tract adjacencies+proximities/'
                  )
prx <- prx.dir %>%
  list.files(pattern = 'adjacencies.rds$'
             ,full.names = T) %>%
  read_rds()


# join sfg to proximities ------------------------------------------------------

sfg <- sfg %>% sfg.seg::cbg.flows2tracts()

sfgw <-
  sfg %>%
  left_join(prx,
            by = c('origin'='geoid'))

# for each tract i, get flow weights based on tracts within distance sphere dst.max
i <- sfgw[1,]

j.ids <- i$below.cutoff[[1]]
js <- sfgw %>%
  filter(origin %in% i &
           dest %in% j.ids) %>%
  filter(origin != dest)

flwws <- js %>%
  mutate(flow.weight =
           n / sum(n)) %>%
  select(dest, flow.weight)

# add geoids as names
flwws$flow.weight %>%
  as.vector() %>%
  setNames(flwws$dest)

# that's % of trips from i to j -- i.e., visited weights for i
# to get visitor weights for i:
i
j.ids
js <- sfgw %>%
  filter(dest %in% i &
           origin %in% j.ids) %>%
  filter(origin != dest)

flwws <- js %>%
  mutate(flow.weight =
           n / sum(n)) %>%
  select(dest, flow.weight)
flwws



#' get.flow.weights.within.distance
#'
#' Get "flow weights," i.e., proportion of visitor/visited from other tracts within
#' distance threshold. Relies on precalculated distance matrix. Rounds to 4 digits
#' for sparseness.
#'
#' @inheritParams get.dist.weighted.composite
#' @param flow.type Whether to get weights based on all visitors to i from other
#'   proximate areas (default); or based on all proximiate areas visited by people in
#'   area i.
#' @param drop.loops Whether to drop loops (where origin==destination)
#' @param prox.col name for list-column in `spws` with all areas within distance
#'   cutoff from each other
#' @param flow.counts data.frame with origin/destination/estimated visits (n).
#' @param weight.floor flow weight floor. Minimum percent of incoming/visited trips
#'   between tracts to be included. 0.1% by default.
#'
#' @export get.flow.weights.within.distance
get.flow.weights.within.distance <- function(i
                                             ,flow.type = c('visitors', 'visiting')
                                             ,drop.loops = T
                                             ,prox.col =  'below.cutoff'
                                             ,flow.counts = sfg
                                             ,spatial.weights = spws
                                             ,weight.floor = .001) {

  flow.type <- flow.type[[1]]

  # get geoids for js within distance threshold
  j.ids <- spatial.weights  %>% filter(geoid %in% i) %>% pull(prox.col) %>% unlist()
  # drop loops if appropriate
  if(drop.loops)
    flow.counts <- flow.counts %>% filter(origin != dest)

  # Get either all trips from Js to I (visitors) or from I to all Js (visiting)
  if(flow.type == 'visitors') {

    js <- flow.counts %>%
      filter(dest %in% i &
               origin %in% j.ids)
  } else if(flow.type == 'visiting') {

    js <- flow.counts %>%
      filter(origin %in% i &
               dest %in% j.ids)
  }


  # get flow weights
  flwws <- js %>%
    mutate(flow.weight =
             n / sum(n)) %>%
    select(dest, flow.weight) %>%
    filter(flow.weight >=
             weight.floor)

  # return as named vector
  flwws$flow.weight %>%
    as.vector() %>%
    round(digits = 4) %>%
    setNames(flwws$dest)
}

sfg

get.flow.weights.within.distance(i = '34001000100'
                                 ,'visitors'
                                 ,flow.counts = sfg
                                 ,spatial.weights = prx
                                 ,weight.floor = 0)

tmp <- prx %>%
  filter(geoid %in% sfg$origin)
test.geoids <- tmp$below.cutoff %>% unlist() %>% unique()

test <- prx %>%
  filter(geoid %in% test.geoids) %>%
  mutate(inc.flow.weights =
           map(geoid
               ,~get.flow.weights.within.distance(i = .x
                                                  ,'visitors'
                                                  ,flow.counts = sfg
                                                  ,spatial.weights = prx
                                                  ,weight.floor = 0.001
                                                  )))
test <- test %>%
  mutate(visited.flow.weights =
           map(geoid
               ,~get.flow.weights.within.distance(i = .x
                                                  ,'visiting'
                                                  ,flow.counts = sfg
                                                  ,spatial.weights = prx
                                                  ,weight.floor = 0.001
               )))

test.geoids %>%
  tibble(state = substr(.,1,2)
         ,co = substr(.,3,5)
           )
# so the approach can be:
# -read in one CZ worth of sfg data
# -get all geoids that are distance eligible (w/in 20 mi)
# -use sqldf::read.csv.sql to read data eligible cbgs outside the CZ
# -map.get.flow.weights.within.distnace
sfg

sfg
test
