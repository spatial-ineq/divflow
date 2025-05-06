# ---------------------------------------------------------------------------
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


# load residence values --------------------------------------------------------



# load existing composite ------------------------------------------------------

# get flow composites
lag.dir <- paste0(ddir,
                  'adjacencies+proximities/spatial-composites/' )
lag.dir %>% list.files()

splags <- lag.dir %>%
  list.files(pattern = 'full-spatial-composites-by-cbsa-pctiles-new.csv'
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
splags %>% filter(var == 'perc_bl') %>% filter(is.na(inc.flww.composite))

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
  list.files(pattern = 'res-chars-long-pctiles-new'
             ,full.names = T) %>%
  vroom::vroom()

resl

# gen res/flow counts by var -------------------------------------------------------------

# these are easy, i just multiply %flow composite by total flows or total pop/hh

# first drop vars you didn't get other counts for
selags <- splags
# %>% filter(var %in%  c('perc_bl', 'perc_wh'
#                               ,'perc_hsp', 'perc_asian'
#                               ,'perc_below.pov', 'perc_above.pov'))

# nas pre-join
selags %>% map_dbl( ~sum(is.na(.x)) )

ccs <- selags %>%
  mutate(rid = geox::fix.geoid(rid, width = 5)) %>%
  left_join(tots)

ccs <- ccs %>%
  mutate(
    n.value = value * weight
    ,n.inc.flw =
      inc.flww.composite * total.incoming
    ,n.vis.flw =
      vis.flww.composite * total.visits.from
  )


# checks ------------------------------------------------------------------
ccs %>%
  geox::geosubset(cz = '24701') %>%
  filter(var == 'perc_bl') %>%
  select(1,var,value, matches('^n|flw'))

# missing values?
ccs %>%  map_dbl( ~sum(is.na(.x)) )
# missing regions (issue fixed):
missing <- ccs %>% anti_join(tots) %>% select(rt, rid) %>% distinct()
missing %>% geox::add.rns()
tots %>% filter(rid %in% missing$rid & rt == 'cbsa')


# trims -------------------------------------------------------------------

# get ~just~ counts from ccs.
ccs <- ccs %>%
  select(geoid,rt,rid,pop,hh,var,value,weight,matches('^n\\.'), everything())

# trim to demographic and pov vars
ccs <- ccs %>%
  filter(var %in% c('perc_bl', 'perc_wh'
                           ,'perc_hsp', 'perc_asian'
                           ,'perc_below.pov', 'perc_above.pov'))

# prx counts -------------------------------------------------------------------

# generated in della
prx.count.dir <- '/scratch/gpfs/km31/adjacencies+proximities/spatial-composites/composite-counts/'

fns <- prx.count.dir %>% list.files(full.names = T)
fns[1]  %>% vroom::vroom()

read.add.region.info <- function(.rt, .rid, .dir) {

  fn <- .dir %>%
    list.files(pattern = glue::glue('{.rt}-{.rid}')
               ,full.names = T)
  if(length(fn) == 0)
    return(glue::glue('missing at {.rt}-{.rid}'))

  x <- fn %>%
    vroom::vroom() %>%
    select(-any_of('...1')) %>%
    mutate(rt = .rt, rid = .rid
           ,.after = geoid)

  return(x)
}

# individual checks
"
read.add.region.info('cbsa', '10100'
                     ,prx.count.dir)

read.add.region.info('cz', '19700'
                     ,prx.count.dir) %>%
  map_dbl(~sum(is.na(.)) / length(.) )


read.add.region.info('cz', '19400'
                     ,prx.count.dir) %>%
  map_dbl(~sum(is.na(.)) / length(.) )
"

# read all
allrs <- geox::rx %>%
  select(cz,cbsa) %>%
  pivot_longer(c(cz, cbsa)
               ,names_to = 'rt'
               ,values_to = 'rid') %>%
  distinct()

tmp <- allrs # %>% tail(40)

prx.ccs <- map2(tmp$rt, tmp$rid,
                    ~read.add.region.info(.x, .y
                                          ,prx.count.dir)
                    )

# check where they read in wrong?
index <- prx.ccs %>% map_dbl( ncol)
prx.ccs[which(index != 6)]

# fix column types and rbind
prx.ccs <- prx.ccs %>%
  map( ~mutate(.x,
              across(matches('^n\\.')
                     ,as.numeric )
              ,across(c(geoid, rt, rid, var)
                      ,as.character )
              )) %>%
  do.call('rbind', .)

prx.ccs

# adj counts --------------------------------------------------------------

adj.count.dir <- '/scratch/gpfs/km31/adjacencies+proximities/spatial-composites/adj-composite-counts/'

adj.ccs <- map2(allrs$rt, allrs$rid,
                    ~read.add.region.info(.x, .y
                                          ,adj.count.dir)
                    )
# check where they read in wrong?
index <- adj.ccs %>% map_dbl( ncol)
adj.ccs[which(index != 5)]

# fix column types and rbind
adj.ccs <- adj.ccs %>%
  map( ~mutate(.x,
               across(matches('^n\\.')
                      ,as.numeric )
               ,across(c(geoid, rt, rid, var)
                       ,as.character )
  )) %>%
  do.call('rbind', .)

adj.ccs

# combine -----------------------------------------------------------------

ccs
prx.ccs
adj.ccs

spcounts <- ccs # purrr::reduce(list(ccs,
              #      prx.ccs,
              #      adj.ccs)
              # ,full_join)

spcounts
spcounts %>% map_dbl( ~sum(is.na(.x)) )

list.dirs('R/creating SpWs/')
# save backup
"
write.csv(spcounts,
          file = 'R/creating SpWs/.intermediate-saves/spcounts-prechecks.csv'
          ,row.names = F)
"


# check -------------------------------------------------------------------

spcounts
spcounts %>% map_dbl(~sum(is.na(.x)) / length(.x))

# remove NA cbsas
spcounts %>%
  filter(is.na(rid)) %>%
  count(rt)

spcounts <- spcounts %>%
  filter(!is.na(rid))

spcounts %>%  map_dbl(~sum(is.na(.x)) / length(.x)  )
spcounts %>% filter_at(vars(everything()), any_vars(is.na(.))) %>% count(geoid)

tmpc <- spcounts %>%
  filter(rid == '00100' &
           rt == 'cz' &
           var %in% 'perc_bl')

tmpp <- splags %>%
  filter(rid == '00100' &
           rt == 'cz' &
           var %in% 'perc_bl')


spcounts %>%
  filter(rid == '19700' &
           rt == 'cz' &
           var %in% 'perc_bl') %>%
  rename(perc = value) %>%
  pivot_longer(matches('^n\\.')
               ,names_prefix = 'n.'
               ,names_to = 'spatial.lag.type'
               ,values_to = 'spatial.lag.value') %>%
  geox::attach.geos() %>%
  st_sf() %>%
  ggplot() +
  geom_sf(aes(fill = spatial.lag.value)
          , color = NA) +
  scale_fill_viridis_c() +
  facet_wrap(vars(spatial.lag.type)) +
  theme_void() +
  theme(legend.position = 'bottom')


# write final copy --------------------------------------------------------

save.path <- '/scratch/gpfs/km31/adjacencies+proximities/spatial-composites/full-spatial-composite-counts-by-cbsa-pctiles-new.csv'
spcounts %>%
  write.csv(
     file = save.path
    ,row.names = F
  )
