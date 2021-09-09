# ws ---------------------------------------------------------------------------
rm(list = ls())
require(tidyverse)
require(sf)
devtools::load_all()

# option setting
sf_use_s2(T)
options(tigris_use_cache = TRUE)
Sys.setenv("VROOM_SHOW_PROGRESS"="false")

wdir <-
  '~/R/all sharkey geoseg work/divflow/R/creating SpWs/spWs/'


# load a sample sfg cz ---------------------------------------------------------
Sys.getenv('drop_dir') %>%

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
                    'adjacencies+proximities/'
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
prx
within.range <- prx %>% filter(geoid %in% sfg$origin) %>% pull(below.cutoff) %>%
  unlist() %>% unique()
cos <- within.range %>% substr(1,5) %>% unique()
czs.in.range <- geox::rx %>%
  filter(countyfp %in% cos) %>%
  pull(cz) %>% unique()

czs.in.range <- czs.in.range %>%
  sfg.seg::read.sfg.CZs(sfg.dir = sfg.dir, year = '2019')

czs.in.range
sfg



query <- paste0('SELECT * from file where SUBSTR(origin_census_block_group,1,5) in (',
                cos[1],
                ')')

# map through czs in range, taking
dover.dir <- paste0(sfg.dir,
                    '19901/2019')
sqldf::read.csv.sql(list.files(dover.dir,full.names = T)
                    ,
                    sql =query)
dover.sfg <- vroom::vroom(list.files(dover.dir,full.names = T))
dover.sfg
sfg.seg::read.sfg.CZs
sfg
sfg
test




# hashing out della fcn for flow-weights ---------------------------------------




# czs <- geox::build.CZs()
# czs
devtools::load_all()
rm(list = ls())
sfg.dir <-
  paste0(Sys.getenv('drop_dir'),
         'sfg-processed/orig_dest_annual/'
  )

prx.dir <- paste0(Sys.getenv('drop_dir'),
                  'adjacencies+proximities/'
                  )

spws <- prx.dir %>%
  list.files(pattern = 'adjacencies.rds$'
             ,full.names = T) %>%
  read_rds()

morristown.flwws <- Della.wrapper_flow.weights(cz = '00200'
                           ,sfg.dir = sfg.dir
                           )
morristown.flwws[1,]$inc.flow.weights

