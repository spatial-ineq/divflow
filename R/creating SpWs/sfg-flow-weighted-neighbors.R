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

prx.dir <- paste0(Sys.getenv('drop_dir'),
                    'tract adjacencies+proximities/'
                  )
prx <- prx.dir %>%
  list.files(pattern = 'adjacencies.rds$'
             ,full.names = T) %>%
  read_rds()

sfg <-
  sfg %>%
  left_
prx
0
