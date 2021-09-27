# ws ---------------------------------------------------------------------------
rm(list = ls())
require(tidyverse)
require(sf)

# option setting
sf_use_s2(T)
options(tigris_use_cache = TRUE)


wdir <-
  '~/R/all sharkey geoseg work/divflow/R/creating SpWs/spWs/'

# start with all census tracts ---------------------------------------------------------------

statefps <- geox::rx$statefp %>% unique()

cts <- map_dfr( statefps,
                ~tigris::tracts(.
                                ,year = 2019
                                ))
cts <- cts %>%
  select(1:4, matches("^A")) %>%
  rename_with(tolower)

# get adjacent neigbhors -------------------------------------------------------


# process illustration / check over test area ---------------------------------
# define test set

cts
tmp <- cts %>% filter(statefp == '01' &
                        countyfp == '055')
"
# use st_relate to get 'queen' contiguity
nbs <- sf::st_relate(tmp
                     ,pattern = 'F***T****' )

nbs[[2]]
# index using neighbor list and get geoids
tmp$nbs <-
  map(nbs,
      ~tmp[., ]$geoid)

tmp$n.nbs <- map_dbl(tmp$nbs, length)
"

# tmp['n.nbs'] %>% mapview::mapview(zcol = 'n.nbs')

# gen over all tracts -----------------------------------------------------------

# get all adjacency
queen.nbs <- sf::st_relate(cts
                           ,pattern = 'F***T****' )

cts$queen.adj.nbs <-
  map(queen.nbs,
      ~cts[., ]$geoid)

# save
tibble(cts) %>%
  select(geoid, queen.adj.nbs) %>%
  saveRDS(file = paste0(wdir,
                        'queen-adj-nbs.rds'))

# rook
rook.nbs <- sf::st_relate(cts
                           ,pattern = 'F***1****' )

cts$rook.adj.nbs <-
  map(rook.nbs,
      ~cts[., ]$geoid)

# save
tibble(cts) %>%
  select(geoid, rook.adj.nbs) %>%
  saveRDS(file = paste0(wdir,
                        'rook-adj-nbs.rds'))


# get distance-weighted neighbors ----------------------------------------------

library(units)
# define cut-off distance
dst.ceiling <- units::set_units(20, 'miles')

# process illustration / check over test area ---------------------------------

# get centroids
ctrs <- tmp %>% st_centroid()
dst.nbs <-
  st_is_within_distance(ctrs
                        , dist =
                          dst.ceiling
                        )


# index using neighbor list and get geoids
tmp$nbs <-
  map(dst.nbs,
      ~tmp[., ]$geoid)

tmp$n.nbs <- map_dbl(tmp$nbs, length)
# units(cts$geometry)

# tmp[,c('geoid', 'n.nbs')] %>% plot()
# for each ct, get distance from all other cts that are within threshold

# function moved to `spW-fncs`
devtools::load_all()
dsts.below.threshold(tmp, i = 1, j.colm = 'nbs')


tmp2 <- tmp %>%
  mutate(dists =
           map(1:nrow(.)
               ,~dsts.below.threshold(tmp, ., 'nbs')))

tmp2[1,]$dists

# do for all cts ---------------------------------------------------------------

# define cut-off distance
dst.ceiling <- units::set_units(20, 'miles')

# get centroids
ctrs <- cts %>% st_centroid()
dst.nbs <-
  st_is_within_distance(ctrs
                        , dist =
                          dst.ceiling
  )


# index using neighbor list and get geoids
cts$below.cutoff <-
  map(dst.nbs,
      ~cts[., ]$geoid)


# save
tibble(cts) %>%
  select(geoid, below.cutoff) %>%
  saveRDS(file = paste0(wdir,
                        'nbs-within-20mi.rds'))

tmp$n.below.cutoff <- map_dbl(cts$below.cutoff, length)



# for each ct, get distance from all other cts that are within threshold

# add below.cutoff neighbors to centroids sf
ctrs$below.cutff <- cts$below.cutoff

cts <- ctrs %>%
  mutate(dists =
           map(1:nrow(ctrs)
               ,~dsts.below.threshold(ctrs, .,
                                      'below.cutff')))

# save
tibble(cts) %>%
  select(geoid, dists) %>%
  saveRDS(file = paste0(wdir,
                        'dists-within-20mi.rds'))



# load aall and turn into one data.frame ---------------------------------------
wdir <-
  '~/R/all sharkey geoseg work/divflow/R/creating SpWs/spWs/'
cgts <- list.files(wdir
                   ,pattern='^dists|^nbs|^queen|^rook')
cgts <-
  cgts %>% map( ~readRDS(paste0(wdir, .)))

cgts


cgts <- purrr::reduce(cgts, left_join)


# check ------------------------------------------------------------------------

cgts$queen.adj.nbs %>% head() %>% map_dbl(length)


cgts
cgts %>%
  mutate(across(matches('adj.nbs')
                ,list(len =
                          ~map_dbl( ., length))
                )) %>%
  #filter(queen.adj.nbs_len != rook.adj.nbs_len)
  filter(queen.adj.nbs_len < rook.adj.nbs_len)

# save
saveRDS(cgts,
        paste0(wdir,
               'tract-adjacencies.rds'))


# some final checks ------------------------------------------------------------

chi <-
  cgts %>%
  geox::geosubset(cz_name = 'Chicago')

chi %>%
  mutate(n.rook.nbs = map_dbl(rook.adj.nbs, length)) %>%
  geox::attach.geos() %>%
  select(1, n.rook.nbs) %>%
  mapview::mapview(zcol = 'n.rook.nbs')

chi <- chi %>%  geox::attach.geos()

near.i <- chi[1,]$dists[[1]]
near.i <- near.i[names(near.i) %in% chi$geoid]

length(near.i)

chi[chi$geoid %in% names(near.i), ] %>%
  ggplot() +
  geom_sf(aes(fill =
                as.numeric(near.i))
          ) +
  geom_sf(data = chi[1, ]
          ,fill='#880033')
