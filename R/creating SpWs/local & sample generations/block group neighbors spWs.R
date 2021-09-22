# ws ---------------------------------------------------------------------------
rm(list = ls())
require(tidyverse)
require(sf)

# option setting
sf_use_s2(T)
options(tigris_use_cache = TRUE)

wdir <-
  #'~/R/all sharkey geoseg work/divflow/R/creating SpWs/spWs/'
  '~/all/divflow/R/creating SpWs/block-group-spWs/'
# start with all block groups ---------------------------------------------------------------

statefps <- geox::rx$statefp %>% unique()

bgs <- map_dfr( statefps,
                ~tigris::block_groups(.
                                      ,year = 2019
                ))
bgs <- bgs %>%
  select(1:4, GEOID, matches("^A")) %>%
  rename_with(tolower)

# get adjacent neigbhors -------------------------------------------------------

# process illustration / check over test area ---------------------------------
# define test set
"
bgs
tmp <- bgs %>% filter(statefp == '01' &
                        countyfp == '055')

# use st_relate to get 'queen' contiguity
nbs <- sf::st_relate(tmp
                     ,pattern = 'F***T****' )

nbs[[2]]
# index using neighbor list and get geoids
tmp$nbs <-
  map(nbs,
      ~tmp[., ]$geoid)
tmp[2,]
tmp$n.nbs <- map_dbl(tmp$nbs, length)
"
# tmp['n.nbs'] %>% mapview::mapview(zcol = 'n.nbs')

# gen over all tracts -----------------------------------------------------------

# get all adjacency
queen.nbs <- sf::st_relate(bgs
                           ,pattern = 'F***T****' )

bgs$queen.adj.nbs <-
  map(queen.nbs,
      ~bgs[., ]$geoid)

# save
"
tibble(bgs) %>%
  select(geoid, queen.adj.nbs) %>%
  saveRDS(file = paste0(wdir,
                        'queen-adj-nbs.rds'))
"
# rook
rook.nbs <- sf::st_relate(bgs
                           ,pattern = 'F***1****' )

bgs$rook.adj.nbs <-
  map(rook.nbs,
      ~bgs[., ]$geoid)

# save
"
tibble(bgs) %>%
  select(geoid, rook.adj.nbs) %>%
  saveRDS(file = paste0(wdir,
                        'rook-adj-nbs.rds'))
"

# get distance-weighted neighbors ----------------------------------------------

library(units)
# define cut-off distance
dst.ceiling <- units::set_units(25
                                , 'miles')

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
# units(bgs$geometry)

# tmp[,c('geoid', 'n.nbs')] %>% plot()
st_crs(ctrs)$wkt
# for each ct, get distance from all other bgs that are within threshold


devtools::load_all()

dsts.below.threshold(ctrs,
                     1)


tmp2 <- ctrs %>%
  mutate(dists =
           map(1:nrow(ctrs)
               ,~dsts.below.threshold(ctrs, ., 'nbs')))

tmp2[1,]$dists

# do for all bgs ---------------------------------------------------------------

# define cut-off distance
dst.ceiling <- units::set_units(25, 'miles')
devtools::load_all()

# get centroids
ctrs <- bgs %>% st_centroid()

dst.nbs <-
  st_is_within_distance(ctrs
                        , dist =
                          dst.ceiling
  )


# index using neighbor list and get geoids
bgs$below.cutoff <-
  map(dst.nbs,
      ~bgs[., ]$geoid)


# save
tibble(bgs) %>%
  select(geoid, below.cutoff) %>%
  saveRDS(file = paste0(wdir,
                        'nbs-within-25mi.rds'))

tmp$n.below.cutoff <- map_dbl(bgs$below.cutoff, length)



# for each ct, get distance from all other bgs that are within threshold

#' dsts.below.threshold
#'
#' Quick fcn with hardcoded items, including assumption of meters from crs length
#' unit
#'
#' @param sfx an sf object to selectively get pairwise distances between, probably
#'   containing point geometries.
#' @param i row # index of given tract/centroid in `sfx`
#' @param js list of geoids in same sf object `sfx` to get distance from
dsts.below.threshold <- function(sfx, i, j.colm = 'nbs') {

  #  browser()
  # get sf corresponding to geoids j
  js <- pull(sfx, j.colm)[[i]]
  geos <- sfx[sfx$geoid %in% js, ]$geometry

  # get pairwise distance between i and j
  as.vector(st_distance(sfx[i, ]
                        , geos)) %>%
    set_units('m') %>%
    set_units('km') %>%
    round(digits = 2) %>%
    setNames(js)

}

# add below.cutoff neighbors to centroids sf
ctrs$below.cutff <- bgs$below.cutoff

bgs <- ctrs %>%
  mutate(dists =
           map(1:nrow(ctrs)
               ,~dsts.below.threshold(ctrs, .,
                                      'below.cutff')))

# save
tibble(bgs) %>%
  select(geoid, dists) %>%
  saveRDS(file = paste0(wdir,
                        'dists-within-20mi.rds'))



# load aall and turn into one data.frame ---------------------------------------
#wdir <- '~/R/all sharkey geoseg work/divflow/R/creating SpWs/spWs/'

cgts <- list.files(wdir
                   ,pattern='^dists|^nbs|^queen|^rook')
cgts <-
  cgts %>% map( ~readRDS(paste0(wdir, .)))

cgts


cgts <- purrr::reduce(cgts, left_join)

# save
#saveRDS(cgts,
#        paste0(wdir,
#               'tract-adjacencies.rds'))


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
