# ws ---------------------------------------------------------------------------
rm(list = ls())
require(tidyverse)
require(sf)

# option setting
sf_use_s2(T)
options(tigris_use_cache = TRUE)


# load spatial weight matrices -------------------------------------------------

spw.dir <- paste0(Sys.getenv('drop_dir'), '/tract adjacencies+proximities/')
spws <- list.files(spw.dir, pattern='adjacencies.rds'
           ,full.names = T) %>%
  readRDS()

spws


# get demovars -----------------------------------------------------------------

# (compiled in data-raw/)
demo.pth <- paste0(Sys.getenv('drop_dir')
                  ,'seg-measures/by tract/broader ineq flows/'
                  ,'res-chars-long.csv')
resl <- vroom::vroom(demo.pth)

# get adjacency weighted avgs --------------------------------------------------

#' get.avg.adjacent
#'
#' @param i a given tract geoid
#' @param x full ct data to subset from
#' @param value.col column name in `x` to calculate weighted avg of adjacent tracts
#'   for.
#' @param weight.col column name in `x`, probably denoted population or households
#' @param nbs string for column name in `spws`
#' @param spws spatial weight matrices with a geoid column and other list columns for
#'   different spatial weights
#'
get.avg.adjacent <- function(i,  x
                             , value.col = 'value'
                             , weight.col = 'weight'
                             ,nb.col
                             ,spatial.weights = spws) {


  #browser()
  #x <- x %>% filter(name %in% !!value)
  j.ids <- unlist(spatial.weights[spatial.weights$geoid == i, nb.col])
  js <- x[x$geoid %in% j.ids, ]

  # spatial mean
  spu <- stats::weighted.mean(js[[value.col]]
                              ,js[[weight.col]]
                              ,na.rm = T)

  return(spu)
}


# split-map-bind
"
# queen-weighted
resl <- resl %>%
  split(.$name) %>%
  map_dfr( ~mutate(.,
                   q.adj.wval =
                     map_dbl(geoid,
                             function(i)
                               get.avg.adjacent(i, .,
                                                nb.col= 'queen.adj.nbs'
                               ))))
spws
# rook-weighted
resl <- resl %>%
  split(.$name) %>%
  map_dfr( ~mutate(.,
                   rook.adj.wval =
                     map_dbl(geoid,
                             function(i)
                               get.avg.adjacent(i, .,
                                                nb.col= 'rook.adj.nbs'
                               ))))
"

# save adjacencies

wdir <-
  '~/R/all sharkey geoseg work/divflow/R/creating SpWs/'

"write.csv(resl
        ,file = paste0(wdir
               ,'spatial-composites/'
               ,'adjacency-weighted-tract-vals.csv')
        ,row.names = F)
"


# map sample area --------------------------------------------------------------
wdir <-
  '~/R/all sharkey geoseg work/divflow/R/creating SpWs/'
resl <- vroom::vroom(paste0(wdir
                        ,'spatial-composites/'
                        ,'adjacency-weighted-tract-vals.csv'))
resl

resl$rook.adj.wval - resl$q

resl %>% count(name)

phl <- resl %>%
  filter(name %in%
           'affluence13_17') %>%
  geox::geosubset(cz_name = "Philadelphia")

phllsf <- phl %>%
  geox::attach.geos() %>%
  select(1, pop,hh, aff13_17 = value, matches('adj.wval')) %>%
  pivot_longer(4:6)

phllsf %>%
  filter(!duplicated(geoid)) %>%
  st_bbox() %>% st_centroid()

phllsf %>%
  filter(name != 'rook.adj.wval') %>%
  ggplot() +
  geom_sf(aes(fill = value)
          ,color = NA) +
  facet_wrap(vars(name)) +
  scale_fill_viridis_c() +
  coord_sf(expand = F) +
  theme_void() +
  theme(legend.position = 'bottom')




# distance-weighted composites -------------------------------------------------

library(units)
# define cut-off distance
dst.ceiling <- units::set_units(20, 'miles')

# see distance decay / SIF fcns in SIF-fcns.R

#' powerlaw.dist2weights
devtools::load_all()
power.law


# get.dist.weighted.composite fcn moved to `spW-fcns`



# split-map-bind ---------------

# 5mi bisq
resl <- resl %>%
  split(.$name) %>%
  map_dfr(
    ~mutate(.,
            bisq.dist.wval_5mi =
              map_dbl(geoid,
                      function(i)
                        get.dist.weighted.composite(
                          i, .
                          ,dist.decay.fcn = bisq.dist2weights
                          ,dist.col= 'dists'
                          ,cutoff = units::set_units(5, 'miles'))
              )
    )
  )

# 2.5mi bisq
resl <- resl %>%
  split(.$name) %>%
  map_dfr(
    ~mutate(.,
            bisq.dist.wval_2.5mi =
              map_dbl(geoid,
                      function(i)
                        get.dist.weighted.composite(
                          i, .
                          ,dist.decay.fcn = bisq.dist2weights
                          ,dist.col= 'dists'
                          ,cutoff = units::set_units(2.5, 'miles'))
              )
    )
  )

# save adjacencies
wdir <-
  '~/R/all sharkey geoseg work/divflow/R/creating SpWs/'
"
write.csv(resl
        ,file = paste0(wdir
               ,'spatial-composites/'
               ,'bisq-dist-weighted-tract-vals.csv')
        ,row.names = F)
"


# map distance-weighted sample area --------------------------------------------

cos <- tigris::counties(state = 42)
cos %>% filter(grepl("Phila", NAME))
sample <- resl %>%
  geox::geosubset(countyfp = '42101') %>%
  filter(name == 'disadvantage2_13_17')

smsf <- sample  %>%
  geox::attach.geos()

smsfl <- smsf %>%
  select(1,value, matches('bisq')) %>%
  tibble() %>%
  pivot_longer(2:4)

smsfl %>%
  st_sf() %>%
  ggplot(aes(fill = value)) +
  geom_sf(color = NA) +
  scale_fill_viridis_c() +
  theme_void() +
  facet_wrap(vars(name)) +
  theme(legend.position = 'bottom')


# scratch ---------------------------------------------------------------------

# check saved
tmp <-
  vroom::vroom(
    file =
      paste0(wdir
             ,'spatial-composites/'
             ,'bisq-dist-weighted-tract-vals.csv')
  )

tmp
spws

# get subsample
phl <- resl %>%
  geox::geosubset(cz='19700')

phl <- phl %>%
  filter(name ==
           'affluence13_17')



phl <- phl %>%
  mutate(bisq.dist.wval_5mi =
            map_dbl(geoid,
                    function(i)
                      get.dist.weighted.composite(
                        i, .
                        ,dist.decay.fcn = bisq.dist2weights
                        ,dist.col= 'dists'
                        ,cutoff = units::set_units(5, 'miles'))
                    )
         )

phl <- phl %>%
  mutate(bisq.dist.wval_15mi =
           map_dbl(geoid,
                   function(i)
                     get.dist.weighted.composite(
                       i, .
                       ,dist.decay.fcn = bisq.dist2weights
                       ,dist.col= 'dists'
                       ,cutoff = units::set_units(15, 'miles'))
           )
  )
phl

phlsf <- phl %>%
  select(1,pop,aff.13_17 = value, matches('bisq')) %>%
  pivot_longer(3:5) %>%
  #geox::geosubset(countyfp = '42101') %>%
  geox::attach.geos()


options(tigris_use_cache = TRUE)
"
philly.layrs <- visaux::add.map.layers(dsafgas
                                       , add.places = 'white'
                                       ,spatial.trim = st_crop)"

phlsf %>%
  ggplot() +
  geom_sf(aes(fill = value)
          ,color = NA) +
  facet_wrap(vars(name)) +
  scale_fill_viridis_c() +
  coord_sf(expand = F) +
  theme_void() +
  theme(legend.position = 'bottom')

# older scratch ----------------------------------------------------------------
i <- phl[1,]$geoid
nbs <- spws %>% filter(geoid %in% i) %>% pull('dists') # %>% unlist()
dists <- tibble(geoid = names(nbs[[1]])
                ,dist.from.i = nbs[[1]] )

js <- phl[phl$geoid %in% dists$geoid, ] %>%
  left_join(dists) %>%
  filter(geoid != i)

js$spatial.weight <-
  bisq.dist2weights(js$dist.from.i
                    ,cutoff = dst.ceiling)
# or: dist.decay.fcn(js$dist.from.i, ...)

# apply both spatial weight and pop or hh weight
spu <- stats::weighted.mean(js$value
                            ,js$weight * js$spatial.weight
                            ,na.rm = T)


js

# hash out code with a single var
tmp <- res %>%
  select(1:5, perc_bl)

chi <- res %>%  #tmp %>%
  geox::geosubset(cz_name="Chicago")

spws

# make long
chil <- chi %>%
  select(-geoid10) %>%
  pivot_longer(6:ncol(.))


# get a weighted column depending on the var
hh.weighted.cols <-
  colnames(res) %>% grep('fam', .,  value = T)
pop.weighted.cols <-
  colnames(res)[6:ncol(res)] %>% grep('^((?!fam).)*$'
                                      , .,  value = T
                                      ,perl = T)

chil <- chil %>%
  mutate(weight = case_when( name %in% hh.weighted.cols ~ hh
                             ,TRUE ~ pop))

chil

#' get.avg.adjacent
#'
#' @param i a given tract geoid
#' @param x full ct data to subset from
#' @param value.col column name in `x` to calculate weighted avg of adjacent tracts
#'   for.
#' @param weight.col column name in `x`, probably denoted population or households
#' @param nbs string for column name in `spws`
#' @param spws spatial weight matrices with a geoid column and other list columns for
#'   different spatial weights
#'
get.avg.adjacent <- function(i,  x
                             , value.col = 'value'
                             , weight.col = 'weight'
                             ,nb.col
                             ,spatial.weights = spws) {


  #browser()
  #x <- x %>% filter(name %in% !!value)
  j.ids <- unlist(spatial.weights[spatial.weights$geoid == i, nb.col])
  js <- x[x$geoid %in% j.ids, ]

  # spatial mean
  spu <- stats::weighted.mean(js[[value.col]]
                              ,js[[weight.col]]
                              ,na.rm = T)

  return(spu)
}

chil <- chil %>%
  split(.$name) %>%
  map_dfr( ~mutate(.,
                   q.adj.wval =
                     map_dbl(geoid,
                             function(i)
                               get.avg.adjacent(i, .,
                                                nb.col= 'queen.adj.nbs'
                               ))))



get.avg.adjacent( '17031010300'
                 ,chil$working_pooled_pooled_mean
                 ,nb.col='queen.adj.nbs')

chil %>%
  group_by(name) %>%
  mutate(q.avg.adj =
           map_dbl( geoid,
                    ~get.avg.adjacent(., )))


#'
#' @param spw Tract ids where spatial weight/adjacency == 1
#' @param x full ct data to subset from
#' @param value.col column to calculate weighted avg of adjacent tracts for.
#' @param weight.col Weights column. Probably population or households
get.avg.adjecent <- function(spw,  x, value.col = 'value', weight.col = 'weight' ) {

  #x <- x %>% filter(name %in% !!value)
  j.ids <- unlist(spw)
  js <- x[x$geoid %in% j.ids, ]

  # spatial mean
  spu <- stats::weighted.mean(js[[value.col]]
                             ,js[[weight.col]])

  return(spu)
}



chi <- res %>%  #tmp %>%
  geox::geosubset(cz_name="Chicago")

chi %>%
  left_join(spws) %>%
  mutate(q_adj_avg.perc_bl =
           map_dbl( queen.adj.nbs,
                    ~get.avg.adjecent(., x = chi,
                                      value.col = 'perc_bl',
                                      weight.col = 'pop')
                    ))

# use for loops because mapping w/in across is confusing

# get a weighted column depending on the var
hh.weighted.cols <-
  colnames(res) %>% grep('fam', .,  value = T)
pop.weighted.cols <-
  colnames(res)[6:ncol(res)] %>% grep('^((?!fam).)*$'
                                      , .,  value = T
                                      ,perl = T)

chi <- chi %>%
  mutate(weight = case_when( name %in% hh.weighted.cols ~ hh
                            ,TRUE ~ pop))


chi %>%
  group_by(name) %>%
  mutate(q.adj.avg =
           map_dbl( queen.adj.nbs,
                    ~get.avg.adjecent(spw = .
                                      ,x=schi$affluence13_17))
  )
# split
schi <- chi %>%
  split(.$name)

schi <- schi %>%
  map( ~left_join(., spws) )


schi$affluence13_17 %>%
  mutate(q.adj.avg =
           map_dbl( queen.adj.nbs,
                    ~get.avg.adjecent(spw = .
                                      ,x=schi$affluence13_17))
  )


schi %>%
  map(
    ~mutate(.,
            q.adj.avg =
             map_dbl( queen.adj.nbs,
                      ~get.avg.adjecent(spw = .
                                        ,x=schi$affluence13_17))
    )
  )



# alternate adj composite apprach
weight.type <- 'queen.adj.nbs'

i <- chi[1,]$geoid
js <- unlist(spws[spws$geoid == i, weight.type])
js <- chi[chi$geoid %in% js, ]

js %>%
  group_by(name) %>%
  mutate(adj.avg =
           stats::weighted.mean(value, weight))


adj.avg.by.var <- function(i, x
                           , weight.col
                           , spatial.weights = spws
                           , var.col = 'name') {

  j.ids <- unlist(spws[spws$geoid == i, weight.type])
  j <- x[x$geoid %in% j.ids, ]

  js %>%
    group


}


