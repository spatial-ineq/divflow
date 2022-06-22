# ws ---------------------------------------------------------------------------
rm(list = ls())
require(tidyverse)
require(sf)

# option setting
sf_use_s2(T)
options(tigris_use_cache = TRUE)

# load spatial weight matrices -------------------------------------------------

# dropbox dir or della base dir
ddir <- #Sys.getenv('drop_dir')
  '/scratch/gpfs/km31/'

spw.dir <- paste0(ddir, 'adjacencies+proximities/')
spws <-
  list.files(spw.dir,
             pattern='tract-adjacencies.rds'
             ,full.names = T) %>%
  readRDS()

spws

# trim cols
spws <- spws %>% select(-matches('adj'))

# get demovars -----------------------------------------------------------------

# (compiled in data-raw/)
demo.pth <- paste0(ddir
                   ,'seg-measures/by tract/broader ineq flows/'
                   ,'res-chars-long.csv')
resl <- vroom::vroom(demo.pth)

# get dst weighted avgs --------------------------------------------------



#' get.dst.weighted.composite
#'
#' Given an area geoid, a set of attribute variables for all areas, and a list-column
#' of spatial weights between areas, calculate a "spatial composite" that
#' incorporates neighbors' attributes, spatial weights, and population or other
#' weights.
#'
#' Note:: I think I can make this be a uniform fcn for both flow- and
#' proximity/distance wegiths. But I generated the proximity distance weights with
#' list-colms in instead of nested tibble-colms. So they are separate fcns for now.
#'
#' @param i a given tract geoid
#' @param x full ct data to subset from
#' @param value.col Attribute column name in `x` to calculate weighted avg of
#'   adjacent tracts for.
#' @param weight.col column name in `x` for non spatial weights. Probably population
#'   or households
#' @param dist.decay.fcn a distance decay or proximity fcn to convert distances to
#'   spatial weights
#' @param dist.col string for column name in `spws` that pairs geoids with distances
#'   from i
#' @param spws spatial weight matrices with a geoid column and other list columns for
#'   spatial weights
#' @param ... passed on to `dst.decay.fcn`
#'
get.dst.weighted.composite <- function(i,  x
                                       , value.col = 'value'
                                       , weight.col = 'weight'
                                       , dist.decay.fcn
                                       ,dist.col = 'dists'
                                       ,spatial.weights = spws
                                       ,...) {

  # all neighbors within dists of i, based on supplied spatial weights
  nbs <- spatial.weights %>% filter(geoid %in% i) %>% pull(dist.col)
  # organize as tibble
  nbs <- tibble(geoid = names(nbs[[1]])
                ,dist.from.i = nbs[[1]] )

  # combine distances and attributes; filter out loops (where i==j)
  js <- x[x$geoid %in% nbs$geoid, ] %>%
    left_join(nbs, by = 'geoid') %>%
    filter(geoid != i)

  # convert distance to weight w/ proximity fcn
  js$spatial.weight <-
    dist.decay.fcn( as.numeric(js$dist.from.i)
                    ,...)

  # apply both spatial weight and pop or hh weight to get spatial mean
  spu <- stats::weighted.mean(js$value
                              ,js$weight * js$spatial.weight
                              ,na.rm = T)
  return(spu)
}


# some simple decay fcns included in pkg:
#divflow::neg.exp
#divflow::bisq.dist2weights

# sample run
"resl %>%
  filter(var == 'perc_hsp') %>%
  get.dst.weighted.composite(
    i = '46013951300'
    ,x = .
    ,dist.decay.fcn = #neg.exp
      bisq.dist2weights
    ,cutoff = 4
  )
resl %>% filter(geoid == '46013951300' & var == 'perc_hsp')
"

# Della wrapper for dsts -------------------------------------------

#' Della.wrapper_dst.composite.by.region
#'
#' Takes in global env: `spws`, with columns for `dists` for each tract `geoid`
#' (distances of other areas js within distance cutoff of 20 miles)
#'
#' and `resl` with tract characteristics, long by variable
#'
#' @param ddir base directory that contains other needed datasets
#' @param save.dir where to save calculatedcomposites
#'
Della.wrapper_dst.composite.by.region <- function(
  region.type
  ,region.id
  ,ddir =
    '/scratch/gpfs/km31/'
  #Sys.getenv('drop_dir')
  ,save.dir =
    paste0(ddir
           ,'adjacencies+proximities/spatial-composites/distance-composites/')
  ,save.subdir = NULL
) {

  require(tidyverse)

  if(! (exists('resl') & exists('spws')) )
    warning('Wrapper expects resl and spws in global env')

  resl <- resl %>%
    filter(
      !!rlang::sym(region.type) == region.id)


  # define dist-decay fcns
  bisq.dist2weights <- function(d, cutoff = dst.ceiling
                                , ... ) {
    if_else( d >= cutoff
             ,0
             ,(1-as.numeric(d/cutoff)^2)^2)
  }

  neg.exp <- function(d, ...) {
    exp(-d)
  }


  # check fcn in env
  if(!exists('get.dst.weighted.composite'))
    warning('composite fcn "get.dst.weighted.composite" missing from env')

  # split-map-bind w incoming and vis
  composites <-
    resl %>%
    split(.$var) %>%
    map_dfr( ~mutate(.,
                     negexp.dst.composite =
                       map_dbl(geoid,
                               function(i)
                                 get.dst.weighted.composite(
                                   i,  x = .
                                   , value.col = 'value'
                                   , weight.col = 'weight'
                                   , dist.decay.fcn = neg.exp
                                   ,dist.col = 'dists'
                                   ,spatial.weights = spws))

                     , bisq.8km.composite =
                       map_dbl(geoid,
                               function(i)
                                 get.dst.weighted.composite(
                                   i,  x = .
                                   , value.col = 'value'
                                   , weight.col = 'weight'
                                   , dist.decay.fcn = bisq.dist2weights
                                   ,dist.col = 'dists'
                                   ,spatial.weights = spws
                                   ,cutoff = 8 # abt 5 mi
                                   ))
    ))

  if(!is.null(save.subdir)) {

    save.dir <- paste0(save.dir,
                       save.subdir)

    if(!exists(save.dir))
      dir.create(save.dir, recursive = T)

    rids <- divflow::region2identifiers( region.type = region.type
                                ,region.id = region.id)

    save.pth <- paste0(save.dir, rids, '.csv')

    write.csv(composites
              ,save.pth)
  }
  return(composites)

}


# test run
"mdbtst <- Della.wrapper_dst.composite.by.region(region.type = 'cz'
                                                ,region.id = '00301')
# looks great:)
mdbtst"



# params ----------------------------------------------

# remove smaller regions
rpops <- geox::rpops %>% filter(pop >= 25e3)
# regions
genczs <- rpops %>% filter(rt == 'cz') %>% pull(rid)
gencbsas <- rpops %>% filter(rt == 'cbsa') %>% pull(rid)

cz.ct.params <-
  tibble( region.type = 'cz'
          ,region.id = genczs
          ,ddir =
            '/scratch/gpfs/km31/'
          #Sys.getenv('drop_dir')
          ,save.dir =
            paste0(ddir
                   ,'adjacencies+proximities/spatial-composites/distance-composites/')
          ,save.subdir = 'tracts/'
  )

cbsa.ct.params <-
  tibble( region.type = 'cbsa'
          ,region.id = gencbsas
          ,ddir =
            '/scratch/gpfs/km31/'
          #Sys.getenv('drop_dir')
          ,save.dir =
            paste0(ddir
                   ,'adjacencies+proximities/spatial-composites/distance-composites/')
          ,save.subdir = 'tracts/'
  )


# send Della jobs ---------------------------------------------------------

# send job
library(rslurm)
cztracts.dstcomposites.dellajob <-
  slurm_apply(f =
                Della.wrapper_dst.composite.by.region,
              params = cz.ct.params,
              jobname = 'cztracts.dstcomposites',
              nodes = 22,
              cpus_per_node = 1,
              slurm_options = list(time = '15:00:00',
                                   'mem-per-cpu' = '10G',
                                   'mail-type' = list('begin', 'end', 'fail'),
                                   'mail-user' = 'km31@princeton.edu'),
              global_objects = c('Della.wrapper_dst.composite.by.region'
                                 ,'get.dst.weighted.composite'
                                 ,'resl'
                                 ,'spws')
  )



cbsatracts.dstcomposites.dellajob <-
  slurm_apply(f =
                Della.wrapper_dst.composite.by.region,
              params = cbsa.ct.params,
              jobname = 'cbsatracts.dstcomposites',
              nodes = 19,
              cpus_per_node = 1,
              slurm_options = list(time = '10:00:00',
                                   'mem-per-cpu' = '10G',
                                   'mail-type' = list('begin', 'end', 'fail'),
                                   'mail-user' = 'km31@princeton.edu'),
              global_objects = c('Della.wrapper_dst.composite.by.region'
                                 ,'get.dst.weighted.composite'
                                 ,'resl'
                                 ,'spws')
  )


# check genr --------------------------------------------------------------
sdir <-
  paste0('/scratch/gpfs/km31/'
         ,'adjacencies+proximities/spatial-composites/distance-composites/'
         ,'tracts/'
         )
fns <- sdir %>% list.files()
fns[grepl('cz', fns)] %>% length()
fns[grepl('cbsa', fns)] %>% length()
geox::rpops %>% filter(pop >= 25e3) %>% count(rt)
fns[grepl('cbsa', fns)][1] %>% paste0(sdir, .) %>% vroom::vroom()

# scratch -----------------------------------------------------------------

# does units matter for this dst decay process?
library(units)
vs <- c(3,2.9,3.5)
xs <- c(5,10,15)
ms <- set_units(xs, 'miles')
kms <- set_units(ms, 'km')
ms <- neg.exp(as.numeric(ms))
kms <- neg.exp(as.numeric(kms))
stats::weighted.mean(vs, ms)
stats::weighted.mean(vs, kms)
