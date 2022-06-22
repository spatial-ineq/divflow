# ws ---------------------------------------------------------------------------
rm(list = ls())
require(tidyverse)
require(sf)

# option setting
sf_use_s2(T)
options(tigris_use_cache = TRUE)


# load spatial weight matrices -------------------------------------------------

# dropbox dir or della base dir
ddir <- Sys.getenv('drop_dir')
  #'/scratch/gpfs/km31/'

spw.dir <- paste0(ddir, 'adjacencies+proximities/')
spws <-
  list.files(spw.dir,
             pattern='tract-adjacencies.rds'
             ,full.names = T) %>%
  readRDS()

spws
spws %>% filter( map_dbl(queen.adj.nbs, length) !=
                    map_dbl(rook.adj.nbs, length))

# get demovars -----------------------------------------------------------------

# (compiled in data-raw/)
demo.pth <- paste0(ddir
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



# Della wrapper for adjacencies -------------------------------------------

#' Della.wrapper_adj.composite.by.region
#'
#' Takes in global env:
#' `spws`, with columns for `queen.adj.nbs`
#' and `rook.adj.nbs` for each tract `geoid`
#'
#' and `resl` with tract characteristics, long by variable
#'
#' @param ddir base directory that contains other needed datasets
#' @param save.dir where to save calculatedcomposites
#'
Della.wrapper_adj.composite.by.region <- function(
  region.type
  ,region.id
  ,ddir =
    '/scratch/gpfs/km31/'
  #Sys.getenv('drop_dir')
  ,save.dir =
    paste0(ddir
           ,'adjacencies+proximities/spatial-composites/adjacency-composites/')
  ,save.subdir = NULL
  ) {

  require(tidyverse)

  if(! (exists('resl') & exists('spws')) )
    warning('Wrapper expects resl and spws in global env')

  resl <- resl %>%
    filter(
      !!rlang::sym(region.type) == region.id)

  # check fcn in env
  if(!exists('get.avg.adjacent'))
    warning('composite fcn missing from env')

  # split-map-bind w incoming and vis
  composites <-
    resl %>%
    split(.$var) %>%
    map_dfr( ~mutate(.,
                     rook.adj.composite =
                       map_dbl(geoid,
                               function(i)
                                 get.avg.adjacent(
                                   i , x = .
                                   , value.col = 'value'
                                   , weight.col = 'weight'
                                   ,nb.col = 'rook.adj.nbs'
                                   ,spatial.weights = spws))
                     , qn.adj.composite =
                       map_dbl(geoid,
                               function(i)
                                 get.avg.adjacent(
                                   i , x = .
                                   , value.col = 'value'
                                   , weight.col = 'weight'
                                   ,nb.col = 'queen.adj.nbs'
                                   ,spatial.weights = spws))
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
mdbtst <- Della.wrapper_adj.composite.by.region(region.type = 'cz'
                                       ,region.id = '00301')

# looks great:)
mdbtst



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
                  ,'adjacencies+proximities/spatial-composites/adjacency-composites/')
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
                   ,'adjacencies+proximities/spatial-composites/adjacency-composites/')
          ,save.subdir = 'tracts/'
  )


# send Della jobs ---------------------------------------------------------

# send job
library(rslurm)
cztracts.adjcomposites.dellajob <-
  slurm_apply(f =
                Della.wrapper_adj.composite.by.region,
              params = cz.ct.params,
              jobname = 'cztracts.adjcomposites',
              nodes = 22,
              cpus_per_node = 1,
              slurm_options = list(time = '15:00:00',
                                   'mem-per-cpu' = '10G',
                                   'mail-type' = list('begin', 'end', 'fail'),
                                   'mail-user' = 'km31@princeton.edu'),
              global_objects = c('Della.wrapper_adj.composite.by.region'
                                 ,'get.avg.adjacent'
                                 ,'resl'
                                 ,'spws')
  )



cbsatracts.adjcomposites.dellajob <-
  slurm_apply(f =
                Della.wrapper_adj.composite.by.region,
              params = cbsa.ct.params,
              jobname = 'cbsatracts.adjcomposites',
              nodes = 19,
              cpus_per_node = 1,
              slurm_options = list(time = '10:00:00',
                                   'mem-per-cpu' = '10G',
                                   'mail-type' = list('begin', 'end', 'fail'),
                                   'mail-user' = 'km31@princeton.edu'),
              global_objects = c('Della.wrapper_adj.composite.by.region'
                                 ,'get.avg.adjacent'
                                 ,'resl'
                                 ,'spws')
  )



# check genr --------------------------------------------------------------

sdir <-
  paste0( '/scratch/gpfs/km31/'
          ,'adjacencies+proximities/spatial-composites/adjacency-composites/'
          ,'tracts/' )

gend.cbsas <- list.files(sdir,pattern = '^cbsa')
gend.czs <- list.files(sdir,pattern = '^cz')

gend.czs %>% length()
gend.cbsas  %>% length()
geox::rpops %>% filter(pop >= 25e3) %>% count(rt)

smplcz <- gend.czs[1] %>% paste0(sdir, .) %>% vroom::vroom()
smplcbsa <- gend.cbsas[1] %>% paste0(sdir, .) %>% vroom::vroom()

adjczs <- list.files(sdir
                     ,pattern = '^cz'
                     ,full.names = T) %>%
  map_dfr( vroom::vroom )


adjcbsas <- list.files(sdir
                     ,pattern = '^cbsa'
                     ,full.names = T) %>%
  head(30) %>%
  map( vroom::vroom )

smplcz
smplcbsa

smplcz %>% count(var)
smplcbsa %>% count(var)
