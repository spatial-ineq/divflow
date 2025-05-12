# boilderplate ---------------------------------------------------------------------------
rm(list = ls())
require(tidyverse)
require(sf)

# option setting
sf_use_s2(T)
options(tigris_use_cache = TRUE)


# ws ---------------------------------------------------------------------------

# dropbox dir or della base dir
ddir <- #Sys.getenv('drop_dir')
  '/scratch/gpfs/km31/'


# flow weights
flww.dir <- paste0(ddir
                   ,'/adjacencies+proximities/flow-weights-include-loops/')

list.files(flww.dir
           ,recursive = T)

# get demovars -----------------------------------------------------------------

# (compiled in data-raw/)
demo.pth <- paste0(ddir
                   ,'seg-measures/by tract/broader ineq flows/'
                   ,'res-chars-long.csv')
resl <- vroom::vroom(demo.pth)

# sample -- generate cz/cz composite for sample area --------------------------

"smpl.id <- '19700'
smpl.rt <- 'cz'
smpl <- resl %>% filter(rid == smpl.id
                        & rt == smpl.rt)

# get flwws for sample
flww.pth <- list.files(flww.dir
                       ,recursive = T
                       ,full.names = T
                       ,pattern =
                         paste0(
                           smpl.rt,'-',smpl.id,'.*')
)
flwws <- flww.pth[grepl('ct',flww.pth)] %>% read_rds()
"

# flow composite fcn --------------------------------------------------


#' get.flow.weighted.composite
#'
#' Given an area geoid i, a set of attribute variables for all areas, and a
#' (tibble)-column of spatial weights between areas, calculate a "spatial composite" that
#' incorporates neighbors' attributes, spatial weights, and population or other
#' weights.
#'
#' @inheritParams get.dist.weighted.composite
#'
#' @export get.dist.weighted.composite
get.flow.weighted.composite <- function(i, x
                                        ,prx.col = c('inc.flwws', 'vis.flwws')
                                        ,spatial.weights = flwws
                                       ) {

  prx.col <- prx.col[1]

  # sometimes visited flwws are NULL for areas with 0 population in sfg sample.
  # Other times there are no visited OR incoming flows and the geoid doesn't
  # exist in flwws at all. Return numeric NA in those cases
  if(! i %in% spatial.weights$geoid) return(as.numeric(NA))

  # all neighbors flow weights
  nbs <- spatial.weights %>% filter(geoid %in% i) %>% pull(prx.col)
  nbs <- nbs[[1]]

  # again, return NA when no visited flows
  if(is.null(nbs)) return(as.numeric(NA))

  # organize as tibble if not already
  if(! 'tbl' %in% class(nbs))
    nbs <- tibble(geoid = names(nbs)
                  ,dist.from.i = nbs )

  # combine distances and attributes; KEEP loops in for attribute calculation (where i==j)
  js <- x[x$geoid %in% nbs$j, ] %>%
    left_join(nbs, by = c('geoid' = 'j'))
  # %>% filter(geoid != i)

  # for flow weight, only use flww, and not population or hhs
  spu <- stats::weighted.mean( js$value
                              ,js$flww
                              ,na.rm = T)
  return(spu)

}



#smpl %>%
#  filter(var == 'perc_bl') %>%
#  get.flow.weighted.composite( i = '34001000200'
#                            ,x = .
#                            ,prx.col = 'vis.flwws'
                            #'inc.flwws'
#                            ,spatial.weights = flwws )

# Della wrapper -----------------------------------------------------


#' Della.wrapper_flow.composite.by.region
#'
#' Takes in global env:
#'
#' fcn to get composites (`get.flow.weighted.composite`)
#'
#' and `resl` with tract characteristics, long by variable
#'
#' @param ddir base directory that contains other needed datasets
#' @param demo.pth data.frame long by variable--As prepped in divflow/data-raw/
#' @param flww.base.dir base directory for flow weights (flwws)
#' @param
#'
Della.wrapper_flow.composite.by.region <- function(
  region.type
  ,region.id
  ,ddir =
    '/scratch/gpfs/km31/'
  #Sys.getenv('drop_dir')
  ,save.dir =
    paste0(ddir
           ,'adjacencies+proximities/spatial-composites/flow-composites/')
  ,save.subdir = NULL
  ,flww.base.dir =
    paste0(ddir
           ,'adjacencies+proximities/flow-weights/')
  ,flww.subdir = 'cts-by-region/'
  ) {

  require(tidyverse)

  if(! (exists('resl') & exists('get.flow.weighted.composite')) )
    warning('Wrapper expects resl and composite fcn "get.flow.weighted.composite" in global env')

  # filter to region
  resl <- resl %>%
    filter(
      !!rlang::sym(region.type) == region.id)

  # load pre-calculated flww weights
  flww.dir <- paste0(flww.base.dir,
                     flww.subdir)
  flww.pth <- list.files(flww.dir
                         ,full.names = T
                         ,pattern =
                           paste0(
                             region.type,'-',
                             region.id,'.*')
                         )

  flwws <- flww.pth %>% read_rds()

  # split-map-bind w incoming and vis
  composites <-
    resl %>%
    split(.$var) %>%
    map_dfr( ~mutate(.,
                     inc.flww.composite =
                       map_dbl(geoid,
                               function(i)
                                 get.flow.weighted.composite(
                                   i
                                   , x = .
                                   ,prx.col = 'inc.flwws'
                                   ,spatial.weights = flwws))
                     , vis.flww.composite =
                       map_dbl(geoid,
                               function(i)
                                 get.flow.weighted.composite(
                                   i
                                   , x = .
                                   ,prx.col = 'vis.flwws'
                                   ,spatial.weights = flwws))
    ))

  if(!is.null(save.subdir)) {

    save.dir <- paste0(save.dir,
                       save.subdir)

    if(!dir.exists(save.dir))
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
# Della.wrapper_flow.composite.by.region('cz', '19700')


# fcns to incl in add_global ----------------------------------------------

get.flow.weighted.composite
Della.wrapper_flow.composite.by.region

# params ------------------------------------------------------------------

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
                   ,'adjacencies+proximities/spatial-composites-include-loops/flow-composites/')
          ,save.subdir = 'tracts/'
          ,flww.base.dir =
            paste0(ddir
                   ,'adjacencies+proximities/flow-weights-include-loops/')
          ,flww.subdir = 'cts-by-region/')


cbsa.ct.params <-
  tibble( region.type = 'cbsa'
          ,region.id = gencbsas
          ,ddir =
            '/scratch/gpfs/km31/'
          #Sys.getenv('drop_dir')
          ,save.dir =
            paste0(ddir
                   ,'adjacencies+proximities/spatial-composites-include-loops/flow-composites/')
          ,save.subdir = 'tracts/'
          ,flww.base.dir =
            paste0(ddir
                   ,'adjacencies+proximities/flow-weights-include-loops/')
          ,flww.subdir = 'cts-by-region/')


# send Della jobs ---------------------------------------------------------

# send job
library(rslurm)
cztracts.flowcomposites.dellajob <-
  slurm_apply(f =
                Della.wrapper_flow.composite.by.region,
              params = cz.ct.params,
              jobname = 'cztracts.flowcomposites include loops correct',
              nodes = 22,
              cpus_per_node = 1,
              slurm_options = list(time = '15:00:00',
                                   'mem-per-cpu' = '20G',
                                   'mail-type' = list('begin', 'end', 'fail'),
                                   'mail-user' = 'km31@princeton.edu'),
              global_objects = c('Della.wrapper_flow.composite.by.region'
                                 ,'get.flow.weighted.composite'
                                 ,'resl')
  )



cbsatracts.flowcomposites.dellajob <-
  slurm_apply(f =
                Della.wrapper_flow.composite.by.region,
              params = cbsa.ct.params,
              jobname = 'cbsatracts.flowcomposites include loops correct',
              nodes = 19,
              cpus_per_node = 1,
              slurm_options = list(time = '10:00:00',
                                   'mem-per-cpu' = '10G',
                                   'mail-type' = list('begin', 'end', 'fail'),
                                   'mail-user' = 'km31@princeton.edu'),
              global_objects = c('Della.wrapper_flow.composite.by.region'
                                 ,'get.flow.weighted.composite'
                                 ,'resl')
  )




# check generation --------------------------------------------------------

sdir <- paste0( '/scratch/gpfs/km31/'
               ,'adjacencies+proximities/spatial-composites/flow-composites/'
               ,'tracts/'
               )

gend.cbsas <- list.files(sdir,pattern = '^cbsa')
gend.czs <- list.files(sdir,pattern = '^cz')

geox::rpops %>% filter(pop > 25e3) %>% count(rt)
gend.czs %>% length()
gend.cbsas  %>% length()

# I Am missing Flwws for some cbsas....

# scratch ----------------------------------------------------------------------



flww.pths <- list.files(flww.dir
                        ,recursive = T
                        ,full.names = T
                        ,pattern='cz-[0-9]{5}'
)
flwws <- flww.pths %>% head() %>% map_dfr(read_rds)
flwws




gend.cbsas <-stringr::str_extract(gend.cbsas, '[0-9]{5}')
ungend <- gencbsas[!gencbsas %in% gend.cbsas]
ungend <- ungend[!is.na(ungend)]
ungend <- geox::rx %>% select(matches('cbsa')) %>% distinct() %>% filter(cbsa %in% ungend)
geox::rpops %>% filter(rid %in% ungend$cbsa & rt == 'cbsa') %>% filter(pop >= 25e3) %>% geox::add.rns()


# check one that didn't generate..?
albtst <- Della.wrapper_flow.composite.by.region(region.type = 'cbsa'
                                                      , region.id = '10700'
)
albtst
# --> Concl for ungenerated CBSAs
# they are missing from resl or flwws; would have to go back

# for CZs
length(gend.czs)
nrow(geox::rpops %>% filter(rt == 'cz' & pop >= 25e3))
gend.czs <-stringr::str_extract(gend.czs, '[0-9]{5}')
ungend <- genczs[!genczs %in% gend.czs]
ungend <- ungend[!is.na(ungend)]
ungend <- geox::rx %>% select(matches('cz')) %>% distinct() %>% filter(cz %in% ungend)
geox::rpops %>% filter(rid %in% ungend$cz & rt == 'cz') %>% filter(pop >= 250e3) %>% geox::add.rns()

resl %>% geox::geosubset(cz = '34901')
ocalatst <- Della.wrapper_flow.composite.by.region(region.type = 'cz'
                                                   , region.id = '07800'
)
ocalaresl <- resl %>% filter(cz == '07800')
noflows.ocala <- ocala %>% filter(! geoid %in% flwws)

ocalaresl %>% get.flow.weighted.composite(
  i
  , x = .
  ,prx.col = 'inc.flwws'
  ,spatial.weights = flwws))
