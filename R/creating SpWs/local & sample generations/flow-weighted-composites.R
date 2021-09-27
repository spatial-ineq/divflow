# boilderplate ---------------------------------------------------------------------------
rm(list = ls())
require(tidyverse)
require(sf)

# option setting
sf_use_s2(T)
options(tigris_use_cache = TRUE)


# ws ---------------------------------------------------------------------------

# dropbox dir
ddir <- Sys.getenv('drop_dir')

# flow-weights
flww.dir <- paste0(ddir
                   ,'/adjacencies+proximities/flow-weights/')
list.dirs(flww.dir)
list.files(flww.dir
           ,recursive = T)

# get demovars -----------------------------------------------------------------

# (compiled in data-raw/)
demo.pth <- paste0(Sys.getenv('drop_dir')
                   ,'seg-measures/by tract/broader ineq flows/'
                   ,'res-chars-long.csv')
resl <- vroom::vroom(demo.pth)

# sample -- generate cz/cz composite for sample area --------------------------

smpl.id <- '19700'
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
flwws <- flww.pth %>% read_rds()


"# plot for kicks
smpl %>%
  filter(var %in% unique(smpl$var)[1:6]) %>%
  geox::attach.geos() %>%
  ggplot() +
  geom_sf(aes(fill = value
              ,alpha = weight)
          ,color = NA) +
  scale_alpha_continuous(range = c(.5,1)) +
  scale_fill_viridis_c() +
  facet_wrap(vars(var)) +
  theme_void()
"




# hash out flow composite fcn --------------------------------------------------


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
get.prx.weighted.composite <- function(i, x
                                        ,prx.col = c('inc.flwws', 'vis.flwws')
                                        ,spatial.weights = flwws
                                        ,...
) {

  prx.col <- prx.col[1]
  # all neighbors within dists of i, based on supplied spatial weights
  nbs <- spatial.weights %>% filter(geoid %in% i) %>% pull(prx.col)
  nbs <- nbs[[1]]

  # sometimes visited flwws are NULL for areas with 0 population in sfg sample;
  # return numeric NA in those cases
  if(is.null(nbs))
    return(as.numeric(NA))

  # organize as tibble if not already
  if(! 'tbl' %in% class(nbs))
    nbs <- tibble(geoid = names(nbs)
                  ,dist.from.i = nbs )

  # combine distances and attributes; filter out loops (where i==j)
  js <- x[x$geoid %in% nbs$j, ] %>%
    left_join(nbs, by = c('geoid' = 'j')) %>%
    filter(geoid != i)

  # for flow weight, only use flww, and not population or hhs
  spu <- stats::weighted.mean(js$value
                              ,js$flww
                              ,na.rm = T)
  return(spu)

}



smpl %>%
  filter(var == 'perc_bl') %>%
  get.prx.weighted.composite( i = '34001000200'
                            ,x = .
                            ,prx.col = 'vis.flwws'
                            #'inc.flwws'
                            ,spatial.weights = flwws
                            )
smpl %>%
  filter(var == 'perc_bl') %>%
  split(.$var) %>%
  map_dfr( ~mutate(.,
                   inc.flww.composite =
                     map_dbl(geoid,
                             function(i)
                               get.prx.weighted.composite(i
                                                          , x = .
                                                          ,prx.col = 'inc.flwws'
                                                          ,spatial.weights = flwws)
                             )))




# hash out for all regions -----------------------------------------------------


#' Della.wrapper_flow.composite.by.region
#'
#'
#' Currently thinking i'll include fcn to get composites
#' (`get.prx.weighted.composite`) in global env rather ttan a pkg
#'
#' @param ddir base directory that contains other needed datasets
#' @param demo.pth data.frame long by variable--As prepped in divflow/data-raw/
#' @param flww.base.dir base directory for flow weights (flwws)
#' @param
#'
Della.wrapper_flow.composite.by.region <- function( region.type
                                                   ,region.id
                                                   ,ddir = Sys.getenv('drop_dir')
                                                   ,save.dir =
                                                     paste0(ddir
                                                            ,'/adjacencies+proximities/spatial-composites/flow-composites/')
                                                   ,save.subdir = NULL
                                                   ,demo.pth =
                                                     paste0(ddir
                                                            ,'seg-measures/by tract/broader ineq flows/'
                                                            ,'res-chars-long.csv')
                                                   ,flww.base.dir =
                                                     paste0(ddir
                                                            ,'/adjacencies+proximities/flow-weights/')
                                                   ,flww.subdir = 'cts-by-region/') {

  require(tidyverse)

  # load long demographics
  resl <- vroom::vroom(demo.pth)
  # filter to region
  resl <- resl %>%
    filter(rid == region.id
           & rt == region.type)

  # load pre-calculated flww weights
  flww.dir <- paste0(flww.base.dir,
                     flww.subdir)
  flww.pth <- list.files(flww.dir
                         ,full.names = T
                         ,pattern =
                           paste0(
                             region.type,'-',region.id,'.*')
  )

  flwws <- flww.pth %>% read_rds()

  # check fcn in env
  if(!exists('get.prx.weighted.composite'))
    warning('composite fcn missing from env')

  # split-map-bind w incoming and vis
  composites <-
    resl %>%
    split(.$var) %>%
    map_dfr( ~mutate(.,
                     inc.flww.composite =
                       map_dbl(geoid,
                               function(i)
                                 get.prx.weighted.composite(
                                   i
                                   , x = .
                                   ,prx.col = 'inc.flwws'
                                   ,spatial.weights = flwws))
                     , vis.flww.composite =
                       map_dbl(geoid,
                               function(i)
                                 get.prx.weighted.composite(
                                   i
                                   , x = .
                                   ,prx.col = 'vis.flwws'
                                   ,spatial.weights = flwws))
    ))

  if(!is.null(save.subdir)) {

    save.dir <- paste0(save.dir,
                       save.subdir)

    if(!exists(save.dir))
      dir.create(save.dir)

    rids <- geox::add.rns(tibble(rt = region.type, rid = region.id ))
    rids <- rids %>% paste0(collapse = '-')

    save.pth <- paste0(save.dir, rids, '.csv')

    write.csv(composites
              ,save.pth)
  }
  return(composites)
}


# test run
Della.wrapper_flow.composite.by.region('cz'
                                       ,'19700')


# scratch ----------------------------------------------------------------------


smpl <- resl  %>% geox::geosubset(cz = '19700')
flwws

tmp <- smpl %>% left_join(flwws)

untmp <- tmp %>%
  unnest(cols = inc.flwws)
untmp <- untmp %>% filter(geoid != j)

# experimenting with faster, tidier way, but I think only solution would be
# tidygraph (because I'm effectively subsetting based on ties--- using js to
# subset from attribute df). But Tidygraph also isn't fast. So i will just keep
# my manual little slow fcn.
tst <- untmp %>%
  filter(var %in% c('perc_hsp'
                    ,'perc_bl')) %>%
  group_by(geoid,
           var) %>%
  summarise(inc.flww.composite =
              stats::weighted.mean(value, flww))


orig <- smpl %>%
  filter(var == 'perc_bl') %>%
  mutate(inc.flww.composite =
           map_dbl(geoid,
          function(i)
            get.flow.weighted.composite(i
                                       , x = .
                                       ,prx.col = 'inc.flwws'
                                       ,spatial.weights = flwws)))
tst
orig
