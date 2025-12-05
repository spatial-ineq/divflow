# ws ---------------------------------------------------------------------------
rm(list = ls())
require(tidyverse)


# dirs --------------------------------------------------------------------

basedir <- '/scratch/gpfs/km31/'

cdir <- paste0(basedir
               ,'/adjacencies+proximities/spatial-composites-urbsub/')

composite.dirs <-
  list.dirs(cdir
          ,recursive = F)
# get directories ending in "-composites"
composite.dirs <- composite.dirs[grepl('-composites$',composite.dirs)]

composite.dirs

composite.dirs[1] %>% list.files(recursive = T) %>% head()

# only tract-level generated rn
composite.dirs <- paste0(composite.dirs, '/tracts/')

composite.dirs
# helper function ---------------------------------------------------------

#' clean.composites
#'
#' I saved them with all a uniform structure. This improves on save structure by
#' making long by region type so czs/cbsas can be saved in single file, dropping
#' rownumber colm, and fixing the geoids as uniform-width character strings.
#'
load.clean.bind <- function(path, r.t = c('cz' , 'cbsa')) {

  r.t <- r.t[1]
  pattern <- paste0(r.t, '.*csv$')

  # load
  fns <- path %>%
    list.files(pattern = pattern, full.names = T)

  # make it fast for testing
  # fns <- head(fns)

  x <- map(fns, vroom::vroom )

  # files often have a xwalk, but sp composites generated based on region
  # type. Only keep identifier for type used in composites
  drop.cols <- c('...1', setdiff(c('cz','cbsa'), r.t))

  # trim cols, make long, fix geoids
  x <- x %>%
    map( ~select( .,
                  -any_of(drop.cols))) %>%
    map( ~geox::region.reorg( ., r.t )) %>%
    map( ~mutate(., across( c(geoid, rid)
                           ,geox::fix.geoid )))

  x <- do.call('rbind', x)
  return(x)
}



# to check sample -------------------------------------------------------------
#composite.dirs[1] %>% list.files(pattern = 'cz', full.names = T) %>%
#  `[`(1) %>% vroom::vroom()


# read clean bind by region/dir -------------------------------------------

spcz <- composite.dirs %>%
  map( ~load.clean.bind(., 'cz'))

spcbsa <- composite.dirs %>%
  map( ~load.clean.bind(., 'cbsa'))

# reduce
spcz <- spcz %>% purrr::reduce(full_join)
spcbsa <- spcbsa %>% purrr::reduce(full_join)

# row bind
spcs <- rbind(spcz
      ,spcbsa)

spcs <- spcbsa


# checks ------------------------------------------------------------------

spcs %>% map( summary )

spcs

spcs %>% count(var)
spcs %>% count(var, rt)
spcs %>% map_dbl( ~sum(is.na(.)) / nrow(spcs) )
spcs %>%
  filter(is.na(inc.flww.composite)) %>%
  count(rt)
# right because I didn't calculate flow weights for many CBSAs yet
spcs %>%
  filter(is.na(vis.flww.composite)) %>%
  count(var) %>% arrange(desc(n))

spcs %>%
  filter(is.na(value)) %>%
  count(rt,rid)


# save --------------------------------------------------------------------

write.csv(spcs,
          file = paste0(cdir
                        ,'full-spatial-composites-by-cbsa-urbsub.csv'
                        )
          ,row.names = F)

