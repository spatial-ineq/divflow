# ws ---------------------------------------------------------------------------

rm(list = ls())

require(tidyverse)
require(sf)

# option setting
sf_use_s2(T)
options(tigris_use_cache = TRUE)

# dropbox dir or della base dir
ddir <- '/scratch/gpfs/km31/' #Sys.getenv('drop_dir')


devtools::load_all()



# workflow ----------------------------------------------------------------

# https://rpubs.com/kmc39/ergmdiv


# wrapper fcn -------------------------------------------------------------


Della.wrapper_region2ergm <- function( cz = NULL
                                       ,cbsa = NULL
                                       ,min.flows = 1
                                       ,tie.str.deciles = 5
                                       ,tracts.or.groups = 'bg' #'ct'
                                       ,trim.loops = T
                                       ,directed = F
                                       ,year = 2019
                                       ,ddir = '/scratch/gpfs/km31/'
                                       ,sfg.dir = 'sfg-processed/orig_dest_annual/'
                                       ,base.save.dir =paste0(ddir, 'ergms-newer/')
                                       ,run.ergms = T
) {
  # ws
  require(tidyverse)
  require(sf)
  devtools::load_all(path = '/home/km31/all/divflow/')

  # option setting
  sf_use_s2(T)
  options(tigris_use_cache = TRUE)
  set.seed(1234) # seed for ergm

  rids <- geox::get.region.identifiers(cz = cz, cbsa = cbsa)
  rids <- paste0(rids, collapse = '-')


  if(!dir.exists(base.save.dir))
    dir.create(base.save.dir, recursive = T)

  # browser()

  ## get sfg/graph objects
  sfg <- rid2sfg(cz = cz, cbsa = cbsa
                 ,min.flows = min.flows
                 ,tracts.or.groups = tracts.or.groups
                 ,trim.loops = trim.loops
                 ,year = year
                 ,ddir = ddir
                 ,sfg.dir = sfg.dir
  )

  gh <- sfg2gh(sfg = sfg
               ,directed = directed)

  ghsf <- spatialize.graph( gh = gh
                            , frame.sf = NULL
                            , tracts.or.groups = tracts.or.groups
                            , directed = directed
                            , year = year)


  ## get dividedeness
  .countyfps <- geox::x2cos(cz=cz, cbsa = cbsa)

  if(!is.null(cz))
    rsf <- geox::build.CZs(cz) %>%
    mutate(rid = cz, rt = 'cz') %>%
    divM::conic.transform()
  else
    rsf <- tigris::core_based_statistical_areas(year = 2019) %>%
    rename_with(tolower) %>%
    mutate(rid = cbsa, rt = 'cbsa') %>%
    filter(geoid %in% rid) %>%
    divM::conic.transform()

  rds <- map_dfr(.countyfps
                 ,~ tigris::roads(state = substr(.x,1,2)
                                  ,county = substr(.x,3,5)
                                  ,year = year)
  ) %>%
    rename_with(tolower) %>%
    divM::conic.transform() %>%
    st_crop(rsf)

  interstates <- rds %>%
    filter(rttyp %in% 'I')

  arterials <- rds %>%
    filter(mtfcc %in%
             Hmisc::Cs(S1100, S1200, S1630))

  rm(rds)

  if(tracts.or.groups[1] == 'ct')
    tigris.call <- tigris::tracts
  else
    tigris.call <- tigris::block_groups

  int.div <-
    divM::gen.cross.tract.dividedness(
      region = rsf
      ,divs = interstates
      ,nbd.query.fcn = tigris.call
      ,year = year
      ,region.id.colm = 'rid'
      ,erase.water = F
    ) %>%
    select(geoid
           , int.poly = poly.id)

  ar.div <-
    divM::gen.cross.tract.dividedness(
      region = rsf
      ,divs = arterials
      ,nbd.query.fcn = tigris.call
      ,year = year
      ,region.id.colm = 'rid'
      ,erase.water = F
    ) %>%
    select(geoid
           , art.poly = poly.id)

  # save div xwalk
  ergm.divs <- full_join( int.div
                          ,ar.div)
  div.dir <- paste0(base.save.dir,
                    'nbhd-divs/')

  if(!dir.exists(div.dir))
    dir.create(div.dir, recursive = T)

  write.csv( ergm.divs
             ,file = paste0(base.save.dir,
                            'nbhd-divs/', rids, '.csv')
             ,row.names = F)

  # merge in w/ nodes
  ghsf <-
    ghsf %>%
    activate('nodes') %>%
    left_join(ergm.divs
              , by = c('name' = 'geoid'))


  ## get other tract-level characteristics
  demos <- sfg.seg::demos.2019 %>%
    select(1,2,matches('wh$|bl$|hsp$'))

  if(tracts.or.groups[1] == 'ct')
    demos <- demos %>% sfg.seg::cbg.demos2tracts()

  ghsf <-
    ghsf %>%
    activate('nodes') %>%
    left_join(demos
              , by = c('name' = 'geoid'))


  # save a visual

  trimmed.ghsf <- apply.flow.filters(ghsf
                                     , frame.sf = NULL
                                     ,tie.str.drop.deciles = 6
                                     ,directed = directed
                                     ,min.flows = min.flows
                                     ,min.tie.str = NULL
                                     ,max.dst = NULL
  ) %>%
    st_transform(4326)

  sttm <- visaux::get.stamen.bkg(st_transform(rsf
                                              ,4326)
                                 ,zoom = 11)

  lonlats <- gh2coords(trimmed.ghsf)

  require(ggraph)
  ghplot <-
    ggmap(sttm
          , base_layer =
            ggraph(trimmed.ghsf # ghsf #
                   , layout = lonlats )
    ) +
    #geom_edge_density(fill = '#00D0D0') +
    geom_edge_fan( aes( edge_alpha = tstr
                        ,edge_width = tstr)
                   ,color = '#008080'
    ) +
    flow.map.base() +
    visaux::bbox2ggcrop(rsf) +
    labs(caption = rids)


  vis.svdir <- paste0(base.save.dir, 'graph-maps/')
  if(!dir.exists(vis.svdir))
    dir.create(vis.svdir, recursive = T)

  #visaux::ragg.wrapper(
  visaux::ggsave.hirez(fn = rids
                       , plot = ghplot
                       , dir = vis.svdir)


  # setup SIF terms
  sfn <- ghsf %>% gh2nodes()
  sfe <- ghsf %>% gh2edges()
  # turn distance into numerics representing km
  sfe <- sfe %>% mutate(dst = as.numeric(dst) / 1000 )

  # browser()

  dst.mat <- build.pairwise.dst.matrix(sfn)
  dst.freq <-
    relative.tie_dst.frequency(
      sfn, sfe
      , dst.mat = dst.mat
      , n.bins = 60
    )

  param.fitting <-
    stats::optim(par= c(1,1.7,1, 10)
                 , fn = power.law.loss.fcn
                 , dst.freq = dst.freq
                 , hessian = T)

  decay.mat <- dst.mat %>%
    power.law(param.fitting$par[1],
              param.fitting$par[2],
              param.fitting$par[3])

  # predicted tie probabilities, based on SIF
  dst.freq$flow.hat <-
    power.law(dst.freq$dst.lvl,
              param.fitting$par[1],
              param.fitting$par[2],
              param.fitting$par[3]
    )

  formated.params <-
    param.fitting$par %>%
    round(digits = 4) %>%
    paste0(collapse = ' - ')

  sif.vis <-
    dst.freq %>%
    select(dst.lvl, actual.flows = avg.flow, fitted.power.law = flow.hat) %>%
    pivot_longer(cols = matches("actual|fitted"),
                 values_to = "average.flow") %>%
    ggplot() +
    geom_path(aes(color = name,
                  y = average.flow,
                  x = dst.lvl)) +
    labs(title = "Relative avg flows",
         subtitle = "actual vs. fitted power law",
         caption = glue::glue('SIF params: {formated.params}\n{rids}'))

  # save SIF terms..
  sif.svdir <- paste0(base.save.dir, 'sifs/')
  if(!dir.exists(sif.svdir))
    dir.create(sif.svdir, recursive = T)
  #visaux::ragg.wrapper(
  visaux::ggsave.hirez(fn = rids
                       , plot = sif.vis
                       , dir = sif.svdir)

  # setup for ERGM

  if(run.ergms) {

    # de-spatialize
    ugh <- ghsf %>%
      as_tbl_graph() %>%
      activate('edges') %>%
      select(-geometry) %>%
      activate('nodes') %>%
      select(-geometry)

    # remove units (should be meters)
    ugh <- ugh %>%
      activate('edges') %>%
      mutate(dst = as.numeric(dst))

    eugh <- ergm.clean.n.convert(ugh)

    require(statnet)
    require(ergm)

    sif.int.fo <-
      formula(
        eugh ~ edges +
          nodematch("int.poly", diff = FALSE) +
          absdiff("perc.wh", pow=1) +
          edgecov(decay.mat, "dst.decay")
      )

    sif.art.fo <-
      formula(
        eugh ~ edges +
          nodematch("art.poly", diff = FALSE) +
          absdiff("perc.wh", pow=1) +
          edgecov(decay.mat, "dst.decay")
      )

    sif.int.ergm <-
      ergm(formula = sif.int.fo
           ,"tstr"
           ,reference = ~Poisson)

    sif.art.ergm <-
      ergm(formula = sif.art.fo
           ,"tstr"
           ,reference = ~Poisson)

    ergm.svdir <- paste0(base.save.dir, 'ergms/')
    if(!dir.exists(ergm.svdir))
      dir.create(ergm.svdir, recursive = T)

    saveRDS(sif.int.ergm
            , file = paste0(ergm.svdir
                            ,'interstate-models/'
                            ,rids,'.rds'))
    saveRDS(sif.art.ergm
            , file = paste0(ergm.svdir
                            ,'arterial-models/'
                            ,rids,'.rds'))
  }

  return(1)
}


?divM::gen.cross.tract.dividedness

# test calls
tmp <- Della.wrapper_region2ergm(cz = '03101'
                                 ,run.ergms = F)


# get areas to run --------------------------------------------------------

library(geox)

arrs <- geox::rpops %>%
  add.rns() %>%
  filter(pop > 5e5)

# set params ---------------------------------------------------------

#'
#' cz = NULL
#' ,cbsa = NULL
#' ,min.flows = 1
#' ,tie.str.deciles = 5
#' ,tracts.or.groups = 'bg' #'ct'
#' ,trim.loops = T
#' ,directed = F
#' ,year = 2019
#' ,ddir = '/scratch/gpfs/km31/'
#' ,sfg.dir = 'sfg-processed/orig_dest_annual/'
#' ,base.save.dir =paste0(ddir, 'ergms-newer/')
#' ,run.ergms = T

base.params <- tibble(
  min.flows = 1
  ,tie.str.deciles = 5
  ,tracts.or.groups = 'bg' #'ct'
  ,trim.loops = T
  ,directed = F
  ,year = 2019
  ,ddir = '/scratch/gpfs/km31/'
  ,sfg.dir = 'sfg-processed/orig_dest_annual/'
  ,base.save.dir = paste0(ddir, 'ergms-block-groups/')
  ,run.ergms = T
)

czs2gen <- arrs %>% filter(rt == 'cz')




# gen CZs -----------------------------------------------------------------

cz.params <-
  tibble(cz = czs2gen$rid
         ,base.params)


library(rslurm)


job <- rslurm::slurm_apply(
  Della.wrapper_region2ergm
  ,params = cz.params
  ,jobname = 'cz-bg-ergms'
  ,nodes = 19,
  cpus_per_node = 1,
  slurm_options = list(time = '10:00:00',
                       'mem-per-cpu' = '20G',
                       'mail-type' = list('begin', 'end', 'fail'),
                       'mail-user' = 'km31@princeton.edu')
  #global_objects = c('spws')
  )

