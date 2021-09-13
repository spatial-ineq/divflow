# ws ---------------------------------------------------------------------------
rm(list = ls())
require(tidyverse)
require(sf)
#devtools::load_all()

# option setting
sf_use_s2(T)
options(tigris_use_cache = TRUE)
Sys.setenv("VROOM_SHOW_PROGRESS"="false")


# set directories ---------------------------------------------------------

#prx.dir <-  '/scratch/gpfs/km31/adjacencies+proximities/'

sfg.dir <-
  '/projects/SHARKEY/safegraph/processed/orig_dest_annual/'


# block groups ------------------------------------------------------------

fns <- sfg.dir %>% list.files(recursive = T, pattern = '2019.*csv$')

## params for job ##

# save dir
save.dir <-
  paste0(prx.dir
         ,'flow-weights/'
         ,'visited-bg-flws/')
# '~/R/all sharkey geoseg work/divflow/R/creating SpWs/spWs/'
save.dir

czs <- geox::rx$cz %>% unique()



## simple della fcn
Della.wrapper_incoming.flww.by.cz <- function(cz
                                              ,agg2tracts = F
                                              ,drop.loops = T
                                              ,weight.floor = .01
                                              ,sfg.dir
                                              ,year = '2019'
                                              ,save.dir = NULL) {

  require(tidyverse)
  sfg <- sfg.seg::read.sfg.CZs(cz,
                               sfg.dir = sfg.dir,
                               year = year)
  if(agg2tracts)
    sfg <- sfg.seg::cbg.flows2tracts(sfg)

  if(drop.loops)
    sfg <- sfg %>% filter(origin != dest)

  flww <- sfg %>%
    group_by(origin) %>%
    mutate(visited.flww =
             n / sum(n) ) %>%
    ungroup() %>%
    filter(visited.flww >= weight.floor)


  if(!is.null(save.dir)) {

    if(!exists(save.dir))
      dir.create(save.dir)

    write.csv(flww
              ,file = paste0(save.dir
                             ,as.character(cz), '-', as.character(year)
                             ,'-flow-weights.csv')
              ,row.names = F)
  }
  return(flww)
}



bg.params <-
  tibble(
    cz = czs
    ,agg2tracts = T
    ,weight.floor = 0.001
    ,drop.loops = T
    ,sfg.dir = sfg.dir
    ,year = '2019'
    ,save.dir = save.dir
  )


# send job
library(rslurm)
bg.flwws.dellajob <-
  slurm_apply(f =
                Della.wrapper_incoming.flww.by.cz,
              params = bg.params,
              jobname = 'cbg visited flow weights no dist filter',
              nodes = 19,
              cpus_per_node = 1,
              slurm_options = list(time = '10:00:00',
                                   'mem-per-cpu' = '20G',
                                   'mail-type' = list('begin', 'end', 'fail'),
                                   'mail-user' = 'km31@princeton.edu')
  )

# scratch -----------------------------------------------------------------

# sample cz
smpl <- sfg.seg::read.sfg.CZs('19700'
                              ,sfg.dir)

smplflw <- smpl %>%
  filter(origin != dest) %>%
  group_by(origin) %>%
  mutate(visited.flww =
           n/sum(n)) %>%
  ungroup() %>%
  filter(visited.flww >= .01)

smplflw
