# ws ---------------------------------------------------------------------------

rm(list = ls())

require(tidyverse)
require(sf)

# option setting
sf_use_s2(T)
options(tigris_use_cache = TRUE)

# dropbox dir or della base dir
ddir <- Sys.getenv('drop_dir')
  # '/scratch/gpfs/km31/'

devtools::load_all()


# sample areas -----------------------------------------------------------------

library(geox)
arrs <- geox::rpops %>%
  filter(pop > 50e3
         ,rt == 'cz') %>%
  add.rns()

arrs %>%
  filter(pop > 500e3) %>%
  arrange(pop) %>%
  head(100) %>%
  pull(rn)

arrs %>%
  filter(grepl('Portland', rn))

czsf <- build.CZs('20100')
plcsf <- geox::places.wrapper(x = czsf)

library(mapview)
czsf %>% mapview() + mapview(plcsf, color = 'red')

# for analysis, i want probably cz/cbsa (or maybe sometimes Place)
# for visuals, i want region network cropped to bbox
# for sample, i'll use Place
smplc <- plcsf %>%
  filter(grepl('^Portland', name)) %>%
  st_transform(4326)
smplc %>% mapview()


# get block groups for place ---------------------------------------------------

bgs <- geox::tracts.from.sf(x = st_bbox(smplc)
                            ,query.fcn = tigris::block_groups)
bgs[3] %>% plot()

# remove water areas (per census naming conventions)
mapview(bgs, zcol = 'awater')
sbgs['awater'] %>% plot()
sbgs <- bgs %>% filter(substr(tractce,1,2) != '99')


# subset to polygon, not bboxx
sbgs <- xwalks::get.spatial.overlap(sbgs, smplc
                                    ,'geoid', 'placefp')

# filter original census query to keep all columns
sbgs <- bgs %>%
  filter(geoid %in% sbgs$geoid)

sbgs %>% mapview()

# just b/c for visuals, remove islands in ad hoc way
sbgs <- sbgs %>% filter(tractce != '002400')

bounds <- st_sf(geometry = st_union(sbgs))

sbgs <- sbgs %>% st_transform(4326)
bounds <- bounds %>% st_transform(4326)

# load safegraph for area ------------------------------------------------------

devtools::load_all()

# flexible function
sfg <- sfx2sfg(bounds
               ,min.flows = 10
               ,tracts.or.groups = 'bg'
               ,trim.loops = T
               ,year=2019)

#trim to just place
sfgt <- sfg %>%
  filter(across(c(origin, dest)
                      , ~.x %in% sbgs$geoid))

# turn to (undirected) graph
ugh <- sfg2gh(sfgt, directed = F)

ugh

ghsf <- spatialize.graph(ugh,
                         frame.sf = NULL
                         ,tracts.or.groups = 'bg'
                         ,directed = F
                         ,year = 2019)
# transform
ghsf <- st_transform(ghsf, 4326)

ghsf
sbgs

# flow summaries -------------------------------------------------------------

edges <- ghsf %>%
  activate(edges) %>%
  as_tibble()

tibble(edges)[Hmisc::Cs(n, tstr, dst)] %>% map(summary)

tibble(edges)[Hmisc::Cs(n, tstr, dst)] %>% map( ~quantile(.x, seq(0,1,.1)))


# trim by flow str -------------------------------------------------------------

devtools::load_all()

bounds
sbgs
# skipped; sample area seems small enough
trimmed.ghsf <- apply.flow.filters(ghsf
                                   ,frame.sf = NULL
                                   ,directed = F
                                   ,min.tie.str = NULL
                                   ,tie.str.deciles = 4 # filter out bottom x deciles
                                   ,min.flows = 10
                                   ,max.dst = units::set_units(20, 'miles'))

trimmed.ghsf
ghsf


# map graph ------------------------------------------------------------

sbgs <- sbgs %>% st_transform(4326)

# get background for mapping
sttm <- visaux::get.stamen.bkg(
  sbgs
  , zoom = 12
  ,maptype = 'toner-background'
  )

# lonlat matrix for notes
lonlats <- ghsf %>%
  activate('nodes') %>%
  as_tibble() %>%
  st_sf() %>%
  st_coordinates()

library(ggraph)

ggmap(sttm
      , base_layer =
        ggraph(trimmed.ghsf # ghsf #
               , layout = lonlats )
        ) +
  geom_edge_density(fill = '#00D0D0') +
  geom_edge_fan( aes( edge_alpha = tstr
                    ,edge_width = tstr)
                 ,color = '#008080'
                ) +
  flow.map.base() +
  visaux::bbox2ggcrop(bounds)


# merge in n'hood level ERGM controls -------------------------------


# to (re)generate nhood divisions ----------------------------------------------



hwys <- geox::get.NHPN(sfx = bounds #%>% st_buffer(1e4)
                 ,dropbox.dir = ddir)
fwhys <- divM::Fix.all.hwys(hwys
                            ,threshold = 500
                            ,return.gap.map = T)

xnbd <- divM::gen.cross.tract.dividedness(
  smplc
  ,hwys
  ,nbd.query.fcn = tigris::block_groups
  ,region.id.colm = 'placefp'
  ,fill.nhpn.gaps = F # this can help in some cases and not others... nhpn data is not perfect
  ,erase.water = F
  ,year = 2019
)


xnbd <- xnbd %>%
  select(geoid, hwy.poly = poly.id)

# merge in w/ nodes
ghsf <- ghsf %>%
  activate('nodes') %>%
  left_join(xnbd
            , by = c('name' = 'geoid'))

# vis
ghsf %>%
  activate('nodes') %>%
  as_tibble() %>% st_sf() %>% mapview(zcol = 'hwy.poly') +
  mapview(hwys)



# alternate with urban arterials/ census road data --------------------------------

#' ( this would mimic Grannis 2005 definition, using census roads and filtering by
#' CFCC (census feature classification codes)). The census codes changed since his
#' paper and new ones referenced here:
#'
#' https://www2.census.gov/geo/pdfs/reference/mtfccs2019.pdf
#' (See MTFCCS: S1100, S1200, S1400)


sbgs
rds <- tigris::roads(state = unique(sbgs$statefp)
              ,county = unique(sbgs$countyfp)) %>%
  rename_with(tolower)

# filter to primary & secondary roads (limited-access and arterials) and hwy ramps
arterials <- rds %>%
  filter(mtfcc %in%
           Hmisc::Cs(S1100, S1200, S1630)) %>%
  st_transform(4326)

arterials <- st_crop(arterials, smplc)

xnbd2 <- divM::gen.cross.tract.dividedness(
  smplc
  ,arterials
  ,nbd.query.fcn = tigris::block_groups
  ,region.id.colm = 'placefp'
  ,fill.nhpn.gaps = F # this can help in some cases and not others... nhpn data is not perfect
  ,erase.water = F
  ,year = 2019
)


xnbd2 <- xnbd2 %>%
  select(geoid, arterial.poly = poly.id)

# merge in w/ nodes
ghsf <- ghsf %>%
  activate('nodes') %>%
  left_join(xnbd2
            , by = c('name' = 'geoid'))

# vis
ghsf %>%
  activate('nodes') %>%
  as_tibble() %>% st_sf() %>% mapview(zcol = 'arterial.poly') +
  mapview(arterials)

# can use this data to treat roads in more differentiated way too!! And census data
# is better than NHPN, if functions accommodate the complexity well enough.


# get other tract-level characteristics ----------------------------------------

demos <- sfg.seg::demos.2019 %>%
  select(1,2,matches('wh$|bl$|hsp$'))

ghsf <-
  ghsf %>%
  activate('nodes') %>%
  left_join(demos
            , by = c('name' = 'geoid'))

#region-level demos
demos %>%
  select(matches('pop|^n')) %>%
  colSums() %>% format(big.mark = ',')

# play with democgraphic flow visualization ------------------------------------

# could be better if not undirected
dghsf <- ghsf %>%
  activate('edges') %>%
  mutate(flow.n.bl =
           n * (.N()$n.bl[from] + .N()$n.bl[to]) / 2
         ,flow.n.hsp =
           n * (.N()$n.hsp[from] + .N()$n.hsp[to]) / 2
         ,flow.n.wh =
           n * (.N()$n.wh[from] + .N()$n.wh[to]) / 2
         )
# .....



# final setup for ergm ---------------------------------------------------------

# ergm can't handle a lot of things: no factors, no geometry, no bundled units, no
# NAs, etc....

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

# convert using statnet packages
# hard to manipulate, but built for the ergm implementation
devtools::load_all()
library(statnet)
eugh <- ergm.clean.n.convert(ugh)

# run ergm ---------------------------------------------------------------------

#install.packages("ergm.count")
#install.packages("ergm.rank")
#install.packages("latentnet")
#update.packages()
library(ergm)
library(ergm.rank)
library(latentnet)

eugh
fo.base <- formula(
  eugh ~ edges +
    nodematch("arterial.poly", diff = FALSE) +
    #nodematch("hwy.poly", diff = FALSE) +
    absdiff("perc.wh", pow=1)
  )
ugh

# ergms take while to run. A good time to use Della
ergm1 <- ergm(formula = fo.base
     ,"n"
     ,reference = ~Poisson)

ergm1 %>% summary()


# save ws ----------------------------------------------------------------------

save.dir <- 'R/ERGM versions/'
save.dir %>% list.files()

#save.image(paste0(save.dir, 'sample ergm ws.rdata'))


# generate spatial interaction fcn ---------------------------------------------

# may little SIF tutorial:
# https://rpubs.com/kmc39/sif

sfn <- ghsf %>%
  activate('nodes') %>%
  as_tibble()

sfe <- ghsf %>%
  activate('edges') %>%
  as_tibble()

# turn distance into numerics representing km
sfe <- sfe %>% mutate(dst = as.numeric(dst) / 1000 )

dst.mat <- build.pairwise.dst.matrix(sfn)
dst.mat %>% head() # values are kilometers

# binned relative frequencies (# of ties over distance, relative to crosswise
# distance distribution)
devtools::load_all()

dst.freq <-
  relative.tie_dst.frequency(
      sfn, sfe
    , dst.mat = dst.mat
    , n.bins = 60
    )

dst.freq

dst.freq %>%
  ggplot(aes(x = dst.lvl,
             y = avg.flow))  +
  geom_path()

# log-log'd
dst.freq %>%
  ggplot(aes(x = log(dst.lvl),
             y = log(avg.flow))) +
  geom_path()

# SIF typically shows decreasing flows with distance, until a certain point where the
# relationship is lost. About 2km here.

# use flow bins to parameterize a power law SIF
# below relies on the `dst.freq` in the general environment
param.fitting <-
  stats::optim(par= c(1,1.7,1, 10),
               fn = power.law.loss.fcn, hessian = T)

# predicted tie probabilities, based on SIF
dst.freq$flow.hat <-
  power.law(dst.freq$dst.lvl,
                    param.fitting$par[1],
                    param.fitting$par[2],
                    param.fitting$par[3]
  )

dst.freq %>%
  select(dst.lvl, avg.flow, flow.hat) %>%
  pivot_longer(cols = contains("flow"),
               values_to = "avg.flow") %>%
  ggplot() +
  geom_path(aes(color = name,
                y = avg.flow,
                x = dst.lvl)) +
  ggtitle("Relative avg flows",
          "actual vs. fitted power law")

decay.mat <- dst.mat %>%
  power.law(param.fitting$par[1],
            param.fitting$par[2],
            param.fitting$par[3])


# ERGM with sif ----------------------------------------------------------------

fo.sif <-
  formula(
  eugh ~ edges +
    nodematch("arterial.poly", diff = FALSE) +
    #nodematch("hwy.poly", diff = FALSE) +
    absdiff("perc.wh", pow=1) +
    edgecov(decay.mat, "dst")
  )
eugh

# ergms take while to run. A good time to use Della
ergm2 <- ergm(formula = fo.sif
              ,"tstr"
              ,reference = ~Poisson)

ergm2 %>% summary()


# using wrapper fcns to replicate similar --------------------------------------
devtools::document()
devtools::load_all()

?setup.gh.wrapper
bounds <- bounds %>% st_transform(4326)
bounds %>% plot()
bxgh <- setup.gh.wrapper(sfx = bounds)
bxgh


# check older saved models -----------------------------------------------------

erdir <- paste0(ddir,'graphs-and-ergms/')

library(tidygraph)


smpgl <- erdir %>%
  list.files(pattern='plc-.*rds'
             ,full.names = T
             ,recursive = T) %>%
  head(1) %>%
  read_rds()

smpgl
