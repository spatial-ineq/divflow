#' purpose of this script is to scratch out where I'm at so far:
#'
#' Show the visualizations; do some models; play with different spatial scales;
#' do the SIF parameter optimization, etc.
#' But without getting too bogged down trying to do everywhere at once.

# setup ws ----------------------------------------------------------------

rm(list=ls())
library(tidyverse)
library(sf)
library(statnet)
library(igraph)
library(tidygraph)

#ddir <- 'scratch/gpfs/km31'
ddir <- Sys.getenv('drop_dir')
ergm.dir <- paste0(ddir, 'graphs-and-ergms/')
cur.dir <- "R/SIF & mapping - mkdown & analysis dev/"

devtools::load_all()

# get hwy data for visuals & transform to node crs
nhpn <- geox::get.NHPN()

# nhpn <- nhpn %>% divM::conic.transform()


# select city -------------------------------------------------------------

city.str <-
  "19700" #philly
  #"24701" # st louis

rids <- geox::get.region.identifiers(cz = city.str)

# bounds
bounds <- geox::build.CZs(rids$rid) #%>% st_transform(4326)

# load prepped graphs -----------------------------------------------------

# at cz & place levels
czdir <- paste0(ergm.dir, "cz-graphs/")
plcdir <- paste0(ergm.dir, "plc-graphs/")

czgh <- list.files(czdir
                   , pattern = city.str, full.names = T) %>% readRDS()
plcgh <- list.files(plcdir
                    , pattern = city.str, full.names = T) %>% readRDS()

# little additions --------------------------------------------------------

# i had formed without trimming by population, just to allow max flexibility
# later
trim.pop <- function(gh, pop.floor = 500){
  gh %>%
    activate("nodes") %>%
    filter(pop > pop.floor)}

# czgh <- czgh %>%  trim.pop
# plcgh <- plcgh %>% trim.pop

czgh <-
  czgh %>% spatialize.graph(bounds)

plcgh <-
  plcgh %>% spatialize.graph(bounds)

# get visuals -------------------------------------------------------------
make.graph.map
czplot <-
  flow.map.wrapper(bounds
                   ,czgh
                   ,tie.str.floor = 5e-3
                   ,flow.colm = "binned.tie.strength"
                   ,trim2plc = F
                   ,n_breaks = 3
                   ,digits = 5)

devtools::load_all()

plcplot <-
  flow.map.wrapper(bounds
                   ,plcgh
                   ,tie.str.floor = 5e-3
                   ,edge.attr = "tie.strength"
                   ,trim2plc = T
                   ,n_breaks = 3
                   ,digits = 5)

czplot
plcplot

# everything for cz or plc  -------------------------------------------------

gh <- czgh
#gh <- plcgh

# some setup --------------------------------------------------------------

# spatial nodes
sfn <- gh %>% get_nodes() %>% st_sf()

# spatial edges
sfe <-  gh %>%
  activate("edges") %>%
  sfnetworks::as_sfnetwork(., directed = T,
                           edges_as_lines = T) %>%
  get_edges() %>%
  st_sf()

# edge distances
sfe$dst <- st_length(sfe$geometry)
sfe$dst <- as.numeric(sfe$dst) / 1000 # to km

# pairwise (tractxtract) distance matrix
dst.mat <- build.pairwise.dst.matrix(sfn)
# max km span
max(dst.mat)

# do spatial interaction function analysis ------------------------------------

# binned relative frequencies (# of ties over distance, relative to crosswise
# distance distribution)
dst.freq <- relative.tie_dst.frequency(
  sfn, sfe, dst.mat = dst.mat,
  n.bins = 60)

dst.freq %>%
  ggplot(aes(x = dst.lvl,
             y = avg.flow)) +
  geom_bar(stat = "identity")

# log-log'd
dst.freq %>%
  ggplot(aes(x = log(dst.lvl),
             y = log(avg.flow))) +
  geom_path()
# this pattern is common for this and shows up in the Daraganova paper and all
# the places I've looked at at any level. It shows the the "power law"
# distribution holds where distance < threshold, and then breaks down.


#' trim to below 50 km, around when power law distribution begins breaking down
ghT <- gh %>%
  activate("edges") %>%
  mutate(dst = sfe$dst) %>%
  filter(dst < 50)

dst.freq <- relative.tie_dst.frequency(
  gh=ghT,
  n.bins = 60)

dst.freq %>%
  ggplot(aes(x = dst.lvl,
             y = avg.flow)) +
  geom_bar(stat = "identity")

# log-log'd, with distant trips trimmed
dst.freq %>%
  ggplot(aes(x = log(dst.lvl),
             y = log(avg.flow))) +
  geom_path()


# estimate parameters for SIF ---------------------------------------------

#' these steps were kinda confusing to me. the `igraph` package has
#' power.law.fit, and another package called `poweRlaw` dedicated to fitting
#' power law distributions, which cites some literature about doing this. Then
#' the Daraganova and Butts and Acton papers do a slightly different thing, with
#' no clear implementation in R. I played around and was kinda confused about
#' this, but I end up (I think) implementing the Daraganova and Butts and Acton
#' strategy, which is to get a general distribution (power law), and then use
#' maximum-likelihood estimation (MLE) to fit the parameters. I use the tutorial
#' here to do the MLE for the power law distribution

# below relies on the `dst.freq` in the general environment
param.fitting <-
  stats::optim(par= c(1,1.7,1, 10),
               fn = power.law.loss.fcn, hessian = T)


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
                x = dst.lvl))

#so there's the general distribution, and we can use those estimates to
#parameterize the SIF in the ERGM
decay.mat <-
  dst.mat %>%
  power.law(
    1, #param.fitting$par[1],
    param.fitting$par[2],
    param.fitting$par[3]
  )



# run an ergm -------------------------------------------------------------
library(statnet)

# to trim ties to save time while hashing out code
sfe$tie.strength %>% quantile(seq(0,1,.05))
ghT <-
  ghT %>%
  activate("edges") %>%
  filter(tie.strength > 1e-3)

# erg model specification
fo <-
  formula(
    ergh ~ edges +
      nodematch("int.poly", diff = FALSE) +
      # nodematch("in.cc") +
      edgecov(decay.mat, "proximity") +
      nodecov("ln.pop.dens", form="sum") +
      absdiff("perc.wh", pow=1) +
      absdiff("perc.hsp", pow=1) +
      absdiff("hh.median.income", pow=1)
  )


# switch to statnet libraries
ergh <<- ergm.clean.n.convert(ghT)

# set seed
set.seed(1234)

# run model
mo <- ergm(
  fo,
  response = "tie.strength",
  reference = ~Poisson
)

mo
summary(mo)

# same ws
save.image(
  paste0(cur.dir,
         "19700",
         "-philadelphia-cz-analysis-env.RData"))


# plotting, facetted by distance ------------------------------------------

# ( i think for most areas, this is only interesting at Place/county level-- CZ
# will be often too zoomed out)

#' Finally, let's facet the plots above by distance bins:
#'
#' 1) Very close; less than a mile (~1.6km)
#'
#' 2) kinda close; less than 5 mi (~8km)
#'
#' 3) below 50 km, around when power law distribution begins breaking down
#'
#' 4) >50km
spatial.brks <-
  c(0, 1.6, 8, 50,1e9)

agh <- ghT %>%
  activate("edges") %>%
  mutate(spatial.scale = cut(dst,
                             spatial.brks,
                             include.lowest = T))

# make graph outside of Wrapper to play with more

# gotta do it outside of wrapper to facet
p <-
  Wrapper_make.graph.map(agh,
                         tie.str.floor = 5e-3,
                         flow.colm = "binned.tie.strength"
                         ,trim2plc = F
                         ,n_breaks = 3
                         ,digits = 5)
p +
  facet_wrap(~spatial.scale)
