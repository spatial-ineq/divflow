# ws ---------------------------------------------------------------------------

rm(list = ls())

require(tidyverse)
require(sf)

# option setting
sf_use_s2(T)
options(tigris_use_cache = TRUE)

# dropbox dir or della base dir
ddir <- # Sys.getenv('drop_dir')
  '/scratch/gpfs/km31/'

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
plcsf <- visaux::places.wrapper(x = czsf) %>% distinct() # (bug  in the helper; i shouldn't need "distinct")

plcsf <- plcsf %>% rename_with(tolower)

library(mapview)
czsf %>% mapview()
plcsf %>% mapview()

# for analysis, i want probably cz/cbsa (or maybe sometimes Place)
# for visuals, i want region network cropped to bbox
# for sample, i'll just use Place
smplc <- plcsf %>%
  filter(grepl('^Portland', name))
smplc %>% mapview()

# load safegraph for area ------------------------------------------------------

# flexible function
sfg <- sfx2sfg(smplc
               ,min.flows = 10
               ,tracts.or.groups = 'bg'
               ,trim.loops = T
               ,year=2019)

sfg

# trim to place largely for speed reasons
bgs <- sfg %>%
  select(geoid = origin) %>%
  distinct() %>%
  geox::attach.geos(query_fcn = tigris::block_groups
                    ,year=2019)


# subset block groups to place
smplc <- smplc %>% st_transform(4326)

sbgs <- xwalks::get.spatial.overlap(bgs, smplc
                                    ,'geoid', 'placefp')

# visualize blocks in sample area
sbgs %>% mapview(zcol = 'perc.area') +
  mapview(st_boundary(smplc)
          ,color = 'red')

# remove islands
sbgs <- sbgs[lengths(st_intersects(sbgs)) > 1, ]
sbgs %>% mapview()

# trimmed sfg
sfgt <- sfg %>%
  filter(origin %in% sbgs$geoid &
           dest %in% sbgs$geoid)

# ugh: Undirected GrapH
ugh <- sfg2gh(sfgt, directed = F)

ugh

#sf_use_s2(F)
ghsf <- spatialize.graph(ugh,
                         frame.sf = NULL
                         ,tracts.or.groups = 'bg'
                         ,directed = F
                         ,year = 2019)

ghsf
sbgs

# flow summaries -------------------------------------------------------------

edges <- ghsf %>%
  activate(edges) %>%
  as_tibble()

tibble(edges)[Hmisc::Cs(n, tstr, dst)] %>% map(summary)

tibble(edges)[Hmisc::Cs(n, tstr, dst)] %>% map( ~quantile(.x, seq(0,1,.1)))


# trim by flow str -------------------------------------------------------------

# skipped; sample area seems small enough
"trimmed.ghsf <- apply.flow.filters(ghsf
                                   ,frame.sf = mn.plc
                                   ,directed = F
                                   ,min.tie.str = NULL
                                   ,tie.str.deciles = 5
                                   ,min.flows = 10
                                   ,max.dst = units::set_units(4, 'miles'))

trimmed.ghsf
"


# setup for mapping ------------------------------------------------------------

sbgs <- sbgs %>% st_transform(4326)

# get background for mapping
sttm <- visaux::get.stamen.bkg(
  sbgs
  , zoom = 12
  ,maptype = 'toner-background'
  )

sbgs <- bgs %>%
  filter(geoid %in% sbgs$geoid)
sbgs

sttm %>%
  ggmap() +
  geom_sf(data = sbgs %>% st_cast('POLYGON') %>% st_boundary()
          ,color = '#009964'
          ,inherit.aes = F) +
  theme_void()


# merge in n'hood level ERGM controls -------------------------------

# dropbox dir or della base dir
ddir <- Sys.getenv('drop_dir')



# to (re)generate nhood divisions ----------------------------------------------
sbgs %>% st_bbox() %>% st_buffer(10) %>% mapview()

boundary <- st_union(sbgs)

nhpn <- visaux::get.NHPN(sfx = boundary #%>% st_buffer(1e4)
                 ,dropbox.dir = ddir)
nhpn <- nhpn %>%
  filter(! st_is_empty(geometry))

# vis:
sttm %>%
  ggmap() +
  geom_sf( data = st_boundary(boundary)
           ,color = 'purple'
           ,linetype = '21'
           ,size = 1.2
           ,inherit.aes = F) +
  geom_sf( data = nhpn
          ,aes(color = signt1)
          ,size = 1.3
          ,inherit.aes = F)

#hwy.divs <-
  boundary %>%
    st_split(nhpn) %>%
    st_cast('POLYGON')

  divM::polygonal.div(boundary, nhpn
                               ,return.sf = T)




