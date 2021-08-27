
# setup ws ---------------------------------------------------------------------

proj.dir <- "~/R/all sharkey geoseg work/divseg/"

source(paste0(proj.dir,
              "R/tidy models/setup to map flows.R"))


# fcns & templates to map ------------------------------------------------------

#' subset.and.spatialize
#'
#' Subsets to an area to map. Downloads tract geos and attaches them.
#'
#' @param ... passed to `sfg.seg::geo.subset.cbgs` to subset to area (cz_name, cz,
#'   etc.)
#' @param flr a flow/residential dataset wide by flow.type, as generated in setup
#'   script
subset.and.spatialize <- function(flr, cutout.water = F, ...) {

  x <- flr %>%
    sfg.seg::geo.subset.cbgs(subset.cols = "geoid",
                             ...)

  # all 5-char state-counties
  counties <- unique(x$countyfp)

  # download tract
  ctsfs <- map_dfr(counties,
                     ~tigris::tracts(state =
                                       substr(., 1, 2),
                                     county =
                                       substr(., 3, 5))
                   )
  # merge in
  x <- left_join(x, ctsfs["GEOID"],
                 by = c("geoid" = "GEOID"))

  x <- x %>% st_sf() %>% st_transform(4326)

  # cut out water
  if(cutout.water)
    x <- download.and.cutout.water(x)

  return(x)
}

# for spatial area x
download.and.cutout.water <- function(x, size.min = 5e6) {

  # all 5-char state-counties
  counties <- unique(x$countyfp)

  # download water
  water <- map_dfr(counties,
                   ~tigris::area_water(state =
                                         substr(., 1, 2),
                                       county =
                                         substr(., 3, 5))
                   )
  # union and explode water
  water <- st_union(water) %>% st_cast("POLYGON") %>% st_sf()
  # filter by size of union'd body
  water <- water %>% filter(as.numeric(st_area(.$geometry)) > size.min )
  water <- water %>% st_transform(st_crs(x))
  #water <- water %>% rmapshaper::ms_simplify()

  xt <- st_difference(x, st_union(water))

  return(xt)
}

phl <- subset.and.spatialize(fl, cz_name = "Philadelphia")
# map.fcns ---------------------------------------------------------------------


# residential/inflow/outflow tiled map

#' make.res.flow.map
#'
#' @param flr a flow/residential dataset wide by flow.type, as generated in setup
#'   script and subsetted w/ `subset.and.spatialize`
#' @param div.sf NULL if no divs included, otherwise an sf object containing
#'   divisions to overlap and a `divtype` column
make.res.flow.map <- function(flr, div.sf = NULL) {

  area.perc.bl <-
    stats::weighted.mean(flr$bl, flr$pop
                         , na.rm= T)

  p <-
    flr %>%
    ggplot() +
    geom_sf(aes(fill =
                  bl)
            ,color = NA) +
    facet_wrap(~ flow.type, ncol = 2) +
    scale_fill_gradientn(
      colours =
        colorspace::diverge_hcl(3 , palette = "Purple-Green"),
      values = c(0,
                 area.perc.bl,
                 1)
    )

  if(!is.null(div.sf))
    p <- p + add.divs(flr, div.sf)

  return(p)
}


# dif from res map tiled inflow/outflow

#' dif.from.res.map
#'
#'
dif.from.res.map <- function(flr, div.sf = NULL) {

  range <- flr$dif.from.res_bl %>% range(na.rm = T)

  bound <- max(abs(r))

  p <-
    flr %>%
    filter(flow.type != "res") %>%
    ggplot() +
    geom_sf(aes(fill =
                  dif.from.res_bl)
            ,color = NA) +
  facet_wrap(~ flow.type,
             ncol = 1) +
    scale_fill_gradientn(
      colours =
        colorspace::diverge_hcl(3 , palette = "Blue-Red")
      ,limits =   c(-bound, bound),
      #values = c(-1*bound, 0, bound)
    ) +
    labs(captions = "dif.from.res = flow - residential
         so when dif.from.res_bl is negative, residents of tract are more likely to be black than visitors
         when it's positive, visitors or more likely to be black than residents.")

  if(!is.null(div.sf))
    p <- p + add.divs(flr, div.sf)

  return(p)
}



# adding divisions fcns --------------------------------------------------------

# get data
hwys <-
  st_read("~/R/shapefiles/National_Highway_Planning_Network-shp/National_Highway_Planning_Network.shp")

hwys <- st_transform(hwys, 4326)


# fcn to trim to mapped area and
add.divs <- function(area, divs, trim = T) {

  if (trim) {
    divs <- st_transform(divs, st_crs(area))
    divs <- st_intersection(divs, st_union(area))
  }

  list(
#    geom_sf(data = divs,
#          color = "#FFFFFF"
#            ,size = .9) ,
      geom_sf(data = divs,
              color = "#2d6b6b"
              ,size = .5)
  )
}

make.res.flow.map(chi, div.sf = ints)

# get other candidate czs ------------------------------------------------------

# where area-wide flow-res dissim metrics are greatest
difs <- mx$cz %>%
  select(-contains("truly")) %>%
  sfg.seg::compare.measures_mvmtVres(measure.ids = c("dissim_wh.bl" , "isol_wh", "isol_bl"))

difs <-
  difs %>%
  select(1:3, pop, contains("perc."), contains("dif")) %>%
  arrange(desc(.$dif_inflow.vs.res_dissim_wh.bl))

difs %>%
  filter(pop > 5e5) %>%
  filter(perc.bl > .1) %>%
  arrange(desc(.$dif_inflow.vs.res_dissim_wh.bl))

# after the population filters above, there are 67 CZs
# (at least .5million residents and 10% Black)

# Atlanta is number10 ~least~ difference res-inflow segregation
# NYC is #1 ~most~ difference between res-inflow seg..

# making maps ------------------------------------------------------------------

co2cz <- xwalks::co2cz


plc.co.cz <- xwalks::plc.co.cz
# czs to map
cz_list <-
  c("Atlanta", "Philadelphia", "St. Louis",
    "Chicago", "New York", "Houston",
    "Charlotte", "Minneapolis", "Detroit"
    )

# also to map all Places w/ >100k population in those CZs (zoom ins!)
plc_list <-
  plc.co.cz %>%
  filter(cz_name %in% cz_list) %>%
  #count(plc2cz.perc.overlap)
  mutate(statefp=substr(countyfp, 1,2)) %>%
  select(1:3, statefp) %>%
  distinct()

plc_list
plcs <-
  map_dfr(plc_list$statefp,
          tigris::places)
plcs <-
  plcs %>%
  filter(NAME %in% cz_list) %>%
  filter(GEOID %in% plc_list$plc) %>%  # to get all placers over 100k pop in CZs
  select(GEOID, NAME, LSAD, CLASSFP) %>%
  tibble() %>% distinct() %>% st_sf()
#
plcs

plc_list <- plcs$GEOID
names(plc_list) <- plcs$NAME


# interpetation notes ----------------------------------------------------------

# dif.from.res = flow - residential
# so when dif.from.res_bl is negative, more % black folks live in a tract than visit.
# When dif.from.res_whi is negative, more % wh folks live in a tract than visit

# Note that in the subtraction above, to get dif.from.res, places like parks,
# airports, golf courses, with 0 residential pop and many visitors are lost. Will
# have to consider how to handle this for visuals.

# sample run -------------------------------------------------------------------
sfg.seg::geo.subset.cbgs

chi <- subset.and.spatialize(fl,
                             cutout.water = F,
                             plc_id = plc_list[["Chicago"]])

make.res.flow.map(chi)
dif.from.res.map(chi) + ggtitle("chicago plc")

# and with divs
ints <- hwys %>% filter(SIGNT1 == "I")
iu

make.res.flow.map(chi, div.sf = ints)
make.res.flow.map(chi, div.sf = hwys)

dif.from.res.map(chi, div.sf = ints)




# iterate through cz/plcs ------------------------------------------------------

# save to output dir
proj.dir <- "~/R/all sharkey geoseg work/divseg/"
save.dir <- paste0(proj.dir,
                   "R/mapping flows/flow maps/")


save.fcn <- function(save.title, p = last_plot()) {
  ggsave(filename = paste0(save.dir, save.title,".png"),
         width = 7, heigh = 8)
}

save.cz.maps <- function(cz, no.div=F) {

  if(!exists("save.dir"))
    stop("specify `save.dir`")

  area <- subset.and.spatialize(fl,
                                cutout.water = T,
                                cz_name = cz
                                )
  # base fn for saved plots
  fn <- paste0(cz, "-cz")

  if(no.div) {
    make.res.flow.map(area) + ggtitle(paste0(cz, " - CZ"))
    save.fcn(paste0(fn,"_flow-maps"))
    dif.from.res.map(area) + ggtitle(paste0(cz, " - CZ"))
    save.fcn(paste0(fn, "_res-flow"))
  }

  # w divs
  make.res.flow.map(area,
                    div.sf = ints) +
    ggtitle(paste0(cz, " - CZ"),
            "with interstates")
  save.fcn(paste0(fn, "_flow-maps_w-divs"))

  dif.from.res.map(area,
                   div.sf = ints) +
    ggtitle(paste0(cz, " - CZ"),
            "with interstates")
  save.fcn(paste0(fn, "_res-flow_w-divs"))

  cat(cz, "done")
}


save.plc.maps <- function(plc, plc_name, no.div=F) {

  if(!exists("save.dir"))
    stop("specify `save.dir`")

  area <- subset.and.spatialize(fl,
                                cutout.water = T,
                                plc_id = plc
                                )

  # base fn for saved plots
  fn <- paste0(plc_name, "-plc")

  if(no.div) {
    make.res.flow.map(area) + ggtitle(paste0(plc_name, " - plc"))
    save.fcn(paste0(fn,"_flow-maps"))
    dif.from.res.map(area) + ggtitle(paste0(plc_name, " - plc"))
    save.fcn(paste0(fn, "_res-flow"))
  }
  # w divs

  make.res.flow.map(area,
                    div.sf = ints) +
    ggtitle(paste0(plc_name, " - plc"),
            "with interstates")

  save.fcn(paste0(fn, "_flow-maps_w-divs"))

  dif.from.res.map(area,
                   div.sf = ints) +
    ggtitle(paste0(plc_name, " - plc"),
            "with interstates")

  save.fcn(paste0(fn, "_res-flow_w-divs"))

  cat(plc_name, "done")
}



# loop and save ----------------------------------------------------------------

ints <- hwys %>% filter(SIGNT1 == "I")

map(cz_list,
      save.cz.maps)

map2(plc_list, names(plc_list),
    save.plc.maps)



map("Minneapolis",
    save.cz.maps)

map("Detroit",
    save.cz.maps)

remaining_places <-
  plc_list[c("Detroit", "Minneapolis")]
map2(remaining_places, names(remaining_places),
     save.plc.maps)



# improving maps scratch -----------------------------------------------------------

# for dif from res, I should make 0 the midpoint and interpolate at equal rates on
# either side -- to make maps represent systematic difs in bl/wh visits

r <- chi$dif.from.res_bl %>% range(na.rm = T)

bound <- max(abs(r))
c(-1*bound, 0, bound)

chi %>%
  filter(flow.type != "res") %>%
  ggplot() +
  geom_sf(aes(fill =
                dif.from.res_bl)
          ,color = NA) +
  facet_wrap(~ flow.type) +
  scale_fill_gradientn(
    colours =
      colorspace::diverge_hcl(3 , palette = "Blue-Red")
    #,values =
    #  c(-bound, 0, bound)
    ,limits =   c(-bound, bound),
    n.breaks = 7
    )


phl["pop"] %>% plot()

hst <- subset.and.spatialize(fl, cutout.water = F,
                             plc = plc_list[["Houston"]])

r <- hst$dif.from.res_bl %>% range(na.rm = T)

bound <- max(abs(r))
c(-1*bound, 0, bound)

hst %>%
  tibble() %>%
  group_by(flow.type) %>%
  summarise(
    stats::weighted.mean(dif.from.res_bl, pop))

hst %>%
  filter(flow.type != "res") %>%
  ggplot() +
  geom_sf(aes(fill =
                dif.from.res_bl)
          ,color = NA) +
  facet_wrap(~ flow.type) +
  scale_fill_gradientn(
    colours =
      colorspace::diverge_hcl(3 , palette = "Blue-Red")
    #,values =
    #  c(-bound, 0, bound)
    ,limits =   c(-bound, bound),
    n.breaks = 7
  )



# checking hwys ----------------------------------------------------------------

# do i have do subset by signt1 or signt2
char <- subset.and.spatialize(fl, cutout.water = F,
                              cz_name = "Charlotte")
char
# charlotte hwys
tmp <- hwys %>%
  st_intersection(st_union(char))
# charlotte interstates
tmp %>% filter(SIGNT1 %in% "I") %>% nrow()
tmp %>%
  filter(if_any(contains("SIGNT"),
                ~{.x %in% c("I")})
         ) %>% nrow()
plot(tmp["SIGNT1"])

