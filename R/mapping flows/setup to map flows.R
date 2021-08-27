# boilerplate ----------------------------------------------------------------------
library(tidyverse)
library(sf)
library(mapview)
devtools::load_all()

rm(list = ls())


# run "tidy models" script -----------------------------------------------------

# gets cz/cbsa data and models, which can be used to identify residuals+outliers
proj.dir <- "~/R/all sharkey geoseg work/divseg/"

source(paste0(proj.dir,
              "R/tidy models/tidy-models.R"))


# workspace review -------------------------------------------------------------

# datasets
mx

# model table
mvs

# remove extra items that were used to create `mvs`
rm(ivs, rts, res.equivs)


# get tract-level data ---------------------------------------------------------

ct <- divseg::get.prepped("ct")
divseg::setup.to.model_all

ct


# get res/flow table -----------------------------------------------------------

flr <-
  ct %>% select(1:6,
                matches("pop$|total."), aland,
                matches("(perc.*)(wh|bl|hsp)"))


flr <- flr %>%
  select(-matches(".*incl$"))


# fix colnames - add res marker to residential values
flr <-
  flr %>%
  rename_with(
    .fn =
      ~{.x %>% gsub("\\.", "_", .) %>% paste0("res.", .)
      },
    # renames where contains "perc" but not "_"
    .cols = matches("(?=.*perc)(?!.*_)(.*+)", perl =T)
  )

colnames(flr)

fl <-
  flr %>%
  #  select(-contains("res.")) %>%
  pivot_longer(cols =
                 #matches("(incoming.*|visited.*)(wh|bl|hsp)"),
                 matches("(res.*|incoming.*|visited.*)(wh|bl|hsp)"),
               names_to = c("flow.type", "demo"),
               names_pattern = "(.*)_(.*)",
               values_to = "perc")


fl <- fl %>%
  filter(demo != "hsp")

fl$flow.type <- gsub("perc|\\.", "", fl$flow.type)

fl <- fl %>%
  pivot_wider(
    names_from = demo,
    values_from = perc
  )

#res <- fl %>% filter(flow.type == "res") %>% select(geoid, wh, bl)
#res <- res %>% rename_with( ~paste0("res_", .), .cols = c("wh", "bl"))

fl

# add flow dif from res
fl <-
  fl %>%
  arrange(geoid, flow.type) %>%
  group_by(geoid) %>%
  mutate(across(c("wh", "bl"),
                list(dif.from.res =
                       ~{.x -.x[2L] }) # 2L indexes to residential row based on arrange call above.
                , .names = "{.fn}_{.col}"
  )
  ) %>%
  ungroup()

fl

# dif.from.res = flow - residential
# so when dif.from.res_bl is negative, more % black folks live in a tract than visit.
# When dif.from.res_whi is negative, more % wh folks live in a tract than visit

# Note that in the subtraction above, to get dif.from.res, places like parks,
# airports, golf courses, with 0 residential pop and many visitors are lost. Will
# have to consider how to handle this for visuals.

# port in divisions ------------------------------------------------------------

