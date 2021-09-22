# ws ---------------------------------------------------------------------------
rm(list = ls())
require(tidyverse)
require(sf)

# option setting
sf_use_s2(T)
options(tigris_use_cache = TRUE)




# get demovars -----------------------------------------------------------------

ddir <- paste0(Sys.getenv('drop_dir'), 'seg-measures/by tract/broader ineq flows/')
flr <- ddir %>%
  list.files(pattern='csv$'
             ,full.names = T) %>%
  vroom::vroom() %>%
  select(-`...1`)

flr
res <- flr %>%
  filter(region.type == 'cz') %>%
  select(1:3,pop,hh,
         matches('^residen'))

res %>% filter(duplicated(geoid)) #(alaska tracts w county splits)

# addl from patrick ------------------------------------------------------------

addlres <- '.local-data/geoid10_segmovement.csv'
addlres <-  vroom::vroom(addlres)

addlres <- select(addlres, geoid10, everything())
addlres$geoid10 <- geox::fix.geoid(addlres$geoid10)

# merge ------------------------------------------------------------------------

# add 2010-19 xtract
xwalks::cts2cts_time.series
res <- res %>%
  left_join( xwalks::cts2cts_time.series[1:2]
             ,by = c('geoid'='geoid15')) %>%
  left_join(addlres)

# other org
res <- res %>%
  rename_with(~gsub('^residential\\.', '', .))

res <- res %>%
  geox::abv.rcols()

colnames(res)
res

# descriptives/other checks -----------------------------------------------------

# most range from abt 0-1
# res %>% map(summary)

res %>% filter(duplicated(geoid))

# make long and organize weights ----------------------------------------------

# make res long
resl <- res %>%
  select(-geoid10) %>%
  pivot_longer(6:ncol(.))


# get a weighted column depending on the var
hh.weighted.cols <-
  colnames(res) %>% grep('fam', .,  value = T)
pop.weighted.cols <-
  colnames(res)[6:ncol(res)] %>% grep('^((?!fam).)*$'
                                      , .,  value = T
                                      ,perl = T)

resl <- resl %>%
  mutate(weight = case_when( name %in% hh.weighted.cols ~ hh
                             ,TRUE ~ pop))


# rename
resl <- resl %>%
  rename(var = name)


# write ------------------------------------------------------------------------

save.pth <- paste0(Sys.getenv('drop_dir')
                   ,'seg-measures/by tract/broader ineq flows/'
                   ,'res-chars-long.csv')
save.pth

# write
write.csv(resl, save.pth
          ,row.names = F)
