# ws ---------------------------------------------------------------------------
rm(list = ls())
require(tidyverse)
require(sf)

# option setting
sf_use_s2(T)
options(tigris_use_cache = TRUE)




# get demovars -----------------------------------------------------------------

ddir <- paste0(Sys.getenv('drop_dir'),
               'seg-measures/by tract/broader ineq flows/')
flr <- ddir %>%
  list.files(pattern='flow-res.csv$'
             ,full.names = T) %>%
  vroom::vroom() %>%
  select(-`...1`)


# flow estimates vary by region type, but not residential. Becuase only flows in a
# given region are used
flr %>% filter(geoid == '46013951900') %>%  select(1:3,pop,hh,
                                                   matches('^residen'))

res <- flr %>%
  select(geoid,pop,hh,
         matches('^residen')) %>% distinct()


# re-add region xwalks
res <- res %>%
  mutate(countyfp = substr(geoid, 1,5)
         ,.after = geoid) %>%
    left_join(geox::rx[c('countyfp','cz','cbsa')])

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

colnames(res)
res


# rearragne --------------------------------------------------------------------

res <- res %>%
  select(1,2,cz,cbsa,pop,hh,
         everything(), -geoid10)

# descriptives/other checks -----------------------------------------------------

# most range from abt 0-1
# res %>% map(summary)

res %>% filter(duplicated(geoid))

# make long and organize weights ----------------------------------------------

colnames(res)
vvars <- setdiff(colnames(res)
                 ,c('geoid', 'countyfp', 'cz','cbsa', 'pop', 'hh') )
# make res long
resl <- res %>%
  pivot_longer(vvars)


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

resl <- resl %>% select(-countyfp)
resl

save.pth <- paste0(Sys.getenv('drop_dir')
                   ,'seg-measures/by tract/broader ineq flows/'
                   ,'res-chars-long.csv')
save.pth

# write
write.csv(resl, save.pth
          ,row.names = F)
