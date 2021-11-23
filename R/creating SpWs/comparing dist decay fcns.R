library(tidyverse)

devtools::load_all()


# to compare: ------------------------------------------------------------------

xlims <- c(0,25)

raw.dists <- seq(xlims[1], xlims[2], by = .25)
bisquares <-
  imap_dfc(c(one = 1, two.5 = 2.5, five=5,twenty.five = 25)
      ,~bisq.dist2weights(raw.dists, cutoff = .x))
bisquares <- bisquares %>%
  rename_with(~paste0('bisq_H.equals.', .x))
neg.exp <-
  exp(-raw.dists)

comps <- tibble(
  raw.dists = raw.dists,
  neg.exp = neg.exp,
  bisquares) %>%
  pivot_longer(2:6
               ,names_to = 'prox_fcn'
               ,values_to = 'spatial.weight')

comps %>%
  ggplot(aes(x = raw.dists
             ,y = spatial.weight
             ,color = prox_fcn)) +
  geom_line(size = 2) +
  scale_color_discrete() +
  ggtitle('Spatial weight by various proximity functions'
          ,'between 0-25 miles')
