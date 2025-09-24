library(ggdist)
library(tidyverse)
library(brms)

library(magrittr)
data(AustinCats, package = "rethinking")

d <-
    d %>% 
    mutate(black = ifelse(color == "Black", "black", "other"))

d <-
    d %>% 
    mutate(adopted  = ifelse(out_event == "Adoption", 1, 0),
           censored = ifelse(out_event != "Adoption", 1, 0))

glimpse(d)

d %>% 
    mutate(censored = factor(censored)) %>% 
    filter(days_to_event < 300) %>% 
    
    ggplot(aes(x = days_to_event, y = censored)) +
    # let's just mark off the 50% intervals
    stat_halfeye(.width = .5,  height = 4) +
    scale_y_discrete(NULL, labels = c("censored == 0", "censored == 1")) +
    coord_cartesian(ylim = c(1.5, 5.1)) + theme_minimal() + 
    theme(axis.ticks.y = element_blank(),
          text = element_text(size=16)) + labs(title = "Cat adoptions")

ggsave("slides/images/survival_ex.png", width=8, h = 3.5)
