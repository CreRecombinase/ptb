library(tidyverse)

full <- read_tsv("ga.full.qtltab")
base <- read_tsv("ga.base.qtltab")

res <- full %>% rename(FDR_full = FDR) %>% 
       select(region, FDR_full) %>%
       full_join(., base, by="region") %>%
       rename(FDR_base = FDR) %>%
       select(region, FDR_base, FDR_full)


plt = ggplot(res) + geom_point(aes(x=FDR_base, y=FDR_full)) +  
      geom_hline(yintercept=0.1, linetype=2) + 
      geom_vline(xintercept=0.1, linetype=2) 
