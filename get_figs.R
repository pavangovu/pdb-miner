setwd('~/Projects/Halide/scripts/')

library(dplyr)
library(tidyr)
library(ggplot2)
library(gghighlight)
library(data.table)
library(ggrepel)
library(ggfortify)

df <- as_data_frame(fread('../full_data_CA.tsv'))
df <- df %>% na.omit()

# freq
df %>% group_by(halide, residue_aa) %>%
  summarise(n = n(), mean_dist=mean(distance), mean_angle=mean(angle)) %>%
  mutate(n=n/first(n)) %>% arrange(desc(n), halide) %>%
  arrange(halide, desc(n)) %>% mutate(id=row_number()) %>%
      ggplot(aes(mean_dist, mean_angle))+
      geom_point(aes(col=residue_aa, size=desc(id)), alpha=0.5, show.legend = F)+
      gghighlight(id<6, label_key = residue_aa, unhighlighted_colour = ggplot2::alpha("grey", 0.1), use_direct_label = F, use_group_by = F)+
      geom_label_repel(aes(label=residue_aa, col=residue_aa), show.legend = F)+
      facet_wrap(~halide)+
      theme_bw()

# distance distrib

df %>% group_by(halide) %>% summarise(mean_dist=mean(distance), mean_angle=mean(angle), sd_dist = sd(distance), sd_angel = sd(angle))%>% 
  ggplot(aes(halide, mean_dist, fill=halide))+
  geom_bar(stat='identity')+
  geom_errorbar(aes(ymin=mean_dist-sd_dist, ymax=mean_dist+sd_dist), width=0.1)+
  theme_bw()
df %>% group_by(halide) %>% summarise(mean_dist=mean(distance), mean_angle=mean(angle), sd_dist = sd(distance), sd_angel = sd(angle))%>% 
  ggplot(aes(halide, mean_angle, fill=halide))+
  geom_bar(stat='identity')+
  geom_errorbar(aes(ymin=mean_angle-sd_angle, ymax=mean_angle+sd_angle), width=0.1)+
  theme_bw()

df %>% 
  ggplot(aes(halide, distance, fill=halide))+
  geom_violin()+
  theme_bw()

# angle distrib
df %>% 
  ggplot(aes(halide, angle, fill=halide))+
  geom_violin()+
  theme_bw()

# coord density
df %>% group_by(model_name, halide, id) %>% summarise(N = n()) %>% ungroup() %>% group_by(halide) %>% 
  ggplot(aes(N, fill=halide))+
  geom_density(col='black')+
  scale_x_continuous('distribution of coordination number')+
  scale_fill_discrete('Halide')+
  facet_wrap(halide~.)+
  theme_bw()

# coord outlayers
df %>% group_by(model_name, halide, id) %>% summarise(N = n()) %>% ungroup() %>% group_by(halide) %>% 
  ggplot(aes(halide, N, fill=halide))+
  geom_boxplot()+
  scale_x_discrete('Halide')+
  scale_y_continuous('coordination number')+
  scale_fill_discrete('Halide')+
  theme_bw()

# corr+scatter
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

df$density <- get_density(df$distance, df$angle, 10)
df %>% ggplot(aes(distance, angle, col=density))+
  geom_point(alpha=0.1, size=0.5)+
  geom_smooth(se = F, col='red')+
  geom_rug(alpha=.2)+
  facet_wrap(~halide)+
  theme_bw()

df %>% ggplot(aes(distance, angle, col=density))+
  geom_point(alpha=0.1, size=0.01)+
  geom_smooth(col='red', size=0.5)+
  geom_rug(alpha=.2)+
  facet_wrap(~residue_aa)+
  theme_bw()

