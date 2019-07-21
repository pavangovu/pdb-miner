setwd('~/Projects/Halide/scripts/')

library(dplyr)
library(tidyr)
library(ggplot2)
library(gghighlight)
library(data.table)
library(ggrepel)
library(ggfortify)
library(wesanderson)
library(plot3D)
library(ggbiplot)


df %>% group_by(model_name, id, halide) %>% summarise(m_dist=mean(distance), m_angle=mean(angle), m_asa=mean(asa)) %>% 
  ggplot(aes(x=m_dist, y=m_angle, z=m_asa, col=halide))+
  geom_point()+
  facet_wrap(~halide)+
  theme_bw()

df <- as_data_frame(fread('../full_data.tsv'))
df <- df %>% na.omit()
df <- df %>% filter(distance>=0.23)

cols <- c('#EFBC7B', '#FE6367', '#5B1A18', '#D97134')
df
df %>% group_by(halide, model_name, id) %>% summarise(asa=mean(asa)) %>% 
  ggplot(aes(x=halide, y=asa, fill=halide))+
  geom_violin(adjust=2, show.legend = F)+
  scale_fill_manual(values = cols)+
  theme_minimal()+
  scale_x_discrete('Halide')+
  scale_y_continuous('Accessible surface area')

df %>% group_by(halide, model_name, id) %>% summarise(asa=mean(asa)) %>% 
  ggplot(aes(asa, fill=halide))+
  geom_histogram(position = 'dodge', bins = 20, show.legend = T)+
  scale_fill_manual('Halide', values = cols)+
  theme_minimal()+
  scale_x_continuous('Accessible surface area',  limits = c(0,100))+
  scale_y_continuous('Count')

df %>% group_by(halide, model_name, id) %>% summarise(asa=mean(asa)) %>% 
  ggplot(aes(asa, fill=halide))+
  geom_density(position = 'stack', alpha=0.5)+
  scale_fill_manual('Halide', values = cols)+
  theme_minimal()+
  scale_x_continuous('Accessible surface area')+
  scale_y_continuous('Density')+
  xlim(0,100)

df %>% group_by(halide, model_name, id) %>% summarise(count=n(), asa=mean(asa)) %>% ungroup() %>% group_by(halide) %>% 
  summarise(mean_asa=mean(asa), se_asa=sd(asa)/sqrt(length(asa))) %>% 
  ggplot(aes(x=halide, y=mean_asa, fill=halide))+
  geom_bar(stat = 'identity', show.legend = F)+
  geom_errorbar(aes(ymin=mean_asa-se_asa, ymax=mean_asa+se_asa), width=.1)+
  theme_minimal()+
  scale_fill_manual('Halide', values = cols)+
  scale_x_discrete('Halide')+
  scale_y_continuous('Average accessible surface area')

df %>% group_by(halide, model_name, id) %>% summarise(count=n(), asa=mean(asa)) %>% 
ggplot(aes(x=asa, y=count, fill=halide, col=halide))+
  geom_point(show.legend = F)+
  geom_smooth(col='firebrick', fill='lightcoral')+
  facet_wrap(~halide)+
  scale_color_manual(values = cols)+
  theme_bw()+
  theme(strip.background =element_rect(fill="indianred4"))+
  theme(strip.text = element_text(colour = 'white'))+
  scale_x_continuous('Accessible surface area')+
  scale_y_continuous('Coordination number', limits=c(0, 20))

  # freq
df %>% group_by(halide, residue_aa) %>%
  summarise(n = n(), mean_dist=mean(distance), mean_angle=mean(angle)) %>%
  mutate(n=n/first(n)) %>% arrange(desc(n), halide) %>%
  arrange(halide, desc(n)) %>% mutate(id=row_number()) %>%
      ggplot(aes(mean_dist, mean_angle))+
      geom_point(aes(col=residue_aa, size=desc(id)), alpha=0.5, show.legend = F)+
      gghighlight(id<6, label_key = residue_aa, unhighlighted_colour = ggplot2::alpha("grey", 0.1), use_direct_label = F, use_group_by = F)+
      geom_label_repel(aes(label=residue_aa, col=residue_aa), show.legend = F)+
      scale_color_manual(values = wes_palette(n = 8, 'GrandBudapest1', type = 'continuous'))+
      facet_wrap(~halide)+
      scale_x_continuous('Average distance of amino acids')+
      scale_y_continuous('Average angle of amino acids')+
      theme_bw()+
      theme(strip.background =element_rect(fill="indianred4"))+
      theme(strip.text = element_text(colour = 'white'))+
      theme(legend.position = "none",
            panel.grid = element_blank())

# distance distrib
df %>% group_by(halide) %>% summarise(mean_dist=mean(distance), sd_dist = sd(distance), se_dist = sd(distance/sqrt(length(distance)))) %>% 
  ggplot(aes(halide, mean_dist, fill=halide))+
  geom_bar(stat='identity', show.legend = F)+
  geom_errorbar(aes(ymin=mean_dist-se_dist, ymax=mean_dist+se_dist), width=0.1)+
  scale_fill_manual(values = wes_palette(n = 4, 'GrandBudapest1'))+
  scale_x_discrete('Halide')+
  scale_y_continuous('Average distance')+
  theme_minimal()

df %>% group_by(halide) %>% dplyr::summarise(mean_angle=mean(angle), sd_angle = sd(angle), se_angle=sd(angle)/sqrt(length(angle)))%>% 
  ggplot(aes(halide, mean_angle, fill=halide))+
  geom_bar(stat='identity', show.legend = F)+
  geom_errorbar(aes(ymin=mean_angle-se_angle, ymax=mean_angle+se_angle), width=0.1)+
  scale_fill_manual(values = wes_palette(n = 4, 'GrandBudapest1'))+
  scale_x_discrete('Halide')+
  scale_y_continuous('Average angle')+
  theme_minimal()

df %>% 
  ggplot(aes(halide, distance, fill=halide))+
  geom_violin(show.legend = F)+
  scale_fill_manual(values = wes_palette(n = 4, 'GrandBudapest1'))+
  scale_x_discrete('Halide')+
  scale_y_continuous('Distance')+
  theme_minimal()
# angle distrib
df %>% 
  ggplot(aes(halide, angle, fill=halide))+
  geom_violin(show.legend = F)+
  scale_fill_manual(values = wes_palette(n = 4, 'GrandBudapest1'))+
  scale_x_discrete('Halide')+
  scale_y_continuous('Angle')+
  theme_minimal()


# coord density
df %>% group_by(model_name, halide, id) %>% summarise(N = n()) %>% ungroup() %>% group_by(halide) %>% 
  ggplot(aes(N, fill=halide))+
  geom_density(col='black', show.legend = F)+
  scale_x_continuous('Distribution of coordination number')+
  scale_y_continuous('Density')+
  scale_fill_manual(values = wes_palette(n = 4, 'GrandBudapest1'))+
  facet_wrap(halide~.)+
  theme_bw()+
  theme(strip.background =element_rect(fill="indianred4"))+
  theme(strip.text = element_text(colour = 'white'))+
  xlim(0,11)

df %>% group_by(model_name, halide, id) %>% summarise(N = n()) %>% ungroup() %>% group_by(halide) %>% 
  ggplot(aes(N, fill=halide))+
  geom_density(col='black', alpha=0.5, position = 'stack')+
  scale_y_continuous('Density')+
  scale_fill_manual('Halide', values = wes_palette(n = 4, 'GrandBudapest1'))+
  theme_bw()+
  theme(strip.background =element_rect(fill="indianred4"))+
  theme(strip.text = element_text(colour = 'white'))+
  scale_x_continuous('Coordination number', limits = c(0,11))

# hist
df %>% group_by(model_name, halide, id) %>% summarise(N = n()) %>% ungroup() %>% group_by(halide) %>% 
  ggplot(aes(N, fill=halide))+
  geom_histogram(col='black', position = 'dodge', bins = 11)+
  scale_y_continuous('Count')+
  scale_fill_manual('Halide', values = wes_palette(n = 4, 'GrandBudapest1'))+
  # facet_wrap(halide~.)+
  theme_bw()+
  theme(strip.background =element_rect(fill="indianred4"))+
  theme(strip.text = element_text(colour = 'white'))+
  scale_x_continuous('Coordination number', limits = c(0,11))

df %>% group_by(model_name, halide, id) %>% summarise(N = n()) %>% ungroup() %>% group_by(halide) %>% 
  ggplot(aes(N, fill=halide))+
  geom_histogram(col='black', bins = 11)+
  scale_y_continuous('Count')+
  scale_fill_manual('Halide', values = wes_palette(n = 4, 'GrandBudapest1'))+
  facet_wrap(halide~.)+
  theme_bw()+
  theme(strip.background =element_rect(fill="indianred4"))+
  theme(strip.text = element_text(colour = 'white'))+
  scale_x_continuous('Coordination number', limits = c(0,11))


# coord outlayers
df %>% group_by(model_name, halide, id) %>% summarise(N = n()) %>% ungroup() %>% group_by(halide) %>% 
  ggplot(aes(halide, N, fill=halide))+
  geom_boxplot(show.legend = F)+
  scale_x_discrete('Halide')+
  scale_y_continuous('Coordination number')+
  scale_fill_manual('Halide', values = wes_palette(n = 4, 'GrandBudapest1'))+
  theme_bw()


c('#FFFEFF', '#34BCD2', '#01A08A', '#FF7D31', '#FF0020')
pal <- wes_palette(n = 500, "Darjeeling1", type = "continuous")

df %>% ggplot(aes(distance, angle))+
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+
  # scale_fill_gradientn('Density', colours = rev(pal))+
  scale_fill_gradientn('Density', colours = c('#FFFEFF', '#34BCD2', '#01A08A', 'goldenrod2', 'tan2', '#FF7D31', '#FF0020', 'indianred4'))+
  scale_x_continuous('Distance', expand = c(0, 0)) +
  scale_y_continuous('Angle', expand = c(0, 0)) +
  facet_wrap(~halide)+
  theme_bw()+
  # theme(legend.position='none')+
  geom_smooth(col='red2', se=F)+
  theme(strip.background =element_rect(fill='indianred4'))+
  theme(strip.text = element_text(colour = 'white'))


  
require(mgcv)
library(scatterplot3d)
library(tidyverse)
library(RNHANES)

df <- as_data_frame(fread('../full_data.tsv'))
df <- df %>% na.omit()
df <- df %>% filter(distance>=0.23)
df$code_halide <- unclass(factor(df$halide))

cols <- c('#EFBC7B', '#FE6367', '#5B1A18', '#D97134')
with(df, 
     scatterplot3d(distance, angle, asa,
                   pch = 16, color = cols[as.numeric(df$code_halide)],
                   xlab = "Distance",
                   ylab = "Angle",
                   zlab = "Accessible surface area", grid = FALSE, box=T))
legend("right", legend = levels(df$code_halide),
       col =  c('#EFBC7B', '#FE6367', '#5B1A18', '#D97134'), pch = 16)


library(gridExtra)
grid.arrange()
devtools::install_github("AckerDWM/gg3D")

library("gg3D")

ggplot(df, aes(x=distance, y=angle, z=asa, col=halide))+
  stat_wireframe(alpha=.5)+
  theme_void()+
  axes_3D()+
  stat_3D()

scatter3D(df$distance, df$angle, df$asa, col = c("#1B9E77", "#D95F02", "#7570B3", "#7570B4"))

plot_ly(data=df, x=df$distance, y=df$angle, z=df$asa, color=df$halide, alpha = 0.01, sizes = 0.01, type="scatter3d")
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
  scale_fill_manual('Halide', values = wes_palette('GrandBudapest1', type = 'continuous'))+
  theme_bw()

df %>% ggplot(aes(distance, angle, col=density))+
  geom_point(alpha=0.1, size=0.01)+
  geom_smooth(col='red', size=0.5)+
  geom_rug(alpha=.2)+
  facet_wrap(~residue_aa)+
  theme_bw()

df_for_PCA_full <- df %>% select(-atom_name, -residue_aa) %>% group_by(model_name, id, halide) %>%
  summarise(n=n(), dist=mean(distance), angle=mean(angle), asa=mean(asa))

df_pca <- df_for_PCA_full %>% ungroup %>% select(-model_name, -halide, -id) %>% 
  dplyr::rename(Distance=dist, Angle=angle, 'Accession surface area'=asa, 'Coord.number'=n)

pca <- prcomp(df_pca, scale = TRUE, center = T)

ggbiplot(pca, obs.scale = 1, var.scale = 1, groups = df_for_PCA_full$halide, ellipse = T, circle = F, alpha = 0.025)+
  theme_minimal()+
  theme(legend.direction = 'horizontal', legend.position = 'top')+
  scale_color_manual('Halide', values = wes_palette("Darjeeling1"))+
  ylim(-3,3)
