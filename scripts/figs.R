library(dplyr)
library(ggplot2)
library(gghighlight)
library(data.table)
library(ggrepel)
library(wesanderson)
library(scatterplot3d)



args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.png"
}


df <- as_data_frame(fread(args[1]))
df <- df %>% rename(mol_type='4', chain='5') %>%
  na.omit() %>%
  filter(mol_type=='PROTEIN')

cols <- c('#EFBC7B', '#FE6367', '#5B1A18', '#D97134')

png(filename = paste0(args[2], 'm_dist_vs_m_angle.png'))
df %>% group_by(model_name, id, halide) %>% filter(angle>0) %>% summarise(m_dist=mean(distance), m_angle=mean(angle), m_asa=mean(asa)) %>% 
  ggplot(aes(x=m_dist, y=m_angle, z=m_asa, col=halide))+
  geom_point(show.legend = F)+
  facet_wrap(~halide)+
  scale_color_manual(values = cols)+
  theme_bw()+
  theme(strip.background =element_rect(fill="indianred4"))+
  theme(strip.text = element_text(colour = 'white'))+
  scale_x_continuous('Average distance')+
  scale_y_continuous('Average angle')
dev.off()

png(filename = paste0(args[2], 'asa_per_halide.png'))
df %>% group_by(halide, model_name, id) %>% summarise(asa=mean(asa)) %>% 
  ggplot(aes(x=halide, y=asa, fill=halide))+
  geom_violin(adjust=2, show.legend = F)+
  scale_fill_manual(values = cols)+
  theme_minimal()+
  scale_x_discrete('Halide')+
  scale_y_continuous('Accessible surface area')
dev.off()

png(filename = paste0(args[2], 'asa_hist.png'))
df %>% group_by(halide, model_name, id) %>% summarise(asa=mean(asa)) %>% 
  ggplot(aes(asa, fill=halide))+
  geom_histogram(position = 'dodge', bins = 20, show.legend = T)+
  scale_fill_manual('Halide', values = cols)+
  theme_minimal()+
  scale_x_continuous('Accessible surface area',  limits = c(0,100))+
  scale_y_continuous('Count')
dev.off()

png(filename = paste0(args[2], 'asa_distriburion.png'))
df %>% group_by(halide, model_name, id) %>% summarise(asa=mean(asa)) %>% 
  ggplot(aes(asa, fill=halide))+
  geom_density(position = 'stack', alpha=0.5)+
  scale_fill_manual('Halide', values = cols)+
  theme_minimal()+
  scale_x_continuous('Accessible surface area')+
  scale_y_continuous('Density')+
  xlim(0,100)
dev.off()

png(filename = paste0(args[2], 'average_asa_per_halide.png'))
df %>% group_by(halide, model_name, id) %>% summarise(count=n(), asa=mean(asa)) %>% ungroup() %>% group_by(halide) %>% 
  summarise(mean_asa=mean(asa), se_asa=sd(asa)/sqrt(length(asa))) %>% 
  ggplot(aes(x=halide, y=mean_asa, fill=halide))+
  geom_bar(stat = 'identity', show.legend = F)+
  geom_errorbar(aes(ymin=mean_asa-se_asa, ymax=mean_asa+se_asa), width=.1)+
  theme_minimal()+
  scale_fill_manual('Halide', values = cols)+
  scale_x_discrete('Halide')+
  scale_y_continuous('Average accessible surface area')
dev.off()

png(filename = paste0(args[2], 'asa_vs_coordination_number.png'))
df %>% group_by(halide, model_name, id) %>% summarise(count=n(), asa=mean(asa)) %>% 
  ggplot(aes(x=asa, y=count, fill=halide, col=halide))+
  geom_point(show.legend = F)+
  geom_smooth(method = 'lm',col='firebrick', fill='lightcoral')+
  facet_wrap(~halide)+
  scale_color_manual(values = cols)+
  theme_bw()+
  theme(strip.background =element_rect(fill="indianred4"))+
  theme(strip.text = element_text(colour = 'white'))+
  scale_x_continuous('Accessible surface area')+
  scale_y_continuous('Coordination number', limits=c(0, 20))
dev.off()

# freq
png(filename = paste0(args[2], 'most_frequent_AA.png'))
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
dev.off()

# distance distrib
png(filename = paste0(args[2], 'average_distance_per_halide.png'))
df %>% group_by(halide) %>% summarise(mean_dist=mean(distance), sd_dist = sd(distance), se_dist = sd(distance/sqrt(length(distance)))) %>% 
  ggplot(aes(halide, mean_dist, fill=halide))+
  geom_bar(stat='identity', show.legend = F)+
  geom_errorbar(aes(ymin=mean_dist-se_dist, ymax=mean_dist+se_dist), width=0.1)+
  scale_fill_manual(values = wes_palette(n = 4, 'GrandBudapest1'))+
  scale_x_discrete('Halide')+
  scale_y_continuous('Average distance')+
  theme_minimal()
dev.off()

png(filename = paste0(args[2], 'average_angle_per_halide.png'))
df %>% filter(angle>0) %>% group_by(halide) %>% dplyr::summarise(mean_angle=mean(angle), sd_angle = sd(angle), se_angle=sd(angle)/sqrt(length(angle)))%>% 
  ggplot(aes(halide, mean_angle, fill=halide))+
  geom_bar(stat='identity', show.legend = F)+
  geom_errorbar(aes(ymin=mean_angle-se_angle, ymax=mean_angle+se_angle), width=0.1)+
  scale_fill_manual(values = wes_palette(n = 4, 'GrandBudapest1'))+
  scale_x_discrete('Halide')+
  scale_y_continuous('Average angle')+
  theme_minimal()
dev.off()

png(filename = paste0(args[2], 'distance_distribution.png'))
df %>% 
  ggplot(aes(halide, distance, fill=halide))+
  geom_violin(show.legend = F)+
  scale_fill_manual(values = wes_palette(n = 4, 'GrandBudapest1'))+
  scale_x_discrete('Halide')+
  scale_y_continuous('Distance')+
  theme_minimal()
dev.off()

# angle distrib
png(filename = paste0(args[2], 'angle_distribution.png'))
df %>% filter(angle>0) %>% 
  ggplot(aes(halide, angle, fill=halide))+
  geom_violin(show.legend = F)+
  scale_fill_manual(values = wes_palette(n = 4, 'GrandBudapest1'))+
  scale_x_discrete('Halide')+
  scale_y_continuous('Angle')+
  theme_minimal()
dev.off()

# coord density
png(filename = paste0(args[2], 'coordination_number_density.png'))
df %>% group_by(model_name, halide, id) %>% summarise(N = n()) %>% ungroup() %>% group_by(halide) %>% 
  ggplot(aes(N, fill=halide))+
  geom_density(col='black', alpha=0.5, position = 'stack')+
  scale_y_continuous('Density')+
  scale_fill_manual('Halide', values = wes_palette(n = 4, 'GrandBudapest1'))+
  theme_bw()+
  theme(strip.background =element_rect(fill="indianred4"))+
  theme(strip.text = element_text(colour = 'white'))+
  scale_x_continuous('Coordination number', limits = c(0,11))
dev.off()

# coord hist
png(filename = paste0(args[2], 'coordination_number_hist.png'))
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
dev.off()

# coord hist 2
png(filename = paste0(args[2], 'coordination_number_hist_2.png'))
df %>% group_by(model_name, halide, id) %>% summarise(N = n()) %>% ungroup() %>% group_by(halide) %>% 
  ggplot(aes(N, fill=halide))+
  geom_histogram(col='black', bins = 11)+
  scale_y_continuous('Count')+
  scale_fill_manual('Halide', values = wes_palette(n = 4, 'GrandBudapest1'))+
  facet_wrap(halide~.)+
  theme_bw()+
  theme(strip.background = element_rect(fill="indianred4"))+
  theme(strip.text = element_text(colour = 'white'))+
  scale_x_continuous('Coordination number', limits = c(0,11))
dev.off()

# coord outlayers
png(filename = paste0(args[2], 'coordination_number_outlayers.png'))
df %>% group_by(model_name, halide, id) %>% summarise(N = n()) %>% ungroup() %>% group_by(halide) %>% 
  ggplot(aes(halide, N, fill=halide))+
  geom_boxplot(show.legend = F)+
  scale_x_discrete('Halide')+
  scale_y_continuous('Coordination number')+
  scale_fill_manual('Halide', values = wes_palette(n = 4, 'GrandBudapest1'))+
  theme_bw()
dev.off()

png(filename = paste0(args[2], 'distance_vs_angle.png'))
df %>% filter(angle>0) %>% ggplot(aes(distance, angle))+
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
dev.off()

df$code_halide <- unclass(factor(df$halide))

png(filename = paste0(args[2], 'distance_vs_angle_vs_asa.png'))
with(df, 
     scatterplot3d(distance, angle, asa,
                   pch = 16, color = cols[as.numeric(df$code_halide)],
                   xlab = "Distance",
                   ylab = "Angle",
                   zlab = "Accessible surface area", grid = FALSE, box=T))
legend("right", legend = levels(df$code_halide),
       col =  c('#EFBC7B', '#FE6367', '#5B1A18', '#D97134'), pch = 16)
dev.off()
