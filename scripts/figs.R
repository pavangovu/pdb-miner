library(dplyr)
library(tidyr)
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
} else if (length(args)<3) {
  # default output file
  args[3] = "out.pdf"
}


df <- as_data_frame(fread(args[1]))
df <- df %>% rename(mol_type='4', chain='5') %>%
  na.omit() %>%
  filter(mol_type=='PROTEIN')

cols <- c('#EFBC7B', '#FE6367', '#5B1A18', '#D97134')

pdf(file = paste0(args[3], 'm_dist_vs_m_angle.pdf'))
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

pdf(file = paste0(args[3], 'asa_per_halide.pdf'))
df %>% group_by(halide, model_name, id) %>% summarise(asa=mean(asa)) %>% 
  ggplot(aes(x=halide, y=asa, fill=halide))+
  geom_violin(adjust=2, show.legend = F)+
  scale_fill_manual(values = cols)+
  theme_minimal()+
  scale_x_discrete('Halide')+
  scale_y_continuous('Accessible surface area')
dev.off()

pdf(file = paste0(args[3], 'asa_hist.pdf'))
df %>% group_by(halide, model_name, id) %>% summarise(asa=mean(asa)) %>% 
  ggplot(aes(asa, fill=halide))+
  geom_histogram(position = 'dodge', bins = 20, show.legend = T)+
  scale_fill_manual('Halide', values = cols)+
  theme_minimal()+
  scale_x_continuous('Accessible surface area',  limits = c(0,100))+
  scale_y_continuous('Count')
dev.off()

pdf(file = paste0(args[3], 'asa_distriburion.pdf'))
df %>% group_by(halide, model_name, id) %>% summarise(asa=mean(asa)) %>% 
  ggplot(aes(asa, fill=halide))+
  geom_density(position = 'stack', alpha=0.5)+
  scale_fill_manual('Halide', values = cols)+
  theme_minimal()+
  scale_x_continuous('Accessible surface area')+
  scale_y_continuous('Density')+
  xlim(0,100)
dev.off()

pdf(file = paste0(args[3], 'average_asa_per_halide.pdf'))
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

pdf(file = paste0(args[3], 'asa_vs_coordination_number.pdf'))
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
d <- df %>% group_by(halide, residue_aa) %>%
  summarise(n = n(), mean_dist=mean(distance), mean_angle=mean(angle)) %>%
  mutate(n=n/first(n)) %>% arrange(desc(n), halide) %>%
  arrange(halide, desc(n)) %>% mutate(id=row_number())

d5 <- d %>% top_n(n=5, wt=n)

pdf(file = paste0(args[3], 'most_frequent_AA.pdf'))
d %>% ggplot(aes(mean_dist, mean_angle))+
  geom_point(col='gray', size=3, alpha=0.5, show.legend = F)+
  facet_wrap(~halide)+
  geom_point(data = d5, aes(col=residue_aa), size=5, alpha=0.8)+
  geom_label_repel(data=d5, aes(label=residue_aa, col=residue_aa), show.legend = F)+
  scale_color_manual('', values = wes_palette(n = 20, 'GrandBudapest1', type = 'continuous'))+
  scale_x_continuous('Average distance of amino acids')+
  scale_y_continuous('Average angle of amino acids')+
  theme_bw()+
  theme(strip.background =element_rect(fill="indianred4"))+
  theme(strip.text = element_text(colour = 'white'))
dev.off()

d <- df %>% group_by(halide, residue_aa, atom_name) %>% summarise(n = n(), mean_dist=mean(distance), mean_angle=mean(angle)) %>% 
  mutate(atom = paste0(atom_name, '(', residue_aa, ')')) %>% ungroup() %>% select(halide, atom, mean_dist, mean_angle, n) %>% 
  mutate(n=n/first(n)) %>% arrange(desc(n), halide) %>% arrange(halide, desc(n)) %>% mutate(id=row_number()) %>% 
  group_by(halide)
d5 <- d %>% top_n(n = 5, wt = n)

pdf(file = paste0(args[3], 'most_frequent_ATOM.pdf'))
d %>% ggplot(aes(mean_dist, mean_angle))+
  geom_point(col='gray', size=3, alpha=0.5, show.legend = F)+
  facet_wrap(~halide)+
  geom_point(data = d5, aes(col=atom), size=5, alpha=0.8)+
  geom_label_repel(data=d5, aes(label=atom, col=atom), show.legend = F)+
  scale_color_manual('', values = wes_palette(n = 20, 'GrandBudapest1', type = 'continuous'))+
  scale_x_continuous('Average distance of amino acids')+
  scale_y_continuous('Average angle of amino acids')+
  theme_bw()+
  theme(strip.background =element_rect(fill="indianred4"))+
  theme(strip.text = element_text(colour = 'white'))
dev.off()

# distance distrib
pdf(file = paste0(args[3], 'average_distance_per_halide.pdf'))
df %>% group_by(halide) %>% summarise(mean_dist=mean(distance), sd_dist = sd(distance), se_dist = sd(distance/sqrt(length(distance)))) %>% 
  ggplot(aes(halide, mean_dist, fill=halide))+
  geom_bar(stat='identity', show.legend = F)+
  geom_errorbar(aes(ymin=mean_dist-se_dist, ymax=mean_dist+se_dist), width=0.1)+
  scale_fill_manual(values = wes_palette(n = 4, 'GrandBudapest1'))+
  scale_x_discrete('Halide')+
  scale_y_continuous('Average distance')+
  theme_minimal()
dev.off()

pdf(file = paste0(args[3], 'average_angle_per_halide.pdf'))
df %>% filter(angle>0) %>% group_by(halide) %>% dplyr::summarise(mean_angle=mean(angle), sd_angle = sd(angle), se_angle=sd(angle)/sqrt(length(angle)))%>% 
  ggplot(aes(halide, mean_angle, fill=halide))+
  geom_bar(stat='identity', show.legend = F)+
  geom_errorbar(aes(ymin=mean_angle-se_angle, ymax=mean_angle+se_angle), width=0.1)+
  scale_fill_manual(values = wes_palette(n = 4, 'GrandBudapest1'))+
  scale_x_discrete('Halide')+
  scale_y_continuous('Average angle')+
  theme_minimal()
dev.off()

pdf(file = paste0(args[3], 'distance_distribution.pdf'))
df %>% 
  ggplot(aes(halide, distance, fill=halide))+
  geom_violin(show.legend = F)+
  scale_fill_manual(values = wes_palette(n = 4, 'GrandBudapest1'))+
  scale_x_discrete('Halide')+
  scale_y_continuous('Distance')+
  theme_minimal()
dev.off()

# angle distrib
pdf(file = paste0(args[3], 'angle_distribution.pdf'))
df %>% filter(angle>0) %>% 
  ggplot(aes(halide, angle, fill=halide))+
  geom_violin(show.legend = F)+
  scale_fill_manual(values = wes_palette(n = 4, 'GrandBudapest1'))+
  scale_x_discrete('Halide')+
  scale_y_continuous('Angle')+
  theme_minimal()
dev.off()

# coord density
pdf(file = paste0(args[3], 'coordination_number_density.pdf'))
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
pdf(file = paste0(args[3], 'coordination_number_hist.pdf'))
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
pdf(file = paste0(args[3], 'coordination_number_hist_2.pdf'))
df %>% group_by(model_name, halide, id) %>% summarise(N = n()) %>% ungroup() %>% group_by(halide) %>% 
  ggplot(aes(N, fill=halide))+
  geom_histogram(col='black', bins = 11)+
  scale_y_continuous('Count')+
  scale_fill_manual('Halide', values = wes_palette(n = 4, 'GrandBudapest1'))+
  facet_wrap(halide~., scales = 'free')+
  theme_bw()+
  theme(strip.background = element_rect(fill="indianred4"))+
  theme(strip.text = element_text(colour = 'white'))+
  scale_x_continuous('Coordination number', limits = c(0,11))
dev.off()

# coord outlayers
pdf(file = paste0(args[3], 'coordination_number_outlayers.pdf'))
df %>% group_by(model_name, halide, id) %>% summarise(N = n()) %>% ungroup() %>% group_by(halide) %>% 
  ggplot(aes(halide, N, fill=halide))+
  geom_boxplot(show.legend = F)+
  scale_x_discrete('Halide')+
  scale_y_continuous('Coordination number')+
  scale_fill_manual('Halide', values = wes_palette(n = 4, 'GrandBudapest1'))+
  theme_bw()
dev.off()

pdf(file = paste0(args[3], 'distance_vs_angle.pdf'))
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

pdf(file = paste0(args[3], 'distance_vs_angle_vs_asa.pdf'))
with(df, 
     scatterplot3d(distance, angle, asa,
                   pch = 16, color = cols[as.numeric(df$code_halide)],
                   xlab = "Distance",
                   ylab = "Angle",
                   zlab = "Accessible surface area", grid = FALSE, box=T))
legend("right", legend = levels(df$code_halide),
       col =  c('#EFBC7B', '#FE6367', '#5B1A18', '#D97134'), pch = 16)
dev.off()


##############################
##### AA compositions ########
##############################

d <- as_data_frame(fread(args[2], header = F))
head(d)
d <- d %>% separate(V1, c('model','halide','ID', 'atoms')) %>% group_by(halide, V2) %>% summarise(n=n()) %>%
  ungroup() %>% group_by(halide) %>% arrange(desc(n)) %>% top_n(n = 20, wt = n)
# %>% filter(halide!='F')
sum_d <- d %>% summarise(sum=sum(n))
d <- as_data_frame(merge(d, sum_d))

sum_d <- d %>% ungroup() %>% group_by(halide) %>% summarise(sum=sum(n), min_n=min(n)/sum)

d <- as_data_frame(merge(d, sum_d))

d_CL <- d %>% dplyr::mutate(N=n/sum) %>% filter(halide=='CL') %>% filter(N>min_n)
d_BR <- d %>% dplyr::mutate(N=n/sum) %>% filter(halide=='BR') %>% filter(N>min_n)
d_I <- d %>% dplyr::mutate(N=n/sum) %>% filter(halide=='I') %>% filter(N>min_n)

pdf(file = paste0(args[3], 'compositions.pdf'))
bind_rows(d_CL, d_BR, d_I) %>% 
  ggplot(aes(x=V2, y=n/sum, fill=halide))+
  geom_bar(stat='identity')+
  # facet_grid(~halide)+
  theme_minimal()+
  scale_fill_manual('Halide', values = cols)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(strip.background =element_rect(fill="indianred4"))+
  theme(strip.text = element_text(colour = 'white'))+
  scale_x_discrete('Compostition')+
  scale_y_continuous('Count')
dev.off()