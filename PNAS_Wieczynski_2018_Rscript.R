


#
##### LOAD PACKAGES -----
library(ggplot2)
library(dplyr)
library(tidyr)
library(raster)
library(stringr)
library(grid)
library(gridExtra)
library(gtable)
library(moments)
library(colorRamps)
library(colorspace)
library(RColorBrewer)
library(corrplot)
library(ggcorrplot)
library(BIEN)
library(scales)
library(BSDA)
library(boot)
library(lemon)
library(vegan)
library(ggfortify)
library(cluster)

#
#####
##### DATA SELECTION -----
#

## Clear Memory
rm(list=ls())


## Set working directory

directory<-"~/"


## Select DATASET from following options:   "All"   "InSitu"   "Bien"
dataset<-"All"


## Select WEIGHTING from following options:   "Abundance"   "Basal Area"   "Species"
weighting<-"Abundance"


#
#####
##### PRIMARY DATA ANALYSIS -----
#####
### Import Data -----
# import raw trait data
trait_data_full<-data.frame(read.csv(paste0(directory, "PNAS_Wieczynski_2018_MainTraitData.csv")))
trait_data<-if(dataset=="All") {trait_data_full} else {filter(trait_data_full, Dataset==dataset)}


# import full plot metadata
metadata<-data.frame(read.csv(paste0(directory, "PNAS_Wieczynski_2018_MainPlotMetadata.csv")))



env_groups<-data.frame(Variable_1=c("Abs_Latitude", "Elevation", "Mean_Annual_Temp", "Isothermality", "Temp_Diurnal_Range", "Temp_Annual_Range", "Temp_Seasonality",
                                    "Annual_Precip", "Precip_Seasonality", "CWD", "vap_yr", "vpd_yr", "pet_yr", "PET_MAP_ratio", "wind_yr", "srad_yr", "Species_Richness"), 
                       Env_grp=c("Geography", "Geography", rep.int("Temperature", 5), rep.int("Precipitation", 3), rep.int("Vapor Pressure", 4), rep.int("Other", 3)),
                       Text=c("|Lat|", "Elev", "MAT", "ISO", "TDR", "TAR", "TS", "MAP", "PS", "CWD", "VAP", "VPD", "PET", "  PET:MAP", "Wind", "SR", "Sp. Richness")) 
env_groups$Env_grp<-factor(env_groups$Env_grp, levels=c("Geography", "Temperature", "Precipitation", "Vapor Pressure", "Other"))

env_list<-env_groups$Variable_1 %>%
  as.character()

t_order<-c("SLA", "Wood_Density", "logHeight", "logLeafArea", "logSeedMass", "percent_C", "percent_N", "percent_P", "CN_ratio", "CP_ratio", "NP_ratio")
trait_groups<-data.frame(Trait=t_order, trait_group=c(rep.int("Morphology", 5), rep.int("Element_composition", 3), rep.int("Element_ratios", 3))) %>%
  mutate(order=row_number())

trait_list<-c("SLA", "Wood_Density", "logHeight", "logLeafArea", "logSeedMass", "percent_C", "percent_N", "percent_P", "NP_ratio")

#
### calculate weighted moments -----

trait_melt_outRemoved<-trait_data[-1] %>%
  gather(Trait, value, -c(1:7)) %>%
  drop_na(value) %>%
  mutate(value=as.numeric(as.character(value))) %>%
  group_by(Trait) %>% 
  filter(value>quantile(value,0.005, na.rm=T), value<quantile(value,0.995, na.rm=T)) %>%
  ungroup() %>%
  rename(Trait_Mean=value)


local_weighted_traits<-if(weighting=="Species")
{trait_melt_outRemoved %>% mutate(abund_wx=Trait_Mean, ba_wx=Trait_Mean, abund_n0=ifelse(individual_count==0,0,1), ba_n0=ifelse(is.na(BA) | BA==0,0,1))} else
    {trait_melt_outRemoved %>% mutate(abund_wx=Trait_Mean*individual_count, ba_wx=Trait_Mean*BA, abund_n0=ifelse(individual_count==0,0,1), ba_n0=ifelse(is.na(BA) | BA==0,0,1))} 


global_weighted_traits<-trait_melt_outRemoved %>%
  group_by(Family,Genus,Species,Trait) %>%
  mutate(Trait_Mean=mean(Trait_Mean, na.rm=T), individual_count=sum(individual_count, na.rm=T), BA=ifelse(all(is.na(BA)), NA, sum(BA, na.rm=T))) %>%
  ungroup %>%
  mutate(abund_wx=Trait_Mean*individual_count, ba_wx=Trait_Mean*BA, abund_n0=ifelse(individual_count==0,0,1), ba_n0=ifelse(is.na(BA) | BA==0,0,1)) %>%
  ungroup()


weighted_moments_func<-function(Data){
  weighted_traits<-eval(parse(text=Data))
  grouping<-if(word(Data, 1, sep="_")=="local"){c("Plot","Trait")} else {c("Trait")}

  weighted_means<-weighted_traits %>%
    group_by_(.dots=grouping) %>%
    summarise(sum_w_stem=sum(individual_count,na.rm=T),sum_w_ba=sum(BA,na.rm=T),sum_wx_stem=sum(abund_wx,na.rm=T),
              sum_wx_ba=sum(ba_wx,na.rm=T),n_prime_stem=sum(abund_n0,na.rm=T),n_prime_ba=sum(ba_n0,na.rm=T)) %>%
    rowwise() %>%
    mutate(abund_weighted_mean=sum_wx_stem/sum_w_stem,ba_weighted_mean=sum_wx_ba/sum_w_ba) %>%
    ungroup()

  weighted_sds<-weighted_traits %>%
    left_join(weighted_means) %>%
    rowwise() %>%
    mutate(SDnum_stem=individual_count*(Trait_Mean-abund_weighted_mean)^2,SDnum_ba=BA*(Trait_Mean-ba_weighted_mean)^2) %>%
    ungroup() %>%
    group_by_(.dots=grouping) %>%
    summarize(sum_wxx_stem=sum(SDnum_stem,na.rm=T),sum_wxx_ba=sum(SDnum_ba,na.rm=T)) %>%
    ungroup() %>%
    left_join(weighted_means) %>%
    rowwise() %>%
    mutate(abund_weighted_sd=sqrt(sum_wxx_stem/sum_w_stem),ba_weighted_sd=sqrt(sum_wxx_ba/sum_w_ba),abund_weighted_variance=sum_wxx_stem/sum_w_stem,ba_weighted_variance=sum_wxx_ba/sum_w_ba) %>%
    ungroup()

  weighted_moments<-weighted_traits %>%
    left_join(weighted_sds) %>%
    rowwise() %>%
    mutate(SKEWnum_stem=individual_count*((Trait_Mean-abund_weighted_mean)/abund_weighted_sd)^3,SKEWnum_ba=BA*((Trait_Mean-ba_weighted_mean)/ba_weighted_sd)^3,
           KURTnum_stem=individual_count*((Trait_Mean-abund_weighted_mean)/abund_weighted_sd)^4,KURTnum_ba=BA*((Trait_Mean-ba_weighted_mean)/ba_weighted_sd)^4) %>%
    ungroup() %>%
    group_by_(.dots=grouping) %>%
    summarize(sum_skew_stem=sum(SKEWnum_stem,na.rm=T),sum_skew_ba=sum(SKEWnum_ba,na.rm=T),
              sum_kurt_stem=sum(KURTnum_stem,na.rm=T),sum_kurt_ba=sum(KURTnum_ba,na.rm=T)) %>%
    ungroup() %>%
    left_join(weighted_means) %>%
    rowwise() %>%
    mutate(abund_weighted_skew=sum_skew_stem/sum_w_stem,ba_weighted_skew=sum_skew_ba/sum_w_ba,
           abund_weighted_kurt=(sum_kurt_stem/sum_w_stem)-3,ba_weighted_kurt=(sum_kurt_ba/sum_w_ba)-3) %>%
    ungroup() %>%
    left_join(weighted_sds) %>%
    dplyr::select(-sum_skew_stem, -sum_skew_ba, -sum_kurt_stem, -sum_kurt_ba, -sum_w_stem, -sum_w_ba, -sum_wx_stem, -sum_wx_ba, -sum_wxx_stem, -sum_wxx_ba, -sum_w_stem, -sum_w_ba, -sum_wx_stem, -n_prime_stem, -n_prime_ba) %>%
    mutate(abund_weighted_sd=ifelse(abund_weighted_sd==0,NA,abund_weighted_sd), ba_weighted_sd=ifelse(ba_weighted_sd==0,NA,ba_weighted_sd),
           abund_weighted_variance=ifelse(abund_weighted_variance==0,NA,abund_weighted_variance), ba_weighted_variance=ifelse(ba_weighted_variance==0,NA,ba_weighted_variance),
           abund_weighted_skew=ifelse(abund_weighted_kurt==-3 | abund_weighted_kurt==-2,NA,abund_weighted_skew), abund_weighted_kurt=ifelse(abund_weighted_kurt==-3 | abund_weighted_kurt==-2,NA,abund_weighted_kurt),
           ba_weighted_skew=ifelse(ba_weighted_kurt==-3 | ba_weighted_kurt==-2,NA, ba_weighted_skew), ba_weighted_kurt=ifelse(ba_weighted_kurt==-3 | ba_weighted_kurt==-2,NA, ba_weighted_kurt))

  return(weighted_moments)
}

gdata<-weighted_moments_func("local_weighted_traits") %>%
  left_join(metadata) %>%
  left_join(trait_groups)
gdata$Trait<-factor(gdata$Trait, levels=trait_groups$Trait)


gdata_outRemoved<-gdata[c(1:12)] %>% 
  gather(Moment, value, -c(1:2)) %>%
  drop_na(value) %>%
  group_by(Trait, Moment) %>%
  filter(value>quantile(value,0.005, na.rm=T), value<quantile(value,0.995, na.rm=T)) %>% 
  ungroup() %>%
  spread(Moment, value) %>%
  left_join(metadata) %>%
  left_join(trait_groups)
gdata_outRemoved$Trait<-factor(gdata_outRemoved$Trait, levels=trait_groups$Trait)

#
### trait correlation analysis -----

data_version<-gdata_outRemoved %>%  
  filter(Trait %in% c("SLA", "Wood_Density", "logHeight", "logLeafArea", "logSeedMass", "percent_C", "percent_N", "percent_P", "NP_ratio"))
trait_list_full<-c("SLA", "Wood_Density", "logHeight", "logLeafArea", "logSeedMass", "percent_C", "percent_N", "percent_P", "NP_ratio")
trait_list<-if(dataset=="InSitu"){c("SLA", "Wood_Density", "logLeafArea", "percent_C", "percent_N", "percent_P", "NP_ratio")} else 
  {c("SLA", "Wood_Density", "logHeight", "logLeafArea", "logSeedMass", "percent_C", "percent_N", "percent_P", "NP_ratio")}  

weighting_pre<-if(weighting=="Abundance"){"abund"} else
  if(weighting=="Basal Area"){"ba"} else
    if(weighting=="Species"){"abund"}

M<-data.frame()
p.mat<-data.frame()
for(i in 1:length(trait_list)){
  df<-filter(data_version, Trait==as.character(trait_list[i])) %>%
    dplyr::select(c(paste(weighting_pre, c("weighted_mean", "weighted_variance", "weighted_skew", "weighted_kurt"), sep="_"), env_list)) %>% 
    drop_na()
  colnames(df)[1:4]<-c("Mean", "Variance", "Skewness", "Kurtosis")
  M<-rbind(M, tibble::rownames_to_column(data.frame(cor(df)[,c(1:4)]), var="Variable_1") %>% mutate(Trait=trait_list[i]))
  p.mat<-rbind(p.mat, tibble::rownames_to_column(data.frame(cor_pmat(df)[,c(1:4)]), var="Variable_1") %>% mutate(Trait=trait_list[i]))
}


corr_data<-gather(M, Variable_2, Corr, c(2:5)) %>%
  full_join(gather(p.mat, Variable_2, p_val, c(2:5)))
corr_data<-filter(corr_data, !Variable_1 %in% c("Mean", "Variance", "Skewness", "Kurtosis"))
corr_data$Variable_1<-factor(corr_data$Variable_1, levels=rev(unique(corr_data$Variable_1)))
corr_data$Variable_2<-factor(corr_data$Variable_2, levels=unique(corr_data$Variable_2))
corr_data$Trait<-factor(corr_data$Trait, levels=trait_list_full)

c.val<-max(abs(c(min(corr_data$Corr),max(corr_data$Corr))))
c.range<-c(-c.val,c.val)


# max correlations
max_corrs<-corr_data %>%
  filter(!Variable_1 %in% c("Species_Richness")) %>%
  group_by(Trait, Variable_2) %>%
  mutate(abs.corr=abs(Corr), max.corr=max(abs.corr), test=abs.corr==max.corr) %>%
  filter(test==TRUE) %>%
  rowwise %>%
  mutate(sig=ifelse(p_val<=0.05, "*", "ns"), p_val=ifelse(p_val<0.001, "< 0.001", as.character(round(p_val, 3)))) %>%
  arrange(Trait, Variable_2)
max_corrs$Variable_1<-factor(max_corrs$Variable_1, levels=rev(unique(env_groups$Variable_1)))


#
#####
##### GRAPHICS -----
#####
### (Figure 1B-E) Community-weighted trait-moment versus Absolute Latitude and Elevation regressions -----

trait_gdata1<-data_version[c("Plot", "Trait", "Abs_Latitude", "Elevation", "abund_weighted_mean", "abund_weighted_variance", "abund_weighted_skew", "abund_weighted_kurt")] %>%
  rename(Mean=abund_weighted_mean, Variance=abund_weighted_variance, Skewness=abund_weighted_skew, Kurtosis=abund_weighted_kurt) %>%
  gather(Moment, moment_val, 5:8)

moment_means<-trait_gdata1 %>%
  group_by(Trait, Moment) %>%
  summarize(moment_mean=mean(moment_val, na.rm=T))

norm_trait_gdata1<-left_join(trait_gdata1, moment_means) %>%
  group_by(Trait, Moment) %>%
  mutate(deviation=moment_val-moment_mean, percent_deviation=(deviation/abs(moment_mean))*100) %>%
  ungroup %>%
  gather(Gradient, gradient_val, c(Abs_Latitude, Elevation))
norm_trait_gdata1$Trait<-factor(norm_trait_gdata1$Trait, levels=trait_list)
norm_trait_gdata1$Moment<-factor(norm_trait_gdata1$Moment, levels=c("Mean", "Variance", "Skewness", "Kurtosis"))

stats_lm1<-norm_trait_gdata1 %>%
  group_by(Trait, Moment, Gradient) %>% 
  do(fit = lm(percent_deviation ~ gradient_val, data = ., na.action=na.omit)) %>%
  mutate(r2=round(summary(fit)$r.squared,5), p_val=round(pf(summary(fit)$fstatistic[1], summary(fit)$fstatistic[2], summary(fit)$fstatistic[3], lower=FALSE)[[1]],5), b=data.frame(fit[1])[1,1], m=data.frame(fit[1])[2,1]) %>%
  dplyr::select(-fit) %>% 
  data.frame()
stats_lm1$Trait<-factor(stats_lm1$Trait, levels=rev(trait_list))
stats_lm1$Moment<-factor(stats_lm1$Moment, levels=c("Mean", "Variance", "Skewness", "Kurtosis"))


stats_lm2<-norm_trait_gdata1 %>%
  group_by(Trait, Moment, Gradient) %>%
  do(fit = lm(moment_val ~ gradient_val, data = ., na.action=na.omit)) %>%
  mutate(r2=round(summary(fit)$r.squared,5), p_val=round(pf(summary(fit)$fstatistic[1], summary(fit)$fstatistic[2], summary(fit)$fstatistic[3], lower=FALSE)[[1]],5), b=data.frame(fit[1])[1,1], m=data.frame(fit[1])[2,1]) %>%
  dplyr::select(-fit) %>% 
  data.frame()
stats_lm2$Trait<-factor(stats_lm2$Trait, levels=rev(trait_list))
stats_lm2$Moment<-factor(stats_lm2$Moment, levels=c("Mean", "Variance", "Skewness", "Kurtosis"))


facet_order<-data.frame(Moment=c("Mean", "Mean", "Variance", "Variance", "Skewness", "Skewness", "Kurtosis", "Kurtosis"), Gradient=rep.int(c("Abs_Latitude", "Elevation"), 4),
                        min_y=c(-120, -25, -105, -80, -1.5, -0.8, -1, 0), 
                        max_y=c(40, 5, 35, 40, 1.5, 1.2, 5, 4), 
                        int_y=c(40, 5, 35, 40, 0.5, 0.4, 1, 1)) %>%
  mutate(Facet_order=row_number())


plot_ranges<-trait_gdata1[-c(1,6)] %>%
  gather(Gradient, gradient_value, 2:3) %>%
  group_by(Moment, Gradient) %>%
  summarize(x_min=0, x_max=max(gradient_value, na.rm=T))


colors<-c("red", "darkorange1", "gold3", "chartreuse3", "forestgreen", "darkturquoise", "royalblue3", "darkorchid3", "maroon1")
even_test<-function(n){n-2*floor(n/2)==0}


## For Figure 1B&C:

g_data<-left_join(stats_lm1, plot_ranges) %>% 
  mutate(y_min=m*x_min+b, y_max=m*x_max+b) %>%
  left_join(facet_order)
g_data$Trait<-factor(g_data$Trait, levels=trait_list)
g_data$Moment<-factor(g_data$Moment, levels=c("Mean", "Variance", "Skewness", "Kurtosis"))
g_data<-arrange(g_data, desc(Trait))


for(facet in 1:4){
  temp_g_data<-filter(g_data, Facet_order==facet)
  temp_lim_data<-filter(facet_order, Facet_order==facet)
  x_break_max<-ifelse(even_test(facet), 3500, 60)
  plot<-
    ggplot()+
    geom_segment(data=temp_g_data, aes(x=x_min, y=y_min, xend=x_max, yend=y_max, color=Trait), linetype="solid", size=1.25, lineend="round", alpha=1) +
    geom_segment(aes(x=0, y=0, xend=ifelse(even_test(facet), 3500, 60), yend=0), color="black", linetype="dotted")+
    scale_color_manual(values=colors)+
    scale_x_continuous(breaks=c(0, x_break_max/2, x_break_max))+ 
    scale_y_continuous(breaks=round(seq(temp_lim_data$min_y, temp_lim_data$max_y, temp_lim_data$int_y), 2),
                       limits=c(min(temp_g_data$y_min, temp_g_data$y_max, temp_lim_data$min_y),
                                max(temp_g_data$y_min, temp_g_data$y_max, temp_lim_data$max_y))) + 
    coord_capped_cart(bottom='both', left='both') +
    theme(panel.background = element_blank(), 
          panel.border = element_blank(),
          axis.line = element_line(),
          axis.title=element_blank(),
          axis.text=element_text(size=10, color="black"),
          #axis.text=element_blank(),
          axis.ticks=element_line(color="black"),
          panel.spacing.y=unit(5,"mm"),
          strip.background=element_blank(),
          strip.text=element_blank(),
          aspect.ratio=1, legend.position="none")
  
  
  #svg(file=paste0(directory, "moments_latElev_rawMoment_", facet, ".svg"), width=1.6, height=1.6)
  print(plot)
  #dev.off()
}





## For Figure 1D&E:

g_data<-left_join(stats_lm2, plot_ranges) %>% 
  mutate(y_min=m*x_min+b, y_max=m*x_max+b) %>%
  left_join(facet_order)
g_data$Trait<-factor(g_data$Trait, levels=trait_list)
g_data$Moment<-factor(g_data$Moment, levels=c("Mean", "Variance", "Skewness", "Kurtosis"))
g_data<-arrange(g_data, desc(Trait))


for(facet in 5:8){
  temp_g_data<-filter(g_data, Facet_order==facet)
  temp_lim_data<-filter(facet_order, Facet_order==facet)
  x_break_max<-ifelse(even_test(facet), 3500, 60)
  plot<-
    ggplot()+
    geom_segment(data=temp_g_data, aes(x=x_min, y=y_min, xend=x_max, yend=y_max, color=Trait), linetype="solid", size=1.25, lineend="round", alpha=1) +
    geom_segment(aes(x=0, y=0, xend=ifelse(even_test(facet), 3500, 60), yend=0), color="black", linetype="dotted")+
    scale_color_manual(values=colors)+
    scale_x_continuous(breaks=c(0, x_break_max/2, x_break_max))+
    scale_y_continuous(breaks=round(seq(temp_lim_data$min_y, temp_lim_data$max_y, temp_lim_data$int_y), 2),
                       limits=c(min(temp_g_data$y_min, temp_g_data$y_max, temp_lim_data$min_y),
                                max(temp_g_data$y_min, temp_g_data$y_max, temp_lim_data$max_y))) + 
    coord_capped_cart(bottom='both', left='both') +
    theme(panel.background = element_blank(), 
          panel.border = element_blank(),
          axis.line = element_line(),
          axis.title=element_blank(),
          axis.text=element_text(size=10, color="black"),
          #axis.text=element_blank(),
          axis.ticks=element_line(color="black"),
          panel.spacing.y=unit(5,"mm"),
          strip.background=element_blank(),
          strip.text=element_blank(),
          aspect.ratio=1, legend.position="none") 
  
  
  #svg(file=paste0(directory, "moments_latElev_rawMoment_", facet,".svg"), width=1.6, height=1.6)
  print(plot)
  #dev.off()
}


#
### (Figure 2) Correlations between individual trait-moments and environmental variables -----

for(i in 1:9){
  trait_list<-c("SLA", "Wood_Density", "logHeight", "logLeafArea", "logSeedMass", "percent_C", "percent_N", "percent_P", "NP_ratio")
  colors<-c("red", "darkorange1", "gold3", "chartreuse3", "forestgreen", "darkturquoise", "royalblue3", "darkorchid3", "maroon1")
  trait<-trait_list[i] 
  color_i<-colors[i]
  
  radar.data<-corr_data %>%
    filter(Trait==trait, !Variable_1 %in% c("Species_Richness")) %>% 
    group_by(Variable_2) %>%
    mutate(x=row_number()) %>%
    ungroup
  link.data<-data.frame(Variable_1=c(rep.int("scale", 4), rep.int(NA, 4)), Trait=trait, Variable_2=rep.int(c("Mean", "Variance", "Skewness", "Kurtosis"), 2), Corr=NA, p_val=0, x=17)
  
  radar.data<-rbind(radar.data, link.data)
  radar.data$Variable_2<-factor(radar.data$Variable_2, levels=c("Mean", "Variance", "Skewness", "Kurtosis")) 
  radar.data$Variable_1<-factor(radar.data$Variable_1, levels=c(as.character(env_groups$Variable_1), "scale"))
  radar.max.data<-max_corrs %>%
    filter(Trait==trait, !Variable_1 %in% c("Species_Richness"))
  radar.max.data$Trait<-factor(radar.max.data$Trait, levels=trait_list)
  radar.data$Variable_2<-factor(radar.data$Variable_2, levels=c("Mean", "Variance", "Skewness", "Kurtosis")) 
  radar.data.ns<-filter(radar.max.data, sig=="ns")
  abs.max.corr<-max(abs(min(radar.max.data$Corr)), abs(max(radar.max.data$Corr)))
  y.max<-min(round(abs.max.corr+abs.max.corr*.3, 1), 1)
  y.min<--y.max
  
  g.env.groups<-filter(env_groups, !Variable_1=="Species_Richness")
  g.env.list<-g.env.groups$Text
  
  radial.axis.text<-data.frame(text=c(as.character(g.env.list), NA), y=y.max+y.max*.5) %>%
    mutate(x=row_number())
  nonsig_data<-filter(radar.data, p_val>.05) %>%
    group_by(Variable_1) %>%
    summarize(min=min(0, Corr), max=max(0, Corr))
  nonsig_ribbon<-left_join(data.frame(Variable_1=unique(radar.data$Variable_1)), nonsig_data) %>%
    mutate(min=ifelse(is.na(min), 0, min), max=ifelse(is.na(max), 0, max)) %>%
    left_join(unique(radar.data[c("Variable_1", "x")]))
  
  
  figure<-
    ggplot() + 
    geom_segment(data=radial.axis.text, aes(x, y=y.min, xend=x, yend=y-y*.34), linetype="solid", size=.5, color="gray") +
    scale_x_discrete(expand = c(0,0)) +
    geom_hline(yintercept=c(y.min, y.max), linetype="solid", size=.5, color="gray") +
    geom_segment(data=filter(radial.axis.text, x %in% c(2,7,10,14,16,17)), aes(x, y=y.max, xend=x+1, yend=y.max), linetype="dashed", size=.5, color="white") +
    geom_segment(data=filter(radial.axis.text, !x %in% c(2,7,10,14,16,17)), aes(x, y=y.max+y.max*1.1, xend=x+1, yend=y.max+y.max*1.1), linetype="solid", size=.7, color="black") +
    geom_point(data=radar.data, aes(Variable_1, y=Corr, group=Variable_2, shape=Variable_2), color=color_i, size=3.5, stroke=1, alpha=0) +
    geom_line(data=filter(radar.data, Variable_2=="Mean"), aes(Variable_1, y=Corr, group=Variable_2), color=color_i, size=1, alpha=1) +
    geom_point(data=filter(radar.data, p_val>0.05), aes(Variable_1, y=Corr, group=Variable_2, shape=Variable_2), color="gray", size=4, stroke=.8, alpha=1) +
    geom_point(data=filter(radar.data, p_val<=0.05), aes(Variable_1, y=Corr, group=Variable_2, shape=Variable_2), color=color_i, size=4, stroke=.8, alpha=1) + 
    geom_point(data=radar.max.data, aes(x=Variable_1, y=Corr, shape=Variable_2), fill=color_i, size=4.5, stroke=1.3, alpha=1) +
    geom_ribbon(data=nonsig_ribbon, aes(x, ymin=min, ymax=max), fill = "gray60", alpha=.7) + 
    geom_hline(yintercept=0, linetype="solid", size=.6, color="black") +
    scale_linetype_manual(values=c("solid", "longdash", "dashed", "dotted"))+
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    scale_shape_manual(values=21:24)+
    scale_alpha_continuous(range = c(0.25,1))+
    scale_size_continuous(range = c(2,6))+
    geom_text(data=radial.axis.text, aes(x, y, label=text), size=7.5) +
    scale_y_continuous(limits=c(y.min+y.min*.5, y.max+y.max*1.1), breaks=c(y.min, y.min/2, 0, y.max/2, y.max)) + 
    coord_polar(start=pi*3/2+(2*pi/17), direction=-1) +
    theme(panel.background = element_blank(), 
          panel.border = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          plot.margin=unit(c(0,0,0,0),"mm"),
          plot.title=element_blank(), legend.position="none")
  
  
  #svg(file=paste0(directory, "corrRadar_", trait,".svg"),width=5,height=5)
  print(figure)
  #dev.off()
  
  print(trait)
  print(y.max)
}




# Average magnitude of trait-moment correlations for each environmental variable

avg_corr_data<-corr_data %>%
  filter(p_val<=0.05, !Variable_1=="Species_Richness") %>%
  mutate(abs_Corr=abs(Corr)) %>%
  group_by(Variable_1) %>%
  summarize(avg_Corr=mean(abs_Corr)) %>%
  arrange(desc(avg_Corr)) %>%
  ungroup %>%
  mutate(order=n()+1-row_number()) %>%
  left_join(env_groups)

avg_corr_data$Env_grp<-factor(avg_corr_data$Env_grp, levels=c("Geography", "Temperature", "Precipitation", "Vapor Pressure", "Other"))

rgb_func<-function(r,g,b){rgb(r,g,b, maxColorValue=255)}

plot<-
  ggplot() +
  geom_bar(data=avg_corr_data, aes(order, avg_Corr, fill=Env_grp), stat="identity", alpha=0.5) +
  geom_text(data=avg_corr_data, aes(order, avg_Corr, label=Text), size=3.25, hjust=1, vjust=0.5, nudge_y=-0.01, color="black") + 
  scale_fill_manual(values=c(rgb_func(77,175,74), rgb_func(228,26,28), rgb_func(55,126,184), rgb_func(152,78,163), rgb_func(255,127,0))) +
  labs(x="Climate variables", y="Mean magnitude of\ncorrelation across all\ntrait moments") +
  guides(fill=guide_legend(title="Climate group")) +
  coord_flip()+
  theme(panel.background = element_blank(), panel.border = element_blank(),
        strip.background = element_blank(), strip.text = element_blank(),
        axis.text.y=element_blank(), axis.ticks.y = element_blank(), 
        #axis.title.x=element_blank(),
        #legend.position="none",
        axis.text.x=element_text(size=14))


#svg(file=paste0(directory, "avg_Corrs.svg"), width=6, height=3.5)
plot
#dev.off()







### (Figure 3 & S4) Principal component analysis -----

## ATTENTION: Select 'dataset<-"Bien"' above and rerun 'PRIMARY DATA ANALYSIS'

## ATTENTION: Select trait list:

# for 9-trait version:
all_traits<-data.frame(Trait=c("SLA", "Wood_Density", "logHeight", "logLeafArea", "logSeedMass", "percent_C", "percent_N", "percent_P", "NP_ratio"),
                       Text=c("SLA", "Wood Density", "Log Height", "Log Leaf Area", "Log Seed Mass", "%C", "%N", "%P", "N:P")) 
trait_list<-all_traits

# for 6-trait version:
# trait_list<-data.frame(Trait=c("SLA", "Wood_Density", "logHeight", "logLeafArea", "logSeedMass", "percent_C", "percent_N", "percent_P", "NP_ratio"),
#                        Text=c("SLA", "Wood Density", "Log Height", "Log Leaf Area", "Log Seed Mass", "%C", "%N", "%P", "N:P")) 
# trait_list<-filter(all_traits, !Trait %in% c("percent_C", "percent_P", "NP_ratio"))



# PCA analysis
env_groups<-data.frame(Variable_1=c("Abs_Latitude", "Elevation", "Mean_Annual_Temp", "Isothermality", "Temp_Diurnal_Range", "Temp_Annual_Range", "Temp_Seasonality",
                                    "Annual_Precip", "Precip_Seasonality", "CWD", "vap_yr", "vpd_yr", "pet_yr", "PET_MAP_ratio", "wind_yr", "srad_yr"), 
                       Env_grp=c("Geography", "Geography", rep.int("Temperature", 5), rep.int("Precipitation", 3), rep.int("Vapor Pressure", 4), rep.int("Other", 2)),
                       Text=c("|Lat|", "Elev", "MAT", "ISO", "TDR", "TAR", "TS", "MAP", "PS", "CWD", "VAP", "VPD", "PET", "PET:MAP", "Wind", "SR"))
env_groups$Env_grp<-factor(env_groups$Env_grp, levels=c("Geography", "Temperature", "Precipitation", "Vapor Pressure", "Other"))
env_list<-env_groups$Variable_1 %>%
  as.character()

all_trait_vals<-trait_melt_outRemoved[c(1:4,8,9)] %>% 
  distinct(.keep_all=T) %>%
  group_by(Family, Genus, Species, Trait) %>% # Plot, 
  summarize(Trait_Mean=mean(Trait_Mean, na.rm=T)) %>%
  rename(abund_weighted_mean=Trait_Mean) 


pca_data<-filter(gdata_outRemoved[c("Plot", "Trait", "abund_weighted_mean")], Trait %in% trait_list$Trait)

pca_df<-pca_data %>%
  spread(Trait, abund_weighted_mean) %>%
  drop_na



## ATTENTION: Insert appropriate PCA equation in 'prcomp' function below:
# for 9-trait version: SLA+Wood_Density+Height+logLeafArea+logSeedMass+percent_C+percent_N+percent_P+NP_ratio
# for 6-trait version: SLA+Wood_Density+Height+logLeafArea+logSeedMass+percent_N

pca<-prcomp(~SLA+Wood_Density+logHeight+logLeafArea+logSeedMass+percent_C+percent_N+percent_P+NP_ratio, data=pca_df, center=T, scale=T) 
summary(pca)
coords<-data.frame(pca$x[,1:2]) %>%
  mutate(PC1=rescale(PC1, to=c(-1,1)), PC2=rescale(PC2, to=c(-1,1)))
stadev<-pca$sdev
loadings<-data.frame(pca$rotation)[1:2] %>%
  mutate(PC1=rescale(PC1, to=c(-1,1)), PC2=rescale(PC2, to=c(-1,1)), variable=rownames(pca$rotation), Text=trait_list$Text)
loadings$variable<-factor(loadings$variable, levels=trait_list$Trait)

pc_latLong<-cbind(pca_df[1], coords) %>%
  left_join(dplyr::select(metadata, c("Plot", env_list))) 


mutate(data.frame(pca$rotation)[1:2], variable=rownames(pca$rotation), Text=trait_list$Text) %>%
  arrange(desc(abs(PC2)))

pc1_var<-round(summary(pca)$importance[2,1]*100, 1)
pc2_var<-round(summary(pca)$importance[2,2]*100, 1)
pc1_var+pc2_var


#colors<-c("red", "darkorange1", "gold3", "chartreuse3", "forestgreen", "royalblue3")
colors<-c("red", "darkorange1", "gold3", "chartreuse3", "forestgreen", "darkturquoise", "royalblue3", "darkorchid3", "maroon1")

pc_plot<-
  ggplot() +
  geom_hline(yintercept=0, linetype="dashed", size=.4, color="black") +
  geom_vline(xintercept=0, linetype="dashed", size=.4, color="black") +
  geom_point(data=pc_latLong, aes(PC1, PC2), shape=20, color="gray50", size=.5, alpha=.8)+
  geom_segment(data=loadings, aes(x=0, y=0, xend=PC1, yend=PC2, color=variable), arrow = arrow(length = unit(0.03, "npc"), type="closed", angle=20), size=1, alpha=1) + 
  scale_color_manual(values=colors)+
  labs(x=paste0("Trait PC1 (",pc1_var,"%)"), y=paste0("Trait PC2 (",pc2_var,"%)")) +
  xlim(-1.1, 1.1)+
  ylim(-1.1, 1.1)+
  theme(panel.background = element_blank(), panel.border = element_rect(color = "black", fill=NA),
        axis.ticks=element_line(color="black"),
        axis.text=element_text(size=14, color="black"), axis.title=element_text(size=18, color="black"), legend.position="none")


#svg(file=paste0(directory, "pc_plot_mean.svg"), width=5.25, height=5.25)
pc_plot
#dev.off()



# correlations b/w PC1, PC2, and env vars
pc_env_data<-cbind(pca_df[1], coords) %>%
  left_join(metadata[c("Plot", as.character(env_groups[[1]]))])

df<-pc_env_data[-1] %>%
  drop_na()
M<-tibble::rownames_to_column(data.frame(cor(df)[,c(1:2)]), var="Variable_1")
p.mat<-tibble::rownames_to_column(data.frame(cor_pmat(df)[,c(1:2)]), var="Variable_1")
corr_data<-gather(M, Variable_2, Corr, c(2:3)) %>%
  full_join(gather(p.mat, Variable_2, p_val, c(2:3)))
corr_data<-filter(corr_data, !Variable_1 %in% c("PC1", "PC2"))


corr_gdata<-corr_data %>%
  group_by(Variable_2) %>%
  arrange(desc(abs(Corr))) %>%
  mutate(order=n()+1-row_number()) %>%
  ungroup %>%
  left_join(env_groups)
corr_gdata$Env_grp<-factor(corr_gdata$Env_grp, levels=c("Geography", "Temperature", "Precipitation", "Vapor Pressure", "Other"))



rgb_func<-function(r,g,b){rgb(r,g,b, maxColorValue=255)}


plot<-
  ggplot() +
  geom_bar(data=corr_gdata, aes(order, abs(Corr), fill=Env_grp), stat="identity", alpha=0.5) +
  geom_text(data=corr_gdata, aes(order, abs(Corr), label=Text), size=3.25, hjust=1, vjust=0.5, nudge_y=-0.01, color="black") + 
  scale_fill_manual(values=c(rgb_func(77,175,74), rgb_func(228,26,28), rgb_func(55,126,184), rgb_func(152,78,163), rgb_func(255,127,0))) +
  labs(y="Climatic variable", x="Magnitude of correlation") +
  guides(fill=guide_legend(title="Climate group")) +
  coord_flip()+
  facet_wrap(~Variable_2, scales="free")+
  theme(aspect.ratio=1, plot.margin=unit(c(1,1,1,1), "mm"), panel.background = element_blank(), panel.border = element_blank(),
        strip.background = element_blank(), strip.text = element_blank(),
        #axis.text.y=element_blank(), axis.ticks.y = element_blank(), axis.title.y=element_blank(),
        #axis.title.x=element_blank(),
        #legend.position="none",
        axis.text.x=element_text(size=14))


#svg(file=paste0(directory, "pcCorrs_mean.svg"), width=6, height=3.5)
plot
#dev.off()



### (Figure 4) Relative moment climate change analysis ----- 

## ATTENTION: Select 'dataset<-"InSitu"' above and rerun 'PRIMARY DATA ANALYSIS'

## ATTENTION: Insert moment to analyze from following options:  "mean"   "skew"  
g_moment<-"mean"


# Relative moment analysis

plot_list_1<-trait_melt_outRemoved %>%
  slice(rep(1:n(), .$individual_count)) %>%
  group_by(Plot, Trait) %>%
  summarize(n_test=n()) %>%
  group_by(Plot) %>%
  summarize(n_test_2=length(unique(n_test))) %>%
  filter(n_test_2==1) %>%
  left_join(metadata[c("Plot", "Dataset")])

test<-filter(trait_melt_outRemoved, Plot %in% plot_list_1$Plot) %>%
  slice(rep(1:n(), .$individual_count)) %>%
  group_by(Plot, Trait) %>%
  summarize(n_test=n())

plot_list_mean<-filter(trait_melt_outRemoved, Plot %in% plot_list_1$Plot) %>%
  group_by(Plot, Trait) %>%
  summarize(n_test=n()) %>%
  spread(Trait, n_test) %>%
  drop_na()

plot_list_skew<-filter(plot_list_mean, SLA>=3)


plot_list<-eval(parse(text=paste0("plot_list_", g_moment)))$Plot

data_version<-gdata_outRemoved

trait_list<-unique(data_version$Trait)

M<-data.frame()
p.mat<-data.frame()
for(i in 1:length(trait_list)){
  df<-filter(gdata_outRemoved, Trait==as.character(trait_list[i])) %>%
    dplyr::select(c("abund_weighted_mean", "abund_weighted_variance", "abund_weighted_skew", "abund_weighted_kurt", env_list)) %>% 
    drop_na() %>%
    rename(Mean=abund_weighted_mean, Variance=abund_weighted_variance, Skewness=abund_weighted_skew, Kurtosis=abund_weighted_kurt)
  M<-rbind(M, tibble::rownames_to_column(data.frame(cor(df)[,c(1:4)]), var="Variable_1") %>% mutate(Trait=trait_list[i]))
  p.mat<-rbind(p.mat, tibble::rownames_to_column(data.frame(cor_pmat(df)[,c(1:4)]), var="Variable_1") %>% mutate(Trait=trait_list[i]))
}


corr_data<-gather(M, Variable_2, Corr, c(2:5)) %>%
  full_join(gather(p.mat, Variable_2, p_val, c(2:5)))
corr_data$Variable_1<-factor(corr_data$Variable_1, levels=rev(unique(corr_data$Variable_1)))
corr_data$Variable_2<-factor(corr_data$Variable_2, levels=unique(corr_data$Variable_2))

corr_data<-filter(corr_data, !Variable_1 %in% c("Mean", "Variance", "Skewness", "Kurtosis"))
stats_corr<-filter(corr_data, Variable_1=="Mean_Annual_Temp", Variable_2=="Mean") %>%
  mutate(corr_sign=ifelse(p_val<0.05, sign(Corr), "ns"))


global_gdata<-weighted_moments_func("global_weighted_traits") %>%
  filter(!Trait %in% c("d13C", "d15N", "Height", "logSeedMass"))
colnames(global_gdata)[-1]<-paste("global", colnames(global_gdata[-1]), sep="_")


relative_moments<-data_version[c("Plot", "Trait", "abund_weighted_mean", "abund_weighted_skew")] %>%
  left_join(stats_corr) %>%
  left_join(global_gdata) %>%
  mutate(abund_relative_mean=abund_weighted_mean-global_abund_weighted_mean, abund_relative_skew=abund_weighted_skew-global_abund_weighted_skew,
         percent_deviation_mean=(abund_relative_mean/abs(global_abund_weighted_mean))*100, percent_deviation_skew=(abund_relative_skew/abs(global_abund_weighted_skew))*100) %>%
  dplyr::select(Plot, Trait, Corr, p_val, corr_sign, abund_relative_mean, abund_relative_skew, percent_deviation_mean, percent_deviation_skew) %>%
  gather(relative_moment, relative_moment_value, c(6:7)) %>%
  gather(percent_deviation, percent_deviation_value, c(6:7)) %>%
  drop_na()


relative_moment_ttest<-relative_moments %>%
  group_by(Trait, Corr, p_val, corr_sign, relative_moment, percent_deviation) %>%
  summarize(average_relative_moment=mean(relative_moment_value, na.rm=T), relative_moment_p_val=round(t.test(relative_moment_value)$p.value,5),
            average_percent_deviation=mean(percent_deviation_value, na.rm=T), percent_deviation_p_val=round(t.test(percent_deviation_value)$p.value,5)) %>%
  ungroup()

trait_text<-data.frame(Trait=c("SLA", "Wood_Density", "Height", "logLeafArea", "logSeedMass", "percent_C", "percent_N", "percent_P", "CN_ratio", "CP_ratio", "NP_ratio"),
                       Text=c("SLA", "Wood\nDensity", "Height", "Log\nLeaf Area", "Log\nSeed Mass", "%C", "%N", "%P", "C:N", "C:P", "N:P"))

trait_order<-filter(relative_moment_ttest, relative_moment=="abund_relative_mean", percent_deviation=="percent_deviation_mean") %>% 
  arrange(desc(average_percent_deviation)) %>%
  left_join(trait_text)

g_data_main<-filter(relative_moments, relative_moment==paste0("abund_relative_",g_moment), percent_deviation==paste0("percent_deviation_",g_moment)) %>%
  left_join(metadata[c("Plot", "t_anom_slope")]) %>%
  left_join(dplyr::select(metadata, Plot, Mean_Annual_Temp))


g_data_main$Trait<-factor(g_data_main$Trait, levels=trait_order$Trait)
sig_data<-filter(relative_moment_ttest, relative_moment==paste0("abund_relative_",g_moment), percent_deviation==paste0("percent_deviation_",g_moment), relative_moment_p_val<=0.05)
ns_data<-filter(relative_moment_ttest, relative_moment==paste0("abund_relative_",g_moment), percent_deviation==paste0("percent_deviation_",g_moment), relative_moment_p_val>0.05)

x_max_abs<-g_data_main$percent_deviation_value %>% abs %>% max %>% ceiling
x_lim<-c(-x_max_abs, x_max_abs)


lm_func<-function(T_sign){
  data<-filter(g_data_main, corr_sign==T_sign)
  fit<-lm(percent_deviation_value~t_anom_slope, data=data)
  reg_summary<-data.frame(corr_sign=T_sign,
                          intercept=data.frame(fit[1])[1,1], 
                          slope=data.frame(fit[1])[2,1], 
                          r2=round(summary(fit)$r.squared,5), 
                          p_val=round(pf(summary(fit)$fstatistic[1], summary(fit)$fstatistic[2], summary(fit)$fstatistic[3], lower=FALSE)[[1]],5))
  return(reg_summary)
}

reg_stats<-rbind(lm_func(-1), lm_func(1))
reg_stats



# relative moment (% deviation)
rm_plot<-
  ggplot(g_data_main) + 
  geom_vline(xintercept=0, size=.5, linetype="dotted", color="black") +
  geom_point(aes(x=percent_deviation_value, y=Trait, color=factor(corr_sign)), size=2.2, shape = 21, stroke = .4, alpha=.6) + 
  scale_colour_manual(values = c(rgb(241, 151, 3, maxColorValue=255), "dodgerblue1", "gray60"), breaks=c("1","ns","-1"), labels=c("Positive", "Not significant", "Negative")) +
  geom_point(data=sig_data, aes(x=average_percent_deviation, y=Trait), color="black", shape = "*", size=10) +
  geom_point(data=ns_data, aes(x=average_percent_deviation, y=Trait), color="black", size=1, shape = "|", stroke = 8) +
  scale_y_discrete(labels=trait_order$Text) +
  theme(axis.text=element_text(size=10),
        panel.background = element_blank(), panel.border = element_rect(color = "black", fill=NA), 
        axis.title.y=element_blank(), legend.position="none")

#svg(file=paste0(directory, "relative_moments_", g_moment,".svg"),width=3.25,height=2.75) # ,width=3.25,height=2.75
rm_plot
#dev.off()



# rate of warming VS relative moment (% deviation)
cc_plot<-ggplot(filter(g_data_main, corr_sign %in% c(-1, 1))) + 
  geom_point(aes(x=t_anom_slope, y=percent_deviation_value, color=factor(corr_sign)), size=2, shape = 21, stroke = .4, alpha=.6) +
  geom_hline(yintercept=0, size=.5, linetype="dotted", color="black")+
  geom_abline(data=reg_stats, aes(intercept=intercept, slope=slope, color=factor(corr_sign)), size=1) +
  scale_colour_manual(values = c(rgb(241, 151, 3, maxColorValue=255), "dodgerblue1"))+
  scale_x_continuous(limits=c(0, 0.03), breaks=c(0, 0.015, 0.03), labels=c(0.0, 0.015, 0.03)) +
  theme(axis.text=element_text(size=10), 
        #axis.text=element_blank(),
        #axis.title=element_blank(),
        panel.background = element_blank(), panel.border = element_rect(color = "black", fill=NA), 
        strip.background = element_blank(), strip.text.x = element_blank(), legend.position="none")

#svg(file=paste0(directory, g_moment, "_vs_TRate.svg"),width=2.25, height=2.5) # ,width=2.25, height=2.5
cc_plot
#dev.off()


### (Figure S2 & S6B) Full correlation matrix diagram -----

g_data<-left_join(corr_data, env_groups)
g_data$Variable_1<-factor(g_data$Variable_1, levels=rev(env_groups$Variable_1))



full_corr_plot<-
  ggplot(g_data) +
  geom_tile(data=g_data, aes(x=Variable_2, y=Variable_1, fill=Corr), color="black", size=.25) +
  scale_fill_gradient2(low="red", mid="white", high="blue", na.value = "grey80", midpoint=0, limits=c.range, breaks=c(-.6, 0, .6), labels=c("-0.6", "0", "0.6"), name="Correlation") +
  geom_text(data=filter(g_data, p_val<=0.05), aes(x=Variable_2, y=Variable_1, label=round(Corr,2)), color="black", size=2.5) +
  geom_point(data=filter(g_data, p_val>0.05), aes(x=Variable_2, y=Variable_1), shape=4, color="black", size=3) +
  facet_grid(Env_grp~Trait, scales="free", space="free_y", labeller=labeller(Trait=c(SLA="SLA", Wood_Density="Wood Density", logHeight="Log Height", logLeafArea="Log Leaf Area", logSeedMass="Log Seed Mass", 
                                                                                     percent_C="%C", percent_N="%N", percent_P="%P", CN_ratio="C:N", CP_ratio="C:P", NP_ratio="N:P"))) +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = .5, size=12), 
        axis.text.y=element_text(size=12), 
        strip.text.x=element_text(size=10), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        strip.background=element_blank(),
        panel.background = element_blank(), 
        panel.spacing.x=unit(0.5,"mm"),
        panel.spacing.y=unit(1,"mm"),
        panel.border = element_blank(), legend.title=element_text(angle=90)) +
  guides(fill=guide_colorbar(title.position = "left", title.hjust=.5))


#svg(file=paste0(directory, "moment_corr_all.svg"),width=13.5,height=6.5)
full_corr_plot
#dev.off()




### (Figure S3) Community-weighted trait-moment versus Absolute Latitude and Elevation regressions -----

elev_to_lat<-function(x){rescale(x, to = range(gdata$Abs_Latitude), from = range(gdata$Elevation))}

data_version<-gdata_outRemoved %>%
  filter(Trait %in% trait_list)


method.id<-"lm"
span.width=NA

colors<-c("red", "darkorange1", "gold2", "chartreuse2", "forestgreen", "darkturquoise", "royalblue3", "darkorchid3", "maroon1")

gradient1<-"Abs_Latitude"
gradient2<-"Elevation"
weighting1<-"abund_weighted_mean"
weighting2<-"abund_weighted_variance"
weighting3<-"abund_weighted_skew"
weighting4<-"abund_weighted_kurt"


plot.function<-function(weighting) {
  stats_lm1<-data_version %>%
    group_by(Trait) %>%
    do(fit = lm(eval(parse(text=weighting)) ~ eval(parse(text=gradient1)), data = ., na.action=na.omit)) %>%
    mutate(r2=round(summary(fit)$r.squared,5),p_val=round(pf(summary(fit)$fstatistic[1], summary(fit)$fstatistic[2], summary(fit)$fstatistic[3], lower=FALSE)[[1]],5), slope=data.frame(fit[1])[2,1], slope_sign=sign(slope)) %>%
    dplyr::select(-fit) %>% 
    rowwise() %>%
    mutate(sig=ifelse(p_val<=0.05, "*", "ns")) %>%
    data.frame()
  
  stats_lm2<-data_version %>%
    group_by(Trait) %>%
    do(fit = lm(eval(parse(text=weighting)) ~ eval(parse(text=gradient2)), data = ., na.action=na.omit)) %>%
    mutate(r2=round(summary(fit)$r.squared,5),p_val=round(pf(summary(fit)$fstatistic[1], summary(fit)$fstatistic[2], summary(fit)$fstatistic[3], lower=FALSE)[[1]],5), slope=data.frame(fit[1])[2,1], slope_sign=sign(slope)) %>%
    dplyr::select(-fit) %>% 
    rowwise() %>%
    mutate(sig=ifelse(p_val<=0.05, "*", "ns")) %>%
    data.frame()
  
  trait_gdata1<-left_join(data_version[c("Trait", "trait_group", weighting, gradient1)], stats_lm1)
  trait_gdata_sig1<-filter(trait_gdata1, sig=="*")
  trait_gdata_nonsig1<-filter(trait_gdata1, sig=="ns")
  trait_gdata1$Trait<-factor(trait_gdata1$Trait, levels=trait_list)
  
  trait_gdata2<-left_join(data_version[c("Trait", "trait_group", weighting, gradient2)], stats_lm2) %>%
    mutate(lat_scaled_elev=elev_to_lat(Elevation))
  trait_gdata_sig2<-filter(trait_gdata2, sig=="*")
  trait_gdata_nonsig2<-filter(trait_gdata2, sig=="ns")
  trait_gdata2$Trait<-factor(trait_gdata2$Trait, levels=trait_list)
  
  y.label<-if(weighting=="abund_weighted_mean") {"Community-weighted mean (CWM)"
  } else if(weighting=="abund_weighted_variance") {"Community-weighted variance (CWV)"
  } else if(weighting=="abund_weighted_skew") {"Community-weighted skewness (CWS)"
  } else if(weighting=="abund_weighted_kurt") {"Community-weighted kurtosis (CWK)"
  }
  
  plot<-ggplot() +
    geom_point(data=trait_gdata1, aes(x=eval(parse(text=gradient1)), y=eval(parse(text=weighting))), color="darkorchid1", size=1, shape = 21, alpha=.7, stroke=.3) + 
    geom_point(data=trait_gdata2, aes(x=lat_scaled_elev, y=eval(parse(text=weighting))), color=rgb(0,.85,0), size=1, shape = 21, alpha=.7, stroke=.3) + 
    stat_smooth(data=trait_gdata_nonsig1, aes(x=eval(parse(text=gradient1)), y=eval(parse(text=weighting))), fullrange=T, color="gray", se = F, size=1.5, span=span.width, method=method.id) +
    stat_smooth(data=trait_gdata_nonsig2, aes(x=lat_scaled_elev, y=eval(parse(text=weighting))), fullrange=T, color="gray", se = F, size=1.5, span=span.width, method=method.id) + 
    stat_smooth(data=trait_gdata_sig1, aes(x=eval(parse(text=gradient1)), y=eval(parse(text=weighting))), fullrange=T, color="darkorchid3", se = F, size=1.5, span=span.width, method=method.id) + 
    stat_smooth(data=trait_gdata_sig2, aes(x=lat_scaled_elev, y=eval(parse(text=weighting))), fullrange=T, color=rgb(0,.7,0), se = F, size=1.5, span=span.width, method=method.id) +
    labs(x="|Latitude| (purple)", y=y.label) +
    scale_x_continuous(sec.axis=dup_axis(name="Elevation (green)", breaks=elev_to_lat(c(0,1500,3000)), labels=c("0","1500","3000"))) +
    facet_wrap(~Trait, nrow=1, scales="free", strip.position="right", labeller=labeller(Trait=c(SLA="SLA", Wood_Density="Wood Density", logHeight="Log Height", logLeafArea="Log Leaf Area", logSeedMass="Log Seed Mass", 
                                                                                                percent_C="%C", percent_N="%N", percent_P="%P", CN_ratio="C:N", CP_ratio="C:P", NP_ratio="N:P"))) +
    theme(panel.background = element_blank(), panel.border = element_rect(color = "black", fill=NA), panel.spacing=unit(0.2, "cm"), 
          strip.background = element_blank(), strip.text = element_blank(), axis.title=element_blank(),
          aspect.ratio=1, legend.position="none")
  return(plot)
  
}

plot.function(weighting1)

g1<-plot.function(weighting1)
g2<-plot.function(weighting2)
g3<-plot.function(weighting3)
g4<-plot.function(weighting4)

g<-rbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), ggplotGrob(g4), size = "max") 
panels <- g$layout$t[grep("panel", g$layout$name)]
g$lengths[panels] <- unit(c(1,1), "null")

# multiple moments
#svg(file=paste0(directory, "moments_latElev_all.svg"), width=10, height=9)
grid.newpage()
grid.draw(g)
#dev.off()

# single moment 
#svg(file=paste0(directory, "moment_latElev.svg"), width=10, height=2)
g1
#dev.off()

#
### (Figure S5) Species Abundance Distribution Analysis -----

SAD_metrics<-trait_melt_outRemoved %>%
  mutate(species=paste(Genus, Species, sep="_")) %>%
  dplyr::select(Plot, species, individual_count) %>%
  unique %>%
  group_by(Plot) %>%
  mutate(n_spp=length(unique(species))) %>%
  ungroup() %>%
  filter(n_spp>=3) %>%
  group_by(Plot, n_spp) %>%
  summarize(fishers_alpha=fisherfit(individual_count)[[1]], H=diversity(individual_count), J=H/log(specnumber(individual_count))) %>%
  ungroup %>%
  left_join(metadata[c("Plot", "Abs_Latitude", "Elevation")]) %>%
  gather(Fisher_stat, fs_value, 2:5) %>%
  gather(Environment, e_value, 2:3)


stats_lm1<-SAD_metrics %>%
  group_by(Fisher_stat, Environment) %>%
  do(fit = lm(fs_value ~ e_value, data = ., na.action=na.omit)) %>%
  mutate(r2=round(summary(fit)$r.squared,5),p_val=pf(summary(fit)$fstatistic[1], summary(fit)$fstatistic[2], summary(fit)$fstatistic[3], lower=FALSE)[[1]], slope=data.frame(fit[1])[2,1], slope_sign=sign(slope)) %>%
  dplyr::select(-fit) %>% 
  rowwise() %>%
  mutate(sig=ifelse(p_val<=0.05, "*", "ns")) %>%
  data.frame()



SAD_plot<-
  ggplot() +
  geom_point(data=SAD_metrics, aes(e_value, fs_value), size=1, stroke=.4, alpha=0.7, shape=1, color="gray30") +
  geom_smooth(data=SAD_metrics, aes(e_value, fs_value), method="lm") +
  facet_grid(Fisher_stat~Environment, scales="free") +
  theme(panel.background = element_blank(), panel.border = element_rect(color = "black", fill=NA),
        #strip.background = element_blank(), strip.text = element_blank(),
        axis.text=element_text(size=8), axis.title=element_text(size=14), panel.spacing=unit(0.2, "cm"))


#svg(file=paste0(directory, "fisher_stats.svg"),width=4,height=6) # ,width=3,height=5
SAD_plot
#dev.off()



### (Figure S6) Bootstrapping trait-moment estimates -----

# weighted moment functions
Mean<-function(x){
  z<-na.omit(x)
  sum(z)/length(z)
}
SD<-function(x){
  z<-na.omit(x)
  n<-length(z)
  sqrt((sum((z-mean(z))^2))/n)
}
Variance<-function(x){
  z<-na.omit(x)
  n<-length(z)
  (sum((z-mean(z))^2))/n
}
Skewness<-function(x){
  z<-na.omit(x)
  n<-length(z)
  stddev<-SD(z)
  sum(((z-mean(z))/stddev)^3)/n
}
Kurtosis<-function(x){
  z<-na.omit(x)
  n<-length(z)
  stddev<-SD(z)
  (sum(((z-mean(z))/stddev)^4)/n)-3
}


# bootstrapping function
bs_func<-function(stat_func, df, i){
  resample=df[i,]
  return(stat_func(resample$Trait_Mean))
}


weighted_traits_sliced<-local_weighted_traits %>% 
  dplyr::select(Plot, individual_count, Trait, Trait_Mean) %>%
  slice(rep(1:n(), .$individual_count))

grouping<-c("Plot", "Trait")
weighted_moments_df<-weighted_traits_sliced %>%
  group_by(.dots=grouping) %>%
  summarize(awm=Mean(Trait_Mean), awv=SD(Trait_Mean)^2, aws=Skewness(Trait_Mean), awk=Kurtosis(Trait_Mean)) %>%
  mutate(awv=ifelse(awv==0, NA, awv))


nValues<-trait_melt_outRemoved %>%
  group_by(Plot, Trait) %>%
  summarize(n_unique=length(unique(Trait_Mean))) %>%
  ungroup

nValues_GreaterThan2<-filter(nValues, n_unique>2) %>%
  dplyr::select(Plot, Trait) %>%
  unique

plot_list<-unique(nValues_GreaterThan2$Plot)

bs_results<-data.frame()
for(i in 1:length(plot_list)){ 
  plot_i<-plot_list[[i]]  
  trait_list<-filter(nValues_GreaterThan2, Plot==plot_i)$Trait
  plot_data<-filter(weighted_traits_sliced, Plot==plot_i, Trait %in% trait_list)
  
  bs_moment_func<-function(Moment){
    bootstrapping_df<-plot_data %>% 
      group_by_(.dots=grouping) %>%
      do(b=boot(., bs_func, R=1000, stat_func=eval(parse(text=Moment)))) %>% 
      mutate(original=b$t0, sim=mean(b$t), bias=sim-original, lower_ci=boot.ci(b, type="bca")$bca[1,4], upper_ci=boot.ci(b, type="bca")$bca[1,5]) %>%
      dplyr::select(-b) %>%
      ungroup
    colnames(bootstrapping_df)[3:7]<-paste(colnames(bootstrapping_df)[3:7], Moment, sep="_")
    return(bootstrapping_df)
  }
  
  bs_Mean<-bs_moment_func("Mean")
  bs_Variance<-bs_moment_func("Variance")
  bs_Skewness<-bs_moment_func("Skewness")
  bs_Kurtosis<-bs_moment_func("Kurtosis")
  
  bs_moments<-bs_Mean %>%
    left_join(bs_Variance) %>%
    left_join(bs_Skewness) %>%
    left_join(bs_Kurtosis)
  
  bs_results<-rbind(bs_results, bs_moments)
}


## Export bootstrapping results
# write.csv(bs_results, paste0(directory, "bs_results.csv"), row.names=F)




## Import bootstrapping results
bs_results_df<-data.frame(read.csv(paste0(directory, "bs_results.csv")))

trait_list<-data.frame(Trait=c("SLA", "Wood_Density", "logHeight", "logLeafArea", "logSeedMass", "percent_C", "percent_N", "percent_P", "NP_ratio"),
                       Trait_text=c("SLA", "Wood Density", "Height", "Log Leaf Area", "Log Seed Mass", "%C", "%N", "%P", "N:P"))

trait_moment_order<-expand.grid(Trait=trait_list$Trait, Moment=c("Mean", "Variance", "Skewness", "Kurtosis")) %>%
  mutate(Trait_Moment=paste(Trait, Moment, sep="_")) %>%
  left_join(trait_list)

outlier_data<-gdata_outRemoved %>%
  dplyr::select(Plot, Trait, abund_weighted_mean, abund_weighted_variance, abund_weighted_skew, abund_weighted_kurt) %>%
  rename(Mean=abund_weighted_mean, Variance=abund_weighted_variance, Skewness=abund_weighted_skew, Kurtosis=abund_weighted_kurt) %>%
  gather(Moment, is_outlier, 3:6 ) %>%
  mutate(is_outlier=ifelse(is.na(is_outlier), 1, 0))

nValues<-trait_melt_outRemoved %>%
  group_by(Plot, Trait) %>%
  summarize(n_unique=length(unique(Trait_Mean))) %>%
  ungroup

nValues_GreaterThan2<-filter(nValues, n_unique>2) %>%
  dplyr::select(Plot, Trait) %>%
  unique

bs_results<-bs_results_df %>%
  left_join(nValues) %>%
  filter(Trait %in% trait_list$Trait) %>%
  mutate(Mean=(upper_ci_Mean-lower_ci_Mean),
         Variance=(upper_ci_Variance-lower_ci_Variance),
         Skewness=(upper_ci_Skewness-lower_ci_Skewness), 
         Kurtosis=(upper_ci_Kurtosis-lower_ci_Kurtosis)) %>%
  dplyr::select(Plot, Trait, n_unique, Mean, Variance, Skewness, Kurtosis) %>%
  gather(Moment, ci_range, 4:7) %>%
  left_join(outlier_data) %>%
  mutate(Trait_Moment=paste(Trait, Moment, sep="_")) %>%
  left_join(trait_list)
bs_results$Trait_Moment<-factor(bs_results$Trait_Moment, levels=trait_moment_order$Trait_Moment)


bs_results_outliers<-bs_results %>%
  drop_na(ci_range) %>%
  group_by(Trait, Moment) %>%
  filter(ci_range>quantile(ci_range, 0.975, na.rm=T)) %>%
  ungroup() %>%
  mutate(outlier_0.01=1)


bs_gdata<-bs_results %>%
  left_join(bs_results_outliers)


# plot

# labels
trait_moment_labels<-c(SLA_Mean="SLA", Wood_Density_Mean="Wood Density", logHeight_Mean="Log Height", logLeafArea_Mean="Log Leaf Area", logSeedMass_Mean="Log Seed Mass", percent_C_Mean="%C", percent_N_Mean="%N", percent_P_Mean="%P", NP_ratio_Mean="N:P",
                       SLA_Variance="SLA", Wood_Density_Variance="Wood Density", logHeight_Mean="Log Height", logLeafArea_Variance="Log Leaf Area", logSeedMass_Variance="Log Seed Mass", percent_C_Variance="%C", percent_N_Variance="%N", percent_P_Variance="%P", NP_ratio_Variance="N:P",
                       SLA_Skewness="SLA", Wood_Density_Skewness="Wood Density", logHeight_Mean="Log Height", logLeafArea_Skewness="Log Leaf Area", logSeedMass_Skewness="Log Seed Mass", percent_C_Skewness="%C", percent_N_Skewness="%N", percent_P_Skewness="%P", NP_ratio_Skewness="N:P",
                       SLA_Kurtosis="SLA", Wood_Density_Kurtosis="Wood Density", logHeight_Mean="Log Height", logLeafArea_Kurtosis="Log Leaf Area", logSeedMass_Kurtosis="Log Seed Mass", percent_C_Kurtosis="%C", percent_N_Kurtosis="%N", percent_P_Kurtosis="%P", NP_ratio_Kurtosis="N:P")


bs_plot<-
  ggplot() +
  geom_point(data=filter(bs_gdata, is.na(outlier_0.01)), aes(n_unique, ci_range), size=1, stroke=.4, alpha=0.7, shape=1, color="gray30") +
  geom_point(data=filter(bs_gdata, outlier_0.01==1), aes(n_unique, ci_range), size=1, stroke=.4, alpha=1, shape=1, color="red") +
  labs(x="Number of unique trait values", y="Confidence interval range") +
  facet_wrap(~Trait_Moment, scales="free", nrow=4, labeller=labeller(Trait_Moment=trait_moment_labels)) +
  theme(panel.background = element_blank(), panel.border = element_rect(color = "black", fill=NA),
        axis.text=element_text(size=8), axis.title=element_text(size=14), plot.title=element_text(size=14, face="bold", hjust=0.5), 
        #strip.background = element_blank(), strip.text = element_blank(),
        panel.spacing=unit(0.2, "cm"))

#svg(file=paste0(directory, "bootstrapping_summary.svg"), width=12, height=5)
bs_plot
#dev.off()


#
### (Figure S7) Comparing the strength of local (in situ) versus database (Bien) trait-environment correlations -----

# local (in situ plots) analysis:
trait_melt_outRemoved<-filter(trait_data, Dataset=="InSitu")[-1] %>%
  gather(Trait, value, -c(1:7)) %>%
  drop_na(value) %>%
  mutate(value=as.numeric(as.character(value))) %>%
  group_by(Trait) %>% 
  filter(value>quantile(value,0.005, na.rm=T), value<quantile(value,0.995, na.rm=T)) %>%
  ungroup() %>%
  rename(Trait_Mean=value)


local_weighted_traits<-if(weighting=="Species")
{trait_melt_outRemoved %>% mutate(abund_wx=Trait_Mean, ba_wx=Trait_Mean, abund_n0=ifelse(individual_count==0,0,1), ba_n0=ifelse(is.na(BA) | BA==0,0,1))} else
{trait_melt_outRemoved %>% mutate(abund_wx=Trait_Mean*individual_count, ba_wx=Trait_Mean*BA, abund_n0=ifelse(individual_count==0,0,1), ba_n0=ifelse(is.na(BA) | BA==0,0,1))} 


gdata<-weighted_moments_func("local_weighted_traits") %>%
  left_join(metadata) %>%
  left_join(trait_groups) 
gdata$Trait<-factor(gdata$Trait, levels=trait_groups$Trait)


gdata_outRemoved<-gdata[c(1:12)] %>% # c(1:12)
  gather(Moment, value, -c(1:2)) %>%
  drop_na(value) %>%
  group_by(Trait, Moment) %>% 
  filter(value>quantile(value,0.005, na.rm=T), value<quantile(value,0.995, na.rm=T)) %>% 
  ungroup() %>%
  spread(Moment, value) %>%
  left_join(metadata) %>%
  left_join(trait_groups)
gdata_outRemoved$Trait<-factor(gdata_outRemoved$Trait, levels=trait_groups$Trait)


data_version<-gdata_outRemoved %>% 
  filter(Trait %in% c("SLA", "Wood_Density", "logHeight", "logLeafArea", "logSeedMass", "percent_C", "percent_N", "percent_P", "NP_ratio"))
trait_list_local<-c("SLA", "Wood_Density", "logLeafArea", "percent_C", "percent_N", "percent_P", "NP_ratio")


weighting_pre<-if(weighting=="Abundance"){"abund"} else
  if(weighting=="Basal Area"){"ba"} else
    if(weighting=="Species"){"abund"}

M<-data.frame()
p.mat<-data.frame()
for(i in 1:length(trait_list_local)){
  df<-filter(data_version, Trait==as.character(trait_list_local[i])) %>%
    dplyr::select(c(paste(weighting_pre, c("weighted_mean", "weighted_variance", "weighted_skew", "weighted_kurt"), sep="_"), env_list)) %>% 
    drop_na()
  colnames(df)[1:4]<-c("Mean", "Variance", "Skewness", "Kurtosis")
  M<-rbind(M, tibble::rownames_to_column(data.frame(cor(df)[,c(1:4)]), var="Variable_1") %>% mutate(Trait=trait_list_local[i]))
  p.mat<-rbind(p.mat, tibble::rownames_to_column(data.frame(cor_pmat(df)[,c(1:4)]), var="Variable_1") %>% mutate(Trait=trait_list_local[i]))
}


corr_data<-gather(M, Variable_2, Corr, c(2:5)) %>%
  full_join(gather(p.mat, Variable_2, p_val, c(2:5)))
corr_data<-filter(corr_data, !Variable_1 %in% c("Mean", "Variance", "Skewness", "Kurtosis"))
corr_data$Variable_1<-factor(corr_data$Variable_1, levels=rev(unique(corr_data$Variable_1)))
corr_data$Variable_2<-factor(corr_data$Variable_2, levels=unique(corr_data$Variable_2))
corr_data$Trait<-factor(corr_data$Trait, levels=trait_list_full)


g_data<-left_join(corr_data, env_groups)
g_data$Variable_1<-factor(g_data$Variable_1, levels=rev(env_groups$Variable_1))

g_data_local<-g_data %>%
  mutate(dataset="local")



# database (Bien plots) analysis:
trait_melt_outRemoved<-filter(trait_data, Dataset=="Bien")[-1] %>%
  gather(Trait, value, -c(1:7)) %>%
  drop_na(value) %>%
  mutate(value=as.numeric(as.character(value))) %>%
  group_by(Trait) %>% 
  filter(value>quantile(value,0.005, na.rm=T), value<quantile(value,0.995, na.rm=T)) %>%
  ungroup() %>%
  rename(Trait_Mean=value)


local_weighted_traits<-if(weighting=="Species")
{trait_melt_outRemoved %>% mutate(abund_wx=Trait_Mean, ba_wx=Trait_Mean, abund_n0=ifelse(individual_count==0,0,1), ba_n0=ifelse(is.na(BA) | BA==0,0,1))} else
{trait_melt_outRemoved %>% mutate(abund_wx=Trait_Mean*individual_count, ba_wx=Trait_Mean*BA, abund_n0=ifelse(individual_count==0,0,1), ba_n0=ifelse(is.na(BA) | BA==0,0,1))} 


gdata<-weighted_moments_func("local_weighted_traits") %>%
  left_join(metadata) %>%
  left_join(trait_groups)
gdata$Trait<-factor(gdata$Trait, levels=trait_groups$Trait)


gdata_outRemoved<-gdata[c(1:12)] %>%
  gather(Moment, value, -c(1:2)) %>%
  drop_na(value) %>%
  group_by(Trait, Moment) %>% 
  filter(value>quantile(value,0.005, na.rm=T), value<quantile(value,0.995, na.rm=T)) %>% 
  ungroup() %>%
  spread(Moment, value) %>%
  left_join(metadata) %>%
  left_join(trait_groups)
gdata_outRemoved$Trait<-factor(gdata_outRemoved$Trait, levels=trait_groups$Trait)


data_version<-gdata_outRemoved %>% 
  filter(Trait %in% c("SLA", "Wood_Density", "logHeight", "logLeafArea", "logSeedMass", "percent_C", "percent_N", "percent_P", "NP_ratio"))
trait_list_database<-c("SLA", "Wood_Density", "logLeafArea", "percent_C", "percent_N", "percent_P", "NP_ratio")


weighting_pre<-if(weighting=="Abundance"){"abund"} else
  if(weighting=="Basal Area"){"ba"} else
    if(weighting=="Species"){"abund"}

M<-data.frame()
p.mat<-data.frame()
for(i in 1:length(trait_list_database)){
  df<-filter(data_version, Trait==as.character(trait_list_database[i])) %>%
    dplyr::select(c(paste(weighting_pre, c("weighted_mean", "weighted_variance", "weighted_skew", "weighted_kurt"), sep="_"), env_list)) %>% 
    drop_na()
  colnames(df)[1:4]<-c("Mean", "Variance", "Skewness", "Kurtosis")
  M<-rbind(M, tibble::rownames_to_column(data.frame(cor(df)[,c(1:4)]), var="Variable_1") %>% mutate(Trait=trait_list_database[i]))
  p.mat<-rbind(p.mat, tibble::rownames_to_column(data.frame(cor_pmat(df)[,c(1:4)]), var="Variable_1") %>% mutate(Trait=trait_list_database[i]))
}


corr_data<-gather(M, Variable_2, Corr, c(2:5)) %>%
  full_join(gather(p.mat, Variable_2, p_val, c(2:5)))
corr_data<-filter(corr_data, !Variable_1 %in% c("Mean", "Variance", "Skewness", "Kurtosis"))
corr_data$Variable_1<-factor(corr_data$Variable_1, levels=rev(unique(corr_data$Variable_1)))
corr_data$Variable_2<-factor(corr_data$Variable_2, levels=unique(corr_data$Variable_2))
corr_data$Trait<-factor(corr_data$Trait, levels=trait_list_full)


g_data<-left_join(corr_data, env_groups)
g_data$Variable_1<-factor(g_data$Variable_1, levels=rev(env_groups$Variable_1))

g_data_database<-g_data %>%
  mutate(dataset="database")





g_data_merge<-rbind(g_data_local, g_data_database) %>%
  mutate(Corr=abs(Corr))



local_db_compare_plot<-
  ggplot()+
  geom_boxplot(data=filter(g_data_merge, p_val<=0.05), aes(Env_grp, Corr, color=dataset))+
  ylim(0 ,1)+
  theme(panel.background = element_blank(), panel.border = element_rect(color = "black", fill=NA),
        axis.text=element_text(size=12),
        axis.title=element_blank())


#svg(file=paste0(directory, "moment_corr_localVSdatabase.svg"),width=7.5,height=5)
local_db_compare_plot
#dev.off()



#####
##### END -----











