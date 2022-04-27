setwd("~/Documents/Research/Small_jawed_project/popgenome")

library(ggplot2)
library(tidyverse)
library(IRanges)
library(GenomicRanges)

exons.genes<-read.table('../exon_gene_locations.txt',header=F, stringsAsFactors = F)
colnames(exons.genes)<-c('Scaffold','Type','Exon_Start','Exon_Stop')

exons<-exons.genes %>% filter(Type!="gene")


ASvP.intro.results<-read.csv('ASvP_unique_snps_intro_stats_summary.csv',header=T, stringsAsFactors = F)%>%filter(!is.na(Scaffold))
ASvP.intro.results$Intro_Stop<-as.integer(ASvP.intro.results$Intro_Stop)

APvS.intro.results<-read.csv('APvS_unique_snps_intro_stats_summary.csv',header=T, stringsAsFactors = F)%>%filter(!is.na(Scaffold))
AvSP.intro.results<-read.csv('SP_shared_snps_intro_stats_summary.csv',header=T, stringsAsFactors = F)%>%filter(!is.na(Scaffold))


gr2<-with(ASvP.intro.results,GRanges(Scaffold,IRanges(start=Intro_Start, end=Intro_Stop)))
gr1<-with(exons,GRanges(Scaffold,IRanges(start=Exon_Start, end=Exon_Stop)))

type1 = findOverlaps(query = gr1, subject = gr2, type = 'within')
type1.df = data.frame(exons[queryHits(type1),], ASvP.intro.results[subjectHits(type1),])

type2 = findOverlaps(query = gr1, subject = gr2, type = 'any')
type2.df = data.frame(exons[queryHits(type2),], ASvP.intro.results[subjectHits(type2),])


ASvP.intro.results.exons<-type2.df %>% dplyr::select(-c(Scaffold.1,Gene,Gene_Start,Gene_Stop))

ASvP.snps<-read.table("~/Documents/Research/Small_jawed_project/Fst/sympatric_ponds/ASvP_95th_P_sweep_snps_20k.txt",header=T,stringsAsFactors = F)

ASvP.intro.results.exons.genes<-left_join(ASvP.intro.results.exons,ASvP.snps)

ASvP.intro.results.exons%>%dplyr::select(Exon_Start)%>%unique()%>%nrow()
ASvP.unique_exons<-ASvP.intro.results.exons%>%dplyr::select(Exon_Start)%>%unique()

gr2<-with(APvS.intro.results,GRanges(Scaffold,IRanges(start=Intro_Start, end=Intro_Stop,names=Gene)))
gr1<-with(exons,GRanges(Scaffold,IRanges(start=Exon_Start, end=Exon_Stop)))

type1 = findOverlaps(query = gr1, subject = gr2, type = 'within')
type1.df = data.frame(exons[queryHits(type1),], APvS.intro.results[subjectHits(type1),])

type2 = findOverlaps(query = gr2, subject = gr1, type = 'any')
type2.df = data.frame(APvS.intro.results[queryHits(type2),], exons[subjectHits(type2),])

APvS.intro.results.exons<-type2.df 

APvS.intro.results.exons%>%dplyr::select(Exon_Start)%>%unique()%>%nrow()

gr2<-with(AvSP.intro.results,GRanges(Scaffold,IRanges(start=Intro_Start, end=Intro_Stop,names=Gene)))
gr1<-with(exons,GRanges(Scaffold,IRanges(start=Exon_Start, end=Exon_Stop)))


type1 = findOverlaps(query = gr1, subject = gr2, type = 'within')
type1.df = data.frame(exons[queryHits(type1),], AvSP.intro.results[subjectHits(type1),])

type2 = findOverlaps(query = gr1, subject = gr2, type = 'any')
type2.df = data.frame(exons[queryHits(type2),], AvSP.intro.results[subjectHits(type2),])

AvSP.intro.results.exons<-type2.df 
#0 - not sure why some are showing up.
AvSP.intro.results.exons%>%dplyr::select(Exon_Start)%>%unique()%>%nrow()


all.intro.windows<-read.csv('10kb_window_OSPPvOSPS_diversity_rndmin_stats.csv',header=T,stringsAsFactors = F)
colnames(all.intro.windows)[1]<-"Scaffold"
all.intro.windows$to.pos<-as.integer(all.intro.windows$to.pos)
all.intro.windows$Scaffold = paste0('HiC_', all.intro.windows$Scaffold )
str(all.intro.windows)
str(exons)

gr2<-with(all.intro.windows,GRanges(Scaffold,IRanges(start=from.pos, end=to.pos)))
gr1<-with(exons,GRanges(Scaffold,IRanges(start=Exon_Start, end=Exon_Stop)))


type1 = findOverlaps(query = gr1, subject = gr2, type = 'within')
type1.df = data.frame(exons[queryHits(type1),], all.intro.windows[subjectHits(type1),])

type2 = findOverlaps(query = gr1, subject = gr2, type = 'any')
type2.df = data.frame(exons[queryHits(type2),], all.intro.windows[subjectHits(type2),])


all.intro.windows.results.exons<-type2.df
all.intro.windows.results.exons %>% dplyr::select(Scaffold,from.pos,to.pos)%>%unique()%>%nrow()
all.intro.windows %>% dplyr::select(Scaffold,from.pos,to.pos)%>%unique()%>%nrow()


all.intro.windows.results.exons2<-all.intro.windows.results.exons %>% group_by(Scaffold,from.pos,to.pos)%>%mutate(exon_count=n())
test<-left_join(all.intro.windows,all.intro.windows.results.exons2)%>%mutate_at('exon_count', ~replace(., is.na(.), 0))

AvSP.intro.results.exons2<-AvSP.intro.results.exons %>% group_by(Scaffold,Intro_Start,Intro_Stop)%>%mutate(count=n())
ASvP.intro.results.exons2<-ASvP.intro.results.exons %>% group_by(Scaffold,Intro_Start,Intro_Stop)%>%mutate(count=n())
APvS.intro.results.exons2<-APvS.intro.results.exons %>% group_by(Scaffold,Intro_Start,Intro_Stop)%>%mutate(count=n())


ecdf_fun <- function(x,perc) ecdf(x)(perc)
ecdf_fun(test$exon_count,32) #four exon in a window is in the 55% percentile of the distribution


ecdf_fun <- function(x,perc) ecdf(x)(perc)
ecdf_fun(test$exon_count,5) #four exon in a window is in the 55% percentile of the distribution
###### all intro regions #######

Intro.df.null.div_popsize_fsc.95<-read.table("ARTout_Intro_fd_D_df_RND_outlier_regions.fsc_95th_ms_sim_threshold.txt",stringsAsFactors = F,header=T)
            
Intro.df.null.div_popsize_msmc.95<-read.table("ARTout_Intro_fd_D_df_RND_outlier_regions.msmc_95th_ms_sim_threshold.txt",stringsAsFactors = F,header=T)



Intro.df.null.div_popsize_fsc.99<-read.table("ARTout_Intro_fd_D_df_RND_outlier_regions.fsc_99th_ms_sim_threshold.txt",stringsAsFactors = F,header=T)

Intro.df.null.div_popsize_msmc.99<-read.table("ARTout_Intro_fd_D_df_RND_outlier_regions.msmc_99th_ms_sim_threshold.txt",stringsAsFactors = F,header=T)
            
  
exon_gene_locations.txt
         
#######

head(exons)
exons2<- exons%>%group_by(Scaffold)%>%
  summarise(min.exon.pos=min(Exon_Start),max.exon.pos=max(Exon_Start))%>%
  mutate(Scaff.Start=plyr::round_any(min.exon.pos,10000,floor))%>%
  mutate(Scaff.Stop=plyr::round_any(max.exon.pos,10000, ceiling))


genome.exon.count.df<-exons%>%group_by(Scaffold)%>%
  mutate(Window.Start=plyr::round_any(Exon_Start,10000,floor))%>%
  mutate(Window.Stop=plyr::round_any(Exon_Stop,10000, ceiling))%>%
  group_by(Scaffold, Window.Start,Window.Stop)%>%
  summarise(exon.number=n())


range(genome.exon.count.df$exon.number)
