library(ggplot2)
library(dplyr)
library(lme4)
library(MASS)
library(pscl)

scale.counts<-read.csv('~/Desktop/Research/Small_jawed_project/gut_contents/APS_scale_counts.csv',header=T,stringsAsFactors = F)

scale.counts$Species<-as.factor(scale.counts$Species)

scale.counts<-scale.counts %>%
  filter(Population=='OSP')


SP.scale.counts<-scale.counts %>%
  filter(Species!='Generalist' )

scale.counts<-scale.counts %>%
mutate(Species = factor(Species, levels=c("Generalist", "Small.Jawed", "Scale.Eater")))
  
group.colors<-c(Scale.Eater="aquamarine3",Generalist="goldenrod1",Small.Jawed="#D95F02")

ggplot(scale.counts,aes(x=Species,y=ScaleCount,fill=Species))+
  #geom_boxplot(aes(col=Species),outlier.shape=NA,lwd=1.5,fill=NA)+
  geom_dotplot(binaxis='y', stackdir='center',method='histodot')+
  #geom_jitter(size=3,aes(alpha = 0.2),position=position_jitter(0.3))+
  #scale_color_manual(values = group.colors)+
  ylab('Scale Count')+
  scale_fill_manual(values = group.colors)+
  theme_classic()+theme(legend.position = 'none')



library(boot)


#pdf(file="~/Desktop/95_CI_morph_measurements_APS_allSL.pdf")

boot.calc.mean.scales<-function(data,i){
  df<-data[i,]
  c(mean(df$ScaleCount[which(df$Species=="Generalist")]),
    mean(df$ScaleCount[which(df$Species=="Small.Jawed")]),
    mean(df$ScaleCount[which(df$Species=="Scale.Eater")])
  )
}


SSI_scalecount<-boot(scale.counts, boot.calc.mean.scales, R=1000)


SSI.scalecount.df<-data.frame(Species=c("Generalist","Small.Jawed","Scale.Eater"),
                          Obs.Mean=c(SSI_scalecount$t0[1],
                                     SSI_scalecount$t0[2],
                                     SSI_scalecount$t0[3]),
                          Mean.LCL=c(boot.ci(SSI_scalecount, type="norm",index=1)$normal[2],
                                     boot.ci(SSI_scalecount, type="norm",index=2)$normal[2],
                                     boot.ci(SSI_scalecount, type="norm",index=3)$normal[2]),
                          Mean.UCL=c(boot.ci(SSI_scalecount, type="norm",index=1)$normal[3],
                                     boot.ci(SSI_scalecount, type="norm",index=2)$normal[3],
                                     boot.ci(SSI_scalecount, type="norm",index=3)$normal[3]))
SSI.scalecount.df$Species<-factor(SSI.scalecount.df$Species,levels=c('Generalist','Small.Jawed','Scale.Eater'))


ggplot(SSI.scalecount.df, aes(x=Species, y=Obs.Mean))+  
  geom_jitter(data=scale.counts,aes(x=Species,y=ScaleCount,alpha=0.2),size=2,position=position_jitter(0.3))+
  geom_point(size=5, aes(col=Species,fill=Species))+
  #scale_color_manual(values=group.colors)+
  ylab("Number of Scales")+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL, col=Species,fill=Species),size = 2)+
  scale_color_manual(values = group.colors)+
  theme_classic()+theme(legend.position = 'none')








ggplot(scale.counts,aes(x=Species,y=ScaleCount))+geom_boxplot(aes(col=Species))+
  theme_classic()



summary(glm(ScaleCount~Species,family=poisson,data=scale.counts))
ggplot(scale.counts, aes(ScaleCount, fill = Species)) +
  geom_histogram(binwidth=1, position="dodge")


summary(glm.nb(ScaleCount~Species,data=scale.counts))
summary(zeroinfl(ScaleCount~Species+StandardLength,data=scale.counts))


summary(glm(ScaleCount~Species,family=poisson,data=SP.scale.counts))
ggplot(SP.scale.counts, aes(ScaleCount, fill = Species)) +
  geom_histogram(binwidth=1, position="dodge")


odTest(glm.nb(ScaleCount~Species,data=SP.scale.counts))

summary(glm.nb(ScaleCount~Species,data=SP.scale.counts))
summary(zeroinfl(ScaleCount~Species,data=SP.scale.counts))

library(dplyr)
scale.counts.binom<-scale.counts %>%
  mutate(Presence = as.integer(ScaleCount >=1 & ScaleCount !=0))

summary(glm(Presence~Species,family = binomial(link="cloglog"),data=scale.counts.binom))



summary(aov(ScaleCount~Species,data=scale.counts))
TukeyHSD((aov(ScaleCount~Species,data=scale.counts)))
summary()
kruskal.test(ScaleCount~Species,data=scale.counts)


wilcox.test(ScaleCount~Species,data=SP.scale.counts, alternative = "two.sided",exact = FALSE)


library(dunn.test)
kw<-kruskal.test(ScaleCount~Species,data=scale.counts)
dunn.test(scale.counts$ScaleCount,scale.counts$Species,method="hs")
dunn.test(scale.counts$ScaleCount,scale.counts$Species,method="bh")
dunn.test(scale.counts$ScaleCount,scale.counts$Species,method="bonferroni")
pairwise.wilcox.test(scale.counts$ScaleCount~scale.counts$Species,p.adjust.method = "BH")
