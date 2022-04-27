library(dplyr)
library(tidyverse)
library(ggplot2)
library(Hmisc)
library(corrplot)
library(factoextra)
library(pcaMethods)
library(boot)
library(MASS)
library(ellipse)
library(boot)

#rm(list=ls())

df<-read.csv('~/Documents/Research/Small_jawed_project/morphology/APS.8measurements_2.08.22.csv',header=T,stringsAsFactors = F)


##### WILD VS LAB POPULATIONS ######


wild.lab.df<-df%>% filter(Round=='1') %>%
  filter(X10.StandardLength<35) %>%
  filter(Population=='OSP.F0'|Population=='OSP.F1'|
           Population=='OSP.F2'|Population=='OSP.F3')%>%
  filter(Notes!='P.like')%>%
  filter(Notes!='S.like')%>%
  filter(Notes!='SP.like')%>%
  filter(Notes!='M.like')



A.wild.test<-wild.lab.df%>%filter(Species=="Generalist" & Population =="OSP.F0")
mean(A.wild.test$X10.StandardLength)
A.lab.test<-wild.lab.df%>%filter(Species=="Generalist") %>% filter(Population =="OSP.F2" | Population == "OSP.F3")
mean(A.lab.test$X10.StandardLength)

P.wild.test<-wild.lab.df%>%filter(Species=="Scale.Eater" & Population =="OSP.F0")
mean(P.wild.test$X10.StandardLength)
P.lab.test<-wild.lab.df%>%filter(Species=="Scale.Eater") %>% filter(Population =="OSP.F2" | Population == "OSP.F3")
mean(P.lab.test$X10.StandardLength)


S.wild.test<-wild.lab.df%>%filter(Species=="Small.Jawed" & Population =="OSP.F0")
mean(S.wild.test$X10.StandardLength)
S.lab.test<-wild.lab.df%>%filter(Species=="Small.Jawed") %>% filter(Population =="OSP.F1" )
mean(S.lab.test$X10.StandardLength)


SP.lab<-rbind(S.lab.test,P.lab.test)

t.test(X14.BuccalWidth~Species,data=SP.lab)

t.test(X1.LowerJawLength~Species,data=SP.lab)
summary(aov(X1.LowerJawLength~Species,data=SP.lab))
wild.lab.resids<-wild.lab.df%>%
  dplyr::select(Species,ID,Population,
                X1.LowerJawLength,
                X10.StandardLength,
                X9.CaudalPeduncleHeight,
                X5.PreopercularHeight,
                X16.Bodywidth,
                X8.DorsaltoAnalDistance,
                X6.EyeDiameter,
                X14.BuccalWidth,
                neurocraniumwidth)


namelist<-c("X1.LowerJawLength.resids",
            "X10.StandardLength.resids",
            "X9.CaudalPeduncleHeight.resids",
            "X5.PreopercularHeight.resids",
            "X16.Bodywidth.resids",
            "X8.DorsaltoAnalDistance.resids",
            "X6.EyeDiameter.resids",
            "X14.BuccalWidth.resids",
            "neurocraniumwidth.resids")

for(i in 4:12){
  fit<-lm(log(wild.lab.resids[,i])~log(wild.lab.resids$X10.StandardLength))
  wild.lab.resids<-modelr::add_residuals(data=wild.lab.resids, fit, var="resid")
  
  names(wild.lab.resids)[names(wild.lab.resids) == "resid"] <- namelist[i-3]
}



fit<-lm(log(wild.lab.df$X5.PreopercularHeight)~log(wild.lab.df$X10.StandardLength))
summary(fit)


wild.lab.df%>%filter(Species== "Scale.Eater" | Species== "Small.Jawed")%>%
  group_by(Species, Population)%>%
  dplyr::summarize(Mean=mean(X5.PreopercularHeight))



wild.lab.df%>%filter(Species== "Scale.Eater" | Species== "Small.Jawed")%>%
  group_by(Species, Population)%>%
  dplyr::summarize(Mean=mean(X14.BuccalWidth))


wild.lab.resids<-wild.lab.resids%>%
  dplyr::select(Species,ID,Population,
                X1.LowerJawLength.resids,
                X10.StandardLength.resids,
                X9.CaudalPeduncleHeight.resids,
                X5.PreopercularHeight.resids,
                X16.Bodywidth.resids,
                X8.DorsaltoAnalDistance.resids,
                X6.EyeDiameter.resids,
                X14.BuccalWidth.resids,
                neurocraniumwidth.resids)




wild.lab.resids<-wild.lab.resids%>%mutate(Population2=str_replace(Population, "OSP.F3", "OSP.F2"))%>%
  filter(ID!="LLLF-060")%>%
  filter(ID!="LLLF-142")

wild.lab.resids%>%filter(Species== "Scale.Eater" | Species== "Small.Jawed")%>%
  group_by(Species, Population)%>%
  dplyr::summarize(Mean=mean(X5.PreopercularHeight.resids))



group.colors<-c(Scale.Eater="aquamarine3",Generalist="goldenrod1",Small.Jawed="#D95F02")


boot.calc.mean.jawlen<-function(data,i){
  df<-data[i,]
  c(mean(df$X1.LowerJawLength.resids[which(df$Species=="Generalist"& df$Population2=="OSP.F0")]),
    mean(df$X1.LowerJawLength.resids[which(df$Species=="Generalist"& df$Population2=="OSP.F2")]),
    mean(df$X1.LowerJawLength.resids[which(df$Species=="Small.Jawed"&df$Population2=="OSP.F0")]),
    mean(df$X1.LowerJawLength.resids[which(df$Species=="Small.Jawed"&df$Population2=="OSP.F1")]),
    mean(df$X1.LowerJawLength.resids[which(df$Species=="Scale.Eater"&df$Population2=="OSP.F0")]),
    mean(df$X1.LowerJawLength.resids[which(df$Species=="Scale.Eater"&df$Population2=="OSP.F2")]))
}

wild.lab.resids$X1.LowerJawLength.resids[which(wild.lab.resids$Species=="Scale.Eater"&wild.lab.resids$Population2=="OSP.F2")]

SSI_jawlen<-boot(wild.lab.resids, boot.calc.mean.jawlen, R=1000)


SSI.jawlen.df<-data.frame(Species=c("Generalist","Generalist","Small.Jawed","Small.Jawed","Scale.Eater","Scale.Eater"),
                          Population=c("Wild","Lab","Wild","Lab","Wild","Lab"),
                          Obs.Mean=c(SSI_jawlen$t0[1],
                                     SSI_jawlen$t0[2],
                                     SSI_jawlen$t0[3],
                                     SSI_jawlen$t0[4],
                                     SSI_jawlen$t0[5],
                                     SSI_jawlen$t0[6]),
                          Mean.LCL=c(boot.ci(SSI_jawlen, type="norm",index=1)$normal[2],
                                     boot.ci(SSI_jawlen, type="norm",index=2)$normal[2],
                                     boot.ci(SSI_jawlen, type="norm",index=3)$normal[2],
                                     boot.ci(SSI_jawlen, type="norm",index=4)$normal[2],
                                     boot.ci(SSI_jawlen, type="norm",index=5)$normal[2],
                                     boot.ci(SSI_jawlen, type="norm",index=6)$normal[2]),
                          Mean.UCL=c(boot.ci(SSI_jawlen, type="norm",index=1)$normal[3],
                                     boot.ci(SSI_jawlen, type="norm",index=2)$normal[3],
                                     boot.ci(SSI_jawlen, type="norm",index=3)$normal[3],
                                     boot.ci(SSI_jawlen, type="norm",index=4)$normal[3],
                                     boot.ci(SSI_jawlen, type="norm",index=5)$normal[3],
                                     boot.ci(SSI_jawlen, type="norm",index=6)$normal[3]))

SSI.jawlen.df$Species<-factor(SSI.jawlen.df$Species,levels=c('Generalist','Small.Jawed','Scale.Eater'))

#relevel(as.factor(SSI.jawlen.df$Population),ref="F0")
SSI.jawlen.df$Population<-factor(SSI.jawlen.df$Population,levels=c('Wild','Lab'))

SSI.jawlen<-ggplot(SSI.jawlen.df, aes(x=Species, y=Obs.Mean, col=Species,
                                      fill=Species,shape=Population))+
  geom_point(size=5,position = position_jitter(width = 0.25, seed = 123))+
  scale_shape_manual(values=c(16,17,17))+
  #geom_point(size=5)+scale_shape_manual(values=c(16,17,17))+
  #geom_jitter(data=wild.lab.resids,col='grey70',aes(x=Species,y=X1.LowerJawLength.resids,shape=Population),size=2)+
  #scale_shape_manual(values=c(16,15,17,16,15,17,18))+
  scale_color_manual(values=group.colors)+
  #geom_text(label=wild.lab.resids$ID))+
  ylab("Lower Jaw Length residuals")+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL,lty=Population),
                 size = 2,position = position_jitter(width = 0.25, seed = 123))+
  theme_classic()+theme(legend.position = 'none',
                        panel.grid.minor = element_blank(),panel.grid.major = element_blank())




boot.calc.mean.buccal<-function(data,i){
  df<-data[i,]
  c(mean(df$X14.BuccalWidth.resids[which(df$Species=="Generalist"& df$Population2=="OSP.F0")]),
    mean(df$X14.BuccalWidth.resids[which(df$Species=="Generalist"& df$Population2=="OSP.F2")]),
    mean(df$X14.BuccalWidth.resids[which(df$Species=="Small.Jawed"&df$Population2=="OSP.F0")]),
    mean(df$X14.BuccalWidth.resids[which(df$Species=="Small.Jawed"&df$Population2=="OSP.F1")]),
    mean(df$X14.BuccalWidth.resids[which(df$Species=="Scale.Eater"&df$Population2=="OSP.F0")]),
    mean(df$X14.BuccalWidth.resids[which(df$Species=="Scale.Eater"&df$Population2=="OSP.F2")]))
}


SSI_buccal<-boot(wild.lab.resids, boot.calc.mean.buccal, R=1000)


SSI.buccal.df<-data.frame(Species=c("Generalist","Generalist","Small.Jawed","Small.Jawed","Scale.Eater","Scale.Eater"),
                          Population=c("Wild","Lab","Wild","Lab","Wild","Lab"),
                          Obs.Mean=c(SSI_buccal$t0[1],
                                     SSI_buccal$t0[2],
                                     SSI_buccal$t0[3],
                                     SSI_buccal$t0[4],
                                     SSI_buccal$t0[5],
                                     SSI_buccal$t0[6]),
                          Mean.LCL=c(boot.ci(SSI_buccal, type="norm",index=1)$normal[2],
                                     boot.ci(SSI_buccal, type="norm",index=2)$normal[2],
                                     boot.ci(SSI_buccal, type="norm",index=3)$normal[2],
                                     boot.ci(SSI_buccal, type="norm",index=4)$normal[2],
                                     boot.ci(SSI_buccal, type="norm",index=5)$normal[2],
                                     boot.ci(SSI_buccal, type="norm",index=6)$normal[2]),
                          Mean.UCL=c(boot.ci(SSI_buccal, type="norm",index=1)$normal[3],
                                     boot.ci(SSI_buccal, type="norm",index=2)$normal[3],
                                     boot.ci(SSI_buccal, type="norm",index=3)$normal[3],
                                     boot.ci(SSI_buccal, type="norm",index=4)$normal[3],
                                     boot.ci(SSI_buccal, type="norm",index=5)$normal[3],
                                     boot.ci(SSI_buccal, type="norm",index=6)$normal[3]))

SSI.buccal.df$Species<-factor(SSI.buccal.df$Species,levels=c('Generalist','Small.Jawed','Scale.Eater'))
SSI.buccal.df$Population<-factor(SSI.buccal.df$Population,levels=c('Wild','Lab'))


SSI.buccal<-ggplot(SSI.buccal.df, aes(x=Species, y=Obs.Mean, col=Species,fill=Species,shape=Population))+
  geom_point(size=5)+scale_shape_manual(values=c(16,17,17))+
  #geom_jitter(data=wild.lab.resids,col='grey70',aes(x=Species,y=X14.BuccalWidth.resids,shape=Population),size=2)+
  #scale_shape_manual(values=c(16,15,17,16,15,17,18))+
  scale_color_manual(values=group.colors)+
  #geom_text(label=wild.lab.resids$ID))+
  ylab("Buccal Width Residuals")+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL,lty=Population),size = 2)+
  theme_classic()+theme(legend.position = 'none',panel.grid.minor = element_blank(),panel.grid.major = element_blank())




boot.calc.mean.poh<-function(data,i){
  df<-data[i,]
  c(mean(df$X5.PreopercularHeight.resids[which(df$Species=="Generalist"& df$Population2=="OSP.F0")]),
    mean(df$X5.PreopercularHeight.resids[which(df$Species=="Generalist"& df$Population2=="OSP.F2")]),
    mean(df$X5.PreopercularHeight.resids[which(df$Species=="Small.Jawed"&df$Population2=="OSP.F0")]),
    mean(df$X5.PreopercularHeight.resids[which(df$Species=="Small.Jawed"&df$Population2=="OSP.F1")]),
    mean(df$X5.PreopercularHeight.resids[which(df$Species=="Scale.Eater"&df$Population2=="OSP.F0")]),
    mean(df$X5.PreopercularHeight.resids[which(df$Species=="Scale.Eater"&df$Population2=="OSP.F2")]))
}


SSI_poh<-boot(wild.lab.resids, boot.calc.mean.poh, R=1000)


SSI.poh.df<-data.frame(Species=c("Generalist","Generalist","Small.Jawed","Small.Jawed","Scale.Eater","Scale.Eater"),
                       Population=c("Wild","Lab","Wild","Lab","Wild","Lab"),
                       Obs.Mean=c(SSI_poh$t0[1],
                                  SSI_poh$t0[2],
                                  SSI_poh$t0[3],
                                  SSI_poh$t0[4],
                                  SSI_poh$t0[5],
                                  SSI_poh$t0[6]),
                       Mean.LCL=c(boot.ci(SSI_poh, type="norm",index=1)$normal[2],
                                  boot.ci(SSI_poh, type="norm",index=2)$normal[2],
                                  boot.ci(SSI_poh, type="norm",index=3)$normal[2],
                                  boot.ci(SSI_poh, type="norm",index=4)$normal[2],
                                  boot.ci(SSI_poh, type="norm",index=5)$normal[2],
                                  boot.ci(SSI_poh, type="norm",index=6)$normal[2]),
                       Mean.UCL=c(boot.ci(SSI_poh, type="norm",index=1)$normal[3],
                                  boot.ci(SSI_poh, type="norm",index=2)$normal[3],
                                  boot.ci(SSI_poh, type="norm",index=3)$normal[3],
                                  boot.ci(SSI_poh, type="norm",index=4)$normal[3],
                                  boot.ci(SSI_poh, type="norm",index=5)$normal[3],
                                  boot.ci(SSI_poh, type="norm",index=6)$normal[3]))

SSI.poh.df$Species<-factor(SSI.poh.df$Species,levels=c('Generalist','Small.Jawed','Scale.Eater'))


SSI.poh<-ggplot(SSI.poh.df, aes(x=Species, y=Obs.Mean, col=Species,fill=Species,shape=Population))+
  geom_point(size=5)+scale_shape_manual(values=c(16,17,17))+
  #geom_jitter(data=wild.lab.resids,col='grey70',aes(x=Species,y=X5.PreopercularHeight.resids,shape=Population),size=2)+
  #scale_shape_manual(values=c(16,15,17,16,15,17,18))+
  scale_color_manual(values=group.colors)+
  ylab("POH residuals")+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL,lty=Population),size = 2)+
  theme_classic()+theme(legend.position = 'none',panel.grid.minor = element_blank(),panel.grid.major = element_blank())


boot.calc.mean.cph<-function(data,i){
  df<-data[i,]
  c(mean(df$X9.CaudalPeduncleHeight.resids[which(df$Species=="Generalist"& df$Population2=="OSP.F0")]),
    mean(df$X9.CaudalPeduncleHeight.resids[which(df$Species=="Generalist"& df$Population2=="OSP.F2")]),
    mean(df$X9.CaudalPeduncleHeight.resids[which(df$Species=="Small.Jawed"&df$Population2=="OSP.F0")]),
    mean(df$X9.CaudalPeduncleHeight.resids[which(df$Species=="Small.Jawed"&df$Population2=="OSP.F1")]),
    mean(df$X9.CaudalPeduncleHeight.resids[which(df$Species=="Scale.Eater"&df$Population2=="OSP.F0")]),
    mean(df$X9.CaudalPeduncleHeight.resids[which(df$Species=="Scale.Eater"&df$Population2=="OSP.F2")]))
}

df$X9.CaudalPeduncleHeight 
SSI_CPH<-boot(wild.lab.resids, boot.calc.mean.cph, R=1000)


SSI.CPH.df<-data.frame(Species=c("Generalist","Generalist","Small.Jawed","Small.Jawed","Scale.Eater","Scale.Eater"),
                          Population=c("Wild","Lab","Wild","Lab","Wild","Lab"),
                          Obs.Mean=c(SSI_CPH$t0[1],
                                     SSI_CPH$t0[2],
                                     SSI_CPH$t0[3],
                                     SSI_CPH$t0[4],
                                     SSI_CPH$t0[5],
                                     SSI_CPH$t0[6]),
                          Mean.LCL=c(boot.ci(SSI_CPH, type="norm",index=1)$normal[2],
                                     boot.ci(SSI_CPH, type="norm",index=2)$normal[2],
                                     boot.ci(SSI_CPH, type="norm",index=3)$normal[2],
                                     boot.ci(SSI_CPH, type="norm",index=4)$normal[2],
                                     boot.ci(SSI_CPH, type="norm",index=5)$normal[2],
                                     boot.ci(SSI_CPH, type="norm",index=6)$normal[2]),
                          Mean.UCL=c(boot.ci(SSI_CPH, type="norm",index=1)$normal[3],
                                     boot.ci(SSI_CPH, type="norm",index=2)$normal[3],
                                     boot.ci(SSI_CPH, type="norm",index=3)$normal[3],
                                     boot.ci(SSI_CPH, type="norm",index=4)$normal[3],
                                     boot.ci(SSI_CPH, type="norm",index=5)$normal[3],
                                     boot.ci(SSI_CPH, type="norm",index=6)$normal[3]))

SSI.CPH.df$Species<-factor(SSI.CPH.df$Species,levels=c('Generalist','Small.Jawed','Scale.Eater'))

#relevel(as.factor(SSI.jawlen.df$Population),ref="F0")
SSI.CPH.df$Population<-factor(SSI.CPH.df$Population,levels=c('Wild','Lab'))

SSI.CPH<-ggplot(SSI.CPH.df, aes(x=Species, y=Obs.Mean, col=Species,
                                      fill=Species,shape=Population))+
  geom_point(size=5,position = position_jitter(width = 0.25, seed = 123))+
  scale_shape_manual(values=c(16,17,17))+
  #geom_point(size=5)+scale_shape_manual(values=c(16,17,17))+
  #geom_jitter(data=wild.lab.resids,col='grey70',aes(x=Species,y=X1.LowerJawLength.resids,shape=Population),size=2)+
  #scale_shape_manual(values=c(16,15,17,16,15,17,18))+
  scale_color_manual(values=group.colors)+
  #geom_text(label=wild.lab.resids$ID))+
  ylab("Caudal Peduncle Height residuals")+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL,lty=Population),
                 size = 2,position = position_jitter(width = 0.25, seed = 123))+
  facet_wrap(~Population)+
  theme_classic()+theme(legend.position = 'none',
                        panel.grid.minor = element_blank(),panel.grid.major = element_blank())




boot.calc.mean.d2a<-function(data,i){
  df<-data[i,]
  c(mean(df$X8.DorsaltoAnalDistance.resids[which(df$Species=="Generalist"& df$Population2=="OSP.F0")]),
    mean(df$X8.DorsaltoAnalDistance.resids[which(df$Species=="Generalist"& df$Population2=="OSP.F2")]),
    mean(df$X8.DorsaltoAnalDistance.resids[which(df$Species=="Small.Jawed"&df$Population2=="OSP.F0")]),
    mean(df$X8.DorsaltoAnalDistance.resids[which(df$Species=="Small.Jawed"&df$Population2=="OSP.F1")]),
    mean(df$X8.DorsaltoAnalDistance.resids[which(df$Species=="Scale.Eater"&df$Population2=="OSP.F0")]),
    mean(df$X8.DorsaltoAnalDistance.resids[which(df$Species=="Scale.Eater"&df$Population2=="OSP.F2")]))
}

df$X8.DorsaltoAnalDistance
SSI_D2A<-boot(wild.lab.resids, boot.calc.mean.d2a, R=1000)


SSI.D2A.df<-data.frame(Species=c("Generalist","Generalist","Small.Jawed","Small.Jawed","Scale.Eater","Scale.Eater"),
                       Population=c("Wild","Lab","Wild","Lab","Wild","Lab"),
                       Obs.Mean=c(SSI_D2A$t0[1],
                                  SSI_D2A$t0[2],
                                  SSI_D2A$t0[3],
                                  SSI_D2A$t0[4],
                                  SSI_D2A$t0[5],
                                  SSI_D2A$t0[6]),
                       Mean.LCL=c(boot.ci(SSI_D2A, type="norm",index=1)$normal[2],
                                  boot.ci(SSI_D2A, type="norm",index=2)$normal[2],
                                  boot.ci(SSI_D2A, type="norm",index=3)$normal[2],
                                  boot.ci(SSI_D2A, type="norm",index=4)$normal[2],
                                  boot.ci(SSI_D2A, type="norm",index=5)$normal[2],
                                  boot.ci(SSI_D2A, type="norm",index=6)$normal[2]),
                       Mean.UCL=c(boot.ci(SSI_D2A, type="norm",index=1)$normal[3],
                                  boot.ci(SSI_D2A, type="norm",index=2)$normal[3],
                                  boot.ci(SSI_D2A, type="norm",index=3)$normal[3],
                                  boot.ci(SSI_D2A, type="norm",index=4)$normal[3],
                                  boot.ci(SSI_D2A, type="norm",index=5)$normal[3],
                                  boot.ci(SSI_D2A, type="norm",index=6)$normal[3]))

SSI.D2A.df$Species<-factor(SSI.D2A.df$Species,levels=c('Generalist','Small.Jawed','Scale.Eater'))

#relevel(as.factor(SSI.jawlen.df$Population),ref="F0")
SSI.D2A.df$Population<-factor(SSI.D2A.df$Population,levels=c('Wild','Lab'))

SSI.D2A<-ggplot(SSI.D2A.df, aes(x=Species, y=Obs.Mean, col=Species,
                                fill=Species,shape=Population))+
  geom_point(size=5,position = position_jitter(width = 0.25, seed = 123))+
  scale_shape_manual(values=c(16,17,17))+
  #geom_point(size=5)+scale_shape_manual(values=c(16,17,17))+
  #geom_jitter(data=wild.lab.resids,col='grey70',aes(x=Species,y=X1.LowerJawLength.resids,shape=Population),size=2)+
  #scale_shape_manual(values=c(16,15,17,16,15,17,18))+
  scale_color_manual(values=group.colors)+
  #geom_text(label=wild.lab.resids$ID))+
  ylab("Dorsal to Anal Fin Distance residuals")+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL,lty=Population),
                 size = 2,position = position_jitter(width = 0.25, seed = 123))+
  theme_classic()+theme(legend.position = 'none',
                        panel.grid.minor = element_blank(),panel.grid.major = element_blank())



boot.calc.mean.eye<-function(data,i){
  df<-data[i,]
  c(mean(df$X6.EyeDiameter.resids[which(df$Species=="Generalist"& df$Population2=="OSP.F0")]),
    mean(df$X6.EyeDiameter.resids[which(df$Species=="Generalist"& df$Population2=="OSP.F2")]),
    mean(df$X6.EyeDiameter.resids[which(df$Species=="Small.Jawed"&df$Population2=="OSP.F0")]),
    mean(df$X6.EyeDiameter.resids[which(df$Species=="Small.Jawed"&df$Population2=="OSP.F1")]),
    mean(df$X6.EyeDiameter.resids[which(df$Species=="Scale.Eater"&df$Population2=="OSP.F0")]),
    mean(df$X6.EyeDiameter.resids[which(df$Species=="Scale.Eater"&df$Population2=="OSP.F2")]))
}

df$X6.EyeDiameter
SSI_eye<-boot(wild.lab.resids, boot.calc.mean.eye, R=1000)


SSI.eye.df<-data.frame(Species=c("Generalist","Generalist","Small.Jawed","Small.Jawed","Scale.Eater","Scale.Eater"),
                       Population=c("Wild","Lab","Wild","Lab","Wild","Lab"),
                       Obs.Mean=c(SSI_eye$t0[1],
                                  SSI_eye$t0[2],
                                  SSI_eye$t0[3],
                                  SSI_eye$t0[4],
                                  SSI_eye$t0[5],
                                  SSI_eye$t0[6]),
                       Mean.LCL=c(boot.ci(SSI_eye, type="norm",index=1)$normal[2],
                                  boot.ci(SSI_eye, type="norm",index=2)$normal[2],
                                  boot.ci(SSI_eye, type="norm",index=3)$normal[2],
                                  boot.ci(SSI_eye, type="norm",index=4)$normal[2],
                                  boot.ci(SSI_eye, type="norm",index=5)$normal[2],
                                  boot.ci(SSI_eye, type="norm",index=6)$normal[2]),
                       Mean.UCL=c(boot.ci(SSI_eye, type="norm",index=1)$normal[3],
                                  boot.ci(SSI_eye, type="norm",index=2)$normal[3],
                                  boot.ci(SSI_eye, type="norm",index=3)$normal[3],
                                  boot.ci(SSI_eye, type="norm",index=4)$normal[3],
                                  boot.ci(SSI_eye, type="norm",index=5)$normal[3],
                                  boot.ci(SSI_eye, type="norm",index=6)$normal[3]))

SSI.eye.df$Species<-factor(SSI.eye.df$Species,levels=c('Generalist','Small.Jawed','Scale.Eater'))

#relevel(as.factor(SSI.jawlen.df$Population),ref="F0")
SSI.eye.df$Population<-factor(SSI.eye.df$Population,levels=c('Wild','Lab'))

SSI.eye<-ggplot(SSI.eye.df, aes(x=Species, y=Obs.Mean, col=Species,
                                fill=Species,shape=Population))+
  geom_point(size=5,position = position_jitter(width = 0.25, seed = 123))+
  scale_shape_manual(values=c(16,17,17))+
  #geom_point(size=5)+scale_shape_manual(values=c(16,17,17))+
  #geom_jitter(data=wild.lab.resids,col='grey70',aes(x=Species,y=X1.LowerJawLength.resids,shape=Population),size=2)+
  #scale_shape_manual(values=c(16,15,17,16,15,17,18))+
  scale_color_manual(values=group.colors)+
  #geom_text(label=wild.lab.resids$ID))+
  ylab("Eye Diameter residuals")+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL,lty=Population),
                 size = 2,position = position_jitter(width = 0.25, seed = 123))+
  theme_classic()+theme(legend.position = 'none',
                        panel.grid.minor = element_blank(),panel.grid.major = element_blank())



boot.calc.mean.neuro<-function(data,i){
  df<-data[i,]
  c(mean(df$neurocraniumwidth.resids[which(df$Species=="Generalist"& df$Population2=="OSP.F0")]),
    mean(df$neurocraniumwidth.resids[which(df$Species=="Generalist"& df$Population2=="OSP.F2")]),
    mean(df$neurocraniumwidth.resids[which(df$Species=="Small.Jawed"&df$Population2=="OSP.F0")]),
    mean(df$neurocraniumwidth.resids[which(df$Species=="Small.Jawed"&df$Population2=="OSP.F1")]),
    mean(df$neurocraniumwidth.resids[which(df$Species=="Scale.Eater"&df$Population2=="OSP.F0")]),
    mean(df$neurocraniumwidth.resids[which(df$Species=="Scale.Eater"&df$Population2=="OSP.F2")]))
}

df$neurocraniumwidth
SSI_neuro<-boot(wild.lab.resids, boot.calc.mean.neuro, R=1000)


SSI.neuro.df<-data.frame(Species=c("Generalist","Generalist","Small.Jawed","Small.Jawed","Scale.Eater","Scale.Eater"),
                       Population=c("Wild","Lab","Wild","Lab","Wild","Lab"),
                       Obs.Mean=c(SSI_neuro$t0[1],
                                  SSI_neuro$t0[2],
                                  SSI_neuro$t0[3],
                                  SSI_neuro$t0[4],
                                  SSI_neuro$t0[5],
                                  SSI_neuro$t0[6]),
                       Mean.LCL=c(boot.ci(SSI_neuro, type="norm",index=1)$normal[2],
                                  boot.ci(SSI_neuro, type="norm",index=2)$normal[2],
                                  boot.ci(SSI_neuro, type="norm",index=3)$normal[2],
                                  boot.ci(SSI_neuro, type="norm",index=4)$normal[2],
                                  boot.ci(SSI_neuro, type="norm",index=5)$normal[2],
                                  boot.ci(SSI_neuro, type="norm",index=6)$normal[2]),
                       Mean.UCL=c(boot.ci(SSI_neuro, type="norm",index=1)$normal[3],
                                  boot.ci(SSI_neuro, type="norm",index=2)$normal[3],
                                  boot.ci(SSI_neuro, type="norm",index=3)$normal[3],
                                  boot.ci(SSI_neuro, type="norm",index=4)$normal[3],
                                  boot.ci(SSI_neuro, type="norm",index=5)$normal[3],
                                  boot.ci(SSI_neuro, type="norm",index=6)$normal[3]))

SSI.neuro.df$Species<-factor(SSI.neuro.df$Species,levels=c('Generalist','Small.Jawed','Scale.Eater'))

#relevel(as.factor(SSI.jawlen.df$Population),ref="F0")
SSI.neuro.df$Population<-factor(SSI.neuro.df$Population,levels=c('Wild','Lab'))

SSI.neuro<-ggplot(SSI.neuro.df, aes(x=Species, y=Obs.Mean, col=Species,
                                fill=Species,shape=Population))+
  geom_point(size=5,position = position_jitter(width = 0.25, seed = 123))+
  scale_shape_manual(values=c(16,17,17))+
  #geom_point(size=5)+scale_shape_manual(values=c(16,17,17))+
  #geom_jitter(data=wild.lab.resids,col='grey70',aes(x=Species,y=X1.LowerJawLength.resids,shape=Population),size=2)+
  #scale_shape_manual(values=c(16,15,17,16,15,17,18))+
  scale_color_manual(values=group.colors)+
  #geom_text(label=wild.lab.resids$ID))+
  ylab("Neurocranium Width residuals")+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL,lty=Population),
                 size = 2,position = position_jitter(width = 0.25, seed = 123))+
  theme_classic()+theme(legend.position = 'none',
                        panel.grid.minor = element_blank(),panel.grid.major = element_blank())





boot.calc.mean.BW<-function(data,i){
  df<-data[i,]
  c(mean(df$X16.Bodywidth.resids[which(df$Species=="Generalist"& df$Population2=="OSP.F0")]),
    mean(df$X16.Bodywidth.resids[which(df$Species=="Generalist"& df$Population2=="OSP.F2")]),
    mean(df$X16.Bodywidth.resids[which(df$Species=="Small.Jawed"&df$Population2=="OSP.F0")]),
    mean(df$X16.Bodywidth.resids[which(df$Species=="Small.Jawed"&df$Population2=="OSP.F1")]),
    mean(df$X16.Bodywidth.resids[which(df$Species=="Scale.Eater"&df$Population2=="OSP.F0")]),
    mean(df$X16.Bodywidth.resids[which(df$Species=="Scale.Eater"&df$Population2=="OSP.F2")]))
}

SSI_BW<-boot(wild.lab.resids, boot.calc.mean.BW, R=1000)


SSI.BW.df<-data.frame(Species=c("Generalist","Generalist","Small.Jawed","Small.Jawed","Scale.Eater","Scale.Eater"),
                         Population=c("Wild","Lab","Wild","Lab","Wild","Lab"),
                         Obs.Mean=c(SSI_BW$t0[1],
                                    SSI_BW$t0[2],
                                    SSI_BW$t0[3],
                                    SSI_BW$t0[4],
                                    SSI_BW$t0[5],
                                    SSI_BW$t0[6]),
                         Mean.LCL=c(boot.ci(SSI_BW, type="norm",index=1)$normal[2],
                                    boot.ci(SSI_BW, type="norm",index=2)$normal[2],
                                    boot.ci(SSI_BW, type="norm",index=3)$normal[2],
                                    boot.ci(SSI_BW, type="norm",index=4)$normal[2],
                                    boot.ci(SSI_BW, type="norm",index=5)$normal[2],
                                    boot.ci(SSI_BW, type="norm",index=6)$normal[2]),
                         Mean.UCL=c(boot.ci(SSI_BW, type="norm",index=1)$normal[3],
                                    boot.ci(SSI_BW, type="norm",index=2)$normal[3],
                                    boot.ci(SSI_BW, type="norm",index=3)$normal[3],
                                    boot.ci(SSI_BW, type="norm",index=4)$normal[3],
                                    boot.ci(SSI_BW, type="norm",index=5)$normal[3],
                                    boot.ci(SSI_BW, type="norm",index=6)$normal[3]))

SSI.BW.df$Species<-factor(SSI.BW.df$Species,levels=c('Generalist','Small.Jawed','Scale.Eater'))

#relevel(as.factor(SSI.jawlen.df$Population),ref="F0")
SSI.BW.df$Population<-factor(SSI.BW.df$Population,levels=c('Wild','Lab'))

SSI.BW<-ggplot(SSI.BW.df, aes(x=Species, y=Obs.Mean, col=Species,
                                    fill=Species,shape=Population))+
  geom_point(size=5,position = position_jitter(width = 0.25, seed = 123))+
  scale_shape_manual(values=c(16,17,17))+
  #geom_point(size=5)+scale_shape_manual(values=c(16,17,17))+
  #geom_jitter(data=wild.lab.resids,col='grey70',aes(x=Species,y=X1.LowerJawLength.resids,shape=Population),size=2)+
  #scale_shape_manual(values=c(16,15,17,16,15,17,18))+
  scale_color_manual(values=group.colors)+
  #geom_text(label=wild.lab.resids$ID))+
  ylab("Body Width residuals")+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL,lty=Population),
                 size = 2,position = position_jitter(width = 0.25, seed = 123))+
  theme_classic()+theme(legend.position = 'none',
                        panel.grid.minor = element_blank(),panel.grid.major = element_blank())






library(gridExtra)
pdf(file="~/Documents/Research/Small_jawed_project/morphology/labvwild_95_CI_all_traits.pdf")

grid.arrange(SSI.jawlen,SSI.buccal,SSI.poh,
             SSI.CPH, SSI.D2A,SSI.BW,
             SSI.neuro, SSI.eye, ncol=3)

dev.off()



#################

SSI.buccal<-ggplot(SSI.buccal.df, aes(x=Species, y=Obs.Mean, col=Species,fill=Species,shape=Population))+
  geom_point(size=5)+scale_shape_manual(values=c(16,17,17))+
  #geom_jitter(data=wild.lab.resids,col='grey70',aes(x=Species,y=X14.BuccalWidth.resids,shape=Population),size=2)+
  #scale_shape_manual(values=c(16,15,17,16,15,17,18))+
  scale_color_manual(values=group.colors)+
  #geom_text(label=wild.lab.resids$ID))+
  ylab("Buccal Width Residuals")+
  facet_wrap(~Population)+
  scale_x_discrete(breaks=c("Generalist", "Small.Jawed", "Scale.Eater"),
                   labels=c("var", "wid", "des"))+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL,lty=Population),size = 2)+
  theme_classic()+theme(legend.position = 'none',panel.grid.minor = element_blank(),panel.grid.major = element_blank())




SSI.poh<-ggplot(SSI.poh.df, aes(x=Species, y=Obs.Mean, col=Species,fill=Species,shape=Population))+
  geom_point(size=5)+scale_shape_manual(values=c(16,17,17))+
  #geom_jitter(data=wild.lab.resids,col='grey70',aes(x=Species,y=X5.PreopercularHeight.resids,shape=Population),size=2)+
  #scale_shape_manual(values=c(16,15,17,16,15,17,18))+
  scale_color_manual(values=group.colors)+
  ylab("POH residuals")+
  facet_wrap(~Population)+
  scale_x_discrete(breaks=c("Generalist", "Small.Jawed", "Scale.Eater"),
                   labels=c("var", "wid", "des"))+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL,lty=Population),size = 2)+
  theme_classic()+theme(legend.position = 'none',panel.grid.minor = element_blank(),panel.grid.major = element_blank())

SSI.D2A<-ggplot(SSI.D2A.df, aes(x=Species, y=Obs.Mean, col=Species,
                                fill=Species,shape=Population))+
  geom_point(size=5,position = position_jitter(width = 0.25, seed = 123))+
  scale_shape_manual(values=c(16,17,17))+
  #geom_point(size=5)+scale_shape_manual(values=c(16,17,17))+
  #geom_jitter(data=wild.lab.resids,col='grey70',aes(x=Species,y=X1.LowerJawLength.resids,shape=Population),size=2)+
  #scale_shape_manual(values=c(16,15,17,16,15,17,18))+
  scale_color_manual(values=group.colors)+
  #geom_text(label=wild.lab.resids$ID))+
  ylab("Dorsal to Anal Fin Distance residuals")+
  facet_wrap(~Population)+
  scale_x_discrete(breaks=c("Generalist", "Small.Jawed", "Scale.Eater"),
                   labels=c("var", "wid", "des"))+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL,lty=Population),
                 size = 2,position = position_jitter(width = 0.25, seed = 123))+
  theme_classic()+theme(legend.position = 'none',
                        panel.grid.minor = element_blank(),panel.grid.major = element_blank())



SSI.eye<-ggplot(SSI.eye.df, aes(x=Species, y=Obs.Mean, col=Species,
                                fill=Species,shape=Population))+
  geom_point(size=5,position = position_jitter(width = 0.25, seed = 123))+
  scale_shape_manual(values=c(16,17,17))+
  #geom_point(size=5)+scale_shape_manual(values=c(16,17,17))+
  #geom_jitter(data=wild.lab.resids,col='grey70',aes(x=Species,y=X1.LowerJawLength.resids,shape=Population),size=2)+
  #scale_shape_manual(values=c(16,15,17,16,15,17,18))+
  scale_color_manual(values=group.colors)+
  #geom_text(label=wild.lab.resids$ID))+
  ylab("Eye Diameter residuals")+
  facet_wrap(~Population)+
  scale_x_discrete(breaks=c("Generalist", "Small.Jawed", "Scale.Eater"),
                   labels=c("var", "wid", "des"))+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL,lty=Population),
                 size = 2,position = position_jitter(width = 0.25, seed = 123))+
  theme_classic()+theme(legend.position = 'none',
                        panel.grid.minor = element_blank(),panel.grid.major = element_blank())


SSI.neuro<-ggplot(SSI.neuro.df, aes(x=Species, y=Obs.Mean, col=Species,
                                    fill=Species,shape=Population))+
  geom_point(size=5,position = position_jitter(width = 0.25, seed = 123))+
  scale_shape_manual(values=c(16,17,17))+
  #geom_point(size=5)+scale_shape_manual(values=c(16,17,17))+
  #geom_jitter(data=wild.lab.resids,col='grey70',aes(x=Species,y=X1.LowerJawLength.resids,shape=Population),size=2)+
  #scale_shape_manual(values=c(16,15,17,16,15,17,18))+
  scale_color_manual(values=group.colors)+
  #geom_text(label=wild.lab.resids$ID))+
  ylab("Neurocranium Width residuals")+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL,lty=Population),
                 size = 2,position = position_jitter(width = 0.25, seed = 123))+
  facet_wrap(~Population)+
  scale_x_discrete(breaks=c("Generalist", "Small.Jawed", "Scale.Eater"),
                   labels=c("var", "wid", "des"))+
  theme_classic()+theme(legend.position = 'none',
                        panel.grid.minor = element_blank(),panel.grid.major = element_blank())

SSI.BW<-ggplot(SSI.BW.df, aes(x=Species, y=Obs.Mean, col=Species,
                              fill=Species,shape=Population))+
  geom_point(size=5,position = position_jitter(width = 0.25, seed = 123))+
  scale_shape_manual(values=c(16,17,17))+
  #geom_point(size=5)+scale_shape_manual(values=c(16,17,17))+
  #geom_jitter(data=wild.lab.resids,col='grey70',aes(x=Species,y=X1.LowerJawLength.resids,shape=Population),size=2)+
  #scale_shape_manual(values=c(16,15,17,16,15,17,18))+
  scale_color_manual(values=group.colors)+
  #geom_text(label=wild.lab.resids$ID))+
  ylab("Body Width residuals")+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL,lty=Population),
                 size = 2,position = position_jitter(width = 0.25, seed = 123))+
  facet_wrap(~Population)+
  scale_x_discrete(breaks=c("Generalist", "Small.Jawed", "Scale.Eater"),
                   labels=c("var", "wid", "des"))+
  theme_classic()+theme(legend.position = 'none',
                        panel.grid.minor = element_blank(),panel.grid.major = element_blank())



SSI.CPH<-ggplot(SSI.CPH.df, aes(x=Species, y=Obs.Mean, col=Species,
                                fill=Species,shape=Population))+
  geom_point(size=5,position = position_jitter(width = 0.25, seed = 123))+
  scale_shape_manual(values=c(16,17,17))+
  #geom_point(size=5)+scale_shape_manual(values=c(16,17,17))+
  #geom_jitter(data=wild.lab.resids,col='grey70',aes(x=Species,y=X1.LowerJawLength.resids,shape=Population),size=2)+
  #scale_shape_manual(values=c(16,15,17,16,15,17,18))+
  scale_color_manual(values=group.colors)+
  #geom_text(label=wild.lab.resids$ID))+
  ylab("Caudal Peduncle Height residuals")+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL,lty=Population),
                 size = 2,position = position_jitter(width = 0.25, seed = 123))+
  facet_wrap(~Population)+
  scale_x_discrete(breaks=c("Generalist", "Small.Jawed", "Scale.Eater"),
                   labels=c("var", "wid", "des"))+
  theme_classic()+theme(legend.position = 'none',
                        panel.grid.minor = element_blank(),panel.grid.major = element_blank())

SSI.jawlen<-ggplot(SSI.jawlen.df, aes(x=Species, y=Obs.Mean, col=Species,
                                      fill=Species,shape=Population))+
  geom_point(size=5,position = position_jitter(width = 0.25, seed = 123))+
  scale_shape_manual(values=c(16,17,17))+
  #geom_point(size=5)+scale_shape_manual(values=c(16,17,17))+
  #geom_jitter(data=wild.lab.resids,col='grey70',aes(x=Species,y=X1.LowerJawLength.resids,shape=Population),size=2)+
  #scale_shape_manual(values=c(16,15,17,16,15,17,18))+
  scale_color_manual(values=group.colors)+
  #geom_text(label=wild.lab.resids$ID))+
  ylab("Lower Jaw Length residuals")+
  facet_wrap(~Population)+
  scale_x_discrete(breaks=c("Generalist", "Small.Jawed", "Scale.Eater"),
                   labels=c("var", "wid", "des"))+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL,lty=Population),
                 size = 2,position = position_jitter(width = 0.25, seed = 123))+
  theme_classic()+theme(legend.position = 'none',
                        panel.grid.minor = element_blank(),panel.grid.major = element_blank())


pdf(file="~/Documents/Research/Small_jawed_project/morphology/labvwild_95_CI_all_traits_WildvLab_panels.pdf")

grid.arrange(SSI.jawlen,SSI.buccal,SSI.poh,
             SSI.CPH, SSI.D2A,SSI.BW,
             SSI.neuro, SSI.eye, ncol=3)
dev.off()


A.lab.test<-wild.lab.df%>%filter(Species=="Generalist") %>% filter(Population =="OSP.F2" | Population == "OSP.F3")

wild.lab.resids<-wild.lab.resids %>% mutate(Population2 = ifelse(Population =="OSP.F0", "Wild", "Lab"))

jawlength.aov<-aov(X1.LowerJawLength.resids ~ Species*Population2,data=wild.lab.resids)
summary(jawlength.aov)
TukeyHSD(jawlength.aov)

poh.aov<-aov(X5.PreopercularHeight.resids ~ Species*Population2,data=wild.lab.resids)
summary(poh.aov)
TukeyHSD(poh.aov)

cph.aov<-aov(X9.CaudalPeduncleHeight.resids ~ Species*Population2,data=wild.lab.resids)
summary(cph.aov)
TukeyHSD(cph.aov)


d2a.aov<-aov(X8.DorsaltoAnalDistance.resids ~ Species*Population2,data=wild.lab.resids)
summary(d2a.aov)
TukeyHSD(d2a.aov)

bw.aov<-aov(X16.Bodywidth.resids ~ Species*Population2,data=wild.lab.resids)
summary(bw.aov)
TukeyHSD(bw.aov)

bw.aov<-aov(X14.BuccalWidth.resids ~ Species,data=wild.lab.resids)
summary(bw.aov)
TukeyHSD(bw.aov)

bw.aov<-aov(X14.BuccalWidth.resids ~ Species*Population2,data=wild.lab.resids)
summary(bw.aov)
TukeyHSD(bw.aov)

neuro.aov<-aov(neurocraniumwidth.resids ~ Species*Population2,data=wild.lab.resids)
summary(neuro.aov)
TukeyHSD(neuro.aov)

eye.aov<-aov(X6.EyeDiameter.resids ~ Species*Population2,data=wild.lab.resids)
summary(eye.aov)
TukeyHSD(eye.aov)




jawlength.aov<-aov(X1.LowerJawLength.resids ~ Species,data=wild.lab.resids)
summary(jawlength.aov)
TukeyHSD(jawlength.aov)

poh.aov<-aov(X5.PreopercularHeight.resids ~ Species,data=wild.lab.resids)
summary(poh.aov)
TukeyHSD(poh.aov)

cph.aov<-aov(X9.CaudalPeduncleHeight.resids ~ Species,data=wild.lab.resids)
summary(cph.aov)
TukeyHSD(cph.aov)


d2a.aov<-aov(X8.DorsaltoAnalDistance.resids ~ Species,data=wild.lab.resids)
summary(d2a.aov)
TukeyHSD(d2a.aov)

bw.aov<-aov(X16.Bodywidth.resids ~ Species,data=wild.lab.resids)
summary(bw.aov)
TukeyHSD(bw.aov)

buccal.aov<-aov(X14.BuccalWidth.resids ~ Species,data=wild.lab.resids)
summary(buccal.aov)
TukeyHSD(buccal.aov)


neuro.aov<-aov(neurocraniumwidth.resids ~ Species,data=wild.lab.resids)
summary(neuro.aov)
TukeyHSD(neuro.aov)

eye.aov<-aov(X6.EyeDiameter.resids ~ Species,data=wild.lab.resids)
summary(eye.aov)
TukeyHSD(eye.aov)

#######################
SSI.poh.df$Trait<-'PreOp Height'
SSI.jawlen.df$Trait<-'Lower Jaw Length'
SSI.buccal.df$Trait<-'Buccal Width'

SSI.traits.df<-rbind(SSI.poh.df,SSI.jawlen.df,SSI.buccal.df)

wild.lab.df2<-wild.lab.df%>%mutate(Population2=str_replace(Population, "OSP.F3", "OSP.F2"))

wild.lab.df2<-wild.lab.df2%>%mutate(Population2=case_when(Species='Generalist' & Population))
wild.lab.df3<-wild.lab.df2%>%mutate(Population2=str_replace(Population2, "OSP.F2", "OSP.F1"))
#wild.lab.df3$Population2<-which[wild.lab.df3$Population2=="OSP.F3"]<- "OSP.F2"                                
pdf(file="~/Documents/Research/Small_jawed_project/morphology/labvwild_SL_vs_jaw_buccal_poh.pdf")



SSI.poh<-ggplot(wild.lab.df3, aes(y=X5.PreopercularHeight, x=X10.StandardLength, group=interaction(Species, Population2)))+
  #geom_point(aes(col=Species,shape=Population2))+
  #geom_text(label=wild.lab.df2$ID,aes(col=Species))+
  ggtitle('Osprey Lake')+  scale_color_manual(values = group.colors)+
  scale_y_continuous(trans = "log10")+scale_x_continuous(trans = "log10")+
  ylab('log(Preopercular Height Insertion) (mm)')+ xlab('log(Standard Length) (mm)')+
  geom_smooth(method="lm",aes(col=Species,lty=Population2))+theme_classic()+
  theme(legend.position = 'none')+facet_wrap(.~Species)


SSI.jaw<-ggplot(wild.lab.df3, aes(y=X1.LowerJawLength, x=X10.StandardLength, group=interaction(Species, Population2)))+
  #geom_point(aes(col=Species,shape=Population2))+
  #geom_text(label=wild.lab.df2$ID,aes(col=Species))+
  ggtitle('Osprey Lake')+  scale_color_manual(values = group.colors)+
  scale_y_continuous(trans = "log10")+scale_x_continuous(trans = "log10")+
  ylab('log(Lower Jaw Length) (mm)')+ xlab('log(Standard Length) (mm)')+
  geom_smooth(method="lm",aes(col=Species,lty=Population2))+theme_classic()+
  theme(legend.position = 'none')+facet_wrap(.~Species)


SSI.buccal<-ggplot(wild.lab.df3, aes(y=X14.BuccalWidth, x=X10.StandardLength, group=interaction(Species, Population2)))+
  #geom_point(aes(col=Species,shape=Population2))+
  #geom_text(label=wild.lab.df2$ID,aes(col=Species))+
  ggtitle('Osprey Lake')+  scale_color_manual(values = group.colors)+
  scale_y_continuous(trans = "log10")+scale_x_continuous(trans = "log10")+
  ylab('log(Buccal Width) (mm)')+ xlab('log(Standard Length) (mm)')+
  geom_smooth(method="lm",aes(col=Species,lty=Population2))+theme_classic()+
  theme(legend.position = 'none')+facet_wrap(.~Species)

grid.arrange(SSI.jaw,SSI.buccal,SSI.poh,nrow=3)


dev.off()




SSI.poh<-ggplot(wild.lab.df3, aes(y=X5.PreopercularHeight, x=X10.StandardLength, group=interaction(Species, Population2)))+
  #geom_point(aes(col=Species,shape=Population2))+
  #geom_text(label=wild.lab.df2$ID,aes(col=Species))+
  ggtitle('Osprey Lake')+  scale_color_manual(values = group.colors)+
  scale_y_continuous(trans = "log10")+scale_x_continuous(trans = "log10")+
  ylab('log(Preopercular Height Insertion) (mm)')+ xlab('log(Standard Length) (mm)')+
  geom_smooth(method="lm",aes(col=Species,lty=Population2))+theme_classic()+
  theme(legend.position = 'none')+facet_wrap(.~Population2)


SSI.jaw<-ggplot(wild.lab.df3, aes(y=X1.LowerJawLength, x=X10.StandardLength, group=interaction(Species, Population2)))+
  #geom_point(aes(col=Species,shape=Population2))+
  #geom_text(label=wild.lab.df2$ID,aes(col=Species))+
  ggtitle('Osprey Lake')+  scale_color_manual(values = group.colors)+
  scale_y_continuous(trans = "log10")+scale_x_continuous(trans = "log10")+
  ylab('log(Lower Jaw Length) (mm)')+ xlab('log(Standard Length) (mm)')+
  geom_smooth(method="lm",aes(col=Species,lty=Population2))+theme_classic()+
  theme(legend.position = 'none')+facet_wrap(.~Population2)


SSI.buccal<-ggplot(wild.lab.df3, aes(y=X14.BuccalWidth, x=X10.StandardLength, group=interaction(Species, Population2)))+
  #geom_point(aes(col=Species,shape=Population2))+
  #geom_text(label=wild.lab.df2$ID,aes(col=Species))+
  ggtitle('Osprey Lake')+  scale_color_manual(values = group.colors)+
  scale_y_continuous(trans = "log10")+scale_x_continuous(trans = "log10")+
  ylab('log(Buccal Width) (mm)')+ xlab('log(Standard Length) (mm)')+
  geom_smooth(method="lm",aes(col=Species,lty=Population2))+theme_classic()+
  theme(legend.position = 'none')+facet_wrap(.~Population2)




A.wild.lab.df<-wild.lab.df%>%filter(Species=="Generalist")
ggplot(A.wild.lab.df, aes(y=X1.LowerJawLength, x=X10.StandardLength, group=Population2))+geom_point(col='goldenrod1',aes(shape=Population))+
  #geom_text(label=OSP.df.2$ID,aes(col=Species))+
  ggtitle('Osprey Lake')+ 
  scale_y_continuous(trans = "log10")+scale_x_continuous(trans = "log10")+
  ylab('log(Lower Jaw Length) (mm)')+ xlab('log(Standard Length) (mm)')+
  geom_smooth(method="lm",aes(lty=Population))+theme_classic()+
  theme(legend.position = 'none')


ggplot(A.wild.lab.df, aes(y=X14.BuccalWidth, x=X10.StandardLength, group=Population))+geom_point(col='goldenrod1',aes(shape=Population))+
  #geom_text(label=OSP.df.2$ID,aes(col=Species))+
  ggtitle('Osprey Lake')+ 
  scale_y_continuous(trans = "log10")+scale_x_continuous(trans = "log10")+
  ylab('log(Lower Jaw Length) (mm)')+ xlab('log(Standard Length) (mm)')+
  geom_smooth(method="lm",aes(lty=Population))+theme_classic()+
  theme(legend.position = 'none')


ggplot(A.wild.lab.df, aes(y=X5.PreopercularHeight, x=X10.StandardLength, group=Population))+geom_point(col='goldenrod1',aes(shape=Population))+
  #geom_text(label=OSP.df.2$ID,aes(col=Species))+
  ggtitle('Osprey Lake')+ 
  scale_y_continuous(trans = "log10")+scale_x_continuous(trans = "log10")+
  ylab('log(Lower Jaw Length) (mm)')+ xlab('log(Standard Length) (mm)')+
  geom_smooth(method="lm",aes(lty=Population))+theme_classic()+
  theme(legend.position = 'none')
