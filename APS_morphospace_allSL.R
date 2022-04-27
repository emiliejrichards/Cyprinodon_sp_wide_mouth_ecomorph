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


df<-read.csv('~/Desktop/Research/Small_jawed_project/morphology/APS.8measurements_8.18.20.csv',header=T,stringsAsFactors = F)
#df.2<-df%>% filter(Round='1') %>% filter(X10.StandardLength<25)%>%filter(Population!='OSP.F1')%>%filter(Population!='OSP.F2')
df.2<-df%>% filter(Round=='1') %>%
  #filter(X10.StandardLength<25) %>%
  filter(Population!='OSP.F1')%>%
  filter(Population!='OSP.F2')%>%
  filter(Notes!='P.like')%>%
  filter(Notes!='S.like')%>%
  filter(Notes!='SP.like')%>%
  filter(Notes!='M.like')

A.resids<-A.df.2%>%
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
  fit<-lm(log(A.resids[,i])~log(A.resids$X10.StandardLength))
  A.resids<-modelr::add_residuals(data=A.resids, fit, var="resid")
  
  names(A.resids)[names(A.resids) == "resid"] <- namelist[i-3]
}



A.df.2<-df.2%>%dplyr::filter(Species=='Generalist')
P.df.2<-df.2%>%dplyr::filter(Species=='Scale.Eater')
S.df.2<-df.2%>%dplyr::filter(Species=='Small.Jawed')




jaw.res<-lm(X1.LowerJawLength ~ X10.StandardLength,data=A.df.2)
jaw.b<-summary(jaw.res)$coefficients[1, 1]
jaw.m<-summary(jaw.res)$coefficients[2, 1]
#pred.res<-A.df.2$X1.LowerJawLength-predict(res)



P.df.2$X1.LowerJawLength.resids<-P.df.2$X1.LowerJawLength-((jaw.m*P.df.2$X10.StandardLength)+jaw.b)
S.df.2$X1.LowerJawLength.resids<-S.df.2$X1.LowerJawLength-((jaw.m*S.df.2$X10.StandardLength)+jaw.b)
A.df.2$X1.LowerJawLength.resids<-A.df.2$X1.LowerJawLength-((jaw.m*A.df.2$X10.StandardLength)+jaw.b)


cph.res<-lm(X9.CaudalPeduncleHeight ~ X10.StandardLength,data=A.df.2)
cph.b<-summary(cph.res)$coefficients[1, 1]
cph.m<-summary(cph.res)$coefficients[2, 1]
P.df.2$X9.CaudalPeduncleHeight.resids<-P.df.2$X9.CaudalPeduncleHeight-((cph.m*P.df.2$X10.StandardLength)+cph.b)
S.df.2$X9.CaudalPeduncleHeight.resids<-S.df.2$X9.CaudalPeduncleHeight-((cph.m*S.df.2$X10.StandardLength)+cph.b)
A.df.2$X9.CaudalPeduncleHeight.resids<-A.df.2$X9.CaudalPeduncleHeight-((cph.m*A.df.2$X10.StandardLength)+cph.b)


poh.res<-lm(X5.PreopercularHeight ~ X10.StandardLength,data=A.df.2)
poh.b<-summary(poh.res)$coefficients[1, 1]
poh.m<-summary(poh.res)$coefficients[2, 1]
P.df.2$X5.PreopercularHeight.resids<-P.df.2$X5.PreopercularHeight-((poh.m*P.df.2$X10.StandardLength)+poh.b)
S.df.2$X5.PreopercularHeight.resids<-S.df.2$X5.PreopercularHeight-((poh.m*S.df.2$X10.StandardLength)+poh.b)
A.df.2$X5.PreopercularHeight.resids<-A.df.2$X5.PreopercularHeight-((poh.m*A.df.2$X10.StandardLength)+poh.b)


bw.res<-lm(X16.Bodywidth ~ X10.StandardLength,data=A.df.2)
bw.b<-summary(bw.res)$coefficients[1, 1]
bw.m<-summary(bw.res)$coefficients[2, 1]
P.df.2$X16.Bodywidth.resids<-P.df.2$X16.Bodywidth-((bw.m*P.df.2$X10.StandardLength)+bw.b)
S.df.2$X16.Bodywidth.resids<-S.df.2$X16.Bodywidth-((bw.m*S.df.2$X10.StandardLength)+bw.b)
A.df.2$X16.Bodywidth.resids<-A.df.2$X16.Bodywidth-((bw.m*A.df.2$X10.StandardLength)+bw.b)

d2a.res<-lm(X8.DorsaltoAnalDistance ~ X10.StandardLength,data=A.df.2)
d2a.b<-summary(d2a.res)$coefficients[1, 1]
d2a.m<-summary(d2a.res)$coefficients[2, 1]
P.df.2$X8.DorsaltoAnalDistance.resids<-P.df.2$X8.DorsaltoAnalDistance-((d2a.m*P.df.2$X10.StandardLength)+d2a.b)
S.df.2$X8.DorsaltoAnalDistance.resids<-S.df.2$X8.DorsaltoAnalDistance-((d2a.m*S.df.2$X10.StandardLength)+d2a.b)
A.df.2$X8.DorsaltoAnalDistance.resids<-A.df.2$X8.DorsaltoAnalDistance-((d2a.m*A.df.2$X10.StandardLength)+d2a.b)

eye.res<-lm(X6.EyeDiameter ~ X10.StandardLength,data=A.df.2)
eye.b<-summary(eye.res)$coefficients[1, 1]
eye.m<-summary(eye.res)$coefficients[2, 1]

P.df.2$X6.EyeDiameter.resids<-P.df.2$X6.EyeDiameter-((eye.m*P.df.2$X10.StandardLength)+eye.b)
S.df.2$X6.EyeDiameter.resids<-S.df.2$X6.EyeDiameter-((eye.m*S.df.2$X10.StandardLength)+eye.b)
A.df.2$X6.EyeDiameter.resids<-A.df.2$X6.EyeDiameter-((eye.m*A.df.2$X10.StandardLength)+eye.b)

buccal.res<-lm(X6.EyeDiameter ~ X10.StandardLength,data=A.df.2)
buccal.b<-summary(buccal.res)$coefficients[1, 1]
buccal.m<-summary(buccal.res)$coefficients[2, 1]

P.df.2$X14.BuccalWidth.resids<-P.df.2$X14.BuccalWidth-((buccal.m*P.df.2$X10.StandardLength)+buccal.b)
S.df.2$X14.BuccalWidth.resids<-S.df.2$X14.BuccalWidth-((buccal.m*S.df.2$X10.StandardLength)+buccal.b)
A.df.2$X14.BuccalWidth.resids<-A.df.2$X14.BuccalWidth-((buccal.m*A.df.2$X10.StandardLength)+buccal.b)

neuro.res<-lm(X6.EyeDiameter ~ X10.StandardLength,data=A.df.2)
neuro.b<-summary(neuro.res)$coefficients[1, 1]
neuro.m<-summary(neuro.res)$coefficients[2, 1]

P.df.2$neurocraniumwidth.resids<-P.df.2$neurocraniumwidth-((neuro.m*P.df.2$X10.StandardLength)+neuro.b)
S.df.2$neurocraniumwidth.resids<-S.df.2$neurocraniumwidth-((neuro.m*S.df.2$X10.StandardLength)+neuro.b)
A.df.2$neurocraniumwidth.resids<-A.df.2$neurocraniumwidth-((neuro.m*A.df.2$X10.StandardLength)+neuro.b)

head(A.df.2)
head(S.df.2)

head(P.df.2)
#y=mx+b

APS.df.Aresids<-rbind(A.df.2,S.df.2,P.df.2)

#ypredicted <- predict(res)
#residuals <- y - ypredicted



OSP.df.2<-df.2%>%dplyr::filter(Population=='OSP.F0')
GRE.df.2<-df.2%>%dplyr::filter(Population=='GRE.F0')
OYS.df.2<-df.2%>%dplyr::filter(Population=='OYS.F0')
LIL.df.2<-df.2%>%dplyr::filter(Population=='LIL.F0')

lab.df.S<-df%>% filter(Round=='1') %>%
  #filter(X10.StandardLength<25) %>%
  filter(Species=='Small.Jawed')%>%
  filter(Population=='OSP.F1' | Population=='OSP.F0')

%>%
  filter(Population!='OSP.F2')%>%
  
  #### first look at linear relationships of traits with SL #######

ggplot(df.2, aes(y=X5.PreopercularHeight, x=X10.StandardLength,group=Species))+geom_point(aes(col=Species))+
  geom_text(label=df.2$ID)+
  geom_smooth(method="lm")+theme_classic()

ggplot(df.2, aes(y=X1.LowerJawLength, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=df.2$ID)+
  geom_smooth(method="lm")+theme_classic()

ggplot(df.2, aes(y=log(X1.LowerJawLength), x=log(X10.StandardLength), group=Species))+geom_point(aes(col=Species))+
  geom_text(label=df.2$ID)+
  geom_smooth(method="lm")+theme_classic()

ggplot(df.2, aes(y=X14.BuccalWidth, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=df.2$ID)+
  geom_smooth(method="lm")+theme_classic()


####### OSP  #######

OSP.jaw<-ggplot(OSP.df.2, aes(y=X1.LowerJawLength, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=OSP.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


OSP.preop<-ggplot(OSP.df.2, aes(y=X5.PreopercularHeight, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=OSP.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


OSP.buccal<-ggplot(OSP.df.2, aes(y=X14.BuccalWidth, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=OSP.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()

OSP.eye<-ggplot(OSP.df.2, aes(y=X6.EyeDiameter, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=OSP.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


OSP.peduncle<-ggplot(OSP.df.2, aes(y=X9.CaudalPeduncleHeight, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=OSP.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


OSP.d2a<-ggplot(OSP.df.2, aes(y=X8.DorsaltoAnalDistance, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=OSP.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()

OSP.neurocran<-ggplot(OSP.df.2, aes(y=neurocraniumwidth, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=OSP.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


OSP.bodywidth<-ggplot(OSP.df.2, aes(y=X16.Bodywidth, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=OSP.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()
####### GRE #######


GRE.jaw<-ggplot(GRE.df.2, aes(y=X1.LowerJawLength, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=GRE.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


GRE.preop<-ggplot(GRE.df.2, aes(y=X5.PreopercularHeight, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=GRE.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


GRE.buccal<-ggplot(GRE.df.2, aes(y=X14.BuccalWidth, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=GRE.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()

GRE.eye<-ggplot(GRE.df.2, aes(y=X6.EyeDiameter, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=GRE.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


GRE.peduncle<-ggplot(GRE.df.2, aes(y=X9.CaudalPeduncleHeight, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=GRE.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


GRE.d2a<-ggplot(GRE.df.2, aes(y=X8.DorsaltoAnalDistance, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=GRE.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()

GRE.neurocran<-ggplot(GRE.df.2, aes(y=neurocraniumwidth, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=GRE.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


GRE.bodywidth<-ggplot(GRE.df.2, aes(y=X16.Bodywidth, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=GRE.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()
###### OYS #########


OYS.jaw<-ggplot(OYS.df.2, aes(y=X1.LowerJawLength, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=OYS.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


OYS.preop<-ggplot(OYS.df.2, aes(y=X5.PreopercularHeight, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=OYS.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


OYS.buccal<-ggplot(OYS.df.2, aes(y=X14.BuccalWidth, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=OYS.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()

OYS.eye<-ggplot(OYS.df.2, aes(y=X6.EyeDiameter, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=OYS.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


OYS.peduncle<-ggplot(OYS.df.2, aes(y=X9.CaudalPeduncleHeight, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=OYS.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


OYS.d2a<-ggplot(OYS.df.2, aes(y=X8.DorsaltoAnalDistance, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=OYS.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()

OYS.neurocran<-ggplot(OYS.df.2, aes(y=neurocraniumwidth, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=OYS.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


OYS.bodywidth<-ggplot(OYS.df.2, aes(y=X16.Bodywidth, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=OYS.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()

####### LIL ########


LIL.jaw<-ggplot(LIL.df.2, aes(y=X1.LowerJawLength, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=LIL.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


LIL.preop<-ggplot(LIL.df.2, aes(y=X5.PreopercularHeight, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=LIL.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


LIL.buccal<-ggplot(LIL.df.2, aes(y=X14.BuccalWidth, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=LIL.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


LIL.eye<-ggplot(LIL.df.2, aes(y=X6.EyeDiameter, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=LIL.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


LIL.peduncle<-ggplot(LIL.df.2, aes(y=X9.CaudalPeduncleHeight, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=LIL.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


LIL.d2a<-ggplot(LIL.df.2, aes(y=X8.DorsaltoAnalDistance, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=LIL.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()

LIL.neurocran<-ggplot(LIL.df.2, aes(y=neurocraniumwidth, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=LIL.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


LIL.bodywidth<-ggplot(LIL.df.2, aes(y=X16.Bodywidth, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=LIL.df.2$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()
######## lab reared #########
lab.jaw<-ggplot(lab.df, aes(y=X1.LowerJawLength, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=lab.df$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


lab.preop<-ggplot(lab.df, aes(y=X5.PreopercularHeight, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=lab.df$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


lab.buccal<-ggplot(lab.df, aes(y=X14.BuccalWidth, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=lab.df$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()

lab.eye<-ggplot(lab.df, aes(y=X6.EyeDiameter, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=lab.df$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


lab.peduncle<-ggplot(lab.df, aes(y=X9.CaudalPeduncleHeight, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=lab.df$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


lab.d2a<-ggplot(lab.df, aes(y=X8.DorsaltoAnalDistance, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=lab.df$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()

lab.neurocran<-ggplot(lab.df, aes(y=neurocraniumwidth, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=lab.df$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()


lab.bodywidth<-ggplot(lab.df, aes(y=X16.Bodywidth, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_text(label=lab.df$ID,aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()



lab.S.jaw<-ggplot(lab.df.S, aes(y=X1.LowerJawLength, x=X10.StandardLength, group=Population))+geom_point(aes(col=Population))+
  geom_text(label=lab.df.S$ID,aes(col=Population))+
  geom_smooth(method="lm",aes(col=Population))+theme_classic()


####### plot traits by SL #######

pdf(file="~/Desktop/LM_morph_measurements_APS.pdf")

grid.arrange(OSP.jaw, GRE.jaw, OYS.jaw,LIL.jaw, ncol=2)
grid.arrange(OSP.preop, GRE.preop, OYS.preop,LIL.preop, ncol=2)
grid.arrange(OSP.buccal, GRE.buccal, OYS.buccal,LIL.buccal, ncol=2)
grid.arrange(OSP.neurocran, GRE.neurocran, OYS.neurocran,LIL.neurocran, ncol=2)
grid.arrange(OSP.d2a, GRE.d2a, OYS.d2a,LIL.d2a, ncol=2)
grid.arrange(OSP.eye, GRE.eye, OYS.eye,LIL.eye, ncol=2)
grid.arrange(OSP.peduncle, GRE.peduncle, OYS.peduncle,LIL.peduncle, ncol=2)
grid.arrange(OSP.bodywidth, GRE.bodywidth, OYS.bodywidth,LIL.bodywidth, ncol=2)

dev.off()
#######################
model<-lm(X1.LowerJawLength ~ X10.StandardLength*Species, data=df.2)
summary(model)
model<-lm(X1.LowerJawLength ~ X10.StandardLength*Species+Population, data=df.2)
summary(model)


model<-lm(X1.LowerJawLength ~ X10.StandardLength*Species, data=OSP.df.2)
summary(model)
model<-lm(X1.LowerJawLength ~ X10.StandardLength*Species, data=OYS.df.2)
summary(model)
model<-lm(X1.LowerJawLength ~ X10.StandardLength*Species, data=GRE.df.2)
summary(model)
model<-lm(X1.LowerJawLength ~ X10.StandardLength*Species, data=LIL.df.2)
summary(model)



model<-lm(X5.PreopercularHeight ~ X10.StandardLength*Species, data=OSP.df.2)
summary(model)
model<-lm(X5.PreopercularHeight ~ X10.StandardLength*Species, data=OYS.df.2)
summary(model)
model<-lm(X5.PreopercularHeight ~ X10.StandardLength*Species, data=GRE.df.2)
summary(model)
model<-lm(X5.PreopercularHeight ~ X10.StandardLength*Species, data=LIL.df.2)
summary(model)


model<-lm(X14.BuccalWidth ~ X10.StandardLength*Species, data=OSP.df.2)
summary(model)
model<-lm(X14.BuccalWidth ~ X10.StandardLength*Species, data=OYS.df.2)
summary(model)
model<-lm(X14.BuccalWidth ~ X10.StandardLength*Species, data=GRE.df.2)
summary(model)
model<-lm(X14.BuccalWidth ~ X10.StandardLength*Species, data=LIL.df.2)
summary(model)







####### Standardize dataset by all ##########
APS.resids<-df.2%>%
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
  fit<-lm(log(APS.resids[,i])~log(APS.resids$X10.StandardLength))
  APS.resids<-modelr::add_residuals(data=APS.resids, fit, var="resid")
  
  names(APS.resids)[names(APS.resids) == "resid"] <- namelist[i-3]
}


APS.resids<-APS.resids%>%
  dplyr::select(Species,ID,Population,
                X1.LowerJawLength.resids,
                #X10.StandardLength.resids,
                X9.CaudalPeduncleHeight.resids,
                X5.PreopercularHeight.resids,
                X16.Bodywidth.resids,
                X8.DorsaltoAnalDistance.resids,
                X6.EyeDiameter.resids,
                X14.BuccalWidth.resids,
                neurocraniumwidth.resids)



###### Standardize by just generalists #########

A.df.2<-df.2%>%dplyr::filter(Species=='Generalist')
P.df.2<-df.2%>%dplyr::filter(Species=='Scale.Eater')
S.df.2<-df.2%>%dplyr::filter(Species=='Small.Jawed')




jaw.res<-lm(X1.LowerJawLength ~ X10.StandardLength,data=A.df.2)
jaw.b<-summary(jaw.res)$coefficients[1, 1]
jaw.m<-summary(jaw.res)$coefficients[2, 1]
#pred.res<-A.df.2$X1.LowerJawLength-predict(res)



P.df.2$X1.LowerJawLength.resids<-P.df.2$X1.LowerJawLength-((jaw.m*P.df.2$X10.StandardLength)+jaw.b)
S.df.2$X1.LowerJawLength.resids<-S.df.2$X1.LowerJawLength-((jaw.m*S.df.2$X10.StandardLength)+jaw.b)
A.df.2$X1.LowerJawLength.resids<-A.df.2$X1.LowerJawLength-((jaw.m*A.df.2$X10.StandardLength)+jaw.b)


cph.res<-lm(X9.CaudalPeduncleHeight ~ X10.StandardLength,data=A.df.2)
cph.b<-summary(cph.res)$coefficients[1, 1]
cph.m<-summary(cph.res)$coefficients[2, 1]
P.df.2$X9.CaudalPeduncleHeight.resids<-P.df.2$X9.CaudalPeduncleHeight-((cph.m*P.df.2$X10.StandardLength)+cph.b)
S.df.2$X9.CaudalPeduncleHeight.resids<-S.df.2$X9.CaudalPeduncleHeight-((cph.m*S.df.2$X10.StandardLength)+cph.b)
A.df.2$X9.CaudalPeduncleHeight.resids<-A.df.2$X9.CaudalPeduncleHeight-((cph.m*A.df.2$X10.StandardLength)+cph.b)


poh.res<-lm(X5.PreopercularHeight ~ X10.StandardLength,data=A.df.2)
poh.b<-summary(poh.res)$coefficients[1, 1]
poh.m<-summary(poh.res)$coefficients[2, 1]
P.df.2$X5.PreopercularHeight.resids<-P.df.2$X5.PreopercularHeight-((poh.m*P.df.2$X10.StandardLength)+poh.b)
S.df.2$X5.PreopercularHeight.resids<-S.df.2$X5.PreopercularHeight-((poh.m*S.df.2$X10.StandardLength)+poh.b)
A.df.2$X5.PreopercularHeight.resids<-A.df.2$X5.PreopercularHeight-((poh.m*A.df.2$X10.StandardLength)+poh.b)


bw.res<-lm(X16.Bodywidth ~ X10.StandardLength,data=A.df.2)
bw.b<-summary(bw.res)$coefficients[1, 1]
bw.m<-summary(bw.res)$coefficients[2, 1]
P.df.2$X16.Bodywidth.resids<-P.df.2$X16.Bodywidth-((bw.m*P.df.2$X10.StandardLength)+bw.b)
S.df.2$X16.Bodywidth.resids<-S.df.2$X16.Bodywidth-((bw.m*S.df.2$X10.StandardLength)+bw.b)
A.df.2$X16.Bodywidth.resids<-A.df.2$X16.Bodywidth-((bw.m*A.df.2$X10.StandardLength)+bw.b)

d2a.res<-lm(X8.DorsaltoAnalDistance ~ X10.StandardLength,data=A.df.2)
d2a.b<-summary(d2a.res)$coefficients[1, 1]
d2a.m<-summary(d2a.res)$coefficients[2, 1]
P.df.2$X8.DorsaltoAnalDistance.resids<-P.df.2$X8.DorsaltoAnalDistance-((d2a.m*P.df.2$X10.StandardLength)+d2a.b)
S.df.2$X8.DorsaltoAnalDistance.resids<-S.df.2$X8.DorsaltoAnalDistance-((d2a.m*S.df.2$X10.StandardLength)+d2a.b)
A.df.2$X8.DorsaltoAnalDistance.resids<-A.df.2$X8.DorsaltoAnalDistance-((d2a.m*A.df.2$X10.StandardLength)+d2a.b)

eye.res<-lm(X6.EyeDiameter ~ X10.StandardLength,data=A.df.2)
eye.b<-summary(eye.res)$coefficients[1, 1]
eye.m<-summary(eye.res)$coefficients[2, 1]

P.df.2$X6.EyeDiameter.resids<-P.df.2$X6.EyeDiameter-((eye.m*P.df.2$X10.StandardLength)+eye.b)
S.df.2$X6.EyeDiameter.resids<-S.df.2$X6.EyeDiameter-((eye.m*S.df.2$X10.StandardLength)+eye.b)
A.df.2$X6.EyeDiameter.resids<-A.df.2$X6.EyeDiameter-((eye.m*A.df.2$X10.StandardLength)+eye.b)

buccal.res<-lm(X6.EyeDiameter ~ X10.StandardLength,data=A.df.2)
buccal.b<-summary(buccal.res)$coefficients[1, 1]
buccal.m<-summary(buccal.res)$coefficients[2, 1]

P.df.2$X14.BuccalWidth.resids<-P.df.2$X14.BuccalWidth-((buccal.m*P.df.2$X10.StandardLength)+buccal.b)
S.df.2$X14.BuccalWidth.resids<-S.df.2$X14.BuccalWidth-((buccal.m*S.df.2$X10.StandardLength)+buccal.b)
A.df.2$X14.BuccalWidth.resids<-A.df.2$X14.BuccalWidth-((buccal.m*A.df.2$X10.StandardLength)+buccal.b)

neuro.res<-lm(X6.EyeDiameter ~ X10.StandardLength,data=A.df.2)
neuro.b<-summary(neuro.res)$coefficients[1, 1]
neuro.m<-summary(neuro.res)$coefficients[2, 1]

P.df.2$neurocraniumwidth.resids<-P.df.2$neurocraniumwidth-((neuro.m*P.df.2$X10.StandardLength)+neuro.b)
S.df.2$neurocraniumwidth.resids<-S.df.2$neurocraniumwidth-((neuro.m*S.df.2$X10.StandardLength)+neuro.b)
A.df.2$neurocraniumwidth.resids<-A.df.2$neurocraniumwidth-((neuro.m*A.df.2$X10.StandardLength)+neuro.b)

#y=mx+b

APS.df.Aresids<-rbind(A.df.2,S.df.2,P.df.2)

APS.Aresids<-APS.df.Aresids%>%
  dplyr::select(Species,ID,Population,
                X1.LowerJawLength.resids,
                #X10.StandardLength.resids,
                X9.CaudalPeduncleHeight.resids,
                X5.PreopercularHeight.resids,
                X16.Bodywidth.resids,
                X8.DorsaltoAnalDistance.resids,
                X6.EyeDiameter.resids,
                X14.BuccalWidth.resids,
                neurocraniumwidth.resids)


###############################





res.pca<- prcomp(APS.resids[,c(4:11)])



OSP.resids<-APS.resids%>%
  dplyr::filter(Population=='OSP.F0')
OSP.Aresids<-APS.Aresids%>%
  dplyr::filter(Population=='OSP.F0')

OYS.APS.resids<-APS.resids%>%
  dplyr::filter(Population=='OYS.F0')

LIL.APS.resids<-APS.resids%>%
  dplyr::filter(Population=='LIL.F0')

MER.APS.resids<-APS.resids%>%
  dplyr::filter(Population=='MER.F0')


GRE.APS.resids<-APS.resids%>%
  dplyr::filter(Population=='GRE.F0')



OSP.APS.Aresids<-APS.Aresids%>%
  dplyr::filter(Population=='OSP.F0')

OYS.APS.Aresids<-APS.Aresids%>%
  dplyr::filter(Population=='OYS.F0')

LIL.APS.Aresids<-APS.Aresids%>%
  dplyr::filter(Population=='LIL.F0')

MER.APS.Aresids<-APS.Aresids%>%
  dplyr::filter(Population=='MER.F0')


GRE.APS.Aresids<-APS.Aresids%>%
  dplyr::filter(Population=='GRE.F0')


'''APS.resids<-APS.resids%>%
dplyr::filter(Population!='OSP.F2') %>% 
dplyr::filter(Population!='OSP.F1')'''

APS.res.pca<- prcomp(APS.resids[,c(4:11)])
scores<-APS.res.pca$x[,1:2]
ggdata<-data.frame(scores,Species=APS.resids$Species,ID=APS.resids$ID)

ggplot(ggdata,aes(x=PC1, y=PC2, color=factor(Species))) +
  geom_point(size=2) +
  geom_text(label=ggdata$ID)+
  stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Species),type="norm"),
               geom="polygon", level=0.90, alpha=0.2)+
  theme_classic()+theme(legend.position = 'none')



APS.Ares.pca<- prcomp(APS.Aresids[,c(4:11)])
scores<-APS.Ares.pca$x[,1:2]
ggdata<-data.frame(scores,Species=APS.Aresids$Species,ID=APS.Aresids$ID)

ggplot(ggdata,aes(x=PC1, y=PC2, color=factor(Species))) +
  geom_point(size=2) +
  geom_text(label=ggdata$ID)+
  stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Species),type="norm"),
               geom="polygon", level=0.90, alpha=0.2)+
  theme_classic()+theme(legend.position = 'none')


OSP.res.pca<- prcomp(OSP.APS.resids[,c(4:11)])
scores<-OSP.res.pca$x[,1:2]

#km<-kmeans(scores,centers=3,nstart=10)
ggdata<-data.frame(scores,Species=OSP.APS.resids$Species,ID=OSP.APS.resids$ID)
ggplot(ggdata,aes(x=PC1, y=PC2, color=factor(Species))) +
  geom_point(size=2) +
  geom_text(label=ggdata$ID)+
  stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Species),type="norm"),
               geom="polygon", level=0.90, alpha=0.2)+
  theme_classic()+theme(legend.position = 'none')


OSP.res.pca<- prcomp(OSP.APS.Aresids[,c(4:11)])
scores<-OSP.res.pca$x[,1:2]

#km<-kmeans(scores,centers=3,nstart=10)
ggdata<-data.frame(scores,Species=OSP.APS.Aresids$Species,ID=OSP.APS.Aresids$ID)
ggplot(ggdata,aes(x=PC1, y=PC2, color=factor(Species))) +
  geom_point(size=2) +
  geom_text(label=ggdata$ID)+
  stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Species),type="norm"),
               geom="polygon", level=0.90, alpha=0.2)+
  theme_classic()+theme(legend.position = 'none')



OYS.res.pca<- prcomp(OYS.APS.resids[,c(4:11)])
scores<-OYS.res.pca$x[,1:2]
#km<-kmeans(scores,centers=3,nstart=10)
ggdata<-data.frame(scores,Species=OYS.APS.resids$Species,ID=OYS.APS.resids$ID)
ggplot(ggdata,aes(x=PC1, y=PC2, color=factor(Species))) +
  geom_point(size=2) +
  geom_text(label=ggdata$ID)+
  stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Species),type="norm"),
               geom="polygon", level=0.90, alpha=0.2)+
  theme_classic()+theme(legend.position = 'none')


MER.res.pca<- prcomp(MER.APS.resids[,c(4:11)])
scores<-MER.res.pca$x[,1:2]
#km<-kmeans(scores,centers=3,nstart=10)
ggdata<-data.frame(scores,Species=MER.APS.resids$Species,ID=MER.APS.resids$ID)
ggplot(ggdata,aes(x=PC1, y=PC2, color=factor(Species))) +
  geom_point(size=2) +
  geom_text(label=ggdata$ID)+
  stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Species),type="norm"),
               geom="polygon", level=0.90, alpha=0.2)+
  theme_classic()+theme(legend.position = 'none')


GRE.res.pca<- prcomp(GRE.APS.resids[,c(4:11)])
scores<-GRE.res.pca$x[,1:2]
#km<-kmeans(scores,centers=3,nstart=10)
ggdata<-data.frame(scores,Species=GRE.APS.resids$Species,ID=GRE.APS.resids$ID)
ggplot(ggdata,aes(x=PC1, y=PC2, color=factor(Species))) +
  geom_point(size=2) +
  geom_text(label=ggdata$ID)+
  stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Species),type="norm"),
               geom="polygon", level=0.90, alpha=0.2)+
  theme_classic()+theme(legend.position = 'none')





#pdf(file="~/Desktop/RD_morph_measurements_APS.pdf")

ggplot(ggdata,aes(x=PC1, y=PC2, color=factor(Species))) +
  geom_point(size=2) +
  geom_text(label=ggdata$ID)
stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Species),type="norm"),
             geom="polygon", level=0.90, alpha=0.2)+
  theme_classic()+theme(legend.position = 'none')
#guides(color=guide_legend(c("Cluster","Lake")),fill=guide_legend(c("Cluster","Lake")))

#scale_shape_manual(values=c(3,16,18))+

####### LDA 
lda.model<-lda(Species~.,na.omit(OSP.APS.resids[c(1,4:11)]))
#just rough plot of individuals labeled by species to see what the LDA looks like
#plot(lda.model)
#get the predicted class values for each individual based on the lda model 
pred.SSI<-predict(lda.model, OSP.APS.resids[c(1,4:11)])
#see what species the individuals have been assigned to based on the predictions from the model - should be very high accuracy

SSI.df <- data.frame(species = OSP.APS.resids$Species, lda = pred.SSI$x)
dat_ell <- data.frame() 
SSI.df$species<-as.factor(SSI.df$species)
for(g in levels(SSI.df$species)){ 
  dat_ell <- rbind(dat_ell, cbind(as.data.frame(with(SSI.df[SSI.df$species==g,], ellipse(cor(lda.LD1, lda.LD2), 
                                                                                         scale=c(sd(lda.LD1),sd(lda.LD2)), 
                                                                                         centre=c(mean(lda.LD1),mean(lda.LD2))))),species=g)) 
} 

#ggplot(SSI.df) + geom_point(aes(lda.LD1, lda.LD2, colour = species), size = 2.5)
ggplot(SSI.df, aes(x=lda.LD1, y=lda.LD2, col=species) ) + 
  geom_point( size = 4, aes(color = species))+
  #geom_text(label=ggdata$ID)+
  geom_path(data=dat_ell,aes(x=x,y=y,color=species),size=1,linetype=2)+
  theme_bw()



lda.model<-lda(Species~.,na.omit(APS.resids[c(1,4:11)]))
#just rough plot of individuals labeled by species to see what the LDA looks like
#plot(lda.model)
#get the predicted class values for each individual based on the lda model 
pred.SSI<-predict(lda.model, APS.resids[c(1,4:11)])
#see what species the individuals have been assigned to based on the predictions from the model - should be very high accuracy

SSI.df <- data.frame(species = APS.resids$Species, lda = pred.SSI$x)



#add in some ellipses spanning mean and 2 sd to the dataframe so we can plot them and more easily see 
dat_ell <- data.frame() 
SSI.df$species<-as.factor(SSI.df$species)
for(g in levels(SSI.df$species)){ 
  dat_ell <- rbind(dat_ell, cbind(as.data.frame(with(SSI.df[SSI.df$species==g,], ellipse(cor(lda.LD1, lda.LD2), 
                                                                                         scale=c(sd(lda.LD1),sd(lda.LD2)), 
                                                                                         centre=c(mean(lda.LD1),mean(lda.LD2))))),species=g)) 
} 

ggplot(SSI.df, aes(x=lda.LD1, y=lda.LD2, col=species) ) + 
  geom_point( size = 4, aes(color = species))+
  #geom_text(label=ggdata$ID)+
  geom_path(data=dat_ell,aes(x=x,y=y,color=species),size=1,linetype=2)+
  theme_bw()



lda.model<-lda(Species~.,na.omit(APS.Aresids[c(1,4:11)]))
#just rough plot of individuals labeled by species to see what the LDA looks like
#plot(lda.model)
#get the predicted class values for each individual based on the lda model 
pred.SSI<-predict(lda.model, APS.Aresids[c(1,4:11)])
#see what species the individuals have been assigned to based on the predictions from the model - should be very high accuracy

SSI.df <- data.frame(species = APS.Aresids$Species, lda = pred.SSI$x)



#add in some ellipses spanning mean and 2 sd to the dataframe so we can plot them and more easily see 
dat_ell <- data.frame() 
SSI.df$species<-as.factor(SSI.df$species)
for(g in levels(SSI.df$species)){ 
  dat_ell <- rbind(dat_ell, cbind(as.data.frame(with(SSI.df[SSI.df$species==g,], ellipse(cor(lda.LD1, lda.LD2), 
                                                                                         scale=c(sd(lda.LD1),sd(lda.LD2)), 
                                                                                         centre=c(mean(lda.LD1),mean(lda.LD2))))),species=g)) 
} 

ggplot(SSI.df, aes(x=lda.LD1, y=lda.LD2, col=species) ) + 
  geom_point( size = 4, aes(color = species))+
  #geom_text(label=ggdata$ID)+
  geom_path(data=dat_ell,aes(x=x,y=y,color=species),size=1,linetype=2)+
  theme_bw()


library(boot)


pdf(file="~/Desktop/95_CI_morph_measurements_APS_allSL.pdf")

boot.calc.mean.jawlen<-function(data,i){
  df<-data[i,]
  c(mean(df$X1.LowerJawLength.resids[which(df$Species=="Generalist")]),
    mean(df$X1.LowerJawLength.resids[which(df$Species=="Small.Jawed")]),
    mean(df$X1.LowerJawLength.resids[which(df$Species=="Scale.Eater")])
  )
}


SSI_jawlen<-boot(APS.resids, boot.calc.mean.jawlen, R=1000)


SSI.jawlen.df<-data.frame(Species=c("Generalist","Small.Jawed","Scale.Eater"),
                          Obs.Mean=c(SSI_jawlen$t0[1],
                                     SSI_jawlen$t0[2],
                                     SSI_jawlen$t0[3]),
                          Mean.LCL=c(boot.ci(SSI_jawlen, type="norm",index=1)$normal[2],
                                     boot.ci(SSI_jawlen, type="norm",index=2)$normal[2],
                                     boot.ci(SSI_jawlen, type="norm",index=3)$normal[2]),
                          Mean.UCL=c(boot.ci(SSI_jawlen, type="norm",index=1)$normal[3],
                                     boot.ci(SSI_jawlen, type="norm",index=2)$normal[3],
                                     boot.ci(SSI_jawlen, type="norm",index=3)$normal[3]))
SSI.jawlen.df$Species<-factor(SSI.jawlen.df$Species,levels=c('Generalist','Small.Jawed','Scale.Eater'))


ggplot(SSI.jawlen.df, aes(x=Species, y=Obs.Mean, col=Species,fill=Species))+
  geom_point(size=5)+geom_jitter(data=APS.resids,aes(x=Species,y=X1.LowerJawLength.resids,shape=Population),size=3)+
  #scale_color_manual(values=group.colors)+
  ylab("Lower Jaw Length residuals")+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL),size = 2)+
  theme_classic()

SSI_jawlen<-boot(APS.Aresids, boot.calc.mean.jawlen, R=1000)


SSI.jawlen.df<-data.frame(Species=c("Generalist","Small.Jawed","Scale.Eater"),
                          Obs.Mean=c(SSI_jawlen$t0[1],
                                     SSI_jawlen$t0[2],
                                     SSI_jawlen$t0[3]),
                          Mean.LCL=c(boot.ci(SSI_jawlen, type="norm",index=1)$normal[2],
                                     boot.ci(SSI_jawlen, type="norm",index=2)$normal[2],
                                     boot.ci(SSI_jawlen, type="norm",index=3)$normal[2]),
                          Mean.UCL=c(boot.ci(SSI_jawlen, type="norm",index=1)$normal[3],
                                     boot.ci(SSI_jawlen, type="norm",index=2)$normal[3],
                                     boot.ci(SSI_jawlen, type="norm",index=3)$normal[3]))
SSI.jawlen.df$Species<-factor(SSI.jawlen.df$Species,levels=c('Generalist','Small.Jawed','Scale.Eater'))

ggplot(SSI.jawlen.df, aes(x=Species, y=Obs.Mean, col=Species,fill=Species))+
  geom_point(size=5)+geom_jitter(data=APS.Aresids,aes(x=Species,y=X1.LowerJawLength.resids,shape=Population),size=3)+
  #scale_color_manual(values=group.colors)+
  ylab("Lower Jaw Length A residuals")+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL),size = 2)+
  theme_classic()


boot.calc.mean.buccal<-function(data,i){
  df<-data[i,]
  c(mean(df$X14.BuccalWidth.resids[which(df$Species=="Generalist")]),
    mean(df$X14.BuccalWidth.resids[which(df$Species=="Small.Jawed")]),
    mean(df$X14.BuccalWidth.resids[which(df$Species=="Scale.Eater")])
  )
}


SSI_buccal<-boot(APS.resids, boot.calc.mean.buccal, R=1000)


SSI.buccal.df<-data.frame(Species=c("Generalist","Small.Jawed","Scale.Eater"),
                          Obs.Mean=c(SSI_buccal$t0[1],
                                     SSI_buccal$t0[2],
                                     SSI_buccal$t0[3]),
                          Mean.LCL=c(boot.ci(SSI_buccal, type="norm",index=1)$normal[2],
                                     boot.ci(SSI_buccal, type="norm",index=2)$normal[2],
                                     boot.ci(SSI_buccal, type="norm",index=3)$normal[2]),
                          Mean.UCL=c(boot.ci(SSI_buccal, type="norm",index=1)$normal[3],
                                     boot.ci(SSI_buccal, type="norm",index=2)$normal[3],
                                     boot.ci(SSI_buccal, type="norm",index=3)$normal[3]))
SSI.buccal.df$Species<-factor(SSI.buccal.df$Species,levels=c('Generalist','Small.Jawed','Scale.Eater'))


ggplot(SSI.buccal.df, aes(x=Species, y=Obs.Mean, col=Species,fill=Species))+
  geom_point(size=5)+geom_jitter(data=APS.resids,aes(x=Species,y=X14.BuccalWidth.resids,shape=Population),size=3)+
  #scale_color_manual(values=group.colors)+
  ylab("Buccal Width residuals")+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL),size = 2)+
  theme_classic()



SSI_buccal<-boot(APS.Aresids, boot.calc.mean.buccal, R=1000)


SSI.buccal.df<-data.frame(Species=c("Generalist","Small.Jawed","Scale.Eater"),
                          Obs.Mean=c(SSI_buccal$t0[1],
                                     SSI_buccal$t0[2],
                                     SSI_buccal$t0[3]),
                          Mean.LCL=c(boot.ci(SSI_buccal, type="norm",index=1)$normal[2],
                                     boot.ci(SSI_buccal, type="norm",index=2)$normal[2],
                                     boot.ci(SSI_buccal, type="norm",index=3)$normal[2]),
                          Mean.UCL=c(boot.ci(SSI_buccal, type="norm",index=1)$normal[3],
                                     boot.ci(SSI_buccal, type="norm",index=2)$normal[3],
                                     boot.ci(SSI_buccal, type="norm",index=3)$normal[3]))
SSI.buccal.df$Species<-factor(SSI.buccal.df$Species,levels=c('Generalist','Small.Jawed','Scale.Eater'))


ggplot(SSI.buccal.df, aes(x=Species, y=Obs.Mean, col=Species,fill=Species))+
  geom_point(size=5)+geom_jitter(data=APS.Aresids,aes(x=Species,y=X14.BuccalWidth.resids,shape=Population),size=3)+
  #scale_color_manual(values=group.colors)+
  ylab("Buccal Width A residuals")+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL),size = 2)+
  theme_classic()




boot.calc.mean.poh<-function(data,i){
  df<-data[i,]
  c(mean(df$X5.PreopercularHeight.resids[which(df$Species=="Generalist")]),
    mean(df$X5.PreopercularHeight.resids[which(df$Species=="Small.Jawed")]),
    mean(df$X5.PreopercularHeight.resids[which(df$Species=="Scale.Eater")])
  )
}


SSI_poh<-boot(APS.resids, boot.calc.mean.poh, R=1000)


SSI.poh.df<-data.frame(Species=c("Generalist","Small.Jawed","Scale.Eater"),
                       Obs.Mean=c(SSI_poh$t0[1],
                                  SSI_poh$t0[2],
                                  SSI_poh$t0[3]),
                       Mean.LCL=c(boot.ci(SSI_poh, type="norm",index=1)$normal[2],
                                  boot.ci(SSI_poh, type="norm",index=2)$normal[2],
                                  boot.ci(SSI_poh, type="norm",index=3)$normal[2]),
                       Mean.UCL=c(boot.ci(SSI_poh, type="norm",index=1)$normal[3],
                                  boot.ci(SSI_poh, type="norm",index=2)$normal[3],
                                  boot.ci(SSI_poh, type="norm",index=3)$normal[3]))
SSI.poh.df$Species<-factor(SSI.poh.df$Species,levels=c('Generalist','Small.Jawed','Scale.Eater'))


ggplot(SSI.poh.df, aes(x=Species, y=Obs.Mean, col=Species,fill=Species))+
  geom_point(size=5)+geom_jitter(data=APS.resids,aes(x=Species,y=X5.PreopercularHeight.resids,shape=Population))+
  #scale_color_manual(values=group.colors)+
  ylab("Preopercular Height residuals")+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL),size = 2)+
  theme_classic()



SSI_poh<-boot(APS.Aresids, boot.calc.mean.poh, R=1000)


SSI.poh.df<-data.frame(Species=c("Generalist","Small.Jawed","Scale.Eater"),
                       Obs.Mean=c(SSI_poh$t0[1],
                                  SSI_poh$t0[2],
                                  SSI_poh$t0[3]),
                       Mean.LCL=c(boot.ci(SSI_poh, type="norm",index=1)$normal[2],
                                  boot.ci(SSI_poh, type="norm",index=2)$normal[2],
                                  boot.ci(SSI_poh, type="norm",index=3)$normal[2]),
                       Mean.UCL=c(boot.ci(SSI_poh, type="norm",index=1)$normal[3],
                                  boot.ci(SSI_poh, type="norm",index=2)$normal[3],
                                  boot.ci(SSI_poh, type="norm",index=3)$normal[3]))
SSI.poh.df$Species<-factor(SSI.poh.df$Species,levels=c('Generalist','Small.Jawed','Scale.Eater'))


ggplot(SSI.poh.df, aes(x=Species, y=Obs.Mean, col=Species,fill=Species))+
  geom_point(size=5)+geom_jitter(data=APS.Aresids,aes(x=Species,y=X5.PreopercularHeight.resids,shape=Population),size=3)+
  #scale_color_manual(values=group.colors)+
  ylab("Preopercular Height A residuals")+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL),size = 2)+
  theme_classic()




dev.off()










df.2<-df%>% filter(Round!=2) %>% filter(X10.StandardLength<25)

par(mfrow=c(2,2))

jaw<-ggplot(df.2, aes(y=X1.LowerJawLength, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()+theme(legend.position = 'none')



POH<-ggplot(df.2, aes(y=X5.PreopercularHeight, x=X10.StandardLength,group=Species))+geom_point(aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()+theme(legend.position = 'none')



buccal<-ggplot(df.2, aes(y=X14.BuccalWidth, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()+theme(legend.position = 'none')



NCH<-ggplot(df.2, aes(y=neurocraniumwidth, x=X10.StandardLength, group=Species,label=ID))+geom_point(aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()+theme(legend.position = 'none')



library(gridExtra)

grid.arrange(jaw, POH, buccal, NCH, ncol=2)

BW<-ggplot(df.2, aes(y=X16.Bodywidth, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()+theme(legend.position = 'none')


D2A<-ggplot(df.2, aes(y=X8.DorsaltoAnalDistance, x=X10.StandardLength,group=Species))+geom_point(aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()+theme(legend.position = 'none')



eye<-ggplot(df.2, aes(y=X6.EyeDiameter, x=X10.StandardLength, group=Species))+geom_point(aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()+theme(legend.position = 'none')



CPH<-ggplot(df.2, aes(y=X9.CaudalPeduncleHeight, x=X10.StandardLength,group=Species))+geom_point(aes(col=Species))+
  geom_smooth(method="lm",aes(col=Species))+theme_classic()+theme(legend.position = 'none')



grid.arrange(BW, D2A, eye, CPH, ncol=2)

dev.off()



remove 2018_0188, 2018_0213, OSP001, OSP002

