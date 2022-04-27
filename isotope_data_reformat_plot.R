library(plyr)
library(dplyr)
library(tidyverse)
library(userfriendlyscience)

isotope.2018<-read.csv("~/Documents/Research/Small_jawed_project/isotope_data/pupfish_cichlid_isotope_data_2018.csv",header=T,stringsAsFactors = F)
isotope.2016<-read.csv("~/Documents/Research/Small_jawed_project/isotope_data/pupfish_isotopes_final_2016.csv",header=T,stringsAsFactors = F)

colnames(isotope.2018)[2]<-"d13C"
colnames(isotope.2018)[5]<-"d15N"
colnames(isotope.2018)[1]<-"ID"

isotope.2018.pups<-isotope.2018 %>%
  filter(str_detect(ID,'CRP|OSP')) %>%
  mutate(population=case_when(str_detect(ID,'CRP') ~ 'CRP',
                              str_detect(ID,'OSP') ~ 'OSP')) %>%
  mutate(species=case_when(str_detect(ID,'CRPA') ~ 'A',
                           str_detect(ID,'CRPP') ~ 'P',
                           str_detect(ID,'CRPM') ~ 'M',
                           str_detect(ID,'OSPA') ~ 'A',
                           str_detect(ID,'OSPP') ~ 'P',
                           str_detect(ID,'OSPS') ~ 'S',
                           str_detect(ID,'OSPM') ~ 'M'))%>%
  mutate(year="2018")%>%
  dplyr::select(species,population,d13C,d15N,year,ID)

#isotope.2018.pups$species[which(isotope.2018.pups$ID=='OSPS0365')]<-'P'


colnames(isotope.2016)[1]<-"species.old"
isotope.2016.pups<-isotope.2016%>%
  filter(species.old=="brontotheroides" | species.old=="variegatus" | species.old=="desquamator") %>%
  filter(population=="CRP" |population=="OSP" |population=="GRE" |population=="LIL" |population=="OYS") %>%
  mutate(species=case_when(str_detect(species.old,'variegatus') ~ 'A',
                           str_detect(species.old,'desquamator') ~ 'P',
                           str_detect(species.old,'brontotheroides') ~ 'M'))%>%
  mutate(species=ifelse(population=='GRE' & species.old=='desquamator', "S", species)) %>%                        
  mutate(year="2016",ID=Sample.ID)%>%
  dplyr::select(species,population,d13C,d15N,year,ID)

  
  

isotope.both.pups<-rbind(isotope.2016.pups,isotope.2018.pups)
isotope.both.pups<-na.omit(isotope.both.pups)

isotope.both.pups.gre.osp.oys<-isotope.both.pups%>%filter(population=="OSP" |population=="GRE"|population=="OYS")

isotope.both.pups.gre.osp<-isotope.both.pups%>%filter(population=="OSP" |population=="GRE" & species!='M')
isotope.both.pups.osp<-isotope.both.pups%>%filter(population=="OSP" & species!='M')
isotope.both.pups.crp<-isotope.both.pups%>%filter(population=="CRP")
isotope.both.pups.lil<-isotope.both.pups%>%filter(population=="LIL")


isotope.both.pups.osp$species<-factor(isotope.both.pups.osp$species,levels=c('A','S','P'))

isotope.both.pups.osp.filtered<-isotope.both.pups.osp%>%filter(ID!="OSPA0220")%>%filter(ID!="OSPP0217" )%>%filter(ID!="OSPP0042" )
isotope.both.pups.gre.osp.filtered<-isotope.both.pups.gre.osp%>%filter(ID!="OSPA0220")%>%filter(ID!="OSPP0217")



group.colors<-c(P="aquamarine3",A="goldenrod1",S="#D95F02",M="mediumorchid3")



ggplot(data=isotope.both.pups.osp.filtered,aes(x=d13C,y=d15N,col=species))+
  #geom_text(label=isotope.both.pups.osp.filtered$ID)+
  geom_point(size=3)+scale_color_manual(values=group.colors)+scale_fill_manual(values=group.colors)+
  stat_ellipse(geom="polygon", alpha=0.10,aes(fill =species),type="t", linetype=1, level=.95)+theme_classic()+theme(legend.position = "none")


ggplot(data=isotope.both.pups.gre.osp.filtered,aes(x=d13C,y=d15N,col=species))+
  geom_point(size=2)+scale_color_manual(values=group.colors)+scale_fill_manual(values=group.colors)+
  stat_ellipse(geom="polygon", alpha=0.10,aes(fill =species),type="t", linetype=1, level=.95)+theme_classic()+theme(legend.position = "none")












group.colors<-c(P="aquamarine3",A="goldenrod1",S="coral1",M="mediumorchid3")
ggplot(data=isotope.both.pups.gre.osp.oys,aes(x=d13C,y=d15N,col=species,shape=population))+geom_point(size=2)+scale_color_manual(values=group.colors)+theme_bw()
ggplot(data=isotope.both.pups.gre.osp,aes(x=d13C,y=d15N,col=species,shape=population))+geom_point(size=2)+scale_color_manual(values=group.colors)+theme_bw()
ggplot(data=isotope.both.pups.osp,aes(x=d13C,y=d15N,col=species))+geom_point(size=2)+scale_color_manual(values=group.colors)+theme_bw()
ggplot(data=isotope.both.pups.crp,aes(x=d13C,y=d15N,col=species))+geom_point(size=2)+scale_color_manual(values=group.colors)+theme_bw()
ggplot(data=isotope.both.pups,aes(x=d13C,y=d15N,col=species,shape=population))+geom_point(size=2)+scale_color_manual(values=group.colors)+theme_bw()
ggplot(data=isotope.both.pups.lil,aes(x=d13C,y=d15N,col=species))+geom_point(size=2)+scale_color_manual(values=group.colors)+theme_bw()



ggplot(data=isotope.both.pups.osp.filtered,aes(x=d13C,y=d15N,col=species))+
  geom_point(size=2)+scale_color_manual(values=group.colors)+scale_fill_manual(values=group.colors)+theme_classic()

ggplot(data=isotope.both.pups.osp,aes(x=d15N,fill=species))+
  geom_density(alpha=0.2)+scale_color_manual(values=group.colors)+
  scale_fill_manual(values=group.colors)+theme_classic()+
  geom_rug(stat="identity",aes(col=species),lwd=2)


ggplot(data=isotope.both.pups.osp,aes(y=d15N,x=species,fill=species))+
  geom_boxplot()+scale_color_manual(values=group.colors)+
  scale_fill_manual(values=group.colors)+theme_classic()

ggplot(data=isotope.both.pups.osp,aes(y=d13C,x=species,fill=species))+
  geom_boxplot()+scale_color_manual(values=group.colors)+
  scale_fill_manual(values=group.colors)+theme_classic()


ggplot(data=isotope.both.pups.osp,aes(x=d13C,fill=species))+
  geom_density(alpha=0.2)+scale_color_manual(values=group.colors)+
  scale_fill_manual(values=group.colors)+theme_classic()+
  geom_rug(stat="identity",aes(col=species),lwd=2)





boot.calc.mean.d15N<-function(data,i){
  df<-data[i,]
  c(mean(df$d15N[which(df$species=="A")]),
    mean(df$d15N[which(df$species=="S")]),
    mean(df$d15N[which(df$species=="P")])
  )
}


OSP_d15N<-boot(isotope.both.pups.osp, boot.calc.mean.d15N, R=1000)


OSP_d15N.df<-data.frame(Species=c("A","S","P"),
                          Obs.Mean=c(OSP_d15N$t0[1],
                                     OSP_d15N$t0[2],
                                     OSP_d15N$t0[3]),
                          Mean.LCL=c(boot.ci(OSP_d15N, type="norm",index=1)$normal[2],
                                     boot.ci(OSP_d15N, type="norm",index=2)$normal[2],
                                     boot.ci(OSP_d15N, type="norm",index=3)$normal[2]),
                          Mean.UCL=c(boot.ci(OSP_d15N, type="norm",index=1)$normal[3],
                                     boot.ci(OSP_d15N, type="norm",index=2)$normal[3],
                                     boot.ci(OSP_d15N, type="norm",index=3)$normal[3]))
OSP_d15N.df$Species<-factor(OSP_d15N.df$Species,levels=c('A','S','P'))


ggplot(OSP_d15N.df, aes(x=Species, y=Obs.Mean, col=Species,fill=Species))+
  geom_point(size=5)+
  scale_color_manual(values=group.colors)+ylab("d15N")+
  geom_linerange(aes(ymin =Mean.LCL, ymax =Mean.UCL),size = 2)+
  theme_classic()






mean(isotope.both.pups.osp$d15N[which(isotope.both.pups.osp$species=='P')]) #10.88562
mean(isotope.both.pups.osp$d15N[which(isotope.both.pups.osp$species=='S')]) #9.529412
mean(isotope.both.pups.osp$d15N[which(isotope.both.pups.osp$species=='A')]) # 8.166634

OSP.AOV<-aov(d15N ~ species,data=isotope.both.pups.osp) #signif
#OSP.AOV<-aov(d13C ~ species,data=isotope.both.pups.osp) #not signif

shapiro.test(isotope.both.pups.osp$d15N)
bartlett.test(d15N~species, data=isotope.both.pups.osp)
summary(OSP.AOV)
'''            Df Sum Sq Mean Sq F value   Pr(>F)    
species      2  90.59   45.30   19.93 1.29e-07 ***
Residuals   72 163.64    2.27   '''

TukeyHSD((OSP.AOV))
'''  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = d15N ~ species, data = isotope.both.pups.osp)

$species
diff        lwr      upr     p adj
S-A 1.362777 0.32566709 2.399888 0.0067503
P-A 2.718991 1.65906008 3.778921 0.0000001
P-S 1.356213 0.09954524 2.612881 0.0313910'''

CRP.AOV<-aov(d15N ~ species,data=isotope.both.pups.crp)
summary(CRP.AOV)
TukeyHSD((CRP.AOV))



OSP.AOV<-oneway.test(d15N ~ species,data=isotope.both.pups.osp)
posthocTGH(isotope.both.pups.osp$d15N, isotope.both.pups.osp$species)
