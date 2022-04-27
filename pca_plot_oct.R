library(ggplot2)

table1<-read.table("~/Documents/Research/Small_jawed_project/pca/AMPS_pruned_SSI_filtered_snps_passed.Q20.MAF0.05.MAXMISS0.9.feb.LDpruned.eigenvec",header=F,sep="")

library(plyr)
library(dplyr)
table1$Group<-c("A","A","A","A","A","M","M","M","M","M","M","M","M","M","M","M","M","P","P","P","P","P","P","P","P","P","P","A","A","A","P","S","S","A","M","M","M","M","M","P","P","P","P","P","S","S","S","S","A","A","A","A","P","P","A","M","M","M","M","M","M","A","A","A","A","A","A","A","A","A","A","A","A","A","M","M","M","M","M","M","M","M","M","M","M","M","M","P","P","P","P","P","P","P","P","P","P","P","S","S","S","S","S","S","S","S","S","S","S","A","A","M","M","M","M","M","M","M","M","P","P","P","P","P","S","S","S","S","S","A","A","A")
table1$Group<-as.factor(table1$Group)

eigenval3<-read.table("~/Documents/Research/Small_jawed_project/pca/AMPS_pruned_SSI_filtered_snps_passed.Q20.MAF0.05.MAXMISS0.9.feb.LDpruned.eigenval",header=F,sep="",stringsAsFactors = F)
sum_eig3<-sum(eigenval3$V1)

group.colors<-c(M="mediumorchid3",P="aquamarine3",A="goldenrod1",S="coral1")
,GEO="gray")

table1$Group[which(table1$V1  == "OYSP1")]<-"S"
table1$Group[which(table1$V1  == "OYSP3")]<-"S"
table1$Group[which(table1$V1  == "OYSP7")]<-"S"
table1$Group[which(table1$V1  == "OYSP4")]<-"S"
table1$Group[which(table1$V1  == "MERP1")]<-"S"
table1$Group[which(table1$V1  == "MERP2")]<-"S"
table1$Group[which(table1$V1  == "OSPP2")]<-"S"

OYSP3<-"S"
OYSP7<-"S"
OYSP4<-"S"
MERP1<-"S"
MERP2<-"S"
OSPP2<-"S"
#rm(OYSS9)

sum_eigs3<-lapply(eigenval3$V1,function(x){
  rt<-(x/sum_eig3)*100
  rt<-round(rt,digits=1)
  return(rt)
})

table1<-table1%>%filter(V1!='CRPA1003')%>%filter(V1!='PAIA1')
ggplot(table1,aes(x=table1$V3,y=table1$V4,label=table1$V1,color=table1$Group))+
  geom_point(size=0)+geom_text(aes(label=table1$V1),size=3,hjust=0, vjust=0)+
  scale_color_manual(values=group.colors)+
  labs(x=paste0("PC1: ",sum_eigs3[[1]],"% variance"),y=paste0("PC2: ",sum_eigs3[[2]],"% variance"),color="locality")+guides(size=F)+
  theme_bw()



ggplot(table1,aes(x=table1$V3,y=table1$V4,color=table1$Group))+
  geom_point(size=3)+
  scale_color_manual(values=group.colors)+
  labs(x=paste0("PC1: ",sum_eigs3[[1]],"% variance"),y=paste0("PC2: ",sum_eigs3[[2]],"% variance"),color="locality")+guides(size=F)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.text = element_text(size=12),
        panel.border = element_rect(fill=NA, colour = "black", size=1.2))







table1<-read.table("~/Documents/Research/Small_jawed_project/pca/APS_GEO_pruned_SSI_filtered_snps_passed.Q20.MAF0.05.MAXMISS0.9.feb.LDpruned.eigenvec",header=F,sep="")


table1$Group<-c("A","A","A","A","A","P","P","P","P","P","P","P","P","P","P","GEO","GEO","GEO","GEO","GEO","GEO","GEO","GEO","GEO","GEO","GEO","GEO","GEO","GEO","GEO","GEO","GEO","A","A","A","P","S","S","A","P","P","P","P","P","S","S","S","S","A","A","A","A","P","P","A","A","A","A","A","A","A","A","A","A","A","A","A","A","P","P","P","P","P","P","P","P","P","P","P","S","S","S","S","S","S","S","S","S","S","S","A","A","P","P","P","P","P","S","S","S","S","S","A","A","A")
#table1$Group<-as.factor(table1$Group)

eigenval3<-read.table("~/Documents/Research/Small_jawed_project/pca/APS_GEO_pruned_SSI_filtered_snps_passed.Q20.MAF0.05.MAXMISS0.9.feb.LDpruned.eigenval",header=F,sep="",stringsAsFactors = F)
sum_eig3<-sum(eigenval3$V1)

table1$Group[which(table1$V1  == "OYSP1")]<-"S"
table1$Group[which(table1$V1  == "OYSP3")]<-"S"
table1$Group[which(table1$V1  == "OYSP7")]<-"S"
table1$Group[which(table1$V1  == "OYSP4")]<-"S"
table1$Group[which(table1$V1  == "MERP1")]<-"S"
table1$Group[which(table1$V1  == "MERP2")]<-"S"
table1$Group[which(table1$V1  == "OSPP2")]<-"S"
#D95F02

#group.colors<-c(P="aquamarine3",A="goldenrod1",S="coral1",GEO="gray")

group.colors<-c(P="aquamarine3",A="goldenrod1",S="#D95F02",GEO="gray")

sum_eigs3<-lapply(eigenval3$V1,function(x){
  rt<-(x/sum_eig3)*100
  rt<-round(rt,digits=1)
  return(rt)
})


ggplot(table1,aes(x=table1$V3,y=table1$V4,label=table1$V1,color=table1$Group))+
  geom_point(size=0)+geom_text(aes(label=table1$V1),size=3,hjust=0, vjust=0)+
  scale_color_manual(values=group.colors)+
  labs(x=paste0("PC1: ",sum_eigs3[[1]],"% variance"),y=paste0("PC2: ",sum_eigs3[[2]],"% variance"),color="locality")+guides(size=F)+
  theme_bw()



ggplot(table1,aes(x=table1$V3,y=table1$V4,color=table1$Group))+
  geom_point(size=3)+
  scale_color_manual(values=group.colors)+
  labs(x=paste0("PC1: ",sum_eigs3[[1]],"% variance"),y=paste0("PC2: ",sum_eigs3[[2]],"% variance"),color="locality")+guides(size=F)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.text = element_text(size=12),
        axis.title =  element_text(size=16),
        legend.position = "none",
        panel.border = element_rect(fill=NA, colour = "black", size=1.2))




table1<-read.table("~/Documents/Research/Small_jawed_project/pca/APS_pruned_SSI_filtered_snps_passed.Q20.MAF0.05.MAXMISS0.9.feb.LDpruned.eigenvec",header=F,sep="")


table1$Group<-c("A","A","A","A","A","P","P","P","P","P","P","P","P","P","P","A","A","A","P","S","S","A","P","P","P","P","P","S","S","S","S","A","A","A","A","P","P","A","A","A","A","A","A","A","A","A","A","A","A","A","A","P","P","P","P","P","P","P","P","P","P","P","S","S","S","S","S","S","S","S","S","S","S","A","A","P","P","P","P","P","S","S","S","S","S","A","A","A")
table1$Group<-as.factor(table1$Group)

eigenval3<-read.table("~/Documents/Research/Small_jawed_project/pca/APS_pruned_SSI_filtered_snps_passed.Q20.MAF0.05.MAXMISS0.9.feb.LDpruned.eigenval",header=F,sep="",stringsAsFactors = F)
sum_eig3<-sum(eigenval3$V1)

table1$Group[which(table1$V1  == "OYSP1")]<-"S"
table1$Group[which(table1$V1  == "OYSP3")]<-"S"
table1$Group[which(table1$V1  == "OYSP7")]<-"S"
table1$Group[which(table1$V1  == "OYSP4")]<-"S"
table1$Group[which(table1$V1  == "MERP1")]<-"S"
table1$Group[which(table1$V1  == "MERP2")]<-"S"
table1$Group[which(table1$V1  == "OSPP2")]<-"S"


group.colors<-c(P="aquamarine3",A="goldenrod1",S="coral1")
group.colors<-c(P="aquamarine3",A="goldenrod1",S="#D95F02")


sum_eigs3<-lapply(eigenval3$V1,function(x){
  rt<-(x/sum_eig3)*100
  rt<-round(rt,digits=1)
  return(rt)
})


ggplot(table1,aes(x=table1$V3,y=table1$V4,label=table1$V1,color=table1$Group))+
  geom_point(size=0)+geom_text(aes(label=table1$V1),size=3,hjust=0, vjust=0)+
  scale_color_manual(values=group.colors)+
  labs(x=paste0("PC1: ",sum_eigs3[[1]],"% variance"),y=paste0("PC2: ",sum_eigs3[[2]],"% variance"),color="locality")+guides(size=F)+
  theme_bw()



ggplot(table1,aes(x=table1$V3,y=table1$V4,color=table1$Group))+
  geom_point(size=3)+
  scale_color_manual(values=group.colors)+
  labs(x=paste0("PC1: ",sum_eigs3[[1]],"% variance"),y=paste0("PC2: ",sum_eigs3[[2]],"% variance"),color="locality")+guides(size=F)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.text = element_text(size=12),
        panel.border = element_rect(fill=NA, colour = "black", size=1.2),
        legend.position = "none")









###### ME2P1 #########

table1<-read.table("~/Documents/Research/Small_jawed_project/pca/AMPS_pruned_SSI_filtered_snps_passed.Q20.MAF0.05.MAXMISS0.9.feb.LDpruned2.eigenvec",header=F,sep="")


#table1$Group<-c("A","A","A","A","A","M","M","M","M","M","M","M","M","M","M","M","M","P","P","P","P","P","P","P","P","P","P","GEO","GEO","GEO","GEO","GEO","GEO","GEO","GEO","GEO","GEO","GEO","GEO","GEO","GEO","GEO","GEO","GEO","A","A","A","P","S","S","A","M","M","M","M","M","P","P","P","P","P","S","S","S","S","A","A","A","A","P","P","A","M","M","M","M","M","M","A","A","A","A","A","A","A","A","A","A","A","A","A","M","M","M","M","M","M","M","M","M","M","M","M","M","P","P","P","P","P","P","P","P","P","P","P","S","S","S","S","S","S","S","S","S","S","S","A","A","M","M","M","M","M","M","M","M","P","P","P","P","P","S","S","S","S","S","A","A","A")
#table1$Group<-as.factor(table1$Group)

eigenval3<-read.table("~/Documents/Research/Small_jawed_project/pca/AMPS_pruned_SSI_filtered_snps_passed.Q20.MAF0.05.MAXMISS0.9.feb.LDpruned2.eigenval",header=F,sep="",stringsAsFactors = F)
sum_eig3<-sum(eigenval3$V1)

group.colors<-c(M="mediumorchid3",P="aquamarine3",A="goldenrod1",S="coral1",GEO="gray")


sum_eigs3<-lapply(eigenval3$V1,function(x){
  rt<-(x/sum_eig3)*100
  rt<-round(rt,digits=1)
  return(rt)
})
table1$color<-'NO'
table1$color[which(table1$V1  == "ME2P1")]<-'YES'
table1$color[which(table1$V1  == "ME2P1-1")]<-'YES'


ggplot(table1,aes(x=table1$V3,y=table1$V4,label=table1$V1))+
  geom_point(size=0)+geom_text(aes(label=table1$V1),size=3,hjust=0, vjust=0)+
  scale_color_manual(values=group.colors)+
  labs(x=paste0("PC1: ",sum_eigs3[[1]],"% variance"),y=paste0("PC2: ",sum_eigs3[[2]],"% variance"),color="locality")+guides(size=F)+
  theme_bw()



ggplot(table1,aes(x=table1$V3,y=table1$V4,col=table1$color))+
  geom_point(size=1)+
  #geom_text(aes(label=table1$V1),subset =.(table1$V1 == "ME2P1"))+
  labs(x=paste0("PC1: ",sum_eigs3[[1]],"% variance"),y=paste0("PC2: ",sum_eigs3[[2]],"% variance"),color="locality")+guides(size=F)+
  theme_bw()
