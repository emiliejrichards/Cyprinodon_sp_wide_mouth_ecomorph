library(tidyverse)
library(RColorBrewer)

######### K CV plot ############

AM.K=read.table("~/Documents/Research/Small_jawed_project/admixture/APS_CV_admixture.txt",sep="")
plot(AM.K$V1,AM.K$V2,xlab="K",ylab="Cross-validation error",main="APS")
lines(AM.K$V1,AM.K$V2,col="red",lty=2)

AM.K=read.table("~/Documents/Research/Small_jawed_project/admixture/AMPS_CV_admixture.txt",sep="")
plot(AM.K$V1,AM.K$V2,xlab="K",ylab="Cross-validation error",main="AMPS")
lines(AM.K$V1,AM.K$V2,col="red",lty=2)
AM.K=read.table("~/Documents/Research/Small_jawed_project/admixture/APS_GEO_CV_admixture.txt",sep="")
plot(AM.K$V1,AM.K$V2,xlab="K",ylab="Cross-validation error",main="APS")
lines(AM.K$V1,AM.K$V2,col="red",lty=2)


#####K of 4 ########
twolabset<-data.frame(species=c("Generalist","Generalist","Generalist","Generalist","Generalist","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Generalist","Generalist","Generalist","Scale-eater","Small jawed","Small jawed","Generalist","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Generalist","Generalist","Generalist","Generalist","Small jawed","Small jawed","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Small jawed","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Generalist","Generalist","Small jawed","Small jawed","Small jawed","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Generalist","Generalist","Generalist"))
twolabset$species<-as.character(twolabset$species)


AMQ<-read.table("~/Documents/Research/Small_jawed_project/admixture/APS_pruned_SSI_filtered_snps_passed.Q20.MAF0.05.MAXMISS0.9.feb.LDpruned.4.Q", header=F)
AMinds<-read.table("~/Documents/Research/Small_jawed_project/admixture/APS_inds.txt")

AMQ$inds<-as.character(AMinds$V1)
AMQ$species<-as.character(twolabset$species)
#AMQ$locality<-as.character(twolabset$locality)

#AMQ.2<-AMQ[which(AMQ$species!="ART"),]
#AMQ.2<-AMQ[which(AMQ$inds!="CRPA1000" & AMQ$inds!="CRPA3" & AMQ$inds!="CRPA1" & AMQ$inds!="OSPH1" & AMQ$inds!="KILA1" & AMQ$inds!="CUNA7" & AMQ$inds!="NCCA11" & AMQ$inds!="OSPP8" & AMQ$inds!="OSPP4" & AMQ$inds!="BAVA10" & AMQ$inds!="CRPA1001" & AMQ$inds!="SPPA1"),]
#AMQ.2$species<-factor(AMQ.2$species,levels=c("M","P","A","GEO","CUN","BAV","VEN","NCC"))
AMQ$species<-factor(AMQ$species,levels=c("Generalist","Small jawed","Scale-eater"))

plot_data <- AMQ %>%
  mutate(id=inds) %>%
  gather('pop','prob',V1:V4) %>%
  group_by(inds) %>%
  mutate(likely_assignment=pop[which.max(prob)],
         assignment_prob=max(prob)) %>%
  arrange(likely_assignment,desc(assignment_prob))%>%
  ungroup()%>%
  mutate(id=forcats::fct_inorder(factor(id)))

AMQ$species[which(AMQ$inds=="GREP1")]<-'Small jawed'


#V1(BAV),V10(CUN),V2(M.1),V4(M.2),V3(VEN),V5(NCC),V6(P.1),V7(A),V8(GEO),V9(P.2),

group.colors<-c(V2="coral1", V1="goldenrod1",V3="aquamarine3",V4="blue")

ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col()+
  facet_grid(~species, scales = 'free', space = 'free')+
  guides(fill=FALSE)+ylab("Ancestry Proportion")+
  scale_fill_manual(values=group.colors)+
  #scale_fill_brewer(palette="Paired")+
  theme_classic()+
 # theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size = 6))


#####K of 5 ########
twolabset<-data.frame(species=c("Generalist","Generalist","Generalist","Generalist","Generalist","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Generalist","Generalist","Generalist","Scale-eater","Small jawed","Small jawed","Generalist","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Generalist","Generalist","Generalist","Generalist","Small jawed","Small jawed","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Small jawed","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Generalist","Generalist","Small jawed","Small jawed","Small jawed","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Generalist","Generalist","Generalist"))
twolabset$species<-as.character(twolabset$species)


AMQ<-read.table("~/Documents/Research/Small_jawed_project/admixture/APS_pruned_SSI_filtered_snps_passed.Q20.MAF0.05.MAXMISS0.9.feb.LDpruned.5.Q", header=F)
AMinds<-read.table("~/Documents/Research/Small_jawed_project/admixture/APS_inds.txt")

AMQ$inds<-as.character(AMinds$V1)
AMQ$species<-as.character(twolabset$species)
#AMQ$locality<-as.character(twolabset$locality)

#AMQ.2<-AMQ[which(AMQ$species!="ART"),]
#AMQ.2<-AMQ[which(AMQ$inds!="CRPA1000" & AMQ$inds!="CRPA3" & AMQ$inds!="CRPA1" & AMQ$inds!="OSPH1" & AMQ$inds!="KILA1" & AMQ$inds!="CUNA7" & AMQ$inds!="NCCA11" & AMQ$inds!="OSPP8" & AMQ$inds!="OSPP4" & AMQ$inds!="BAVA10" & AMQ$inds!="CRPA1001" & AMQ$inds!="SPPA1"),]
#AMQ.2$species<-factor(AMQ.2$species,levels=c("M","P","A","GEO","CUN","BAV","VEN","NCC"))
AMQ$species<-factor(AMQ$species,levels=c("Generalist","Small jawed","Scale-eater"))

AMQ$species[which(AMQ$inds=="GREP1")]<-'Small jawed'
plot_data <- AMQ %>%
  mutate(id=inds) %>%
  gather('pop','prob',V1:V5) %>%
  group_by(inds) %>%
  mutate(likely_assignment=pop[which.max(prob)],
         assignment_prob=max(prob)) %>%
  arrange(likely_assignment,desc(assignment_prob))%>%
  ungroup()%>%
  mutate(id=forcats::fct_inorder(factor(id)))



#V1(BAV),V10(CUN),V2(M.1),V4(M.2),V3(VEN),V5(NCC),V6(P.1),V7(A),V8(GEO),V9(P.2),

group.colors<-c(V2="coral1", V1="goldenrod1",V4="aquamarine3",V3="blue", V5="darkorange")

ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col()+
  facet_grid(~species, scales = 'free', space = 'free')+
  guides(fill=FALSE)+ylab("Ancestry Proportion")+
  scale_fill_manual(values=group.colors)+
  #scale_fill_brewer(palette="Paired")+
  theme_classic()+
  # theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size = 6))


###### K of 3 #########


AMQ<-read.table("~/Documents/Research/Small_jawed_project/admixture/APS_pruned_SSI_filtered_snps_passed.Q20.MAF0.05.MAXMISS0.9.feb.LDpruned.3.Q", header=F)
AMinds<-read.table("~/Documents/Research/Small_jawed_project/admixture/APS_inds.txt")

AMQ$inds<-as.character(AMinds$V1)
AMQ$species<-as.character(twolabset$species)
#AMQ$locality<-as.character(twolabset$locality)

#AMQ.2<-AMQ[which(AMQ$species!="ART"),]
#AMQ.2<-AMQ[which(AMQ$inds!="CRPA1000" & AMQ$inds!="CRPA3" & AMQ$inds!="CRPA1" & AMQ$inds!="OSPH1" & AMQ$inds!="KILA1" & AMQ$inds!="CUNA7" & AMQ$inds!="NCCA11" & AMQ$inds!="OSPP8" & AMQ$inds!="OSPP4" & AMQ$inds!="BAVA10" & AMQ$inds!="CRPA1001" & AMQ$inds!="SPPA1"),]
#AMQ.2$species<-factor(AMQ.2$species,levels=c("M","P","A","GEO","CUN","BAV","VEN","NCC"))
AMQ$species<-factor(AMQ$species,levels=c("Generalist","Small jawed","Scale-eater"))
AMQ$species[which(AMQ$inds=="GREP1")]<-'Small jawed'

plot_data <- AMQ %>%
  mutate(id=inds) %>%
  gather('pop','prob',V1:V3) %>%
  group_by(inds) %>%
  mutate(likely_assignment=pop[which.max(prob)],
         assignment_prob=max(prob)) %>%
  arrange(likely_assignment,desc(assignment_prob))%>%
  ungroup()%>%
  mutate(id=forcats::fct_inorder(factor(id)))



#V1(BAV),V10(CUN),V2(M.1),V4(M.2),V3(VEN),V5(NCC),V6(P.1),V7(A),V8(GEO),V9(P.2),

group.colors<-c(V2="coral1", V3="goldenrod1",V1="aquamarine3")

ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col()+
  facet_grid(~species, scales = 'free', space = 'free')+
  guides(fill=FALSE)+ylab("Ancestry Proportion")+
  scale_fill_manual(values=group.colors)+
  #scale_fill_brewer(palette="Paired")+
  theme_classic()+
  # theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size = 6))

######## APS GEO K of 4 ##################


twolabset<-data.frame(species=c("Generalists","Generalists","Generalists","Generalists","Generalists","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Generalists","Generalists","Generalists","Small jawed","Small jawed","Small jawed","Generalists","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Generalists","Generalists","Generalists","Generalists","Small jawed","Small jawed","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Small jawed","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Generalists","Generalists","Small jawed","Small jawed","Small jawed","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Generalists","Generalists","Generalists"))
twolabset$species<-as.character(twolabset$species)


AMQ<-read.table("~/Documents/Research/Small_jawed_project/admixture/APS_GEO_pruned_SSI_filtered_snps_passed.Q20.MAF0.05.MAXMISS0.9.feb.LDpruned.4.Q", header=F)
AMinds<-read.table("~/Documents/Research/Small_jawed_project/admixture/APS_GEO_inds.txt")

AMQ$inds<-as.character(AMinds$V1)
AMQ$species<-as.character(twolabset$species)
#AMQ$locality<-as.character(twolabset$locality)

#AMQ.2<-AMQ[which(AMQ$species!="ART"),]
#AMQ.2<-AMQ[which(AMQ$inds!="CRPA1000" & AMQ$inds!="CRPA3" & AMQ$inds!="CRPA1" & AMQ$inds!="OSPH1" & AMQ$inds!="KILA1" & AMQ$inds!="CUNA7" & AMQ$inds!="NCCA11" & AMQ$inds!="OSPP8" & AMQ$inds!="OSPP4" & AMQ$inds!="BAVA10" & AMQ$inds!="CRPA1001" & AMQ$inds!="SPPA1"),]
#AMQ.2$species<-factor(AMQ.2$species,levels=c("M","P","A","GEO","CUN","BAV","VEN","NCC"))
AMQ$species<-factor(AMQ$species,levels=c("Rum Cay","Generalists","Small jawed","Scale-eater"))

plot_data <- AMQ %>%
  mutate(id=inds) %>%
  gather('pop','prob',V1:V4) %>%
  group_by(inds) %>%
  mutate(likely_assignment=pop[which.max(prob)],
         assignment_prob=max(prob)) %>%
  arrange(likely_assignment,desc(assignment_prob))%>%
  ungroup()%>%
  mutate(id=forcats::fct_inorder(factor(id)))

#AMQ$species[which(AMQ$inds=="GREP1")]<-'Small jawed'


#V1(BAV),V10(CUN),V2(M.1),V4(M.2),V3(VEN),V5(NCC),V6(P.1),V7(A),V8(GEO),V9(P.2),

group.colors<-c(V2="goldenrod1", V1="aquamarine3",V3="gray",V4="coral1")

ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col()+
  facet_grid(~species, scales = 'free', space = 'free')+
  guides(fill=FALSE)+ylab("Ancestry Proportion")+
  scale_fill_manual(values=group.colors)+
  #scale_fill_brewer(palette="Paired")+
  theme_classic()+
  # theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size = 6))


####### APS GEO K of 5 #########

#twolabset<-data.frame(species=c("Generalists","Generalists","Generalists","Generalists","Generalists","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Generalists","Generalists","Generalists","Small jawed","Small jawed","Small jawed","Generalists","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Generalists","Generalists","Generalists","Generalists","Small jawed","Small jawed","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Small jawed","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Generalists","Generalists","Small jawed","Small jawed","Small jawed","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Generalists","Generalists","Generalists"))
#twolabset$species<-as.character(twolabset$species)


AMQ<-read.table("~/Documents/Research/Small_jawed_project/admixture/APS_GEO_pruned_SSI_filtered_snps_passed.Q20.MAF0.05.MAXMISS0.9.feb.LDpruned.5.Q", header=F)
AMinds<-read.table("~/Documents/Research/Small_jawed_project/admixture/APS_GEO_inds.txt")

AMQ$inds<-as.character(AMinds$V1)
AMQ$species<-as.character(twolabset$species)
#AMQ$locality<-as.character(twolabset$locality)

#AMQ.2<-AMQ[which(AMQ$species!="ART"),]
#AMQ.2<-AMQ[which(AMQ$inds!="CRPA1000" & AMQ$inds!="CRPA3" & AMQ$inds!="CRPA1" & AMQ$inds!="OSPH1" & AMQ$inds!="KILA1" & AMQ$inds!="CUNA7" & AMQ$inds!="NCCA11" & AMQ$inds!="OSPP8" & AMQ$inds!="OSPP4" & AMQ$inds!="BAVA10" & AMQ$inds!="CRPA1001" & AMQ$inds!="SPPA1"),]
#AMQ.2$species<-factor(AMQ.2$species,levels=c("M","P","A","GEO","CUN","BAV","VEN","NCC"))
AMQ$species<-factor(AMQ$species,levels=c("Rum Cay","Generalists","Small jawed","Scale-eater"))

plot_data <- AMQ %>%
  mutate(id=inds) %>%
  gather('pop','prob',V1:V5) %>%
  group_by(inds) %>%
  mutate(likely_assignment=pop[which.max(prob)],
         assignment_prob=max(prob)) %>%
  arrange(likely_assignment,desc(assignment_prob))%>%
  ungroup()%>%
  mutate(id=forcats::fct_inorder(factor(id)))

#AMQ$species[which(AMQ$inds=="GREP1")]<-'Small jawed'


#V1(BAV),V10(CUN),V2(M.1),V4(M.2),V3(VEN),V5(NCC),V6(P.1),V7(A),V8(GEO),V9(P.2),

group.colors<-c(V2="gray", V1="goldenrod1",V3="aquamarine3",V4="blue",V5="coral1")

ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col()+
  facet_grid(~species, scales = 'free', space = 'free')+
  guides(fill=FALSE)+ylab("Ancestry Proportion")+
  scale_fill_manual(values=group.colors)+
  #scale_fill_brewer(palette="Paired")+
  theme_classic()+
  # theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size = 6))





####### AMPS K of 6 #########

twolabset<-data.frame(species=c("Generalists","Generalists","Generalists","Generalists","Generalists","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Generalists","Generalists","Generalists","Small jawed","Small jawed","Small jawed","Generalists","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Generalists","Generalists","Generalists","Generalists","Scale-eater","Scale-eater","Generalists","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Generalists","Generalists","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Small jawed","Small jawed","Small jawed","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Generalists","Generalists","Generalists"))
twolabset$species<-as.character(twolabset$species)


AMQ<-read.table("~/Documents/Research/Small_jawed_project/admixture/AMPS_pruned_SSI_filtered_snps_passed.Q20.MAF0.05.MAXMISS0.9.feb.LDpruned.6.Q", header=F)
AMinds<-read.table("~/Documents/Research/Small_jawed_project/admixture/AMPS_inds.txt")

AMQ$inds<-as.character(AMinds$V1)
AMQ$species<-as.character(twolabset$species)
#AMQ$locality<-as.character(twolabset$locality)

#AMQ.2<-AMQ[which(AMQ$species!="ART"),]
#AMQ.2<-AMQ[which(AMQ$inds!="CRPA1000" & AMQ$inds!="CRPA3" & AMQ$inds!="CRPA1" & AMQ$inds!="OSPH1" & AMQ$inds!="KILA1" & AMQ$inds!="CUNA7" & AMQ$inds!="NCCA11" & AMQ$inds!="OSPP8" & AMQ$inds!="OSPP4" & AMQ$inds!="BAVA10" & AMQ$inds!="CRPA1001" & AMQ$inds!="SPPA1"),]
#AMQ.2$species<-factor(AMQ.2$species,levels=c("M","P","A","GEO","CUN","BAV","VEN","NCC"))
AMQ$species<-factor(AMQ$species,levels=c("Rum Cay","Generalists","Small jawed","Scale-eater","Molluscivore"))


AMQ$species[which(AMQ$inds=="MERP1")]<-'Small jawed'
AMQ$species[which(AMQ$inds=="MERP2")]<-'Small jawed'
AMQ$species[which(AMQ$inds=="OSPP2")]<-'Small jawed'


plot_data <- AMQ %>%
  mutate(id=inds) %>%
  gather('pop','prob',V1:V6) %>%
  group_by(inds) %>%
  mutate(likely_assignment=pop[which.max(prob)],
         assignment_prob=max(prob)) %>%
  arrange(likely_assignment,desc(assignment_prob))%>%
  ungroup()%>%
  mutate(id=forcats::fct_inorder(factor(id)))



#V1(BAV),V10(CUN),V2(M.1),V4(M.2),V3(VEN),V5(NCC),V6(P.1),V7(A),V8(GEO),V9(P.2),

group.colors<-c(V2="plum2", V1="purple",V3="goldenrod1",V4="coral1",V5="mediumorchid3",V6="aquamarine3")

ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col()+
  facet_grid(~species, scales = 'free', space = 'free')+
  guides(fill=FALSE)+ylab("Ancestry Proportion")+
  scale_fill_manual(values=group.colors)+
  #scale_fill_brewer(palette="Paired")+
  theme_classic()+
  # theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size = 6))



####### AMPS K of 7 #########

twolabset<-data.frame(species=c("Generalists","Generalists","Generalists","Generalists","Generalists","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Generalists","Generalists","Generalists","Small jawed","Small jawed","Small jawed","Generalists","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Generalists","Generalists","Generalists","Generalists","Scale-eater","Scale-eater","Generalists","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Generalists","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Generalists","Generalists","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Small jawed","Small jawed","Small jawed","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Generalists","Generalists","Generalists"))
twolabset$species<-as.character(twolabset$species)


AMQ<-read.table("~/Documents/Research/Small_jawed_project/admixture/AMPS_pruned_SSI_filtered_snps_passed.Q20.MAF0.05.MAXMISS0.9.feb.LDpruned.7.Q", header=F)
AMinds<-read.table("~/Documents/Research/Small_jawed_project/admixture/AMPS_inds.txt")

AMQ$inds<-as.character(AMinds$V1)
AMQ$species<-as.character(twolabset$species)
#AMQ$locality<-as.character(twolabset$locality)

#AMQ.2<-AMQ[which(AMQ$species!="ART"),]
#AMQ.2<-AMQ[which(AMQ$inds!="CRPA1000" & AMQ$inds!="CRPA3" & AMQ$inds!="CRPA1" & AMQ$inds!="OSPH1" & AMQ$inds!="KILA1" & AMQ$inds!="CUNA7" & AMQ$inds!="NCCA11" & AMQ$inds!="OSPP8" & AMQ$inds!="OSPP4" & AMQ$inds!="BAVA10" & AMQ$inds!="CRPA1001" & AMQ$inds!="SPPA1"),]
#AMQ.2$species<-factor(AMQ.2$species,levels=c("M","P","A","GEO","CUN","BAV","VEN","NCC"))
AMQ$species<-factor(AMQ$species,levels=c("Generalists","Molluscivore","Small jawed","Scale-eater"))


AMQ$species[which(AMQ$inds=="MERP1")]<-'Small jawed'
AMQ$species[which(AMQ$inds=="MERP2")]<-'Small jawed'
AMQ$species[which(AMQ$inds=="OSPP2")]<-'Small jawed'
AMQ <- AMQ %>%filter(inds!='CRPA1003')%>%filter(inds!='PAIA1')


plot_data <- AMQ %>%
  mutate(id=inds) %>%
  gather('pop','prob',V1:V7) %>%
  group_by(inds) %>%
  mutate(likely_assignment=pop[which.max(prob)],
         assignment_prob=max(prob)) %>%
  arrange(likely_assignment,desc(assignment_prob))%>%
  ungroup()%>%
  mutate(id=forcats::fct_inorder(factor(id)))

#plot_data <- plot_data %>%filter(inds=='CRPA1003')%>%filter(inds=='PAIA1')

#V1(BAV),V10(CUN),V2(M.1),V4(M.2),V3(VEN),V5(NCC),V6(P.1),V7(A),V8(GEO),V9(P.2),

group.colors<-c(V2="aquamarine3", V1="purple",V3="goldenrod1",V4="mediumorchid3",V5="coral1",V6="plum2",V7="blue")

ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col()+
  facet_grid(~species, scales = 'free', space = 'free')+
  guides(fill=FALSE)+ylab("Ancestry Proportion")+
  scale_fill_manual(values=group.colors)+
  #scale_fill_brewer(palette="Paired")+
  theme_classic()+
  # theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size = 6))






##### AMPS GEO 

twolabset<-data.frame(species=c("Generalist","Generalist","Generalist","Generalist","Generalist","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Rum Cay","Generalist","Generalist","Generalist","Small jawed","Small jawed","Small jawed","Generalist","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Generalist","Generalist","Generalist","Generalist","Small jawed","Small jawed","Generalist","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Generalist","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Small jawed","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Generalist","Generalist","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Molluscivore","Small jawed","Small jawed","Small jawed","Scale-eater","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Small jawed","Generalist","Generalist","Generalist"))
twolabset$species<-as.character(twolabset$species)


AMQ<-read.table("~/Documents/Research/Small_jawed_project/admixture/AMPS_GEO_pruned_SSI_filtered_snps_passed.Q20.MAF0.05.MAXMISS0.9.feb.LDpruned.8.Q", header=F)
AMinds<-read.table("~/Documents/Research/Small_jawed_project/admixture/AMPS_GEO_inds.txt")

AMQ$inds<-as.character(AMinds$V1)
AMQ$species<-as.character(twolabset$species)
#AMQ$locality<-as.character(twolabset$locality)

#AMQ.2<-AMQ[which(AMQ$species!="ART"),]
#AMQ.2<-AMQ[which(AMQ$inds!="CRPA1000" & AMQ$inds!="CRPA3" & AMQ$inds!="CRPA1" & AMQ$inds!="OSPH1" & AMQ$inds!="KILA1" & AMQ$inds!="CUNA7" & AMQ$inds!="NCCA11" & AMQ$inds!="OSPP8" & AMQ$inds!="OSPP4" & AMQ$inds!="BAVA10" & AMQ$inds!="CRPA1001" & AMQ$inds!="SPPA1"),]
#AMQ.2$species<-factor(AMQ.2$species,levels=c("M","P","A","GEO","CUN","BAV","VEN","NCC"))
AMQ$species<-factor(AMQ$species,levels=c("Rum Cay","Generalist","Molluscivore","Small jawed","Scale-eater"))


AMQ$species[which(AMQ$inds=="MERP1")]<-'Small jawed'
AMQ$species[which(AMQ$inds=="MERP2")]<-'Small jawed'
AMQ$species[which(AMQ$inds=="OSPP2")]<-'Small jawed'


plot_data <- AMQ %>%
  mutate(id=inds) %>%
  gather('pop','prob',V1:V8) %>%
  group_by(inds) %>%
  mutate(likely_assignment=pop[which.max(prob)],
         assignment_prob=max(prob)) %>%
  arrange(likely_assignment,desc(assignment_prob))%>%
  ungroup()%>%
  mutate(id=forcats::fct_inorder(factor(id)))



#V1(BAV),V10(CUN),V2(M.1),V4(M.2),V3(VEN),V5(NCC),V6(P.1),V7(A),V8(GEO),V9(P.2),

group.colors<-c(V2="goldenrod1", V1="purple",V3="coral1",V4="blue",V5="plum1",V6="gray",V7="aquamarine3",V8="mediumorchid3")

ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col()+
  facet_grid(~species, scales = 'free', space = 'free')+
  guides(fill=FALSE)+ylab("Ancestry Proportion")+
  scale_fill_manual(values=group.colors)+
  #scale_fill_brewer(palette="Paired")+
  theme_classic()+
  # theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size = 6))


