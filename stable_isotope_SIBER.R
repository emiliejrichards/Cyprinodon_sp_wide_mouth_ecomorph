install.packages("SIBER")
library(SIBER)


set.seed(1)

# load in the included demonstration dataset
data("demo.siber.data")
#
# create the siber object
siber.example <- createSiberObject(demo.siber.data)


# Or if working with your own data read in from a *.csv file, you would use
# This *.csv file is included with this package. To find its location
# type
# fname <- system.file("extdata", "demo.siber.data.csv", package = "SIBER")
# in your command window. You could load it directly by using the
# returned path, or perhaps better, you could navigate to this folder
# and copy this file to a folder of your own choice, and create a 
# script from this vingette to analyse it. This *.csv file provides
# a template for how your own files should be formatted.

mydata <- read.csv(fname, header=T)


isotope.both.pups.osp<-read.csv("~/Documents/Research/Small_jawed_project/isotope_data/pupfish_OSP_isotope_data_reformatted.csv",
                                 header=T)


isotope.both.pups.osp.simple<-
  isotope.both.pups.osp%>%
  filter(year=="2018")%>%
  filter(species!="Silverside")%>%
  filter(ID!="OSPA0220")%>%filter(ID!="OSPP0217" )%>%filter(ID!="OSPP0042" )%>%
  dplyr::select(species,population,d13C,d15N)

m1<-lm(d15N ~ species,data=isotope.both.pups.osp.simple)
summary(m1)
summary(aov(d15N ~ species,data=isotope.both.pups.osp.simple))
TukeyHSD(aov(d15N ~ species,data=isotope.both.pups.osp.simple))
isotope.both.pups.osp.simple%>%group_by(species)%>%dplyr::summarise(mean.d15N=mean(d15N),
                                                                    min.d15N=min(d15N),
                                                                    max.d15N=max(d15N))


isotope.both.pups.osp.simple<-
  isotope.both.pups.osp%>%
  filter(year=="2018")%>%
  filter(species!="Silverside")%>%
  filter(ID!="OSPA0220")%>%filter(ID!="OSPP0217" )%>%filter(ID!="OSPP0042" )%>%
  dplyr::select(species,population,d13C,d15N)


isotope.both.pups.osp.simple %>% filter(species=="P") %>% nrow()
isotope.both.pups.osp.simple %>% filter(species=="S") %>% nrow()
isotope.both.pups.osp.simple %>% filter(species=="A") %>% nrow()


m1<-lm(d15N ~ species,data=isotope.both.pups.osp.simple)
summary(m1)
summary(aov(d15N ~ species,data=isotope.both.pups.osp.simple))
TukeyHSD(aov(d15N ~ species,data=isotope.both.pups.osp.simple))
isotope.both.pups.osp.simple%>%group_by(species)%>%dplyr::summarise(mean.d15N=mean(d15N),
                                                                    min.d15N=min(d15N),
                                                                    max.d15N=max(d15N))

summary(aov(d13C~ species,data=isotope.both.pups.osp.simple))
TukeyHSD(aov(d13C~ species,data=isotope.both.pups.osp.simple))

m2<-lm(d13C ~ species,data=isotope.both.pups.osp.simple)

anova_stats(m2)


df<-isotope.both.pups.osp.simple %>% filter(species=="S") 
quantile(df$d15N,c(0.01,0.20)) 


library(sjstats)
effectsize::eta_squared(m1)
#So if you end up with η² = 0.45, you can assume the effect size is very large. 
#It also means that 45% of the change in the DV can be accounted for by the IV.
effectsize::omega_squared(m1)
effectsize::cohens_f(m1)

anova_stats(m1)

df<-isotope.both.pups.osp.simple %>% filter(species!="A") 
t.test(iso2 ~ species,data=df)
var.test(d15N ~ species,data=df)


df2<-isotope.both.pups.osp.simple %>% filter(species!="S") 

var.test(d15N ~ species,data=df2)


df3<-isotope.both.pups.osp.simple %>% filter(species!="P") 
var.test(d15N ~ species,data=df3)


isotope.both.pups.osp.simple %>% filter(species=="P") %>% nrow()
isotope.both.pups.osp.simple %>% filter(species=="S") %>% nrow()
isotope.both.pups.osp.simple %>% filter(species=="A") %>% nrow()

new_df <- df %>% group_by(species) %>% slice_sample(n=14)
var.test(d15N ~ species,data=new_df)
new_df2 <- df2 %>% group_by(species) %>% slice_sample(n=14)
var.test(d15N ~ species,data=new_df2)
new_df3 <- df3 %>% group_by(species) %>% slice_sample(n=14)
test<-var.test(d15N ~ species,data=new_df3)


AvS.pvalue<-vector()
for(i in 1:1000){
  new_df3 <- df3 %>% group_by(species) %>% slice_sample(n=14)
  test<-var.test(d15N ~ species,data=new_df3)
  AvS.pvalue[i]<-test$p.value
}
quantile(AvS.pvalue,c(0.025,0.05,0.95))
hist(AvS.pvalue)
sum(AvS.pvalue<0.05)/1000
var.test(d15N ~ species,data=df3)

shapiro.test(df2$d15N)
bartlett.test(d15N~species,data=isotope.both.pups.osp.simple)
library(car)
leveneTest(d15N~species,data=isotope.both.pups.osp.simple)
shapiro.test(df3$d15N)
fligner.test(d15N~species,data=isotope.both.pups.osp.simple)

ggplot(data=isotope.both.pups.osp.simple,aes(x=species,y=d15N))+
  geom_jitter()+
  geom_boxplot()+theme_classic()

group.colors2<-c(A="goldenrod1",P="aquamarine3",S="#D95F02")
isotope.both.pups.osp.simple$species<-factor(isotope.both.pups.osp.simple$species,levels=c('A','S','P'))

d15N.boxplot<-ggplot(data=isotope.both.pups.osp.simple,aes(x=species,y=d15N,fill=species))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(col="dimgrey",width = 0.25)+
  scale_color_manual(values=group.colors2)+
  scale_fill_manual(values=group.colors2)+
  scale_x_discrete(breaks=c("A", "S", "P"),
                  labels=c("var", "wid", "des"))+
  theme_classic()+theme(legend.position = 'none')


d13C.boxplot<-ggplot(data=isotope.both.pups.osp.simple,aes(x=species,y=d13C,fill=species))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(col="dimgrey",width = 0.25)+
  scale_color_manual(values=group.colors2)+
  scale_fill_manual(values=group.colors2)+
  scale_x_discrete(breaks=c("A", "S", "P"),
                   labels=c("var", "wid", "des"))+
  theme_classic()+theme(legend.position = 'none')


pdf(file="~/Documents/Research/Small_jawed_project/isotope_data/d15N_d13C_panels.pdf")

grid.arrange(d15N.boxplot,d13C.boxplot ,ncol=2)
dev.off()

group.colors


TukeyHSD(aov(d15N ~ species,data=isotope.both.pups.osp))

colnames(isotope.both.pups.osp.simple)<-c("group","community","iso1","iso2")
isotope.both.pups.osp.simple <- isotope.both.pups.osp.simple[, c(3, 4, 1,2)]


siber.example <- createSiberObject(isotope.both.pups.osp.simple)


# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args     <- list(lty = 2, col = "grey20")


par(mfrow=c(1,1))
plotSiberObject(siber.example,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)


par(mfrow=c(1,1))

community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 10000, p.interval = 0.95, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

# this time we will make the points a bit smaller by 
# cex = 0.5
plotSiberObject(siber.example,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = T, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab=expression({delta}^13*C~'\u2030'),
                ylab=expression({delta}^15*N~'\u2030'),
                #cex = 0.5,
                points.order = c(21,23)
)

palette(c("orchid2","oldlace","mistyrose1"))
z=as.factor(c('blue','orange','purple'))


# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(siber.example)
print(group.ML)
#>            1.1       1.2       1.3       2.1       2.2       2.3
#> TA   21.924922 10.917715 17.945127 3.0714363 11.476354 1.4818061
#> SEA   5.783417  3.254484  5.131601 0.8623300  3.458824 0.4430053
#> SEAc  5.989967  3.370715  5.314872 0.8931275  3.582354 0.4588269


# You can add more ellipses by directly calling plot.group.ellipses()
# Add an additional p.interval % prediction ellilpse
plotGroupEllipses(siber.example, n = 10000, p.interval = 0.4,
                  lty = 1, lwd = 2)

# or you can add the XX% confidence interval around the bivariate means
# by specifying ci.mean = T along with whatever p.interval you want.
plotGroupEllipses(siber.example, n = 1000, p.interval = 0.95, ci.mean = T,
                  lty = 1, lwd = 2)



c.id <- 1 # specify the community ID
g.id <- 1 # specify the group ID within the community
g.id2<- 2 # specify the group ID within the community
g.id3 <- 3 # specify the group ID within the community

# see help file for addEllipse for more information
# NB i am using the group identifier g.id to select the colour
# of the ellipse line so that it matches the one created by 
# plotSiberObject(), but you could override this if you wish.
# The function addEllipse returns the coordinates it used for plotting, 
# but more than likely you dont need this information. Here I store these in
# a new variable coords for clarity, but you could just as easily call this tmp.
# See help file for addEllipse for more details on the options, but in short:
# the first two entries look up the means and covariance matrix of the data you
# specified using the group and commmunity indices above.
# m = NULL is used as we are not plotting an ellipse around the mean.
# n = 100 just determines how many points are used to draw a smooth ellipse.
# p.interval = 0.95 for a 95% ellipse of the data
# ci.mean = FALSE as we are not plotting an ellipse around the mean.
# col = your choice of colour.
# lty = your choice of line type.
# lwd = your choice of line width.

ggplot(siber.example$original.data,aes(x=iso1,y=iso2,col=group))+geom_point() 



palette(c("goldenrod1","aquamarine3","#D95F02"))

# base scatter plot
plot(siber.example$original.data$iso1,siber.example$original.data$iso2,xlim=c(-24,-8),ylim=c(3,15),
     xlab=expression({delta}^13*C~'\u2030'),
     ylab=expression({delta}^15*N~'\u2030'),
     col=siber.example$original.data$group,pch=20
)
# add 95% prediction ellipse
addEllipse(siber.example$ML.mu[[c.id]][ , , g.id],
                     siber.example$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = g.id,
                     lty = 3,
                     lwd = 2)
#add standard SEAc ellipse
addEllipse(siber.example$ML.mu[[c.id]][ , , g.id],
           siber.example$ML.cov[[c.id]][ , , g.id],
           m = 17,
           small.sample = T,
           n = 100,
           p.interval = NULL,
           ci.mean = FALSE,
           col = g.id,
           lty = 1,
           lwd = 2)

#add 95% CI around mean ellipse
addEllipse(siber.example$ML.mu[[c.id]][ , , g.id],
           siber.example$ML.cov[[c.id]][ , , g.id],
           m=17,
           n = 1000,
           p.interval = 0.95,
           ci.mean = T,
           col = g.id,
           lty = 1,
           lwd = 2)

#add 95% predicition ellipse
addEllipse(siber.example$ML.mu[[c.id]][ , , g.id2],
           siber.example$ML.cov[[c.id]][ , , g.id2],
           m = 17,
           n = 100,
           p.interval = 0.95,
           ci.mean = FALSE,
           small.sample = T,
           col = g.id3,
           lty = 3,
           lwd = 2)
#add standard SEAc ellipse
addEllipse(siber.example$ML.mu[[c.id]][ , , g.id2],
           siber.example$ML.cov[[c.id]][ , , g.id2],
           m=17,
           small.sample = T,
           n = 1000,
           p.interval = NULL,
           ci.mean = F,
           col = g.id3,
           lty = 1,
           lwd = 2)
#add 95% CI around mean ellipse
addEllipse(siber.example$ML.mu[[c.id]][ , , g.id2],
           siber.example$ML.cov[[c.id]][ , , g.id2],
           m=17,
           n = 1000,
           p.interval = 0.95,
           ci.mean = T,
           col = g.id3,
           lty = 1,
           lwd = 2)

addEllipse(siber.example$ML.mu[[c.id]][ , , g.id3],
           siber.example$ML.cov[[c.id]][ , , g.id3],
           m = NULL,
           n = 100,
           p.interval = 0.95,
           ci.mean = FALSE,
           col = g.id2,
           lty = 3,
           lwd = 2)
#add 95% prediction ellipse
addEllipse(siber.example$ML.mu[[c.id]][ , , g.id3],
           siber.example$ML.cov[[c.id]][ , , g.id3],
           m=14,
           n = 1000,
           p.interval = 0.95,
           ci.mean = T,
           col = g.id2,
           lty = 1,
           lwd = 2)

#add standard SEAc ellipse
addEllipse(siber.example$ML.mu[[c.id]][ , , g.id3],
           siber.example$ML.cov[[c.id]][ , , g.id3],
           m=14,
           n = 1000,
           p.interval = NULL,
           ci.mean = F,
           small.sample = T,
           col = g.id2,
           lty = 1,
           lwd = 2)
#add 95% mean CI ellipse
addEllipse(siber.example$ML.mu[[c.id]][ , , g.id3],
           siber.example$ML.cov[[c.id]][ , , g.id3],
           m=14,
           n = 1000,
           p.interval = 0.95,
           ci.mean = T,
           col = g.id2,
           lty = 1,
           lwd = 2)


# A second plot provides information more suitable to comparing
# the two communities based on the community-level Layman metrics

# this time we will make the points a bit smaller by 
# cex = 0.5
plotSiberObject(siber.example,
                ax.pad = 2, 
                hulls = T, community.hulls.args, 
                ellipses = F, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab=expression({delta}^13*C~'\u2030'),
                ylab=expression({delta}^15*N~'\u2030'),
                cex = 0.5
)

# or you can add the XX% confidence interval around the bivariate means
# by specifying ci.mean = T along with whatever p.interval you want.
plotGroupEllipses(siber.example, n = 100, p.interval = 0.95,
                  ci.mean = T, lty = 1, lwd = 2)


# A second plot provides information more suitable to comparing
# the two communities based on the community-level Layman metrics

# this time we will make the points a bit smaller by 
# cex = 0.5
plotSiberObject(siber.example,
                ax.pad = 2, 
                hulls = T, community.hulls.args, 
                ellipses = F, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab=expression({delta}^13*C~'\u2030'),
                ylab=expression({delta}^15*N~'\u2030'),
                cex = 0.5
)

# or you can add the XX% confidence interval around the bivariate means
# by specifying ci.mean = T along with whatever p.interval you want.
plotGroupEllipses(siber.example, n = 100, p.interval = 0.95,
                  ci.mean = T, lty = 1, lwd = 2)


#### Area and overlap comparisons ########
# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber.example, parms, priors)


SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)


# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)



siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                xlab = c("Community | Group"),
                ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                bty = "L",
                las = 1,
                main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)

# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)






# extract the posterior means
mu.post <- extractPosteriorMeans(siber.example, ellipses.posterior)



# calculate the corresponding distribution of layman metrics
A.P.B <- bayesianLayman(mu.post)


overlap.G3.G4 <- maxLikOverlap("OSP.P", "OSP.S", siber.example, p = 0.95, n =100)
prop.of.first <- as.numeric(overlap.G3.G4["overlap"] / overlap.G3.G4["area.1"])
prop.of.second <- as.numeric(overlap.G3.G4["overlap"] / overlap.G3.G4["area.2"])
print(prop.of.second)
prop.of.both <- as.numeric(overlap.G3.G4["overlap"] / (overlap.G3.G4["area.1"] + overlap.G3.G4["area.2"]))
print(prop.of.both)
prop.of.both.less.overlap <- as.numeric(overlap.G3.G4["overlap"] / (overlap.G3.G4["area.1"] + overlap.G3.G4["area.2"] - overlap.G3.G4["overlap"]))
print(prop.of.both.less.overlap)


print(prop.of.first)
bayesianOverlap(
  1.1,
  1.2,
  ellipses.posterior,
  draws = 10,
  p.interval = 0.95,
  n = 100,
  do.plot = FALSE
)


bayes.overlap.G2.G3 <- bayesianOverlap(
  "OSP.P",
  "OSP.S",
  ellipses.posterior,
  draws = 10,
  p.interval = 0.95,
  n = 100,
  do.plot = TRUE
)

overlap.credibles <- lapply(
  as.data.frame(bayes.overlap.G2.G3), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(overlap.credibles)


overlap.S.A <- maxLikOverlap("OSP.A", "OSP.S", siber.example, p = 0.95, n =100)
prop.of.first <- as.numeric(overlap.S.A["overlap"] / overlap.S.A["area.1"])
prop.of.second <- as.numeric(overlap.S.A["overlap"] / overlap.S.A["area.2"])
print(prop.of.second)
prop.of.both <- as.numeric(overlap.S.A["overlap"] / (overlap.S.A["area.1"] + overlap.S.A["area.2"]))
print(prop.of.both)
prop.of.both.less.overlap <- as.numeric(overlap.S.A["overlap"] / (overlap.S.A["area.1"] + overlap.S.A["area.2"] - overlap.S.A["overlap"]))
print(prop.of.both.less.overlap)


print(prop.of.first)
bayesianOverlap(
  1.1,
  1.2,
  ellipses.posterior,
  draws = 10,
  p.interval = 0.95,
  n = 100,
  do.plot = FALSE
)

hist(bayes.overlap.A.P[,3],10)

bayes.overlap.S.A <- bayesianOverlap(
  "OSP.A",
  "OSP.S",
  ellipses.posterior,
  draws = 10,
  p.interval = 0.95,
  n = 100,
  do.plot = FALSE
)

overlap.credibles.S.A <- lapply(
  as.data.frame(bayes.overlap.S.A), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(overlap.credibles.S.A)


overlap.S.P.over<-(bayes.overlap.A.P[,3]/bayes.overlap.A.P[,1]-bayes.overlap.A.P[,3])
hist(overlap.S.P.over,10)

overlap.S.P <- maxLikOverlap("OSP.S", "OSP.P", siber.example, p = 0.95, n =100)
prop.of.first <- as.numeric(overlap.S.P["overlap"] / overlap.S.P["area.1"])
prop.of.second <- as.numeric(overlap.S.P["overlap"] / overlap.S.P["area.2"])
print(prop.of.second)
prop.of.both <- as.numeric(overlap.S.P["overlap"] / (overlap.S.P["area.1"] + overlap.S.P["area.2"]))
print(prop.of.both)
prop.of.both.less.overlap <- as.numeric(overlap.S.P["overlap"] / (overlap.S.P["area.1"] + overlap.S.P["area.2"] - overlap.S.P["overlap"]))
print(prop.of.both.less.overlap)


print(prop.of.first)


bayes.overlap.S.P <- bayesianOverlap(
  "OSP.P",
  "OSP.S",
  ellipses.posterior,
  draws = 10,
  p.interval = 0.95,
  n = 100,
  do.plot = FALSE
)

overlap.credibles.S.P <- lapply(
  as.data.frame(bayes.overlap.S.P), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(overlap.credibles.S.P)




overlap.A.P <- maxLikOverlap("OSP.A", "OSP.P", siber.example, p = 0.95, n =100)
prop.of.first <- as.numeric(overlap.A.P["overlap"] / overlap.A.P["area.1"])
prop.of.second <- as.numeric(overlap.A.P["overlap"] / overlap.A.P["area.2"])
print(prop.of.second)
prop.of.both <- as.numeric(overlap.A.P["overlap"] / (overlap.A.P["area.1"] + overlap.A.P["area.2"]))
print(prop.of.both)
prop.of.both.less.overlap <- as.numeric(overlap.A.P["overlap"] / (overlap.A.P["area.1"] + overlap.A.P["area.2"] - overlap.A.P["overlap"]))
print(prop.of.both.less.overlap)


print(prop.of.first)


bayes.overlap.A.P2 <- bayesianOverlap(
  "OSP.A",
  "OSP.P",
  ellipses.posterior,
  draws = 100,
  p.interval = 0.95,
  n = 100,
  do.plot = FALSE
)

overlap.credibles.A.P <- lapply(
  as.data.frame(bayes.overlap.A.P), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)


overlap.credibles.A.P2 <- lapply(
  as.data.frame(bayes.overlap.A.P2), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(overlap.credibles.A.P)
print(overlap.credibles.A.P2)


print(overlap.credibles.A.P$overlap)
print(overlap.credibles.S.P$overlap)
print(overlap.credibles.S.A$overlap)










#####graveyared #######
# Calculate the various Layman metrics on each of the communities.
community.ML <- communityMetricsML(siber.example) 
print(community.ML)












# --------------------------------------
# Visualise the first community
# --------------------------------------
siberDensityPlot(layman.B[[1]], xticklabels = colnames(layman.B[[1]]), 
                 bty="L", ylim = c(0,20))

# add the ML estimates (if you want). Extract the correct means 
# from the appropriate array held within the overall array of means.
comm1.layman.ml <- laymanMetrics(siber.example$ML.mu[[1]][1,1,],
                                 siber.example$ML.mu[[1]][1,2,]
)
points(1:6, comm1.layman.ml$metrics, col = "red", pch = "x", lwd = 2)


# --------------------------------------
# Visualise the second community
# --------------------------------------
siberDensityPlot(layman.B[[2]], xticklabels = colnames(layman.B[[2]]), 
                 bty="L", ylim = c(0,20))

# add the ML estimates. (if you want) Extract the correct means 
# from the appropriate array held within the overall array of means.
comm2.layman.ml <- laymanMetrics(siber.example$ML.mu[[2]][1,1,],
                                 siber.example$ML.mu[[2]][1,2,]
)
points(1:6, comm2.layman.ml$metrics, col = "red", pch = "x", lwd = 2)


# --------------------------------------
# Alternatively, pull out TA from both and aggregate them into a 
# single matrix using cbind() and plot them together on one graph.
# --------------------------------------

# go back to a 1x1 panel plot
par(mfrow=c(1,1))

siberDensityPlot(cbind(layman.B[[1]][,"TA"], layman.B[[2]][,"TA"]),
                 xticklabels = c("Community 1", "Community 2"), 
                 bty="L", ylim = c(0,20),
                 las = 1,
                 ylab = "TA - Convex Hull Area",
                 xlab = "")



