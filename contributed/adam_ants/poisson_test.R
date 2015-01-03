#Test for Poisson?

int_cat<-with(dat, interaction(spCode, Subplot, Plant_Species, Month))
avg_cat<-c(by(log(dat$Count), int_cat, mean))
var_cat<-c(by(log(dat$Count), int_cat, var))

plot(avg_cat, var_cat)
abline(a=0, b=1)

###
dat_mean<-mean(dat$Count)
dat_var<-var(dat$Count)

qqplot(dat$Count, rpois(10000,dat_mean))
qqplot(log(dat$Count), rpois(10000,dat_mean))

#######################################
# Bolker cookbook
#######################################
setwd("~/Dropbox/Active Work/BAC ants/multinom_analysis/")

#Load data
dat.tf<-read.csv("bacants.csv")
for(i in 5:7) {
  dat.tf[,i]<-as.factor(as.character(dat.tf[,i]))
}
dat.tf$PlotID<-as.factor(dat.tf$PlotID)
dat.tf$Larson<-as.factor(dat.tf$Larson)
dat.tf$Month<-as.factor(dat.tf$Month)
dat.tf$Year<-as.factor(dat.tf$Year)

spmat<-with(dat.tf, table(rep(paste(PlotID, Subplot, Month, Year, Plant_Species), Count), rep(reduced_spCode, Count)))
envmat<-data.frame(t(matrix(nrow=5, data=unlist(strsplit(rownames(spmat), " ")))))
colnames(envmat)<-c("Plot", "Subplot", "Month", "Year", "PlantSp")
for(i in 1:4) {
  envmat[,i]<-as.factor(envmat[,i])
}
envmat[,5]<-as.numeric(as.character(envmat[,5]))
#dat.tf$Month<-as.numeric(dat.tf$Month)
#dat.tf$Month[dat.tf$smpDate=="08/06/11"]<-8
#dat.tf$Month[dat.tf$smpDate=="09/11/11"]<-9
#dat.tf$Month[dat.tf$smpDate=="8/19/2012"]<-8
#dat.tf$Month[dat.tf$smpDate=="9/23/2012"]<-9
#dat.tf$Month<-as.factor(dat.tf$Month)

#Load packages
library(lme4)
library(coefplot2) ## for coefplot2
library(reshape)
library(plyr)
library(ggplot2)
library(gridExtra)
library(emdbook)
source("glmm_funs.R")

## use within() to make multiple changes within a data frame
dat.tf_save<-dat.tf
dat.tf<-dat.red

dat.tf <- within(dat.tf,
{
  spmc <- interaction(Subplot,Plant_Species,Month,spGroup)
  spmc <- reorder(spmc, Count, mean)
})

## plot data
print(bwplot(log(Count + 1) ~ spmc,
             data=dat.tf,
             scales=list(x=list(rot=90)) ## rotate x-axis labels
))

grpMeans <- with(dat.tf,tapply(Count, list(spmc), mean))
summary(grpMeans)

grpVars <- with(dat.tf,tapply(Count, list(spmc), var)) ## variance of UNTRANSFORMED data

lm1 <- lm(grpVars~grpMeans-1) #'quasi-Poisson' fit
phi.fit <- coef(lm1)

lm2 <- lm(grpVars~I(grpMeans^2)+offset(grpMeans)-1) #negative binomial fit
k.fit <- 1/coef(lm2)



plot( grpVars ~ grpMeans, xlab="group means", ylab="group variances" )
abline(c(0,1), lty=2)
text(40,500,"Poisson")
curve(phi.fit*x, col=2,add=TRUE)
  ## bquote() is used to substitute numeric values
  ## in equations with symbols
text(40,10000,
     bquote(paste("QP: ",sigma^2==.(round(phi.fit,1))*mu)),
     col=2)
curve(x*(1+x/k.fit),col=4,add=TRUE)
text(40,20000,paste("NB: k=",round(k.fit,1),sep=""),col=4)
## loess fit
Lfit <- loess(grpVars~grpMeans)
mvec <- 0:120
lines(mvec,predict(Lfit,mvec),col=5)
text(48,4000,"loess",col=5)
#Negative binomial fits very well
dat.tf<-dat.tf_save


#Plot responses vs. treatments
print(stripplot(log(Count+1) ~ Subplot|Plant_Species,
                data=dat.tf,
                groups=Month,
                auto.key=list(lines=TRUE, columns=2),
                strip=strip.custom(strip.names=c(TRUE,TRUE)),
                type=c('p','a'), layout=c(4,1,1),
                scales=list(x=list(rot=90)), main="Panel: Subplot"))

print(stripplot(log(Count+1) ~ Subplot|Month, dat.tf,
                groups=reduced_spCode,
                strip=strip.custom(strip.names=c(TRUE,TRUE)),
                type=c('p', 'a'), ## points and panel-average value## see ?panel.xyplot
                scales=list(x=list(rot=90)),
                main="Panel: Month, Color: Species"))

#Fit poisson by group
glm.lis <- lmList(Count~Month*Subplot|spCode,data=dat.tf,
                  family="poisson")
print(plot(glm.lis,scale=list(x=list(relation="free"))))
print(qqmath(glm.lis))



#Fit full glmm
#mp1 <- glmer(Count ~ Month*Subplot +
#               Plant_Species +
#               (Month*Subplot|PlotID)+
#               (Month*Subplot|spCode),
#             data=dat.tf, family="poisson")
mp1 <- glmer(Count ~ (Month+Subplot+Plant_Species)^2 +
               (Month+Subplot+Plant_Species|spCode),
             data=dat.tf, family="poisson")

#Check for overdispersion (which exists)
overdisp_fun(mp1)
deviance(mp1)
summary(mp1)

#Add observation level random effect to account for overdispersion
mp2 <- update(mp1,.~.+(1|X))
summary(mp2)
#Check out correlations
attr(VarCorr(mp2)$spCode,"correlation")

#Print var/covar matrix, with all perfect correlations marked
printvc(mp2)

overdisp_fun(mp2) #No longer overdispersed
deviance(mp2)
#Bolker has problems here with perfect correlations - we don't have them in this dataset, but the fit is still very bad.


## Inference
# Plot varcov plot
coefplot2(mp2,ptype="vcov",intercept=TRUE)

mp2_v1<-update(mp2, .~.-(Month+Subplot+Plant_Species|spCode))
anova(mp2, mp2_v1)

## Plot random effects?
## squash margins a bit ...
pp <- list(layout.widths=list(left.padding=0, right.padding=0),
             layout.heights=list(top.padding=0, bottom.padding=0))
r2 <- ranef(mp2,postVar=TRUE)
d2 <- dotplot(r2, par.settings=pp)
library(gridExtra)
grid.arrange(d2$spCode,d2$X,nrow=1)


##Try negative binomial model?
detach("package:coefplot2")
detach("package:lme4")
library(glmmADMB)
library(lme4) ## load lme4 AFTERWARDS
library(coefplot2)

gnb2 <- glmmadmb(Count ~ Subplot*Plant_Species +
                   Month,
                 random=~(Subplot*Plant_Species+Month|spCode),
                 data=dat.tf, family="nbinom")
summary(gnb2)
gnb2ZI <- glmmadmb(Count ~ Subplot*Plant_Species +
                   Month,
                 random=~(Subplot*Plant_Species+Month|spCode),
                 data=dat.tf, family="nbinom1", zeroInflation=TRUE)
summary(gnb2ZI)

gnb3 <- glmmadmb(Count ~ Subplot+Plant_Species+Month+Subplot:Plant_Species+Plant_Species:Month+Subplot:Month,
                 random=~(Subplot+Plant_Species+Month|spCode),
                 data=dat.tf, family="nbinom2")
summary(gnb3)
coefficients(gnb3)

#KISS models
dat.tf$Year<-as.factor(dat.tf$Year)
dat.tf$PlotID<-as.factor(dat.tf$PlotID)
dat.tf$Month<-as.factor(dat.tf$Month)
dat.tf$Larson<-as.factor(dat.tf$Larson)

gnb4 <- glmmadmb(Count ~ Subplot,
                 data=dat.tf, family="nbinom2")
summary(gnb4)
gnb4ZI <- glmmadmb(Count ~ Subplot,
                 data=dat.tf, family="nbinom2", zeroInflation=T)
summary(gnb4ZI)
anova(gnb4, gnb4ZI)

#Add random effect?
gnb5ZI <- glmmadmb(Count ~ Subplot,
                   random=~(1|reduced_spCode),
                   zeroInflation=TRUE,
                   data=dat.tf, family="nbinom2")
summary(gnb5ZI)

gnb6ZI <- glmmadmb(Count ~ Subplot+Plant_Species,
                   random=~(1|reduced_spCode),
                   zeroInflation=TRUE,
                   data=dat.tf, family="nbinom2")
summary(gnb6ZI)

gnb7ZI <- glmmadmb(Count ~ Subplot+Plant_Species+Month,
                   random=~(1|reduced_spCode),
                   zeroInflation=TRUE,
                   data=dat.tf, family="nbinom2")
summary(gnb7ZI)
#Woohoo!

gnb8ZI <- glmmadmb(Count ~ Subplot+Plant_Species+Month,
                   random=~(1+Month|reduced_spCode),
                   zeroInflation=TRUE,
                   data=dat.tf, family="nbinom2")
summary(gnb8ZI)
gnb9ZI <- glmmadmb(Count ~ Subplot+Plant_Species+Month,
                   random=~(1+Month+Subplot|reduced_spCode),
                   zeroInflation=TRUE,
                   data=dat.tf, family="nbinom2")
summary(gnb9ZI)


gnb10ZI <- glmmadmb(Count ~ Subplot+Plant_Species+Month,
                   random=~(1+Month+Subplot+Plant_Species|reduced_spCode),
                   zeroInflation=TRUE,
                   data=dat.tf, family="nbinom2")
summary(gnb10ZI)

gnb10 <- glmmadmb(Count ~ Subplot+Plant_Species+Month,
                    random=~(1+Month+Subplot+Plant_Species|reduced_spCode),
                    data=dat.tf, family="nbinom2")
summary(gnb10)

gnb11 <- glmmadmb(Count ~ Subplot+Plant_Species+Month,
                  random=~(1|reduced_spCode/Year/PlotID),
                  data=dat.tf, family="nbinom2")
summary(gnb11)


gnb12 <- glmmadmb(Count ~ Subplot*Plant_Species+Month,
                  random=~(1|reduced_spCode/Year/PlotID),
                  data=dat.tf, family="nbinom2")
summary(gnb12)

gnb12ZI <- glmmadmb(Count ~ Subplot+Plant_Species+Month,
                  random=~(1|reduced_spCode/Year/PlotID),
                  zeroInflation=TRUE,
                  data=dat.tf, family="nbinom2")
summary(gnb12ZI)

gnb13 <- glmmadmb(Count ~ Subplot+Plant_Species+Month,
                    random=~(Subplot|reduced_spCode)+(1|Year/PlotID),
                    data=dat.tf, family="nbinom2")
summary(gnb13)

gnb14 <- glmmadmb(Count ~ Subplot+Plant_Species+Month,
                  random=~(Subplot|reduced_spCode)+(Subplot|Year/PlotID),
                  data=dat.tf, family="nbinom2")
summary(gnb14)
anova(gnb13, gnb14)

gnb15 <- glmmadmb(Count ~ Subplot+Plant_Species+Month,
                  random=~(Subplot|reduced_spCode)+(Subplot|Year/PlotID)+(Month/Year),
                  data=dat.tf, family="nbinom2")
summary(gnb15)
anova(gnb15, gnb15) #Not worth it

gnb16 <- glmmadmb(Count ~ Subplot*(Plant_Species+Month),
                  random=~(Subplot|reduced_spCode)+(Subplot|Year/PlotID),
                  data=dat.tf, family="nbinom2")
summary(gnb16)
ranef(gnb16)
hist(ranef(gnb16)$reduced_spCode[,2])
plot(density(ranef(gnb16)$reduced_spCode[,2]))

#mcmc?
gnb16mcmc <- glmmadmb(Count ~ Subplot*(Plant_Species+Month),
                  random=~(Subplot|reduced_spCode)+(Subplot|Year/PlotID),
                  data=dat.tf, family="nbinom2", mcmc=TRUE)
summary(gnb16mcmc)
library(coda)
library(scapeMCMC)
library(coefplot2)
m <- as.mcmc(gnb16mcmc$mcmc)
plotTrace(m)

(gg <- geweke.diag(m))
head(t(apply(m,2,quantile,c(0.025,0.975))))
head(t(apply(m,2,quantile,c(0.5))))

#29Dec2013
#Make data frame for reduced_Species, with zeros included
if(FALSE) { #Make 
dat.red<-with(dat.tf, aggregate(Count, list(PlotID, Subplot, reduced_spCode, Month, Year), FUN=sum))
for(i in 1:6) {
  dat.red[,i]<-as.character(dat.red[,i])
}
colnames(dat.red)<-c("PlotID", "Subplot", "spGroup", "Month", "Year", "Count")
datarr<-with(dat.red, table(spGroup, PlotID, Subplot, Year, Month))
datarrnames<-dimnames(datarr)
dat.red_old<-dat.red

n<-dim(dat.red)[1]+1
tmp<-matrix(nrow=length(unlist(datarr))-dim(dat.red)[1], ncol=dim(dat.red)[2], data=NA)
colnames(tmp)<-c("PlotID", "Subplot", "spGroup", "Month", "Year", "Count")
dat.red<-rbind(dat.red, tmp)
rm(tmp)



for(i in 1:12) {
  for(j in 1:38) {
    for(k in 1:2) {
      for(l in 1:2) {
        for(m in 1:2) {
          subs<-sum(dat.red_old$spGroup==datarrnames$spGroup[i]&
            dat.red_old$PlotID==datarrnames$PlotID[j]&
            dat.red_old$Subplot==datarrnames$Subplot[k]&
            dat.red_old$Year==datarrnames$Year[l]&
            dat.red_old$Month==datarrnames$Month[m])
          
          if(subs==0) {
            dat.red[n,]<-c(datarrnames$PlotID[j], datarrnames$Subplot[k], datarrnames$spGroup[i], datarrnames$Month[m], datarrnames$Year[l], 0)
            n<-n+1
          }
          print(n/nrow(dat.red))
        }
      }
    }
  }
}

rm(dat.red_old)
head(dat.red)
for(i in 1:5) {
  dat.red[,i]<-as.factor(dat.red[,i])
}
dat.red$Plant_Species<-dat.tf[match(dat.red$PlotID, dat.tf$PlotID),"Plant_Species"]
dat.red$Count<-as.numeric(as.character(dat.red$Count))

write.csv(dat.red, "datred.csv", row.names=FALSE)
}
dat.red<-read.csv("datred.csv")
for(i in c(1,4,5)) {
  dat.red[,i]<-as.factor(dat.red[,i])
}

ffit_nbin1 <- glmmadmb(Count ~ Subplot+(Plant_Species+Month),
                      #random=~(1|reduced_spCode),#+(Subplot|Year/PlotID),
                      data=dat.red, family="nbinom2",
                      zeroInflation=TRUE,
                      verbose=TRUE)
summary(ffit_nbin1)
ffit_nbin1
coefplot2(ffit_nbin1)
coefplot2(gnb16mcmc)

ranef(ffit_nbin1)

ffit_nbin2 <- glmmadmb(Count ~ Subplot*(Plant_Species+Month),
                       #random=~(1|reduced_spCode),#+(Subplot|Year/PlotID),
                       data=dat.red, family="nbinom2",
                       zeroInflation=TRUE,
                       verbose=TRUE)
anova(ffit_nbin1, ffit_nbin2) #Interaction does not add significant improvement

ffit_nbin3 <- glmmadmb(Count ~ Subplot+(Plant_Species+Month),
                       random=~(1|spGroup),#+(Subplot|Year/PlotID),
                       data=dat.red, family="nbinom2",
                       zeroInflation=TRUE,
                       verbose=TRUE)
anova(ffit_nbin1, ffit_nbin3)
summary(ffit_nbin3)
coefplot2(ffit_nbin3)

ranef(ffit_nbin3)
plot(ranef(ffit_nbin3)$spGroup+coef(ffit_nbin3)[2])
abline(h=0)


mod <- ffit_nbin4 <- glmmadmb(Count ~ Subplot+(Plant_Species+Month),
                       random=~(1|spGroup)+(1|Year/PlotID),
                       data=dat.red, family="nbinom2",
                       verbose=TRUE)

anova(ffit_nbin1, ffit_nbin3, mod)
summary(mod)
coefplot2(mod)

#Plot output?
plot(1:12, ranef(ffit_nbin4, sd=FALSE)$spGroup[,1]+coefficients(ffit_nbin4)[1],
     axes=FALSE,
     xlab="Species", ylab="Intercept")
axis(2)


mod <- ffit_nbin5 <- glmmadmb(Count ~ Subplot+(Plant_Species+Month),
                              random=~(1+Subplot|spGroup)+(1|Year/PlotID),
                              data=dat.red, family="nbinom2",
                              verbose=TRUE)

anova(ffit_nbin1, ffit_nbin3, ffit_nbin4, mod)
summary(mod)
coefplot2(mod)
#No improvement

mod <- ffit_nbin6 <- glmmadmb(Count ~ Subplot+(Plant_Species+Month),
                              random=~(1+Subplot+Month|spGroup)+(1+Subplot|Year/PlotID),
                              data=dat.red, family="nbinom2",
                              zeroInflation=TRUE,
                              verbose=TRUE)

system("espeak 'simulation complete'")
anova(ffit_nbin1, ffit_nbin3, ffit_nbin4, mod)
summary(mod)
coefplot2(mod)
ranef(mod, sd=TRUE)
ranef(mod, sd=FALSE)

mod <- ffit_nbin7 <- glmmadmb(Count ~ Subplot+(Plant_Species+Month),
                              random=~(1+Subplot+Month|spGroup)+(1+Subplot|Year/PlotID),
                              data=dat.red, family="nbinom2",
                              zeroInflation=TRUE,
                              verbose=TRUE)

system("espeak 'simulation complete'")
anova(ffit_nbin1, ffit_nbin3, ffit_nbin4, ffit_nbin6, mod)
summary(mod)
coefplot2(mod)
ranef(mod, sd=TRUE)
ranef(mod, sd=FALSE)


mod <- ffit_nbin7 <- glmmadmb(Count ~ (Subplot+Plant_Species+Month)*spGroup,
                              data=dat.red,
                              family="nbinom2",
                              zeroInflation=TRUE,
                              random=~(1+Subplot+Month|Year/PlotID),
                              verbose=TRUE)
summary(mod)
coefplot2(mod, intercept=TRUE)
ranef(mod, sd=TRUE)
ranef(mod, sd=FALSE)

tmp<-ranef(mod, sd=FALSE)$spGroup[,2]+coef(mod)[2]
plot(1:12, tmp,
     ylim=c(-2.2, 2.8),
     xlab="Species", ylab="Subplot")
segments(1:12,
         tmp+2*ranef(mod, sd=TRUE)[[1]][,2],
         1:12,
         tmp-2*ranef(mod, sd=TRUE)[[1]][,2],
         lwd=1)
segments(1:12,
         tmp+ranef(mod, sd=TRUE)[[1]][,2],
         1:12,
         tmp-ranef(mod, sd=TRUE)[[1]][,2],
         lwd=2)
abline(h=0)
########################################################################################################################################################################################
############ Marker
########################################################################################################################################################################################

# AIC comparisons among models
library(bbmle)
AICtab(fit_zipoiss,fit_zinbinom,fit_zinbinom1,fit_zinbinom1_bs)


#Look out outputs
ranef(gnb13) #random effects estimates
predict(..., type="response", se.fit=T, interval="confidence", level=0.95)
#NOTE can add mcmc=TRUE to regressions...

vcov(); VarCorr() #var/cov for main and mixed effects

#Show coefplot
library(coefplot2)
ff <- fits[c("gnb2","gnb3")]
names(ff) <- c("NB1, full","NB2, full")
coefplot2(gnb8ZI,legend.x="right",legend=TRUE)



#Check out location scale plot
op <- par(mfrow=c(1,2))
locscaleplot(mp4)
locscaleplot(gnb8ZI)
par(op)

