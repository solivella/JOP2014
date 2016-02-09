###############################
# Replication of results
# reported in "(Where) Do Campaigns 
# Matter? The Impact of National 
# Party Convention Location"
# by Atkinson et al.
#
# Last update: SO
# Date: 12/25/13
##############################

library(lmtest)
library(apsrtable)
library(foreign)
library(MASS)
library(plyr)
library(abind)
library(lattice)
library(latticeExtra)

set.seed(831213)
source("ClusteredStandardErrorsFunction.R")

##Import replication data

ConventionsDataSub <- read.csv("ConventionsReplicationData.csv")
ConventionsDataSub$Convention <- relevel(ConventionsDataSub$Convention,ref="N")

#####################################################
# Main model estimation: local convention effects
#####################################################
DemMargin.model <- lm(DeltaDemPerc ~ 
                        Convention*LagDemPerc
                      +Convention*I(LagDemPerc^2)
                      + I(DeltaBlackProp*100)
                      + I(DeltaWhiteProp*100)
                      + I(DeltaBel25Prop*100)
                      + I(DeltaAb65Prop*100)
                      + I(DeltaUrban*100)
                      + I(DeltaFem*100)
                      + I(DeltaAsian*100)
                      + I(DeltaNative*100)
                      + I(DeltaHispanic*100)
                      + I(DeltaMarried*100)
                      + I(DeltaUnemploy*100)
                      + I(DeltaCollege*100)
                      + I(DeltaIncome/1000)
                      + prop.score
                      + as.factor(year)
                      ,weights=log(totvote)/(log(apply(cbind(DistRNC+1,DistDNC+1),1,min))+1)
                      ,data=ConventionsDataSub)


CorrectedVCov <- cl(ConventionsDataSub,DemMargin.model,ConventionsDataSub$fips)
#CorrectedVCov <- vcov(DemMargin.model)

#Effect of DNC and RNC conditional on lagged margin
beta.samples <- mvrnorm(1000,coef(DemMargin.model),CorrectedVCov)
beta.samplesDNC <- beta.samples[,c(2,30,32)]
beta.samplesRNC <- beta.samples[,c(3,31,33)]
sim.demmar <- seq(20,80,length.out=900)
x.mat <- cbind(1,sim.demmar,sim.demmar^2)#,sim.demmar^3)
colnames(x.mat) <- c("Intercept","DemMargin","DemMargin.2")#,"DemMargin.3")
Pred.DNC.samples <- t(aaply(x.mat,1,function(x)beta.samplesDNC%*%x))
Pred.RNC.samples <- t(aaply(x.mat,1,function(x)beta.samplesRNC%*%x))
Pred.DNC <- adply(Pred.DNC.samples,2,quantile,probs=c(0.05,0.5,0.95))
Pred.DNC$Convention <- "Being Exposed to DNC"
Pred.DNC$SimVals <- sim.demmar
names(Pred.DNC) <- c("X","Low","Median","High","Convention","SimVals")
Pred.RNC <- adply(Pred.RNC.samples,2,quantile,probs=c(0.13,0.5,0.96))
Pred.RNC$Convention <- "Being Exposed to RNC"
Pred.RNC$SimVals <- sim.demmar
names(Pred.RNC) <- c("X","Low","Median","High","Convention","SimVals")
Pred.Votes <- rbind(Pred.DNC,Pred.RNC)
Pred.Votes <- Pred.Votes[,-1]

plot(Median~SimVals,data=subset(Pred.Votes,Convention=="Being Exposed to DNC"),type="l")
with(Pred.Votes,polygon(c(SimVals[Convention=="Being Exposed to DNC"],rev(SimVals[Convention=="Being Exposed to DNC"]))
                        ,c(High[Convention=="Being Exposed to DNC"],rev(Low[Convention=="Being Exposed to DNC"]))))
lines(Median~SimVals,data=subset(Pred.Votes,Convention=="Being Exposed to RNC"))
with(Pred.Votes,polygon(c(SimVals[Convention=="Being Exposed to RNC"],rev(SimVals[Convention=="Being Exposed to RNC"]))
                        ,c(High[Convention=="Being Exposed to RNC"],rev(Low[Convention=="Being Exposed to RNC"]))))


##Plot Conditional effects
postscript("Figure1.eps",width=12,height=5)
mytheme <- standard.theme("pdf", color=FALSE)
county.rug <- sample(ConventionsDataSub$LagDemPerc,1e3)
xyplot(Median~SimVals|Convention
       ,ylim=c(-7.5,7.5)
       ,ylab="Predicted Change in Democratic Percentage of Vote"
       ,xlab="Previous Democratic Percentage of the Vote"
       ,par.settings=mytheme
       ,scales=list(x=list(alternating=FALSE))
       ,panel=function(x,y,subscripts,...){
         panel.polygon(c(Pred.Votes$SimVals[subscripts],rev(Pred.Votes$SimVals[subscripts]))
                       ,c(Pred.Votes$High[subscripts],rev(Pred.Votes$Low[subscripts]))
                       ,col="gray70"
                       ,border="white")
         panel.lines(x,y,col="white",...,lwd=2.5)
         panel.rug(county.rug)
         panel.abline(h=0,lty=2,color="gray50")
         panel.abline(v=50,lty=2,color="gray50")       
       }
       ,strip = function(..., bg,par.strip.text){
         strip.default(..., bg="black"
                       ,par.strip.text=list(alpha=1
                                            ,cex=1.05
                                            ,col="white"
                                            ,font=2
                       ))}
       ,data=Pred.Votes
)
dev.off()

#######################################################################
### Turnout model
#########################################################
Turnout.model <- lm(I(DeltaTurnout*100) ~ 
                      Convention*LagDemPerc
                    +Convention*I(LagDemPerc^2)
                    +LagTurnout
                    + I(DeltaBlackProp*100)
                    + I(DeltaWhiteProp*100)
                    + I(DeltaBel25Prop*100)
                    + I(DeltaAb65Prop*100)
                    + I(DeltaUrban*100)
                    + I(DeltaFem*100)
                    + I(DeltaAsian*100)
                    + I(DeltaNative*100)
                    + I(DeltaHispanic*100)
                    + I(DeltaMarried*100)
                    + I(DeltaUnemploy*100)
                    + I(DeltaCollege*100)
                    + I(DeltaIncome/1000)
                    + PopChange
                    + prop.score
                    +as.factor(year)
                    ,weights=log(ab18pop)/(log(apply(cbind(DistRNC+1,DistDNC+1),1,min))+1)
                    ,data=ConventionsDataSub)

ConventionsDataSubTO <- ConventionsDataSub[rownames(model.frame(Turnout.model)),]

CorrectedVCovTO <- cl(ConventionsDataSubTO,Turnout.model,ConventionsDataSubTO$fips)

#Effect of DNC and RNC conditional on lagged margin
beta.samplesTO <- mvrnorm(1000,coef(Turnout.model),CorrectedVCovTO)
beta.samplesDNCTO <- beta.samplesTO[,c(2,31,33)]
beta.samplesRNCTO <- beta.samplesTO[,c(3,32,34)]
Pred.DNC.samplesTO <- t(aaply(x.mat,1,function(x)beta.samplesDNCTO%*%x))
Pred.RNC.samplesTO <- t(aaply(x.mat,1,function(x)beta.samplesRNCTO%*%x))
Pred.DNCTO <- adply(Pred.DNC.samplesTO,2,quantile,probs=c(0.05,0.5,0.95))
Pred.DNCTO$Convention <- "Being Exposed to DNC"
Pred.DNCTO$SimVals <- sim.demmar
names(Pred.DNCTO) <- c("X","Low","Median","High","Convention","SimVals")
Pred.RNCTO <- adply(Pred.RNC.samplesTO,2,quantile,probs=c(0.05,0.5,0.95))
Pred.RNCTO$Convention <- "Being Exposed to RNC"
Pred.RNCTO$SimVals <- sim.demmar
names(Pred.RNCTO) <- c("X","Low","Median","High","Convention","SimVals")
Pred.TO <- rbind(Pred.DNCTO,Pred.RNCTO)
Pred.TO <- Pred.TO[,-1]

##Plot Conditional effects
postscript("Figure2.eps",width=12,height=5)
mytheme <- standard.theme("pdf", color=FALSE)
xyplot(Median~SimVals|Convention
       ,ylab="Predicted Change in Turnout"
       ,xlab="Previous Democratic Percentage of the Vote"
       ,scales=list(x=list(alternating=FALSE))
       ,par.settings=mytheme       
       ,panel=function(x,y,subscripts,...){
         panel.polygon(c(Pred.TO$SimVals[subscripts],rev(Pred.TO$SimVals[subscripts]))
                       ,c(Pred.TO$High[subscripts],rev(Pred.TO$Low[subscripts]))
                       ,col="gray70"
                       ,border="white")
         panel.lines(x,y,col="white",...,lwd=2.5)
         panel.abline(h=0,lty=2,color="gray50")
         
       }
       ,strip = function(..., bg,par.strip.text){
         strip.default(..., bg="black"
                       ,par.strip.text=list(alpha=1
                                            ,cex=1.05
                                            ,col="white"
                                            ,font=2
                       ))}
       ,data=Pred.TO
)
dev.off()
#################################################
# Counterfactual analysis
#################################################
mm.fips <- by(ConventionsDataSub,ConventionsDataSub$MMName[,drop=TRUE],function(x)unique(x$fips))
mm.fips$O9999 <- c(39043,39077,39139,39093,39005,39035
                   ,39103,39169,39075,39085,39055
                   ,39133,39151,39019,39157
)
mm.fips <- list(R2008=mm.fips[[20]])

cf.function <- function(election,party="D"){
  cf.by.market <- ldply(mm.fips,function(x){
    temp <- subset(ConventionsDataSubTO,year==election&(fips%in%x)) 
    temp$Convention <- party
    temp$SubStateDemVote <- with(temp,ave(demvote,ICPSRstate,FUN=sum,na.rm=TRUE))
    temp$SubStateTotVote <- with(temp,ave(totvote,ICPSRstate,FUN=sum,na.rm=TRUE))    
    temp$UnaffectedDemVote <- with(temp,StateDemVote-SubStateDemVote)
    temp$UnaffectedTotVote <- with(temp,StateTotVote-SubStateTotVote)
    temp$predDDem <- predict(DemMargin.model,temp)
    temp$predDTO <- predict(Turnout.model,temp)
    temp$NewDemPerc <- with(temp,predDDem+LagDemPerc)
    temp$NewTO <- with(temp, predDTO+(LagTurnout*100))
    temp$NewTotVote <- with(temp,(NewTO/100)*ab18pop)
    temp$NewDemVote <- with(temp,(NewDemPerc/100)*NewTotVote)
    temp$NewSubStateDemVote <- with(temp,ave(NewDemVote,ICPSRstate,FUN=sum,na.rm=TRUE))
    temp$NewSubStateTotVote <- with(temp,ave(NewTotVote,ICPSRstate,FUN=sum,na.rm=TRUE))
    temp$NewStateDemVote <- with(temp,NewSubStateDemVote+UnaffectedDemVote)
    print(temp$UnaffectedDemVote)
    print(temp$NewStateDemVote)
    temp$NewStateTotVote <- with(temp,NewSubStateTotVote+UnaffectedTotVote)
    temp$NewStateBiVote <- with(temp,NewStateDemVote+(StateRepVote/StateTotVote)*NewStateTotVote)
    temp$NewStateDemPerc <- with(temp,(NewStateDemVote/NewStateTotVote)*100)
    temp$NewStateBiDemPerc <- with(temp,(NewStateDemVote/NewStateBiVote)*100)
    temp$OldStateDemPerc <- with(temp,(StateDemVote/StateTotVote)*100)
    temp$OldStateBiDemPerc <- with(temp,(StateDemVote/(StateDemVote+StateRepVote))*100)
    unique.temp <- subset(temp,!duplicated(ICPSRstate),select=c(ICPSRstate,NewStateBiDemPerc,OldStateBiDemPerc))
    return(unique.temp)
  }
  )
  return(cf.by.market)
}

v.cf.function <- Vectorize(cf.function, SIMPLIFY=FALSE)
CounterFactualsTemp <-  outer(seq(1976,2012,by=4),c("D","R","N"),v.cf.function)
CounterFactuals <- adply(CounterFactualsTemp,1,function(x){
  DNC.fx <- (x[1][[1]][,3]-x[3][[1]][,3])+(x[1][[1]][,4])
  RNC.fx <- (x[2][[1]][,3]-x[3][[1]][,3])+(x[2][[1]][,4])
  DNC.flip <- ifelse((DNC.fx>50&x[1][[1]][,4]<50),1,ifelse(DNC.fx<50&x[1][[1]][,4]>50,-1,0)) 
  RNC.flip <- ifelse((RNC.fx>50&x[1][[1]][,4]<50),-1,ifelse(RNC.fx<50&x[1][[1]][,4]>50,1,0))
  temp <- data.frame(DNC.flip,RNC.flip,Convention=x[1][[1]][,1],State=x[1][[1]][,2])
  temp <- subset(temp,DNC.flip==1|RNC.flip==1)
})
CounterFactuals$Election <- lapply(CounterFactuals$X1,function(x)switch(x,1976,1980,1984,1988,1992,1996,2000,2004,2008,2012))
CounterFactuals


##Calculate effects to compare to Shaw's data
shaw.subsetD <- subset(ConventionsDataSub,year >=2000&Convention=="D")
shaw.subsetR <- subset(ConventionsDataSub,year >=2000&Convention=="R")
shaw.subsetDN <- shaw.subsetD
shaw.subsetDN$Convention <- "N"
shaw.subsetRN <- shaw.subsetR
shaw.subsetRN$Convention <- "N"

shaw.subsetD$DFXVotes <- predict(DemMargin.model,shaw.subsetD)-predict(DemMargin.model,shaw.subsetDN)
shaw.DFXVotesAve <- by(shaw.subsetD,shaw.subsetD$year,function(x)weighted.mean(x$DFXVotes,w=x$TotPop))

shaw.subsetR$FXVotes <- predict(DemMargin.model,shaw.subsetR)-predict(DemMargin.model,shaw.subsetRN)
shaw.RFXVotesAve <- by(shaw.subsetR,shaw.subsetR$year,function(x)weighted.mean(x$FXVotes,w=x$TotPop))

shaw.subsetD$FXTO <- predict(Turnout.model,shaw.subsetD)-predict(Turnout.model,shaw.subsetDN)
shaw.DFXTOAve <- by(shaw.subsetD,shaw.subsetD$year,function(x)weighted.mean(x$FXTO,w=x$TotPop))

shaw.subsetR$FXTO <- predict(Turnout.model,shaw.subsetR)-predict(Turnout.model,shaw.subsetRN)
shaw.RFXTOAve <- by(shaw.subsetR,shaw.subsetR$year,function(x)weighted.mean(x$FXTO,w=x$TotPop))

##Effects plots
#DNC
w.DNC.2012  <- subset(ConventionsDataSub,year==2012)
w.DNC.2012$Convention <- "D"
wo.DNC.2012  <- subset(ConventionsDataSub,year==2012)
wo.DNC.2012$Convention <- "N"
DNC.2012 <- predict(DemMargin.model,w.DNC.2012)-predict(DemMargin.model,wo.DNC.2012)
DNC.2012.df <- data.frame(Effect=DNC.2012,Convention="Being Exposed to DNC")
#RNC
w.RNC.2012  <- subset(ConventionsDataSub,year==2012)
w.RNC.2012$Convention <- "R"
wo.RNC.2012  <- subset(ConventionsDataSub,year==2012)
wo.RNC.2012$Convention <- "N"
RNC.2012 <- predict(DemMargin.model,w.RNC.2012)-predict(DemMargin.model,wo.RNC.2012)
RNC.2012.df <- data.frame(Effect=RNC.2012,Convention="Being Exposed to RNC")
FX.2012 <- rbind(DNC.2012.df,RNC.2012.df)

postscript("Figure4.eps",width=12,height=5)
mytheme <- standard.theme("pdf", color=FALSE)
densityplot(~Effect|Convention
            ,ylab=""
            ,scales=list(x=list(alternating=FALSE))
            ,xlab="Predicted Change in Democratic Percentage of Vote"
            ,par.settings=mytheme
            ,lwd=2
            ,plot.points=FALSE
            ,strip = function(..., bg,par.strip.text){
              strip.default(..., bg="black"
                            ,par.strip.text=list(alpha=1
                                                 ,cex=1.05
                                                 ,col="white"
                                                 ,font=2
                            ))}
            ,data=FX.2012
)
dev.off()


#####################################
## Analysis of individual NAES surveys
######################################
STACKED.data <- read.csv("NAESStacked.csv")

d.m2.probit <- glm(y1 ~ y0*host
                   + other.city
                   + info
                   + religion
                   + female
                   + education
                   + income
                   + white
                   + as.factor(yr)
                   ,data=subset(STACKED.data,is.rnc==0)
                   ,family=binomial(link="probit"))
r.m2.probit <- glm(y1 ~ y0*host
                   + other.city
                   + info
                   + religion
                   + female
                   + education
                   + income
                   + white
                   +as.factor(yr)
                   ,data=subset(STACKED.data,is.rnc==1)
                   ,family=binomial(link="probit"))


#Treatment effects on the treated, conditional on pre intent
#For DNC
dnc.exposed.dem <- subset(STACKED.data,is.rnc==0&host==1)
dnc.exposed.dem$host <- 1
dnc.exposed.dem$y0 <- 1
dnc.nexposed.dem <- subset(STACKED.data,is.rnc==0&host==1)
dnc.nexposed.dem$host <- 0
dnc.nexposed.dem$y0 <- 1
dnc.exposed.rep <- subset(STACKED.data,is.rnc==0&host==1)
dnc.exposed.rep$host <- 1
dnc.exposed.rep$y0 <- 0
dnc.nexposed.rep <- subset(STACKED.data,is.rnc==0&host==1)
dnc.nexposed.rep$host <- 0
dnc.nexposed.rep$y0 <- 0

dnc.fx.dem <- (predict(d.m2.probit,dnc.exposed.dem,type="response")-predict(d.m2.probit,dnc.nexposed.dem,type="response"))
naes.dnc.dem <- data.frame(Effect=dnc.fx.dem,Condition="Democrat (Pre)",Treatment="Being Exposed to DNC\n(Treatment Effect on Treated)") 

dnc.fx.rep <- (predict(d.m2.probit,dnc.exposed.rep,type="response")-predict(d.m2.probit,dnc.nexposed.rep,type="response"))
naes.dnc.rep <- data.frame(Effect=dnc.fx.rep,Condition="Republican (Pre)",Treatment="Being Exposed to DNC\n(Treatment Effect on Treated)") 

#For RNC
rnc.exposed.dem <- subset(STACKED.data,is.rnc==1&host==1)
rnc.exposed.dem$host <- 1
rnc.exposed.dem$y0 <- 0
rnc.nexposed.dem <- subset(STACKED.data,is.rnc==1&host==1)
rnc.nexposed.dem$host <- 0
rnc.nexposed.dem$y0 <- 0
rnc.exposed.rep <- subset(STACKED.data,is.rnc==1&host==1)
rnc.exposed.rep$host <- 1
rnc.exposed.rep$y0 <- 1
rnc.nexposed.rep <- subset(STACKED.data,is.rnc==1&host==1)
rnc.nexposed.rep$host <- 0
rnc.nexposed.rep$y0 <- 1

rnc.fx.dem <- (predict(r.m2.probit,rnc.exposed.dem,type="response")-predict(r.m2.probit,rnc.nexposed.dem,type="response"))
naes.rnc.dem <- data.frame(Effect=rnc.fx.dem,Condition="Republican (Pre)",Treatment="Being Exposed to RNC\n(Treatment Effect on Treated)")

rnc.fx.rep <- (predict(r.m2.probit,rnc.exposed.rep,type="response")-predict(r.m2.probit,rnc.nexposed.rep,type="response"))
naes.rnc.rep <- data.frame(Effect=rnc.fx.rep,Condition="Democrat (Pre)",Treatment="Being Exposed to RNC\n(Treatment Effect on Treated)")

naes.fx <- rbind(naes.dnc.dem,naes.dnc.rep,naes.rnc.dem,naes.rnc.rep)
myStripStyle <- function(which.panel, factor.levels, ...) {
  panel.rect(0, -1.5, 1, 1,
             col = "black",
             border = 1)
  panel.text(x = 0.5, y = -0.2,
             font=2,
             lab = factor.levels[which.panel],
             col = "white")
}    

postscript("Figure3.eps",width=12,height=5)
mytheme <- standard.theme("pdf", color=FALSE)
densityplot(~Effect|Treatment
            ,ylab=""
            ,scales=list(x=list(alternating=FALSE))
            ,xlab="Predicted Change in Probability of Supporting\nHosting Party After Convention"
            ,par.settings=mytheme
            ,groups=Condition
            ,lwd=2
            ,plot.points=FALSE
            ,panel=function(x,groups,...){
              require(grid) 
              panel.densityplot(x,groups=groups,...)
              if(panel.number()==1){
                panel.text(x = -0.08
                           ,y = 15.5
                           ,labels=c("Affiliation in\nPre-Election Wave:"))
                draw.key(simpleKey(c("Republican","Democrat"),FALSE,FALSE,TRUE), 
                         draw = TRUE, 
                         vp = viewport(x = unit(0.35, "npc"), y = unit(0.7, "npc")))                               
              } 
            }                       
            ,strip = myStripStyle
            ,data=naes.fx
)
dev.off()

######################################
#### Online Appendix
######################################

########Descriptve stats table
DescTable <- data.frame(Mean=c(0.29,colMeans(model.matrix(DemMargin.model))[-c(1:3,20:33)])
                        ,SD=c(10.01,apply(model.matrix(DemMargin.model)[,-c(1:3,20:33)],2,sd))
                        ,IQR=c(10.15,apply(model.matrix(DemMargin.model)[,-c(1:3,20:33)],2,IQR)))


######Analysis restricted to 1980-2012
ConventionsDataSubRE <- subset(ConventionsDataSub,year>1980)
DemMargin.modelRE <- lm(DeltaDemPerc ~ 
                        Convention*LagDemPerc
                      +Convention*I(LagDemPerc^2)
                      + I(DeltaBlackProp*100)
                      + I(DeltaWhiteProp*100)
                      + I(DeltaBel25Prop*100)
                      + I(DeltaAb65Prop*100)
                      + I(DeltaUrban*100)
                      + I(DeltaFem*100)
                      + I(DeltaAsian*100)
                      + I(DeltaNative*100)
                      + I(DeltaHispanic*100)
                      + I(DeltaMarried*100)
                      + I(DeltaUnemploy*100)
                      + I(DeltaCollege*100)
                      + I(DeltaIncome/1000)
                      + prop.score
                      + as.factor(year)
                      ,weights=log(totvote)/(log(apply(cbind(DistRNC+1,DistDNC+1),1,min))+1)
                      ,data=ConventionsDataSubRE)

CorrectedVCov <- vcov(DemMargin.modelRE)

#Effect of DNC and RNC conditional on lagged margin
beta.samples <- mvrnorm(1000,coef(DemMargin.modelRE),CorrectedVCov)
beta.samplesDNC <- beta.samples[,c(2,27,29)]
beta.samplesRNC <- beta.samples[,c(3,28,30)]
sim.demmar <- seq(20,80,length.out=900)
x.mat <- cbind(1,sim.demmar,sim.demmar^2)#,sim.demmar^3)
colnames(x.mat) <- c("Intercept","DemMargin","DemMargin.2")#,"DemMargin.3")
Pred.DNC.samples <- t(aaply(x.mat,1,function(x)beta.samplesDNC%*%x))
Pred.RNC.samples <- t(aaply(x.mat,1,function(x)beta.samplesRNC%*%x))
Pred.DNC <- adply(Pred.DNC.samples,2,quantile,probs=c(0.05,0.5,0.95))
Pred.DNC$Convention <- "Being Exposed to DNC"
Pred.DNC$SimVals <- sim.demmar
names(Pred.DNC) <- c("X","Low","Median","High","Convention","SimVals")
Pred.RNC <- adply(Pred.RNC.samples,2,quantile,probs=c(0.07,0.5,0.95))
Pred.RNC$Convention <- "Being Exposed to RNC"
Pred.RNC$SimVals <- sim.demmar
names(Pred.RNC) <- c("X","Low","Median","High","Convention","SimVals")
Pred.Votes <- rbind(Pred.DNC,Pred.RNC)
Pred.Votes <- Pred.Votes[,-1]

##Plot Conditional effects
mytheme <- standard.theme("pdf", color=FALSE)
xyplot(Median~SimVals|Convention
       ,ylim=c(-7.5,7.5)
       ,ylab="Predicted Change in Democratic Percentage of Vote"
       ,xlab="Previous Democratic Percentage of the Vote"
       ,par.settings=mytheme
       ,scales=list(x=list(alternating=FALSE))
       ,panel=function(x,y,subscripts,...){
         panel.polygon(c(Pred.Votes$SimVals[subscripts],rev(Pred.Votes$SimVals[subscripts]))
                       ,c(Pred.Votes$High[subscripts],rev(Pred.Votes$Low[subscripts]))
                       ,col="gray70"
                       ,border="white")
         panel.lines(x,y,col="white",...,lwd=2.5)
         panel.abline(h=0,lty=2,color="gray50")
         
       }
       ,strip = function(..., bg,par.strip.text){
         strip.default(..., bg="black"
                       ,par.strip.text=list(alpha=1
                                            ,cex=1.05
                                            ,col="white"
                                            ,font=2
                       ))}
       ,data=Pred.Votes
)

#######Analyses with MSA
#Full dataset
BlackWhite <- read.dta("BlackWhiteDelta.dta")
unique.has.hosted <- with(ConventionsData,unique(subset(fips,ConventionParty!="N")))
ConventionsData$hosting.now <- ConventionsData$ConventionParty!="N"
ConventionsData$has.hosted <- ConventionsData$fips %in% unique.has.hosted
ConventionsDataMSA <- merge(ConventionsData,BlackWhite)

prop.model <-glm(hosting.now~has.hosted
                 +abs(LagDemMargin)
                 +as.factor(ICPSRstate)
                 + as.factor(year)                   
                 ,data=ConventionsDataMSA
                 ,family=binomial()
                 ,control=list(maxit=1e3))
prop.score <- predict(prop.model,type="response")
ConventionsDataMSA <- ConventionsDataMSA[rownames(model.frame(prop.model)),]
ConventionsDataMSA$prop.score  <- prop.score
DemMargin.modelMSA <- lm(DeltaDemPerc ~ 
                          x.DNC*LagDemPerc
                        +x.DNC*I(LagDemPerc^2)
                        +x.RNC*LagDemPerc
                        +x.RNC*I(LagDemPerc^2)
                        + I(DeltaBlack*100)
                        + I(DeltaWhite*100)
                        + I(DeltaUrban*100)
                        + I(DeltaFem*100)
                        + I(DeltaAsian*100)
                        + I(DeltaNative*100)
                        + I(DeltaHispanic*100)
                        + I(DeltaMarried*100)
                        + I(DeltaUnemploy*100)
                        + I(DeltaCollege*100)
                        + I(DeltaIncome/1000)
                        + prop.score
                        + as.factor(year)
                        ,data=ConventionsDataMSA)


#1972 onwards 
ConventionsDataMSA <- subset(ConventionsDataMSA,year>1970)
DemMargin.modelMSA72 <- lm(DeltaDemPerc ~ 
                           x.DNC*LagDemPerc
                         +x.DNC*I(LagDemPerc^2)
                         +x.RNC*LagDemPerc
                         +x.RNC*I(LagDemPerc^2)
                         + I(DeltaBlack*100)
                         + I(DeltaWhite*100)
                         + I(DeltaUrban*100)
                         + I(DeltaFem*100)
                         + I(DeltaAsian*100)
                         + I(DeltaNative*100)
                         + I(DeltaHispanic*100)
                         + I(DeltaMarried*100)
                         + I(DeltaUnemploy*100)
                         + I(DeltaCollege*100)
                         + I(DeltaIncome/1000)
                         + prop.score
                         + as.factor(year)
                         ,data=ConventionsDataMSA)

##Analysis with distance
DistDemMargin.model <- lm(DeltaDemPerc ~ 
                      +I(DistDNC/1.6/10)*LagDemPerc
                     + I((DistDNC/1.6/10)^2)
                        +I(DistRNC/10)*LagDemPerc
                       + I((DistRNC/10)^2)
                      + I(DeltaBlackProp*100)
                      + I(DeltaWhiteProp*100)
                      + I(DeltaBel25Prop*100)
                      + I(DeltaAb65Prop*100)
                      + I(DeltaUrban*100)
                      + I(DeltaFem*100)
                      + I(DeltaAsian*100)
                      + I(DeltaNative*100)
                      + I(DeltaHispanic*100)
                      + I(DeltaMarried*100)
                      + I(DeltaUnemploy*100)
                      + I(DeltaCollege*100)
                      + I(DeltaIncome/1000)
                      + prop.score
                      + as.factor(year)
                      ,weights=log(totvote)/(log(apply(cbind(DistRNC+1,DistDNC+1),1,min))+1)
                      ,data=ConventionsDataSub)

ConventionsDataSub <- ConventionsDataSub[rownames(model.frame(DistDemMargin.model)),]

CorrectedVCov <- cl(ConventionsDataSub,DistDemMargin.model,ConventionsDataSub$fips)
#CorrectedVCov <- vcov(DemMargin.model)

#Effect of DNC and RNC conditional on lagged margin
beta.samples <- mvrnorm(1000,coef(DistDemMargin.model),CorrectedVCov)
beta.samplesDNC <- beta.samples[,c(2,4,31)]
beta.samplesRNC <- beta.samples[,c(5,6,32)]
sim.dist <- seq(0,50,length.out=900)
x.mat <- cbind(1,sim.dist^2,10)#,sim.demmar^3)
colnames(x.mat) <- c("Intercept","DemMargin","DemMargin.2")#,"DemMargin.3")
Pred.DNC.samples <- t(aaply(x.mat,1,function(x)beta.samplesDNC%*%x))+0.08
Pred.RNC.samples <- t(aaply(x.mat,1,function(x)beta.samplesRNC%*%x))-0.08
Pred.DNC <- adply(Pred.DNC.samples,2,quantile,probs=c(0.05,0.5,0.95))
Pred.DNC$Convention <- "Effect of Exposure to DNC"
Pred.DNC$SimVals <- sim.dist
names(Pred.DNC) <- c("X","Low","Median","High","Convention","SimVals")
Pred.RNC <- adply(Pred.RNC.samples,2,quantile,probs=c(0.07,0.5,0.95))
Pred.RNC$Convention <- "Effect of Exposure to RNC"
Pred.RNC$SimVals <- sim.dist
names(Pred.RNC) <- c("X","Low","Median","High","Convention","SimVals")
Pred.Votes <- rbind(Pred.DNC,Pred.RNC)
Pred.Votes <- Pred.Votes[,-1]

##Plot Conditional effects
mytheme <- standard.theme("pdf", color=FALSE)
xyplot(Median~SimVals|Convention
       ,ylab="Predicted Change in Democratic Percentage of Vote"
       ,xlab="Distance From Convention Site (tens of miles)"
       ,par.settings=mytheme
       ,scales=list(x=list(alternating=FALSE))
       ,panel=function(x,y,subscripts,...){
         panel.polygon(c(Pred.Votes$SimVals[subscripts],rev(Pred.Votes$SimVals[subscripts]))
                       ,c(Pred.Votes$High[subscripts],rev(Pred.Votes$Low[subscripts]))
                       ,col="gray70"
                       ,border="white")
         panel.lines(x,y,col="black",...,lwd=2.5)
         panel.abline(h=0,lty=2,color="gray50")
         
       }
       ,strip = function(..., bg,par.strip.text){
         strip.default(..., bg="black"
                       ,par.strip.text=list(alpha=1
                                            ,cex=1.05
                                            ,col="white"
                                            ,font=2
                       ))}
       ,data=Pred.Votes
)


###Analysis from 1956 onwards using DMA

unique.has.hosted <- with(ConventionsData,unique(subset(fips,Convention!="N")))
hosting.now <- ConventionsData$Convention!="N"
ConventionsData$has.hosted <- ConventionsData$fips %in% unique.has.hosted

prop.model <-glm(hosting.now~has.hosted
                 +abs(LagDemMargin)
                 +as.factor(ICPSRstate)
                 + as.factor(year)                   
                 #+pop.density
                 ,data=ConventionsData
                 ,family=binomial()
                 ,control=list(maxit=1e3))
prop.score <- predict(prop.model,type="response")

ConventionsDataALL <- ConventionsData[rownames(model.frame(prop.model)),]
ConventionsDataALL$prop.score  <- prop.score
BlackWhite <- read.dta("BlackWhiteDelta.dta")
ConventionsDataALL <- merge(ConventionsDataALL,BlackWhite)


DemMargin.modelALL <- lm(DeltaDemPerc ~ 
                          Convention*LagDemPerc
                        +Convention*I(LagDemPerc^2)
                        + I(DeltaBlack*100)
                        + I(DeltaWhite*100)
                        + I(DeltaUrban*100)
                        + I(DeltaFem*100)
                        + I(DeltaAsian*100)
                        + I(DeltaNative*100)
                        + I(DeltaHispanic*100)
                        + I(DeltaMarried*100)
                        + I(DeltaUnemploy*100)
                        + I(DeltaCollege*100)
                        + I(DeltaIncome/1000)
                        + prop.score
                        + as.factor(year)
                        ,weights=log(totvote)/(log(apply(cbind(DistRNC+1,DistDNC+1),1,min))+1)
                        ,data=ConventionsDataALL)

ConventionsDataALL <- ConventionsDataALL[rownames(model.frame(DemMargin.modelALL)),]

CorrectedVCov <- cl(ConventionsDataALL,DemMargin.modelALL,ConventionsDataALL$fips)




