
# ------------------------------------------------------------------------------
# Purcell GxE model simulation
## ------------------------------------------------------------------------------						
# Revision History
# ------------------------------------------------------------------------------		
#  
# James Clifford -- 11 21 2013
#
# ------------------------------------------------------------------------------

rm(list = ls(all = TRUE))   				# Remove old objects

# Require packages
require(OpenMx)
require(MASS)
source('~/Desktop/Super wicked important R files/GenEpiHelperFunctions.R', chdir = TRUE)

## Assign pathway coefficients for ACE, set C = 0 for an AE model
asim <- sqrt(.8) 
csim <- sqrt(0)
esim <- sqrt(.2)

Asim <- asim*asim
Csim <- csim*csim
Esim <- esim*esim

## Generating the covariance matrices
covMZ<- rbind(cbind( Asim+Csim+Esim, Asim+Csim),
					cbind(Asim+Csim, Asim+Csim+Esim))
					
covDZ<- rbind(cbind( Asim+Csim+Esim, (0.5*Asim)+Csim),
					cbind((0.5*Asim)+Csim, Asim+Csim+Esim))
								
								
## Simulate data for twins with MASS

mzdata<-mvrnorm(1000, c(10,10), covMZ)
dzdata<-mvrnorm(1000, c(10,10), covDZ)
summary(mzdata)
summary(dzdata)
head(mzdata)
head(dzdata)


## Simulate environmental data								
envimz<- as.matrix(rnorm(1000,3,1))
envidz<- as.matrix(rnorm(1000,3,1))
head(envimz)					
mzyg<-rep(1,1000)
dzyg<-rep(2,1000)
famno<-rep(1:500,2)
famno2<-rep(501:1000, 2)
id<-rep(1:2,500)
id2<-rep(1:2,500)



## Add in environmental data to twin datat
datamz<-cbind(mzdata, envimz, mzyg, famno, id)
datadz<-cbind(dzdata, envidz,dzyg, famno2, id2)
colnames(datamz) <- c("Trait 1", "Trait 2", "Environment", "Zygosity", "FamID", "ID")
colnames(datadz) <-  c("Trait 1", "Trait 2", "Environment", "Zygosity", "FamID", "ID")
summary(datamz)
summary(datadz)
data<-rbind(datamz,datadz)
summary(data)

# Histograms of frequency distributions

par(mfrow=c(1,2))
hist(datamz[,3], main="Environment for MZ", xlab="Environmental Measure")
hist(datadz[,3], main="Environment for DZ", xlab="Environmental Measure")

cor(datamz)
cor(datadz)

## Prove means aren't different via eyeball test
mean(datamz[,1])
mean(datamz[,2])
mean(datadz[,1])
mean(datadz[,2])


## Regressions for abline based on means
reg1 <- glm(datamz[,1]~datamz[,2])
reg2 <- glm(datadz[,1]~datadz[,2])


## Plot data for correlations
par(mfrow=c(1,2))
plot(datamz[,1], datamz[,2], main="MZ Twins on Simulated Trait", xlab="Twin 1", ylab="Twin 2")
abline(reg1, col = "red")
plot(datadz[,1], datadz[,2], main="DZ Twins on Simulated Trait", xlab="Twin 1", ylab="Twin 2")
abline (reg2, col="red")

## Put data into OpenMx framework

mzdata<-mxData(datamz, type="raw")
dzdata<-mxData(datadz, type="raw")
summary(datamz)
## Create mx matricies

pathA <-mxMatrix(type = "Full", nrow = 1, ncol = 1, free = T, values = .6, label = "a11", name = "a")
pathC <-mxMatrix(type = "Full", nrow = 1, ncol = 1, free = T, values = .1, label = "c11",name = "c")
pathE <-mxMatrix(type = "Full", nrow = 1, ncol = 1, free = T, values = .3, label = "e11", name = "e")
# OK to leave them all at .6 because the square times 3 is roughly around 1, the variance of the data

covA <- mxAlgebra(a %*% t(a), name = "A")
covC <- mxAlgebra(c %*% t(c), name = "C")
covE <- mxAlgebra(e %*% t(e), name = "E")

expMeans <-mxMatrix(type = "Full", nrow = 1, ncol= 2, free=T, values = 10, labels = "mean", name = "expMeans")

expCovMZ <- mxAlgebra(rbind(cbind( A+C+E, A+C),
					cbind(A+C, A+C+E)), name="expCovMZ")
					
expCovDZ <- mxAlgebra(rbind(cbind( A+C+E, (0.5*A)+C),
					cbind((0.5*A)+C, A+C+E)),name="expCovDZ")

selvars<-c("Trait 1" , "Trait 2")
## Create Mx Objective
objMZ <-mxFIMLObjective(covariance ="expCovMZ" , means="expMeans" ,dimnames=selvars) 
objDZ <- mxFIMLObjective(covariance = "expCovDZ", means = "expMeans", dimnames =selvars)


## Combine groups

pars <- list (pathA, pathC, pathE, covA, covC, covE, expMeans)


## Put into a model

modelMZ<- mxModel(pars, expCovMZ, objMZ, mzdata, name = "MZ")
modelDZ <- mxModel(pars, expCovDZ, objDZ, dzdata, name= "DZ")
minus2ll<-mxAlgebra(expression=MZ.objective + DZ.objective, name="m2LL" )
obj	<-mxAlgebraObjective( "m2LL" )
ACEmodel <-mxModel("ACEmodel", pars, modelMZ, modelDZ, minus2ll, obj)

## Run the model

ACEfit<- mxRun(ACEmodel)
ACEsumm <- summary(ACEfit)
ACEsumm
tableFitStatistics(ACEfit)

## Drop C to double check

ACEmodelnoc <- mxModel(ACEfit, name="No C")
ACEmodelnoc<-omxSetParameters(ACEmodelnoc, labels=c("c11"), free=F, values=0)
ACEmodelnocfit<-mxRun(ACEmodelnoc)
ACEmodelnocfitsumm<-summary(ACEmodelnocfit)

tableFitStatistics(ACEfit, ACEmodelnocfit)

## Drop A...should worsen model fit

ACEmodelnoa <- mxModel(ACEfit, name="No A")
ACEmodelnoa<-omxSetParameters(ACEmodelnoa, labels=c("a11"), free=F, values=0)
ACEmodelnoafit<-mxRun(ACEmodelnoa)
ACEmodelnoafitsumm<-summary(ACEmodelnoafit)

satnested<-list(ACEmodelnocfit, ACEmodelnoafit)

tableFitStatistics(ACEfit, satnested)


# ------------------------------------------------------------------------------ #
# PREPARE Bivariate ACE MODEL
# ------------------------------------------------------------------------------ #

# Split data by twins

twindat <- function(dat, famid, twinid, zygosity) {
  datA <- dat[dat[,twinid]==min(dat[,twinid]),]    #twin1
  datB <- dat[dat[,twinid]==max(dat[,twinid]),]    #twin2 
  DAT <- merge(datA, datB, by=famid, all.x=TRUE, all.y=TRUE, suffixes=c("_T1","_T2"))
  DAT[,paste(twinid,"_T1",sep="")] <- NULL
  DAT[,paste(twinid,"_T2",sep="")] <- NULL
  DAT[,zygosity] <- ifelse(is.na(DAT[,paste(zygosity,"_T1",sep="")]),DAT[,paste(zygosity,"_T2",sep="")],DAT[,paste(zygosity,"_T1",sep="")])
  DAT[,paste(zygosity,"_T1",sep="")] <- NULL  
  DAT[,paste(zygosity,"_T2",sep="")] <- NULL  
  return(DAT)
}

dataR12 <- twindat(dat=data, famid= "FamID", twinid= "ID", zygosity= "Zygosity")
varnames <- colnames(dataR12)
summary(dataR12)
selvars<-c("Trait 1_T1", "Trait 2_T1", "Trait 1_T2", "Trait 2_T2")

datamz<-subset(dataR12, Zygosity==1, selvars)
datadz<-subset(dataR12, Zygosity==2, selvars)

nv<-2 # number of variables per twin
ntv<-nv*2 # Number of variables per pair

mzdata<-mxData(datamz, type="raw")
dzdata<-mxData(datadz, type="raw")

pathA <-mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = T, values = .6, label = c("a11","a21","a22"), name = "a")
pathC <-mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = T, values = .1, label = c("c11","c21","c22") ,name = "c")
pathE <-mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = T, values = .3, label = c("e11","e21","e22") , name = "e")


covA <- mxAlgebra(a %*% t(a), name = "A")
covC <- mxAlgebra(c %*% t(c), name = "C")
covE <- mxAlgebra(e %*% t(e), name = "E")

expMeans <-mxMatrix(type = "Full", nrow = 1, ncol= ntv, free=T, values = 10, labels = c("mean Trait 1 MZ", "mean trait 2 MZ", "mean trait 1 DZ", "mean trait 2 DZ"), name = "expMeans")

expCovMZ <- mxAlgebra(rbind(cbind( A+C+E, A+C),
					cbind(A+C, A+C+E)), name="expCovMZ")
					
expCovDZ <- mxAlgebra(rbind(cbind( A+C+E, (0.5%x%A)+C),
					cbind((0.5%x%A)+C, A+C+E)),name="expCovDZ")

## Create Mx Objective
objMZ <-mxFIMLObjective(covariance ="expCovMZ" , means="expMeans" ,dimnames=selvars) 
objDZ <- mxFIMLObjective(covariance = "expCovDZ", means = "expMeans", dimnames =selvars)


## Combine groups

pars <- list (pathA, pathC, pathE, covA, covC, covE, expMeans)

## Put into a model

modelMZ<- mxModel(pars, expCovMZ, objMZ, mzdata, name = "MZ")
modelDZ <- mxModel(pars, expCovDZ, objDZ, dzdata, name= "DZ")
minus2ll<-mxAlgebra(expression=MZ.objective + DZ.objective, name="m2LL" )
obj	<-mxAlgebraObjective( "m2LL" )
ACEbimodel <-mxModel("ACEbimodel", pars, modelMZ, modelDZ, minus2ll, obj)

## Run the model

ACEfitbi<- mxRun(ACEbimodel)
ACEsummbi <- summary(ACEfitbi)
ACEsummbi
tableFitStatistics(ACEfitbi)

## Drop C from the bivariate

bimodelnoc <- mxModel(ACEfitbi, name="No C")
bimodelnoc<-omxSetParameters(bimodelnoc, labels=c("c11", "c21", "c22"), free=F, values=0)
bimodelnocfit<-mxRun(bimodelnoc)
bimodelnocfitsumm<-summary(bimodelnocfit)

tableFitStatistics(ACEfitbi, bimodelnocfit)

# ------------------------------------------------------------------------------		
#  
# Prepare Moderation model
#
# ------------------------------------------------------------------------------



# Algebra to compute total variances and standard deviations (diagonal only)
covP      <- mxAlgebra( expression=A+C+E, name="V" )
matI      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I")
invSD     <- mxAlgebra( expression=solve(sqrt(I*V)), name="iSD")

# Algebras generated to create summary Table of Derived Variance Components
rowVars   <- rep('vars',nv)
colVars   <- rep(c('A','C','E','SA','SC','SE'),each=nv)
estVars   <- mxAlgebra( expression=cbind(A,C,E,A/V,C/V,E/V), name="Vars", dimnames=list(rowVars,colVars) )

# Matrix for expected Mean 
meanG     <- mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values= 10, label="mean", name="Mean" )

# Matrix for moderating/interacting variable
def    <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.Trait 1"), name="Trait1") # Moderation of the variance components
def2   <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.Trait 1"), name="Trait1") # for moderation of the means 

# Matrices declared to store moderated a, c, and e Path Coefficients
pathAI    <- mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.6, label=c("aI11","aI21", "aI22"), name="aI" ) 
pathCI    <- mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.1, label=c("cI11", "cI21", "cI22"),name="cI" )
pathEI    <- mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.3, label=c("eI11","eI21", "eI22"), name="eI" )

# Matrices declared to store linear and quadratic Regression Coefficients for covariate
pathB     <- mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values= .01, label=c("l11","q11"), name="b" )

# Matrices generated to hold moderated A, C, and E computed Variance Components
covAI     <- mxAlgebra( expression=(a+ def%*%aI) %*% t(a+ def%*%aI), name="AI" )
covCI     <- mxAlgebra( expression=(c+ def%*%cI) %*% t(c+ def%*%cI), name="CI" )
covEI     <- mxAlgebra( expression=(e+ def%*%eI) %*% t(e+ def%*%eI), name="EI" )

# Algebra for expected Mean and Variance/Covariance Matrices in MZ & DZ twins
meanAge   <- mxAlgebra( expression= b%*%rbind(def,def2), name="Mod")
meanGI    <- mxAlgebra( expression= cbind((Mean + Mod),(Mean + Mod)), name="expMean")
covMZI    <- mxAlgebra( expression= rbind ( cbind(AI+CI+EI , AI+CI),
                                            cbind(AI+CI   , AI+CI+EI)), name="expCovMZ" )
covDZI    <- mxAlgebra( expression= rbind( cbind(AI+CI+EI     , 0.5%x%AI+CI),
                                           cbind(0.5%x%AI+CI , AI+CI+EI)), name="expCovDZ" )

# Objective objects for Multiple Groups
objMZ     <- mxFIMLObjective( covariance="expCovMZ", means="expMean", dimnames=selvars )
objDZ     <- mxFIMLObjective( covariance="expCovDZ", means="expMean", dimnames=selvars )

# Matrices & Algebra to plot Means and Variance Components by age (not required for model fitting)
ages      <- mxMatrix( type="Full", nrow=5, ncol=1, values=c(3,6,8,10,12), name="Mod")
age       <- mxMatrix( type="Full", nrow=2, ncol=5, values=c(.3,.09,.4,.16,.5,.25,.6,.36,.7,.49), name="Agelq")
unit      <- mxMatrix( type="Unit", nrow=5, ncol=1, name="unit")
meanI     <- mxAlgebra( expression=unit%x%Mean+ t(b%*%Agelq), name="Mi")
compAI    <- mxAlgebra( expression=diag2vec((unit%x%a+ Mod%x%aI) %*% t(unit%x%a+ AgeMods%x%aI)), name="Ai" )
compCI    <- mxAlgebra( expression=diag2vec((unit%x%c+ Mod%x%cI) %*% t(unit%x%c+ Mod%*%cI)), name="Ci" )
compEI    <- mxAlgebra( expression=diag2vec((unit%x%e+ Mod%x%eI) %*% t(unit%x%e+ Mod%*%eI)), name="Ei" )
compPI    <- mxAlgebra( expression=Ai+Ci+Ei, name="Vi" )

# Algebras generated to create summary Table of Derived Variance Components by Age
rowsAge   <- c("Mod3","Mod6","Mod8","Mod10","Mod12")
colsAge   <- c("MI","AI","CI","EI","VI")
estAge    <- mxAlgebra( expression=cbind(Mi,Ai,Ci,Ei,Vi), name="byAge", dimnames=list(rowsAge,colsAge) )
estPrAge  <- mxAlgebra( expression=cbind(Mi,Ai/Vi,Ci/Vi,Ei/Vi,Vi), name="byPrAge", dimnames=list(rowsAge,colsAge) )
propAge   <- list( ages, age, unit, meanI, compAI, compCI, compEI, compPI, estAge, estPrAge)

# Combine Groups
pars      <- list( pathA, pathC, pathE, pathAI, pathCI, pathEI, pathB, meanG, covA, covC, covE, covP, matI, invSD, estVars )
defs      <- list( def,def2, covAI, covCI, covEI, meanAge)
modelMZ   <- mxModel( pars, defs, meanGI, covMZI, datamz, objMZ, name="MZ" )
modelDZ   <- mxModel( pars, defs, meanGI, covDZI, datadz, objDZ, name="DZ" )
minus2ll  <- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj       <- mxAlgebraObjective( "m2LL" )
ModAceLqModel  <- mxModel( "ModACElq", pars, modelMZ, modelDZ, minus2ll, obj, propAge )

# ------------------------------------------------------------------------------
# RUN MODELS

# Run Moderated ACE model + Linear & Quadratic Moderated Means
ModAceLqFit    <- mxRun(ModAceLqModel)
ModAceLqSumm   <- summary(ModAceLqFit)
ModAceLqSumm
round(ModAceLqFit@output$estimate,4)   ## Gives us the path coefficients and means, rounded to 4 decimals
round(ModAceLqFit$Vars@result,4)	   ## Gives us the square path coefficients & unmoderated standardized coefficients

round(ModAceLqFit$byAge@result,4) 	   ## Gives us the moderated path coefficients by Age (specific values of interest), linear and quadratic components included: mean changes are significiant, path estimates are not (see below)

# Generate Output with Functions
source("GenEpiHelperFunctions.R")
parameterSpecifications(ModAceLqFit)
expectedMeansCovariances(ModAceLqFit)
tableFitStatistics(ModAceLqFit)

# ------------------------------------------------------------------------------
# FIT SUBMODELS

# Run non-Moderated ACE model
OneAceLqModel  <- mxModel( ModAceLqFit, name="OneACElq" )
OneAceLqModel  <- omxSetParameters( OneAceLqModel, labels=c("aI11","cI11","eI11"), free=FALSE, values=0 )
OneAceLqFit    <- mxRun(OneAceLqModel)
round(OneAceLqFit@output$estimate,4)
round(OneAceLqFit$Vars@result,4)
round(OneAceLqFit$byAge@result,4)

LRTmod <- mxCompare(ModAceLqFit,list(OneAceLqFit, ModAceLinFit), all=T)
LRT2 <-mxCompare(OneAceLqFit, ModAceLinFit)



# Fit Moderated ACE model + Linear Moderated Means
ModAceLinModel <- mxModel( OneAceLqFit, name="ModACEl" )
ModAceLinModel <- omxSetParameters( ModAceLinModel, labels="q11", free=FALSE, values=0 )
ModAceLinFit   <- mxRun(ModAceLinModel)
round(ModAceLinFit@output$estimate,4)
round(ModAceLinFit$Vars@result,4)
round(ModAceLinFit$byAge@result,4)

# Fit Moderated ACE model + no Moderated Means
ModAceModel    <- mxModel( ModAceLinFit, name="ModACE" )
ModAceModel    <- omxSetParameters( ModAceModel, labels="l11", free=FALSE, values=0 )
ModAceFit      <- mxRun(ModAceModel)
round(ModAceFit@output$estimate,4)
round(ModAceFit$Vars@result,4)
round(ModAceFit$byAge@result,4)

# ------------------------------------------------------------------------------

# Print Comparative Fit Statistics
AceNested <- list(OneAceLqFit, ModAceLinFit, ModAceFit)
tableFitStatistics(ModAceLqFit,AceNested)

round(rbind(ModAceLqFit@output$estimate,OneAceLqFit@output$estimate,ModAceLinFit@output$estimate,ModAceFit@output$estimate),4)




