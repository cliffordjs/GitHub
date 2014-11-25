rm(list = ls())

# ------------------------------------------------------------------------------
# Written from scratch with Brad using Hermine's scripts form AGES workshop off vipbg.vcu.edu
# -------|---------|---------|---------|---------|---------|---------|

# Load Library
source('http://openmx.psyc.virginia.edu/getOpenMx.R')
source('http://openmx.psyc.virginia.edu/getOpenMxBeta.R')
source("http://www.vipbg.vcu.edu/tc2012/GenEpiHelperFunctions.R") # Update regularly from VIPBG website
require(OpenMx)
require(psych)
mxOption(NULL, "Default optimizer", "NPSOL")
#mxOption(NULL, "Default optimizer", "CSOLNP")


# PREPARE DATA
#setwd("/Users/ngillespie/Documents/work/papers/report_GxE")
data1<-read.csv("~/Desktop/ResearchProjects/Hettema/GxE/0. Initial Data/anx_ff_pairs_jc.csv", na="-99", header = TRUE, stringsAsFactors = FALSE)
data1<-data1[order(data1$famno),]
data1$famno<-data1$famno-1000000
data2<- read.csv("~/Desktop/ResearchProjects/Hettema/GxE/0. Initial Data/anx_pheno_jc.csv", na="-99", header = TRUE, stringsAsFactors = FALSE)
data2<-data2[order(data2$famno),]

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

data3<-twindat(dat=data2, famid="famno", twinid="id", zygosity="zygosity")
data3<-data3[order(data3$famno),]
data1<-data1[order(data1$famno),]	
data3<-data3[order(data3$famno),]	
findat<-merge(data1,data3,by=c("famno"))
findat<-findat[order(findat$famno),]


#----------------------------------------------------------##
#
# Parental loss											   ##	
#
#----------------------------------------------------------##


# Select variables
selVars<- c("t1_par_loss","mdd_T1", "t2_par_loss", "mdd_T2")
mzdata<-subset(findat, zygosity==1, selVars)
dzdata<-subset(findat, zygosity==2, selVars)

mzdata$mdd_T1<-mxFactor(mzdata$mdd_T1, levels=c(0:1))
mzdata$mdd_T2<-mxFactor(mzdata$mdd_T2, levels=c(0:1))
dzdata$mdd_T1<-mxFactor(dzdata$mdd_T1, levels=c(0:1))
dzdata$mdd_T2<-mxFactor(dzdata$mdd_T2, levels=c(0:1))
mzdata$t1_par_loss<-mxFactor(mzdata$t1_par_loss, levels=c(0:1))
mzdata$t2_par_loss<-mxFactor(mzdata$t2_par_loss, levels=c(0:1))
dzdata$t1_par_loss<-mxFactor(dzdata$t1_par_loss, levels=c(0:1))
dzdata$t2_par_loss<-mxFactor(dzdata$t2_par_loss, levels=c(0:1))

# Starting Values
nv        <- 2       # number of variables
ntv       <- nv*2    # number of total variables
thVals 	  <- .5
svPa      <- .6
svPas     <- diag(svPa,nv,nv)
lth       <- -1.5    # start value for first threshold
ith       <- 1       # start value for increments
nth 	  	<- 1	
thVal     <- matrix(rep(c(lth,(rep(ith,nth-1)))),nrow=nth,ncol=nv)
thLabMZ   <- c(paste("t",1:nth,"MZ1",sep=""),paste("t",1:nth,"MZ2",sep=""))
thLabDZ   <- c(paste("t",1:nth,"DZ1",sep=""),paste("t",1:nth,"DZ2",sep=""))
thLB      <- matrix(rep(c(-3,(rep(0.001,nth-1))),nv),nrow=nth,ncol=nv)
lbrVal    <- -0.99   # start value for lower bounds
ubrVal    <- 0.99    # start value for upper bounds

# Cholesky 

# Set Starting Values
paVal     <- .6                        # start value for path coefficient
paValD    <- vech(diag(paVal,nv,nv))   # start values for diagonal of covariance matrix
peVal     <- .8                        # start value for path coefficient for e
peValD    <- vech(diag(paVal,nv,nv))   # start values for diagonal of covariance matrix
paLBo     <- .0001                     # start value for lower bounds
paLBoD    <- diag(paLBo,nv,nv)         # lower bounds for diagonal of covariance matrix
paLBoD[lower.tri(paLBoD)] <- -10       # lower bounds for below diagonal elements
paLBoD[upper.tri(paLBoD)] <- NA        # lower bounds for above diagonal elements

# Create Labels for Lower Triangular Matrices
aLabs     <- paste("a",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")
cLabs     <- paste("c",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")
eLabs     <- paste("e",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")

# Matrices declared to store a, c, and e Path Coefficients
pathA     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(0.1, 0.5,0.3), labels=aLabs, lbound=-.99, ubound=.99, name="a" )
pathC     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(0.9,0.1,0.1), labels=cLabs, lbound=-.99, ubound=2, name="c" )
pathE     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(0.1,0.4,0.6), labels=eLabs, lbound=-.99, ubound=.99, name="e" )
	
# Matrices generated to hold A, C, and E computed Variance Components
covA      <- mxAlgebra( expression=a %*% t(a), name="A" )
covC      <- mxAlgebra( expression=c %*% t(c), name="C" ) 
covE      <- mxAlgebra( expression=e %*% t(e), name="E" )

# Algebra to compute total variances and standard deviations (diagonal only)
covP      <- mxAlgebra( expression=A+C+E, name="V" )
matI      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I")
invSD     <- mxAlgebra( expression=solve(sqrt(I*V)), name="iSD")

# Constraint on variance of Ordinal variables
matUnv   <- mxMatrix( type="Unit", nrow=nv, ncol=1, name="Unv1" )
var1     <- mxConstraint( expression=diag2vec(V)==Unv1, name="Var1" )

# Matrix & Algebra for expected means vector and expected thresholds
meanG    <- mxMatrix( type="Zero", nrow=1, ncol=nv, name="Mean" )
meanT    <- mxAlgebra( expression=cbind(Mean,Mean), name="expMean" )
threG    <- mxMatrix( type="Full", nrow=nth, ncol=nv, free=T, values = c(0.9,0.2),lbound=-10, ubound=10, labels=c("par_loss","mdd"), name="Thre" )
Inc      <- mxMatrix( type="Lower", nrow=nth, ncol=nth, free=FALSE, values=1, name="Inc" )
threT    <- mxAlgebra( expression=cbind(Inc %*% Thre,Inc %*% Thre), name="expThre" )

# Algebra for expected Mean and Variance/Covariance Matrices in MZ & DZ twins
covMZ    <- mxAlgebra( expression= rbind( cbind(A+C+E , A+C),       cbind(A+C, A+C+E)), name="expCovMZ" )
covDZ    <- mxAlgebra( expression= rbind( cbind(A+C+E , 0.5%x%A+C), cbind(0.5%x%A+C , A+C+E)), name="expCovDZ" )

# Data objects for Multiple Groups
dataMZ    <- mxData( observed=mzdata, type="raw" )
dataDZ    <- mxData( observed=dzdata, type="raw" )

# Objective objects for Multiple Groups
#objMZ    <- mxFIMLObjective( covariance="expCovMZ", means="expMean", dimnames=selVars, thresholds="expThre" )
#objDZ    <- mxFIMLObjective( covariance="expCovDZ", means="expMean", dimnames=selVars, thresholds="expThre" )
objMZ    <- mxExpectationNormal( covariance="expCovMZ", means="expMean", thresholds="expThre", dimnames=selVars )
objDZ    <- mxExpectationNormal( covariance="expCovDZ", means="expMean", thresholds="expThre", dimnames=selVars )
fiML     <- mxFitFunctionML()

# Combine Groups
pars      <- list( pathA, pathC, pathE, covA, covC, covE, covP, matI, invSD, matUnv, var1, meanG, meanT, threG, Inc, threT )
modelMZ   <- mxModel( dataMZ, pars, covMZ, objMZ, fiML, name="MZ" )
modelDZ   <- mxModel( dataDZ, pars, covDZ, objDZ, fiML, name="DZ" )
minus2ll  <- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
#obj      <- mxAlgebraObjective( "m2LL" )
obj       <- mxFitFunctionAlgebra( "m2LL" )
CholAceModelpar  <- mxModel( "CholACE_Par_loss", modelMZ, modelDZ, minus2ll, obj)

CholAceFitpar   <- mxTryHard(CholAceModelpar,greenOK=FALSE,scale=0.0,extraTries=20)
CholAceFitpar    <- mxRun(CholAceFitpar)
CholAceSummpar   <- summary(CholAceFitpar)
CholAceSummpar
tableFitStatistics(CholAceFitpar)


# Cholesky with GbyE

mzdata <- na.omit(mzdata)
dzdata <- na.omit(dzdata)

# Declared pathway coefficients (a,c,e) & regression weights for the moderated pathways (aB,cB,eB)
a     	<- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=c(0.3,0.4,0.1), labels=aLabs, lbound=-.99, ubound=.99, name="a" )
c   	<- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=c(0.7,0.2,0.1), labels=cLabs, lbound=-.99, ubound=.99, name="c" )
e	    <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=c(0.2,0.1,0.1), labels=eLabs, lbound=-.99, ubound=.99, name="e" )
aB     	<- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F, T, T), values=c(0,0.9,0.9), label=c("aMod11", "aMod21", "aMod22"), lbound = -3, ubound = 3, name="aB" ) 
cB     	<- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F, T, T), values=c(0,0.9,0.2), label=c("cMod11", "cMod21", "cMod22"), lbound = -3, ubound = 3, name="cB" )
eB     	<- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F, T, T), values=c(0,0.9,0.9), label=c("eMod11", "eMod21", "eMod22"), lbound = -3, ubound = 3, name="eB" )

# Declare means & thresholds: http://openmx.psyc.virginia.edu/wiki/matrix-operators-and-functions
means    	<- mxMatrix( type="Zero", nrow=1, ncol=nv, name="means" )
beta    	<- mxMatrix( type="Full", nrow = 1, ncol = 2, free = c(T,T), labels=c("b","b"), values=1.7, lbound = -2.5, ubound = 4.5, name="beta" )  
defvarT1 	<- mxMatrix( type="Full", nrow = 1, ncol = 1, free = F, labels = "data.t1_par_loss", name = "defvarT1")
defvarT2 	<- mxMatrix( type="Full", nrow = 1, ncol = 1, free = F, labels = "data.t2_par_loss", name = "defvarT2")
modthrT1	<- mxAlgebra( expression = defvarT1 %x% beta, name = "modthrT1")	
modthrT2	<- mxAlgebra( expression = defvarT2 %x% beta, name = "modthrT2")  	
expM  	 	<- mxAlgebra( expression= cbind(means + modthrT1,means + modthrT2), name="expM" ) 


thresh    	<- mxMatrix( type="Full", nrow=nth, ncol=nv, free=T,  values=c(4.9,4.8), lbound=-5, ubound=5, labels=c("par_loss","mdd"), name="thresh" )
expThre 	<- mxAlgebra( expression = cbind(thresh, thresh), name="expThre")

# Declare vVariance components
VarT1	<- mxAlgebra( expression = (a+defvarT1%x%aB)%*%t(a+defvarT1%x%aB) + (c+defvarT1%x%cB)%*%t(c+defvarT1%x%cB) + (e+defvarT1%x%eB)%*%t(e+defvarT1%x%eB), name="VarT1" )
VarT2	<- mxAlgebra( expression = (a+defvarT2%x%aB)%*%t(a+defvarT2%x%aB) + (c+defvarT2%x%cB)%*%t(c+defvarT2%x%cB) + (e+defvarT2%x%eB)%*%t(e+defvarT2%x%eB), name="VarT2" )
covA	<- mxAlgebra( expression = (a+defvarT1%x%aB)%*%t(a+defvarT2%x%aB) , name="covA" )
covC	<- mxAlgebra( expression = (c+defvarT1%x%cB)%*%t(c+defvarT2%x%cB) , name="covC" )

# Constraint on variance of Ordinal variables
unitM  		<- mxMatrix( type="Unit", nrow=4, ncol=1, name="unitM" )
Var1     	<- mxConstraint( expression = rbind(diag2vec(VarT1),diag2vec(VarT2))==unitM, name="Var1" )
covMZ   	<- mxAlgebra( expression = rbind( cbind(VarT1, covA+covC), cbind(covA+covC , VarT2)), name="expCovMZ" )
covDZ    <- mxAlgebra( expression = rbind( cbind(VarT1,0.5%x%covA+covC), cbind(0.5%x%covA+covC, VarT2)), name="expCovDZ" )

# Data objects for Multiple Groups
dataMZ    <- mxData( observed=mzdata, type="raw" )
dataDZ    <- mxData( observed=dzdata, type="raw" )

# Objective objects for Multiple Groups
objMZ    <- mxExpectationNormal( covariance="expCovMZ", means="expM", thresholds="expThre", dimnames=selVars )
objDZ    <- mxExpectationNormal( covariance="expCovDZ", means="expM", thresholds="expThre", dimnames=selVars )
fiML     <- mxFitFunctionML()

# Combine Groups
pars      	<- list( a,c,e,aB,cB,eB, thresh,expThre  )
defs		<- list( means,beta,defvarT1,defvarT2,modthrT1,modthrT2,expM, VarT1,VarT2,covA,covC,unitM,Var1 )
modelMZ   	<- mxModel( pars, defs, dataMZ, covMZ, objMZ, fiML, name="MZ" )
modelDZ   	<- mxModel( pars, defs, dataDZ, covDZ, objDZ, fiML, name="DZ" )
minus2ll  	<- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj       	<- mxFitFunctionAlgebra( "m2LL" )
ModAceModelpar <- mxModel( "ModACE_Par_Loss", pars, modelMZ, modelDZ, minus2ll, obj )

ModAceFitpar    <- mxTryHard(ModAceModelpar,greenOK=FALSE)
ModAceSummpar   <- summary(ModAceFitpar)
ModAceSummpar
tableFitStatistics(ModAceFitpar,CholAceFitpar)

model3<-mxModel(ModAceModelpar)
model3<-omxSetParameters(model3, labels=c("b", "aMod21","aMod22", "cMod21", "cMod22","eMod21","eMod22"), values=0, free=F) 
model3<-omxSetParameters(model3, labels="c11", values=0.8, free=T)

model3fit<-mxTryHard(model3, greenOK=F, extraTries=20)
summary(model3fit)
?mxTryHard
eigen(model3fit$MZ$expCovMZ$result)
eigen(model3fit$DZ$expCovDZ$result)
tableFitStatistics(ModAceFitpar, model3fit)
#----------------------------------------------------------##
#
# Parental Warmth										   ##	
#
#----------------------------------------------------------##



# Select variables
selVars<- c("t1_par_warm_ord","mdd_T1", "t2_par_warm_ord", "mdd_T2")
mzdata<-subset(findat, zygosity==1, selVars)
dzdata<-subset(findat, zygosity==2, selVars)
table(findat$t1_par_warm_ord, findat$t2_par_warm_ord)

mzdata$mdd_T1<-mxFactor(mzdata$mdd_T1, levels=c(0:1))
mzdata$mdd_T2<-mxFactor(mzdata$mdd_T2, levels=c(0:1))
dzdata$mdd_T1<-mxFactor(dzdata$mdd_T1, levels=c(0:1))
dzdata$mdd_T2<-mxFactor(dzdata$mdd_T2, levels=c(0:1))
mzdata$t1_par_warm_ord<-mxFactor(mzdata$t1_par_warm_ord, levels=c(0:4))
mzdata$t2_par_warm_ord<-mxFactor(mzdata$t2_par_warm_ord, levels=c(0:4))
dzdata$t1_par_warm_ord<-mxFactor(dzdata$t1_par_warm_ord, levels=c(0:4))
dzdata$t2_par_warm_ord<-mxFactor(dzdata$t2_par_warm_ord, levels=c(0:4))

# Starting Values
nv        <- 2       # number of variables
ntv       <- nv*2    # number of total variables
thVals 	  <- .5
svPa      <- .6
svPas     <- diag(svPa,nv,nv)
lth       <- -1.5    # start value for first threshold
ith       <- 1       # start value for increments
nth 	  	<- 4	
thVal     <- matrix(rep(c(lth,(rep(ith,nth-1)))),nrow=nth,ncol=nv)
thLabMZ   <- c(paste("t",1:nth,"MZ1",sep=""),paste("t",1:nth,"MZ2",sep=""))
thLabDZ   <- c(paste("t",1:nth,"DZ1",sep=""),paste("t",1:nth,"DZ2",sep=""))
thLB      <- matrix(rep(c(-3,(rep(0.001,nth-1))),nv),nrow=nth,ncol=nv)
lbrVal    <- -0.99   # start value for lower bounds
ubrVal    <- 0.99    # start value for upper bounds

# Cholesky 

# Set Starting Values
paVal     <- .6                        # start value for path coefficient
paValD    <- vech(diag(paVal,nv,nv))   # start values for diagonal of covariance matrix
peVal     <- .8                        # start value for path coefficient for e
peValD    <- vech(diag(paVal,nv,nv))   # start values for diagonal of covariance matrix
paLBo     <- .0001                     # start value for lower bounds
paLBoD    <- diag(paLBo,nv,nv)         # lower bounds for diagonal of covariance matrix
paLBoD[lower.tri(paLBoD)] <- -10       # lower bounds for below diagonal elements
paLBoD[upper.tri(paLBoD)] <- NA        # lower bounds for above diagonal elements

# Create Labels for Lower Triangular Matrices
aLabs     <- paste("a",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")
cLabs     <- paste("c",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")
eLabs     <- paste("e",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")

# Matrices declared to store a, c, and e Path Coefficients
pathA     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(0.3, 0.5,0.3), labels=aLabs, lbound=-.99, ubound=.99, name="a" )
pathC     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(0.2,0.1,0.1), labels=cLabs, lbound=-.99, ubound=2, name="c" )
pathE     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(0.5,0.4,0.6), labels=eLabs, lbound=-.99, ubound=.99, name="e" )
	
# Matrices generated to hold A, C, and E computed Variance Components
covA      <- mxAlgebra( expression=a %*% t(a), name="A" )
covC      <- mxAlgebra( expression=c %*% t(c), name="C" ) 
covE      <- mxAlgebra( expression=e %*% t(e), name="E" )

# Algebra to compute total variances and standard deviations (diagonal only)
covP      <- mxAlgebra( expression=A+C+E, name="V" )
matI      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I")
invSD     <- mxAlgebra( expression=solve(sqrt(I*V)), name="iSD")

# Constraint on variance of Ordinal variables
matUnv   <- mxMatrix( type="Unit", nrow=nv, ncol=1, name="Unv1" )
var1     <- mxConstraint( expression=diag2vec(V)==Unv1, name="Var1" )

# Matrix & Algebra for expected means vector and expected thresholds
meanG    <- mxMatrix( type="Zero", nrow=1, ncol=nv, name="Mean" )
meanT    <- mxAlgebra( expression=cbind(Mean,Mean), name="expMean" )
threG    <- mxMatrix( type="Full", nrow=nth, ncol=nv, free=c(T,T,T,T, T,F,F,F),  values=c(0.1, 0.2, 0.3, 0.99, 4.8, 0,0,0),lbound=-10, ubound=10, labels=c("par_warm1","par_warm2","par_warm3","par_warm4","mdd","mdd2","mdd3","mdd4"), name="Thre" )
Inc      <- mxMatrix( type="Lower", nrow=nth, ncol=nth, free=FALSE, values=1, name="Inc" )
threT    <- mxAlgebra( expression=cbind(Inc %*% Thre,Inc %*% Thre), name="expThre" )

# Algebra for expected Mean and Variance/Covariance Matrices in MZ & DZ twins
covMZ    <- mxAlgebra( expression= rbind( cbind(A+C+E , A+C),       cbind(A+C, A+C+E)), name="expCovMZ" )
covDZ    <- mxAlgebra( expression= rbind( cbind(A+C+E , 0.5%x%A+C), cbind(0.5%x%A+C , A+C+E)), name="expCovDZ" )

# Data objects for Multiple Groups
dataMZ    <- mxData( observed=mzdata, type="raw" )
dataDZ    <- mxData( observed=dzdata, type="raw" )

# Objective objects for Multiple Groups
#objMZ    <- mxFIMLObjective( covariance="expCovMZ", means="expMean", dimnames=selVars, thresholds="expThre" )
#objDZ    <- mxFIMLObjective( covariance="expCovDZ", means="expMean", dimnames=selVars, thresholds="expThre" )
objMZ    <- mxExpectationNormal( covariance="expCovMZ", means="expMean", thresholds="expThre", dimnames=selVars )
objDZ    <- mxExpectationNormal( covariance="expCovDZ", means="expMean", thresholds="expThre", dimnames=selVars )
fiML     <- mxFitFunctionML()

# Combine Groups
pars      <- list( pathA, pathC, pathE, covA, covC, covE, covP, matI, invSD, matUnv, var1, meanG, meanT, threG, Inc, threT )
modelMZ   <- mxModel( dataMZ, pars, covMZ, objMZ, fiML, name="MZ" )
modelDZ   <- mxModel( dataDZ, pars, covDZ, objDZ, fiML, name="DZ" )
minus2ll  <- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
#obj      <- mxAlgebraObjective( "m2LL" )
obj       <- mxFitFunctionAlgebra( "m2LL" )
CholAceModelwarm  <- mxModel( "CholACE_Par_warm", modelMZ, modelDZ, minus2ll, obj)

CholAceFitwarm   <- mxTryHard(CholAceModelwarm,greenOK=FALSE,scale=0.0,extraTries=20)
CholAceFitwarm    <- mxRun(CholAceFitwarm)
CholAceSummwarm   <- summary(CholAceFitwarm)
CholAceSummwarm
tableFitStatistics(CholAceFitwarm)

# Cholesky with GbyE

mzdata <- na.omit(mzdata)
dzdata <- na.omit(dzdata)

# Declared pathway coefficients (a,c,e) & regression weights for the moderated pathways (aB,cB,eB)
a     	<- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=c(0.1,0.3,0.5), labels=aLabs, lbound=-.99, ubound=.99, name="a" )
c   	<- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=c(0.9,0.1,0.1), labels=cLabs, lbound=-.99, ubound=.99, name="c" )
e	    <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=c(0.2,0.1,0.1), labels=eLabs, lbound=-.99, ubound=.99, name="e" )
aB     	<- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F, T, T), values=c(0,0.9,0.9), label=c("aMod11", "aMod21", "aMod22"), lbound = -3, ubound = 3, name="aB" ) 
cB     	<- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F, T, T), values=c(0,0.9,0.2), label=c("cMod11", "cMod21", "cMod22"), lbound = -3, ubound = 3, name="cB" )
eB     	<- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F, T, T), values=c(0,0.9,0.9), label=c("eMod11", "eMod21", "eMod22"), lbound = -3, ubound = 3, name="eB" )

# Declare means & thresholds: http://openmx.psyc.virginia.edu/wiki/matrix-operators-and-functions
means    	<- mxMatrix( type="Zero", nrow=1, ncol=nv, name="means" )
beta    	<- mxMatrix( type="Full", nrow = 1, ncol = 2, free = c(T,T), labels=c("b","b"), values=1.7, lbound = -2.5, ubound = 4.5, name="beta" )  
defvarT1 	<- mxMatrix( type="Full", nrow = 1, ncol = 1, free = F, labels = "data.t1_par_warm_ord", name = "defvarT1")
defvarT2 	<- mxMatrix( type="Full", nrow = 1, ncol = 1, free = F, labels = "data.t2_par_warm_ord", name = "defvarT2")
modthrT1	<- mxAlgebra( expression = defvarT1 %x% beta, name = "modthrT1")	
modthrT2	<- mxAlgebra( expression = defvarT2 %x% beta, name = "modthrT2")  	
expM  	 	<- mxAlgebra( expression= cbind(means + modthrT1,means + modthrT2), name="expM" ) 

thresh    	<- mxMatrix( type="Full", nrow=nth, ncol=nv, free=c(T,T,T,T, T,F,F,F),  values=c(-0.7, -0.2, 0.5, 0.8 , 0.2,0,0,0), lbound=-5, ubound=5, labels=c("par_warm1","par_warm2","par_warm3","par_warm4","mdd","mdd2","mdd3","mdd4"), name="thresh" )
Inc      <- mxMatrix( type="Full", nrow=nth, ncol=nth, free=FALSE, values=1, name="Inc" )
expThre 	<- mxAlgebra( expression = cbind(Inc%*%thresh, Inc%*%thresh), name="expThre")

# Declare vVariance components
VarT1	<- mxAlgebra( expression = (a+defvarT1%x%aB)%*%t(a+defvarT1%x%aB) + (c+defvarT1%x%cB)%*%t(c+defvarT1%x%cB) + (e+defvarT1%x%eB)%*%t(e+defvarT1%x%eB), name="VarT1" )
VarT2	<- mxAlgebra( expression = (a+defvarT2%x%aB)%*%t(a+defvarT2%x%aB) + (c+defvarT2%x%cB)%*%t(c+defvarT2%x%cB) + (e+defvarT2%x%eB)%*%t(e+defvarT2%x%eB), name="VarT2" )
covA	<- mxAlgebra( expression = (a+defvarT1%x%aB)%*%t(a+defvarT2%x%aB) , name="covA" )
covC	<- mxAlgebra( expression = (c+defvarT1%x%cB)%*%t(c+defvarT2%x%cB) , name="covC" )

# Constraint on variance of Ordinal variables
unitM  		<- mxMatrix( type="Unit", nrow=4, ncol=1, name="unitM" )
Var1     	<- mxConstraint( expression = rbind(diag2vec(VarT1),diag2vec(VarT2))==unitM, name="Var1" )
covMZ   	<- mxAlgebra( expression = rbind( cbind(VarT1, covA+covC), cbind(covA+covC , VarT2)), name="expCovMZ" )
covDZ    <- mxAlgebra( expression = rbind( cbind(VarT1,0.5%x%covA+covC), cbind(0.5%x%covA+covC, VarT2)), name="expCovDZ" )

# Data objects for Multiple Groups
dataMZ    <- mxData( observed=mzdata, type="raw" )
dataDZ    <- mxData( observed=dzdata, type="raw" )

# Objective objects for Multiple Groups
objMZ    <- mxExpectationNormal( covariance="expCovMZ", means="expM", thresholds="expThre", dimnames=selVars )
objDZ    <- mxExpectationNormal( covariance="expCovDZ", means="expM", thresholds="expThre", dimnames=selVars )
fiML     <- mxFitFunctionML()

# Combine Groups
pars      	<- list( a,c,e,aB,cB,eB, thresh,expThre , Inc )
defs		<- list( means,beta,defvarT1,defvarT2,modthrT1,modthrT2,expM, VarT1,VarT2,covA,covC,unitM,Var1 )
modelMZ   	<- mxModel( pars, defs, dataMZ, covMZ, objMZ, fiML, name="MZ" )
modelDZ   	<- mxModel( pars, defs, dataDZ, covDZ, objDZ, fiML, name="DZ" )
minus2ll  	<- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj       	<- mxFitFunctionAlgebra( "m2LL" )
ModAceModelwarm <- mxModel( "ModACE_Par_Warm", pars, modelMZ, modelDZ, minus2ll, obj )
ModAceFitwarm    <- mxTryHard(ModAceModelwarm,greenOK=FALSE)
ModAceFitwarm    <- mxRun(ModAceModelwarm)
ModAceSummwarm   <- summary(ModAceFitwarm)
ModAceSummwarm
tableFitStatistics(ModAceFitwarm,CholAceFitwarm)

ciVC	<-mxCI(c('b'), interval = 0.95)
ModAceModelwarm2 <- mxModel( "ModACE_Par_warm2", pars, modelMZ, modelDZ, minus2ll, obj, ciVC )

#ModAceFitwarm2    <- mxTryHard(ModAceModelwarm2, unsafe=T, intervals=T)
ModAceFitwarm2    <- mxRun(ModAceModelwarm2, unsafe=T, intervals=T)
ModAceSummwarm2   <- summary(ModAceFitwarm2)
ModAceSummwarm2

