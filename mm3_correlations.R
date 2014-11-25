# 
# 
# Author: Nathan Gillespie
# Library: http://www.vipbg.vcu.edu/tc2012/tc2012_OpenMx.shtml
# Date: 02 20 2010 

# ------------------------------------------------------------------------------						
# File Description
# ------------------------------------------------------------------------------						
# Matrix style model inputs - Raw data input
# Data created here:
# SPS: /Users/ngillespie/Documents/ngillespie/papers/report_k99/data_prep1/MF2_multi_drug_use_v3(new_stem_all_males).sps
# Section: COMMENT Data preparation for Emily Olivares: Merge PARENTAL MONITORING & MF2 drug diagnoses
# Data: /Users/ngillespie/Documents/ngillespie/papers/report_k99/data_prep1/MF2MM3_monitor_dsm4.dat
# 
#
# ------------------------------------------------------------------------------						
# Revision History
# ------------------------------------------------------------------------------						

# Nathan Gillespie 	-- 03 16 2012 

# ------------------------------------------------------------------------------						
# Clear worksace & name folder location
# ------------------------------------------------------------------------------						

# 1. Change and set a new working directory
rm(list = ls(all = TRUE))
setwd("~/Desktop/FYP")
getwd()					# Checks the folder location

# ------------------------------------------------------------------------------						
# Update / load required libraries
# ------------------------------------------------------------------------------						

source('http://openmx.psyc.virginia.edu/getOpenMx.R')
require(OpenMx)
#source("http://www.vipbg.vcu.edu/tc2012/GenEpiHelperFunctions.R") # Update regularly from VIPBG website

# 2. Load in the 'psych' library. THis will be used to run some basic correlations
# Instruction for installing 'psych' library package:
# Step 1. Go to Packages & Data>Package Installer> 
# Step 2. Select CRAN (sources) in drop down menu
# Step 3. Type in psych, hit 'Get List', highlight 'psych' in list & hit Install Selected 
require(psych)


# ------------------------------------------------------------------------------						
# Prepare / Read in Data
# ------------------------------------------------------------------------------	

allvars <- c("FamNo", "ID", "h1", "h51", "h63",  "f1", "f2", "f3", "f4", "zyg", "age", "temp1", "temp2", "temp3", "temp4", "epoch2", "epoch3")
# 3. Change the location of the data file 			  
dataR1 <- read.table("~/Desktop/MM3/Data files/mm3_final.dat",header=F, na.strings="99", col.names=allvars)
summary(dataR1)
dim(dataR1)
dimnames(dataR1)
			
# ------------------------------------------------------------------------------						
# Variable description
# ------------------------------------------------------------------------------	


# ------------------------------------------------------------------------------						
# Select variables for analyses
# ------------------------------------------------------------------------------

selvars <- c("temp1", "epoch2", "epoch3") 

# ------------------------------------------------------------------------------
# Subset, extract data sets, recode & specify as ordinal data
# ------------------------------------------------------------------------------

#dataR1[,c('age_int')] <- dataR1[,c('age_int')]			# Select covariates
data <- cbind(dataR1[,3:5], dataR1[,12:15])					# Re-combines variables & age 
summary(data) 


# ------------------------------------------------------------------------------
# Esimate Spearman & Kendall correlations
# ------------------------------------------------------------------------------

# If use is "complete.obs" then missing values are handled by casewise deletion
# If use is "pairwise.complete.obs" then the correlation or covariance between each pair of variables is computed using 
# all complete pairs of observations on those variables.
# Spearman rho & Kendall's tau-b for rank correlations
# See http://blog.codalism.com/?p=59
# Pearson's r for quantiative variables

# 4. Read up on the differences between Pearson vs Spearman & Kendall and their correlations
# 5. Note the difference between 'complete' and 'pairwise deletion'
# 6. Run "pearson" correlations (complete and pairwise) and then do the same for "kendall" and "pearson

cor(data[,c(1:6)], method="pearson", use="complete.obs")
cor(data[,c(1:6)], method="pearson", use="pairwise.complete.obs")




# ------------------------------------------------------------------------------
# Convert / declare ordinal data using mx Factor
# ------------------------------------------------------------------------------

# 'Order' STEM variables in columns 1, 2, 3 in data with 2 levels
 Part1 <- mxFactor( x=data[,(3:5)], levels=c(0:1))
# 'Order' DIAGNOSES in columns 4,5, 6, 7 in data with 3 levels		
 Part2 <- mxFactor( x=data[,(6:9)], levels=c(0:2),ordered=T)	
head(data)
DataOrdF <- cbind(Part1, Part2, data[,10:11])

summary(DataOrdF)

#OpenMx requires ordinal data to be ordered. R's factor function doesn't enforce this. Relying on the data to specify the data is foolhardy for the following reasons: The factor function will skip levels missing from the data: Specifying these in levels leaves the list of levels complete. Data will often not explore the min and max level that the user knows are possible. For these reasons this function forces you to write out all possible levels explicitly.

# ------------------------------------------------------------------------------
# Print Descriptive Statistics
# ------------------------------------------------------------------------------

summary(DataOrdF)
dim(DataOrdF)									# Size of data set

# ------------------------------------------------------------------------------						
# Update / load 'polycor' library
# ------------------------------------------------------------------------------						
# Install library / packages instructions:
# Step 1. Go to Packages & Data>Package Installer> 
# Step 2. Select CRAN (sources) in drop down menu
# Step 3. Type in polycor, hit 'Get List', highlight 'polycor' in list & hit Install Selected 

# ------------------------------------------------------------------------------
# Polychoric correlations
# ------------------------------------------------------------------------------

# polychoric correlation
library(polycor)				# THis loads up the 'polycor' library
hetcor(DataOrdF[,c(1:6)]) 		# TASK: Run correlations between MONITOR (1) & STEM items (2:4) 


# ------------------------------------------------------------------------------
# Run Pairwise Polychoric wrapper
# ------------------------------------------------------------------------------

result <- polychoricMatrix(DataOrdF, useDeviations=TRUE)
summary(result)
warnings()
result

result <- polypairwise(DataOrdF, useDeviations=TRUE)
summary(result)
warnings()
result
polycor_result<-result

save.image("pgd.RData")

#Write out the R object
write.table(result$SE,"~/Desktop/SEmatrix.txt")
write.table(result$R,"~/Desktop/correlations1.txt")
commas <- read.table("~/Desktop/correlations1.txt")
SEmatrix <- read.table("~/Desktop/SEmatrix.txt")

#Read in the R object
corMatrix 			<- as.matrix(result$R)
seMatrix 			<- as.matrix(result$SE)
corMatrix 			<- as.matrix(result$R)
corMatrix 			<- as.matrix(commas)
dimnames(corMatrix)


# ------------------------------------------------------------------------------						
# Prepare / Read in Data
# ------------------------------------------------------------------------------	

# 1. Change and set a new working directory
rm(list = ls(all = TRUE))
setwd("~/Desktop/mm3")
source('http://openmx.psyc.virginia.edu/getOpenMx.R')
2
require(OpenMx)

source("http://www.vipbg.vcu.edu/tc2012/GenEpiHelperFunctions.R") # Update regularly from VIPBG website
#source("GenEpiHelperFunctions.R") # Update regularly from VIPBG website


allvars <- c( "FamNo", "ID", "h1", "h51", "h63",  "f1", "f2", "f3", "f4", "zyg", "temp1", "temp2", "temp3", "temp4")
# 3. Change the location of the data file 			  
dataR1 <- read.table("~/Desktop/mm3_substance_b.dat",header=F, na.strings="99", col.names=allvars)
summary(dataR1)

# 'Order' STEM variables in columns 1, 2, 3 4 in DataOrd with 4 levels
 Part1 <- mxFactor( x=dataR1[,c(3:5)], levels=c(0:1))
# 'Order' DIAGNOSES in columns 5, 6, 7 in DataOrd with 2 levels		
 Part2 <- mxFactor( x=dataR1[,c(11:14)], levels=c(0:2),ordered=T)	
#DataOrdF <- cbind(Part1,Part2, DataOrd[,8])					# Re-combines variables & age into 1 file
DataOrdF <- cbind(Part1, Part2, dataR1[,1:2], dataR1[,6:10])			
summary(DataOrdF)


# ------------------------------------------------------------------------------
# Convert to twin format: 2 cases per line
# ------------------------------------------------------------------------------


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

dataR12 <- twindat(dat=DataOrdF, famid= "FamNo", twinid= "ID", zygosity= "zyg")
head(dataR12)
dimnames(dataR12)
summary(dataR12)
varnames <- colnames(dataR12)


# ------------------------------------------------------------------------------						
# Univariate models with thresholds
# ------------------------------------------------------------------------------						


# Program: UnivSat&ACE_Bin.R  
# Fruhling Rijsdijk: http://ibg.colorado.edu/cdrom2012/rijsdijk/
# http://ibg.colorado.edu/cdrom2012/rijsdijk/ThresholdLiabilityModels/Practical/UnivSat&ACE_Bin.R

require(psych)
require(OpenMx)
# ------------------------------------------------------------------------------
# Subset data and split into zygosity groups
# ------------------------------------------------------------------------------
summary(dataR12)
str(dataR12)

selVars <- c("temp3_T1" , "h1_T1", "h51_T1", "h63_T1", "temp4_T1","temp3_T2" , "h1_T2", "h51_T2", "h63_T2", "temp4_T2")
mzData <- subset(dataR12, zyg==3, selVars)
dzData <- subset(dataR12, zyg==4, selVars)
head(dzData)


summary(mzData)

# ------------------------------------------------------------------------------
# Specify starting values & create labels used in the script
# ------------------------------------------------------------------------------

vars	<- 'PGD'
vars2	<- 'Drug_Initiation'

nv 		<- 3		# number of variables per twin
ntv 	<- nv*2		# number of variables per pair

nth 	<- 2		# Number of thresholds
lth     <- -1.5    	# start values for first threshold
ith     <-  0.7    	# start values for threshold increments

th	    <- matrix(rep( c(lth,(rep(ith,nth-1))) ),nrow=nth,ncol=nv)			# Start values. Increments foreced to be +ve
thVal	<- cbind(th,th)
thLB      <- matrix(rep( c(-4,(rep(0.001,nth-1)) ),nv),nrow=nth,ncol=nv)	# Lower bounds for increments are +ve
thUB      <-  4.0															# Upper bounds threshold parameters	

corVals   <-  0.5			# Start value for correlations
lbrVal    <- -0.99   		# Lower bounds for a, c & e path coefficients
ubrVal    <-  0.99    		# Upper bounds for a, c & e path coefficients
thLab     <- paste("t",1:nth,"Z",sep="")
paVal     <- 0.65

# Fancy way of labelling threshold labels
THlabs   <- c(paste("MZ1_t1",1:nv,sep=""),paste("MZ1_t2",1:nv,sep=""),paste("MZ1_t3",1:nv,sep=""))

mz1th     <- paste("MZ1","_th",1:nth, sep="")
mz2th     <- paste("MZ2","_th",1:nth, sep="")
dz1th     <- paste("DZ1","_th",1:nth, sep="")
dz2th     <- paste("DZ2","_th",1:nth, sep="")


# Starting values for correaltions
corValsM  <-.8    	    # start value for MZ correlations
corValsD  <-.4     	# start value for DZ correlations


# ---------------------------------------------------------------------
# Print Descriptive Statistics
# ---------------------------------------------------------------------

summary(mzData)
summary(dzData)
table(mzData$temp3_T2, mzData$temp3_T2 )
table(dzData$temp3_T2, dzData$temp3_T2)

# ___________________________________
# PREPARE SATURATED MODEL
# ___________________________________

# Matrices for expected Means & Thresholds (on liabilities) in MZ & DZ twin pairs
meanG	<-mxMatrix( type="Zero", nrow=1, ncol=ntv, name="expMean" )

threMZ	<-mxMatrix(type="Full", nrow=nth, ncol=ntv, free=T, values=thVal, lbound=thLB, ubound=thUB, labels=cbind(mz1th,mz2th), name="threMZ" )
threDZ	<-mxMatrix(type="Full", nrow=nth, ncol=ntv, free=T, values=thVal, lbound=thLB, ubound=thUB, labels=cbind(dz1th,dz2th), name="threDZ" )
corMZ	<-mxMatrix(type="Stand", nrow=ntv, ncol=ntv, free=T, values=corValsM, lbound=-.99, ubound=.99, labels=c("rMZ"), name="corMZ") 
corDZ	<-mxMatrix(type="Stand", nrow=ntv, ncol=ntv, free=T, values=corValsD, lbound=-.99, ubound=.99, labels=c("rDZ"), name="corDZ") 

Inc      	<- mxMatrix( type="Lower", nrow=nth, ncol=nth, free=F, values=1, name="Inc" )
expThreMZ   <- mxAlgebra( expression= cbind( Inc %*% threMZ ), name="expThreMZ" )
expThreDZ   <- mxAlgebra( expression= cbind( Inc %*% threDZ ), name="expThreDZ" )

# Data objects for Multiple Groups
dataMZ	<-mxData(mzData, type="raw")
dataDZ	<-mxData(dzData, type="raw")

# Objective objects for Multiple Groups
objMZ	<-mxFIMLObjective( covariance="corMZ", means="expMean", dimnames=selVars, thresholds="expThreMZ" )
objDZ	<-mxFIMLObjective( covariance="corDZ", means="expMean", dimnames=selVars, thresholds="expThreDZ" )
 
# Combine Groups	
groupMZ		<-mxModel("MZ", meanG, threMZ, corMZ, Inc, expThreMZ, dataMZ, objMZ )
groupDZ		<-mxModel("DZ", meanG, threDZ, corDZ, Inc, expThreDZ, dataDZ, objDZ )
minus2ll	<-mxAlgebra( MZ.objective + DZ.objective, name="minus2loglikelihood" )
obj			<-mxAlgebraObjective("minus2loglikelihood") 
ciCor		<-mxCI(c('MZ.corMZ[2,1]', 'DZ.corDZ[2,1]'), interval=0.95 )
ciThre		<-mxCI( c('MZ.expThreMZ','DZ.expThreDZ'), interval=0.95 )
twinSatModel   <-mxModel( "twinSat", minus2ll, obj, groupMZ, groupDZ, ciCor, ciThre ) 
      
      
# -----------------------------------------------------------------------
#  RUN SATURATED MODEL (Tetrachoric correlations) 
# -----------------------------------------------------------------------

twinSatFit	<- mxRun(twinSatModel, intervals=T)
twinSatSumm 	<- summary(twinSatFit)
round(twinSatFit@output$estimate,4)
twinSatSumm



# -----------------------------------------------------------------------
# Generate Saturated Model Output
#--------------------------------------------------------------------------

rMZ<- twinSatFit$MZ.corMZ@values[2,1]
rDZ<- twinSatFit$DZ.corDZ@values[2,1]
rMZ 
rDZ

tMZ<- twinSatFit$MZ.expThreMZ@result
tDZ<- twinSatFit$DZ.expThreDZ@result
tMZ
tDZ

twinSatFit@output$confidenceIntervals

twinSatLLL     <-twinSatFit@output$Minus2LogLikelihood
twinSatAIC     <-twinSatSumm$AIC
twinSatOS      <-twinSatSumm$observedStatistics
twinSatDF      <-twinSatSumm$degreesOfFreedom
twinSatNP      <-length(twinSatSumm$parameters[[1]])

#twinSatMatrices <- c("expThreMZ", "expThreDZ")
#twinSatLabels <- c("MZthresholds","DZthresholds")
#formatOutputMatrices(twinSatFit,twinSatMatrices,twinSatLabels,selVars,4)

# Use Helper Functions
source("GenEpiHelperFunctions.R")
#source("http://www.vipbg.vcu.edu/tc2012/GenEpiHelperFunctions.R") # Update regularly from VIPBG website
expectedMeansCovariances(twinSatFit)
tableFitStatistics(twinSatFit)


# -------------------------------------------------------------------------
# RUN SUBMODELS - Equating thresholds
# -------------------------------------------------------------------------

# SubModel 1: constraining Thresholds across Twin 1 and Twin 2 within zyg group to be equal

eqThresholdsTwinModel    <- twinSatFit
eqThresholdsTwinModel    <- omxSetParameters( eqThresholdsTwinModel, label=cbind(mz1th,mz2th), free=TRUE, values=thVal, newlabels=cbind(mz1th,mz1th) )
eqThresholdsTwinModel    <- omxSetParameters( eqThresholdsTwinModel, label=cbind(dz1th,dz2th), free=TRUE, values=thVal, newlabels=cbind(dz1th,dz1th) )

eqThresholdsTwinModel@submodels$MZ@matrices$threMZ@labels
eqThresholdsTwinModel@submodels$DZ@matrices$threDZ@labels

eqThresholdsTwinFit     <- mxRun( eqThresholdsTwinModel, intervals=T )
eqThresholdsTwinSumm	<- summary( eqThresholdsTwinFit )
eqThresholdsTwinLLL	<- eqThresholdsTwinFit@output$Minus2LogLikelihood
tableFitStatistics(twinSatFit, eqThresholdsTwinFit)
eqThresholdsTwinSumm

expectedMeansCovariances(eqThresholdsTwinFit)
tableFitStatistics(eqThresholdsTwinFit)


# SubModel 2: constraining Thresholds across Twin 1 and Twin 2 and across zyg group to be equal

eqThresholdsZygModel     <- eqThresholdsTwinModel
eqThresholdsZygModel     <- omxSetParameters( eqThresholdsZygModel, label=cbind(dz1th,dz1th), free=TRUE, values=thVal, newlabels=cbind(mz1th,mz1th) )

eqThresholdsZygModel@submodels$MZ@matrices$threMZ@labels
eqThresholdsZygModel@submodels$DZ@matrices$threDZ@labels

eqThresholdsZygFit      <- mxRun( eqThresholdsZygModel, intervals=F )
eqThresholdsZygSumm	<- summary( eqThresholdsZygFit )
eqThresholdsZygLLL	<- eqThresholdsZygFit@output$Minus2LogLikelihood
eqThresholdsZygSumm
eqThresholdsZygLLL
submodels <- c(eqThresholdsTwinFit,eqThresholdsZygFit)
tableFitStatistics(twinSatFit, submodels)


# Print Comparative Fit Statistics for Saturated Models
# -----------------------------------------------------------------------------

submodels <- c(eqThresholdsTwinFit,eqThresholdsZygFit)
tableFitStatistics(twinSatFit, submodels)


# -------------------------------------------------------------------------
# RUN SUBMODELS - Equating correlations
# -------------------------------------------------------------------------

# SubModel 3: constraining MZ and DZ correlations to be equal

eqCorrelationsModel     <- eqThresholdsZygFit
eqCorrelationsModel     <- omxSetParameters( eqCorrelationsModel, labels="rDZ", values=corValsM, free=TRUE, newlabels="rMZ" )
eqCorrelationsModel     <- omxSetParameters( eqCorrelationsModel, labels="rMZ", values=corValsM)

eqCorrelationsModel@submodels$MZ@matrices$corMZ@labels
eqCorrelationsModel@submodels$DZ@matrices$corDZ@labels

eqCorrelationsFit      	<- mxRun( eqCorrelationsModel, intervals=T )
eqCorrelationsFitSumm	<- summary( eqCorrelationsFit )
eqCorrelationsFitLL		<- eqCorrelationsFit@output$Minus2LogLikelihood
tableFitStatistics(eqThresholdsZygFit, eqCorrelationsFit)
eqCorrelationsFitSumm

# SubModel 4: constraining MZ and DZ correlations to equal ZERO

eq0CorrelationsModel     <- eqThresholdsZygFit
eq0CorrelationsModel     <- omxSetParameters( eq0CorrelationsModel, labels="rMZ", values=0, free=F)
eq0CorrelationsModel     <- omxSetParameters( eq0CorrelationsModel, labels="rDZ", values=0, free=F)

eq0CorrelationsModel@submodels$MZ@matrices$corMZ@values
eq0CorrelationsModel@submodels$DZ@matrices$corDZ@values

eq0CorrelationsFit      	<- mxRun( eq0CorrelationsModel, intervals=T )
eq0CorrelationsFitSumm	<- summary( eq0CorrelationsFit )
eq0CorrelationsFitLL		<- eq0CorrelationsFit@output$Minus2LogLikelihood
tableFitStatistics(eqThresholdsZygFit, eq0CorrelationsFit)
eq0CorrelationsFitSumm


# Print Comparative Fit Statistics for Saturated Models
# -----------------------------------------------------------------------------
SatNested <- list(eqCorrelationsFit, eq0CorrelationsFit)
tableFitStatistics(eqThresholdsZygFit, SatNested)


# ___________________________________
# PREPARE GENETIC MODEL
# ___________________________________


# ACE Model with one overall Threshold
# Matrices to store a, c, and e Path Coefficients
pathA    <-mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=.6, label="a11", name="a" ) 
pathC    <-mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=.6, label="c11", name="c" )
pathE    <-mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=.6, label="e11", name="e" )
	
# Algebra to generate Matrices to hold A, C, and E computed Variance Components
covA     <-mxAlgebra( expression=a %*% t(a), name="A" )
covC     <-mxAlgebra( expression=c %*% t(c), name="C" ) 
covE     <-mxAlgebra( expression=e %*% t(e), name="E" )

# Algebra to compute Total Variance
covP     <-mxAlgebra( expression=A+C+E, name="V" )

# Matrices for expected Means & Thresholds (on liabilities) 
meanG    <-mxMatrix( type="Zero", nrow=1, ncol=ntv, name="expMean" )
Inc      <- mxMatrix( type="Lower", nrow=nth, ncol=nth, free=F, values=1, name="Inc" )
threMZ	 <-mxMatrix(type="Full", nrow=nth, ncol=ntv, free=T, values=thVal, lbound=thLB, ubound=thUB, labels=cbind(mz1th,mz1th), name="threMZ" )
expThreMZ   <- mxAlgebra( expression= cbind( Inc %*% threMZ ), name="expThreMZ" )

# Algebra to generate Matrices to hold Parameter Estimates and Derived Variance Components
rowVars  <-rep('vars',nv)
colVars  <-rep(c('A','C','E','SA','SC','SE'),each=nv)
estVars  <-mxAlgebra( expression=cbind(A,C,E,A/V,C/V,E/V), name="Vars", dimnames=list(rowVars,colVars))

# Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covMZ    <-mxAlgebra( expression= rbind( cbind(A+C+E , A+C),
                                          cbind(A+C   , A+C+E)), name="expCovMZ" )
covDZ    <-mxAlgebra( expression= rbind( cbind(A+C+E     , 0.5%x%A+C),
                                          cbind(0.5%x%A+C , A+C+E)), name="expCovDZ" )

# Constraint on variance of the liability of Binary variables (assumed to have a SND) 
matUnv	<-mxMatrix( type="Unit", nrow=nv, ncol=1, name="Unv1" )
var1	<-mxConstraint( expression=diag2vec(V)==Unv1, name="Var1" )

# Data objects for Multiple Groups
dataMZ   <-mxData( observed=mzData, type="raw" )
dataDZ   <-mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ	<-mxFIMLObjective( covariance="expCovMZ", means="expMean", dimnames=selVars, thresholds="expThreMZ" )
objDZ	<-mxFIMLObjective( covariance="expCovDZ", means="expMean", dimnames=selVars, thresholds="expThreMZ" )

# Combine Groups
pars	<-list( pathA, pathC, pathE, covA, covC, covE, covP, meanG, Inc, threMZ, expThreMZ, matUnv, estVars )
modelMZ	<-mxModel( pars, covMZ, dataMZ, objMZ, name="MZ" )
modelDZ	<-mxModel( pars, covDZ, dataDZ, objDZ, name="DZ" )
minus2ll<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj	<-mxAlgebraObjective( "m2LL" )
ciVC	<-mxCI(c('A', 'C', 'E'), interval=0.95 )
AceModel <-mxModel( "ACE", pars, var1, modelMZ, modelDZ, minus2ll, obj, ciVC)


# -------------------------------------------------
# RUN GENETIC MODEL
# -------------------------------------------------

# Run ACE Model
AceFit  <- mxRun(AceModel, intervals=T )
AceSumm <- summary(AceFit)
AceSumm
round(AceFit@output$estimate,2)
round(AceFit$Vars@result,2)

AceFit@matrices$threMZ@values
AceFit@algebras$expThreMZ

# Generate ACE Model Output (using helper function GenEpiHelperFunctions.R)
parameterSpecifications(AceFit)
expectedMeansCovariances(AceFit)
tableFitStatistics(AceFit)

ACEcovMatrices	<- c("A","C","E","V","A/V","C/V","E/V")
ACEcovLabels	<- c("covs_A","covs_C","covs_E","Var","prop_A","prop_C","prop_E")
formatOutputMatrices(AceFit,ACEcovMatrices,ACEcovLabels,vars,4)

# -----------------------------------------------------------
# RUN SUBMODELS
# -----------------------------------------------------------

# Fit AE model
# -----------------------------------------------------------------
AeModel	<- mxModel( AceFit, name="AE")
AeModel	<- omxSetParameters( AeModel, labels="c11", free=FALSE, values=0 )
AeFit	<- mxRun(AeModel, intervals=T)
AeSumm  <- summary(AeFit)
AeSumm
round(AeFit@output$estimate,4)
round(AeFit$Vars@result,4)


# Run CE model
# ------------------------------------------------------------------
CeModel	<- mxModel( AceFit, name="CE")
CeModel	<- omxSetParameters( CeModel, labels="a11", free=FALSE, values=0 )
CeFit	<- mxRun(CeModel)
round(CeFit@output$estimate,4)
round(CeFit$Vars@result,4)

# ------------------------------------------------------------------------------

# Run E model
# ------------------------------------------------------------------
eModel	<- mxModel( AceFit, name="E")
eModel	<- omxSetParameters( eModel, labels="a11", free=FALSE, values=0 )
eModel	<- omxSetParameters( eModel, labels="c11", free=FALSE, values=0 )
eFit	<- mxRun(eModel)
round(eFit@output$estimate,4)
round(eFit$Vars@result,4)

# ------------------------------------------------------------------------------

# Print Comparative Fit Statistics
AceNested	<- list(AeFit, CeFit, eFit)
tableFitStatistics(AceFit,AceNested)

#round(rbind(AceFit@output$estimate,AeFit@output$estimate,CeFit@output$estimate,eFit@output$estimate),4)
round(rbind(AceFit$Vars@result,AeFit$Vars@result,CeFit$Vars@result,eFit$Vars@result),4)












# ------------------------------------------------------------------------------
# Subset data and split into zygosity groups: Cannabis stem item
# ------------------------------------------------------------------------------
summary(dataR12)
str(dataR12)

selVars2 <- c("stem_ca_T1" , "stem_ca_T2")
mzData2 <- subset(dataR12, zygosity==3, selVars2)
dzData2 <- subset(dataR12, zygosity==4, selVars2)
head(dzData2)


# ------------------------------------------------------------------------------
# Specify starting values & create labels used in the script
# ------------------------------------------------------------------------------

vars2	<- 'cannabis'
nv 		<- 1		# number of variables per twin
ntv 	<- nv*2		# number of variables per pair

nth 	<- 3		# Number of thresholds
lth     <- -1.5    	# start values for first threshold
ith     <-  0.7    	# start values for threshold increments

th	    <- matrix(rep( c(lth,(rep(ith,nth-1))) ),nrow=nth,ncol=nv)			# Start values. Increments foreced to be +ve
thVal	<- cbind(th,th)
thLB      <- matrix(rep( c(-4,(rep(0.001,nth-1)) ),nv),nrow=nth,ncol=nv)	# Lower bounds for increments are +ve
thUB      <-  4.0															# Upper bounds threshold parameters	

corVals   <-  0.5			# Start value for correlations
lbrVal    <- -0.99   		# Lower bounds for a, c & e path coefficients
ubrVal    <-  0.99    		# Upper bounds for a, c & e path coefficients
thLab     <- paste("t",1:nth,"Z",sep="")
paVal     <- 0.65

# Fancy way of labelling threshold labels
THlabs   <- c(paste("MZ1_t1",1:nv,sep=""),paste("MZ1_t2",1:nv,sep=""),paste("MZ1_t3",1:nv,sep=""))

mz1th     <- paste("MZ1","_th",1:nth, sep="")
mz2th     <- paste("MZ2","_th",1:nth, sep="")
dz1th     <- paste("DZ1","_th",1:nth, sep="")
dz2th     <- paste("DZ2","_th",1:nth, sep="")


# Starting values for correaltions
corValsM  <-.8    	    # start value for MZ correlations
corValsD  <-.4     	# start value for DZ correlations


# ---------------------------------------------------------------------
# Print Descriptive Statistics
# ---------------------------------------------------------------------

summary(mzData2)
summary(dzData2)
table(mzData2$stem_ca_T1, mzData2$stem_ca_T2)
table(dzData2$stem_ca_T1, dzData2$stem_ca_T2)

# ___________________________________
# PREPARE SATURATED MODEL for thresholds
# ___________________________________

# Matrices for expected Means & Thresholds (on liabilities) in MZ & DZ twin pairs
meanG	<-mxMatrix( type="Zero", nrow=1, ncol=ntv, name="expMean" )

threMZ	<-mxMatrix(type="Full", nrow=nth, ncol=ntv, free=T, values=thVal, lbound=thLB, ubound=thUB, labels=cbind(mz1th,mz2th), name="threMZ" )
threDZ	<-mxMatrix(type="Full", nrow=nth, ncol=ntv, free=T, values=thVal, lbound=thLB, ubound=thUB, labels=cbind(dz1th,dz2th), name="threDZ" )
corMZ	<-mxMatrix(type="Stand", nrow=ntv, ncol=ntv, free=T, values=corValsM, lbound=-.99, ubound=.99, labels=c("rMZ"), name="corMZ") 
corDZ	<-mxMatrix(type="Stand", nrow=ntv, ncol=ntv, free=T, values=corValsD, lbound=-.99, ubound=.99, labels=c("rDZ"), name="corDZ") 

Inc      	<- mxMatrix( type="Lower", nrow=nth, ncol=nth, free=F, values=1, name="Inc" )
expThreMZ   <- mxAlgebra( expression= cbind( Inc %*% threMZ ), name="expThreMZ" )
expThreDZ   <- mxAlgebra( expression= cbind( Inc %*% threDZ ), name="expThreDZ" )

# Data objects for Multiple Groups
dataMZ2	<-mxData(mzData2, type="raw")
dataDZ2	<-mxData(dzData2, type="raw")

# Objective objects for Multiple Groups
objMZ	<-mxFIMLObjective( covariance="corMZ", means="expMean", dimnames=selVars2, thresholds="expThreMZ" )
objDZ	<-mxFIMLObjective( covariance="corDZ", means="expMean", dimnames=selVars2, thresholds="expThreDZ" )
 
# Combine Groups	
groupMZ		<-mxModel("MZ", meanG, threMZ, corMZ, Inc, expThreMZ, dataMZ2, objMZ )
groupDZ		<-mxModel("DZ", meanG, threDZ, corDZ, Inc, expThreDZ, dataDZ2, objDZ )
minus2ll	<-mxAlgebra( MZ.objective + DZ.objective, name="minus2loglikelihood" )
obj			<-mxAlgebraObjective("minus2loglikelihood") 
ciCor		<-mxCI(c('MZ.corMZ[2,1]', 'DZ.corDZ[2,1]'), interval=0.95 )
ciThre		<-mxCI( c('MZ.expThreMZ','DZ.expThreDZ'), interval=0.95 )
Sat_CA_Stem_Model   <-mxModel( "Sat_CA_Stem", minus2ll, obj, groupMZ, groupDZ, ciCor, ciThre ) 
      
      
# -----------------------------------------------------------------------
#  RUN SATURATED MODEL (Tetrachoric correlations) 
# -----------------------------------------------------------------------

Sat_CA_Stem_ModelFit	<- mxRun(Sat_CA_Stem_Model, intervals=T)
Sat_CA_Stem_ModelSumm 	<- summary(Sat_CA_Stem_ModelFit)
round(Sat_CA_Stem_ModelFit@output$estimate,4)
Sat_CA_Stem_ModelSumm


# -----------------------------------------------------------------------
# Generate Saturated Model Output
#--------------------------------------------------------------------------

rMZ<- twinSatFit$MZ.corMZ@values[2,1]
rDZ<- twinSatFit$DZ.corDZ@values[2,1]
rMZ 
rDZ

tMZ<- twinSatFit$MZ.expThreMZ@result
tDZ<- twinSatFit$DZ.expThreDZ@result
tMZ
tDZ

twinSatFit@output$confidenceIntervals

twinSatLLL     <-twinSatFit@output$Minus2LogLikelihood
twinSatAIC     <-twinSatSumm$AIC
twinSatOS      <-twinSatSumm$observedStatistics
twinSatDF      <-twinSatSumm$degreesOfFreedom
twinSatNP      <-length(twinSatSumm$parameters[[1]])

# Use Helper Functions
source("GenEpiHelperFunctions.R")
#source("http://www.vipbg.vcu.edu/tc2012/GenEpiHelperFunctions.R") # Update regularly from VIPBG website
expectedMeansCovariances(twinSatFit)
tableFitStatistics(twinSatFit)


# -------------------------------------------------------------------------
# RUN SUBMODELS - Equating thresholds
# -------------------------------------------------------------------------

# SubModel 1: constraining Thresholds across Twin 1 and Twin 2 within zyg group to be equal

eqThresholdsTwinModel    <- twinSatFit
eqThresholdsTwinModel    <- omxSetParameters( eqThresholdsTwinModel, label=cbind(mz1th,mz2th), free=TRUE, values=thVal, newlabels=cbind(mz1th,mz1th) )
eqThresholdsTwinModel    <- omxSetParameters( eqThresholdsTwinModel, label=cbind(dz1th,dz2th), free=TRUE, values=thVal, newlabels=cbind(dz1th,dz1th) )

eqThresholdsTwinModel@submodels$MZ@matrices$threMZ@labels
eqThresholdsTwinModel@submodels$DZ@matrices$threDZ@labels

eqThresholdsTwinFit     <- mxRun( eqThresholdsTwinModel, intervals=T )
eqThresholdsTwinSumm	<- summary( eqThresholdsTwinFit )
eqThresholdsTwinLLL	<- eqThresholdsTwinFit@output$Minus2LogLikelihood
tableFitStatistics(twinSatFit, eqThresholdsTwinFit)
eqThresholdsTwinSumm

# SubModel 2: constraining Thresholds across Twin 1 and Twin 2 and across zyg group to be equal

eqThresholdsZygModel     <- eqThresholdsTwinModel
eqThresholdsZygModel     <- omxSetParameters( eqThresholdsZygModel, label=cbind(dz1th,dz1th), free=TRUE, values=thVal, newlabels=cbind(mz1th,mz1th) )

eqThresholdsZygModel@submodels$MZ@matrices$threMZ@labels
eqThresholdsZygModel@submodels$DZ@matrices$threDZ@labels

eqThresholdsZygFit      <- mxRun( eqThresholdsZygModel, intervals=F )
eqThresholdsZygSumm	<- summary( eqThresholdsZygFit )
eqThresholdsZygLLL	<- eqThresholdsZygFit@output$Minus2LogLikelihood
eqThresholdsZygSumm
eqThresholdsZygLLL
submodels <- c(eqThresholdsTwinFit,eqThresholdsZygFit)
tableFitStatistics(twinSatFit, submodels)


# Print Comparative Fit Statistics for Saturated Models
# -----------------------------------------------------------------------------

submodels <- c(eqThresholdsTwinFit,eqThresholdsZygFit)
tableFitStatistics(twinSatFit, submodels)


# -------------------------------------------------------------------------
# RUN SUBMODELS - Equating correlations
# -------------------------------------------------------------------------

# SubModel 3: constraining MZ and DZ correlations to be equal

eqCorrelationsModel     <- eqThresholdsZygFit
eqCorrelationsModel     <- omxSetParameters( eqCorrelationsModel, labels="rDZ", values=corValsM, free=TRUE, newlabels="rMZ" )
eqCorrelationsModel     <- omxSetParameters( eqCorrelationsModel, labels="rMZ", values=corValsM)

eqCorrelationsModel@submodels$MZ@matrices$corMZ@labels
eqCorrelationsModel@submodels$DZ@matrices$corDZ@labels

eqCorrelationsFit      	<- mxRun( eqCorrelationsModel, intervals=T )
eqCorrelationsFitSumm	<- summary( eqCorrelationsFit )
eqCorrelationsFitLL		<- eqCorrelationsFit@output$Minus2LogLikelihood
tableFitStatistics(eqThresholdsZygFit, eqCorrelationsFit)
eqCorrelationsFitSumm

# SubModel 4: constraining MZ and DZ correlations to equal ZERO

eq0CorrelationsModel     <- eqThresholdsZygFit
eq0CorrelationsModel     <- omxSetParameters( eq0CorrelationsModel, labels="rMZ", values=0, free=F)
eq0CorrelationsModel     <- omxSetParameters( eq0CorrelationsModel, labels="rDZ", values=0, free=F)

eq0CorrelationsModel@submodels$MZ@matrices$corMZ@values
eq0CorrelationsModel@submodels$DZ@matrices$corDZ@values

eq0CorrelationsFit      	<- mxRun( eq0CorrelationsModel, intervals=T )
eq0CorrelationsFitSumm	<- summary( eq0CorrelationsFit )
eq0CorrelationsFitLL		<- eq0CorrelationsFit@output$Minus2LogLikelihood
tableFitStatistics(eqThresholdsZygFit, eq0CorrelationsFit)
eq0CorrelationsFitSumm


# Print Comparative Fit Statistics for Saturated Models
# -----------------------------------------------------------------------------
SatNested <- list(eqCorrelationsFit, eq0CorrelationsFit)
tableFitStatistics(eqThresholdsZygFit, SatNested)


# ___________________________________
# PREPARE GENETIC MODEL
# ___________________________________


# ACE Model with one overall Threshold
# Matrices to store a, c, and e Path Coefficients
pathA    <-mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=.6, label="a11", name="a" ) 
pathC    <-mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=.6, label="c11", name="c" )
pathE    <-mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=.6, label="e11", name="e" )
	
# Algebra to generate Matrices to hold A, C, and E computed Variance Components
covA     <-mxAlgebra( expression=a %*% t(a), name="A" )
covC     <-mxAlgebra( expression=c %*% t(c), name="C" ) 
covE     <-mxAlgebra( expression=e %*% t(e), name="E" )

# Algebra to compute Total Variance
covP     <-mxAlgebra( expression=A+C+E, name="V" )

# Matrices for expected Means & Thresholds (on liabilities) 
meanG    <-mxMatrix( type="Zero", nrow=1, ncol=ntv, name="expMean" )
Inc      <- mxMatrix( type="Lower", nrow=nth, ncol=nth, free=F, values=1, name="Inc" )
threMZ	 <-mxMatrix(type="Full", nrow=nth, ncol=ntv, free=T, values=thVal, lbound=thLB, ubound=thUB, labels=cbind(mz1th,mz1th), name="threMZ" )
expThreMZ   <- mxAlgebra( expression= cbind( Inc %*% threMZ ), name="expThreMZ" )

# Algebra to generate Matrices to hold Parameter Estimates and Derived Variance Components
rowVars  <-rep('vars',nv)
colVars  <-rep(c('A','C','E','SA','SC','SE'),each=nv)
estVars  <-mxAlgebra( expression=cbind(A,C,E,A/V,C/V,E/V), name="Vars", dimnames=list(rowVars,colVars))

# Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covMZ    <-mxAlgebra( expression= rbind( cbind(A+C+E , A+C),
                                          cbind(A+C   , A+C+E)), name="expCovMZ" )
covDZ    <-mxAlgebra( expression= rbind( cbind(A+C+E     , 0.5%x%A+C),
                                          cbind(0.5%x%A+C , A+C+E)), name="expCovDZ" )

# Constraint on variance of the liability of Binary variables (assumed to have a SND) 
matUnv	<-mxMatrix( type="Unit", nrow=nv, ncol=1, name="Unv1" )
var1	<-mxConstraint( expression=diag2vec(V)==Unv1, name="Var1" )

# Data objects for Multiple Groups
dataMZ   <-mxData( observed=mzData, type="raw" )
dataDZ   <-mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ	<-mxFIMLObjective( covariance="expCovMZ", means="expMean", dimnames=selVars, thresholds="expThreMZ" )
objDZ	<-mxFIMLObjective( covariance="expCovDZ", means="expMean", dimnames=selVars, thresholds="expThreMZ" )

# Combine Groups
pars	<-list( pathA, pathC, pathE, covA, covC, covE, covP, meanG, Inc, threMZ, expThreMZ, matUnv, estVars )
modelMZ	<-mxModel( pars, covMZ, dataMZ, objMZ, name="MZ" )
modelDZ	<-mxModel( pars, covDZ, dataDZ, objDZ, name="DZ" )
minus2ll<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj	<-mxAlgebraObjective( "m2LL" )
ciVC	<-mxCI(c('A', 'C', 'E'), interval=0.95 )
AceModel <-mxModel( "ACE", pars, var1, modelMZ, modelDZ, minus2ll, obj, ciVC)


# -------------------------------------------------
# RUN GENETIC MODEL
# -------------------------------------------------

# Run ACE Model
AceFit  <- mxRun(AceModel, intervals=T )
AceSumm <- summary(AceFit)
AceSumm
round(AceFit@output$estimate,2)
round(AceFit$Vars@result,2)

AceFit@matrices$threMZ@values
AceFit@algebras$expThreMZ

# Generate ACE Model Output (using helper function GenEpiHelperFunctions.R)
parameterSpecifications(AceFit)
expectedMeansCovariances(AceFit)
tableFitStatistics(AceFit)

ACEcovMatrices	<- c("A","C","E","V","A/V","C/V","E/V")
ACEcovLabels	<- c("covs_A","covs_C","covs_E","Var","prop_A","prop_C","prop_E")
formatOutputMatrices(AceFit,ACEcovMatrices,ACEcovLabels,vars,4)

# -----------------------------------------------------------
# RUN SUBMODELS
# -----------------------------------------------------------

# Fit AE model
# -----------------------------------------------------------------
AeModel	<- mxModel( AceFit, name="AE")
AeModel	<- omxSetParameters( AeModel, labels="c11", free=FALSE, values=0 )
AeFit	<- mxRun(AeModel, intervals=T)
AeSumm  <- summary(AeFit)
AeSumm
round(AeFit@output$estimate,4)
round(AeFit$Vars@result,4)


# Run CE model
# ------------------------------------------------------------------
CeModel	<- mxModel( AceFit, name="CE")
CeModel	<- omxSetParameters( CeModel, labels="a11", free=FALSE, values=0 )
CeFit	<- mxRun(CeModel)
round(CeFit@output$estimate,4)
round(CeFit$Vars@result,4)

# ------------------------------------------------------------------------------

# Run E model
# ------------------------------------------------------------------
eModel	<- mxModel( AceFit, name="E")
eModel	<- omxSetParameters( eModel, labels="a11", free=FALSE, values=0 )
eModel	<- omxSetParameters( eModel, labels="c11", free=FALSE, values=0 )
eFit	<- mxRun(eModel)
round(eFit@output$estimate,4)
round(eFit$Vars@result,4)

# ------------------------------------------------------------------------------

# Print Comparative Fit Statistics
AceNested	<- list(AeFit, CeFit, eFit)
tableFitStatistics(AceFit,AceNested)

#round(rbind(AceFit@output$estimate,AeFit@output$estimate,CeFit@output$estimate,eFit@output$estimate),4)
round(rbind(AceFit$Vars@result,AeFit$Vars@result,CeFit$Vars@result,eFit$Vars@result),4)


