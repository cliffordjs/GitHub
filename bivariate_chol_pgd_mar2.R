
# ------------------------------------------------------------------------------
# CCC MODEL for MM3 PGD and Marijuana
## ------------------------------------------------------------------------------						
# Revision History
# ------------------------------------------------------------------------------		
#  
# James Clifford -- 12 15 2012
# James Clifford -- 01 25 2013: added Model 2b
# James Clifford -- 01 31 2013: Added models 2b.1/2, 2c.1/2
#
# ------------------------------------------------------------------------------

rm(list = ls(all = TRUE))
setwd("~/Desktop/ResearchProjects/MM3")
source('http://openmx.psyc.virginia.edu/getOpenMx.R')
1
require(OpenMx)
require(psych)
require(polycor)

source('~/Desktop/Super wicked important R files/GenEpiHelperFunctions.R', chdir = TRUE)

allvars <- c("FamNo", "ID", "h1", "h51", "h63",  "f1", "f2", "f3", "f4", "zyg", "age", "temp1", "temp2", "temp3", "temp4", "epoch2", "epoch3", "marsx", "fst1", "fst2")

# Specify data file location 			  
dataR1 <- read.table("~/Desktop/ResearchProjects/MM3/Peers_substance_DOC/1. Initial Data/Data files/mm3_final_12.15.2012_cannabis.dat",header=F, na.strings="99", col.names=allvars)
summary(dataR1)


# Combines Ordinal data prior to factor with IDs and Zyg	
data <- cbind(dataR1[,1:5], dataR1[,12:20], dataR1[,10:11])	
summary(data)
describe(data)


# ------------------------------------------------------------------------------
# Convert / declare ordinal data using mx Factor
# ------------------------------------------------------------------------------

 Part1 <- mxFactor( x=data[,c(3:5)], levels=c(0:1)) 			# Sets drug measures to ORD data with 2 levels
 Part2 <- mxFactor( x=data[,c(6:11)], levels=c(0:2),ordered=T)  # Sets levels for factor scores and epoch ordinal data	
# Part3 <- mxFactor( x=data[,c(12)], levels=c(3:4),ordered=T)  
# Part3 <- mxFactor( x=cut(data[,13]/100, breaks=c(.24,.32,.39,.47,.62), ordered_result=T, labels=c(0,1,2,3)), levels=c(0:3),ordered=T)
 Part4<-mxFactor(x=data[,12:14],levels=c(0:2),ordered=T)  # Declares sx counts ordinal

summary(Part1)
summary(Part2)
summary(Part4)

DataOrdF <- cbind(data[,1:2], data[,15:15], Part1, Part2, data[,16]/100, Part4) # Recombines ordinalized data with identifiers and covariate
colnames <- c( "FamNo", "ID", "zyg", "h1", "h51", "h63", "temp1", "temp2", "temp3", "temp4", "epoch2", "epoch3","age", "marsx", "fst1", "fst2")
colnames(DataOrdF) <- colnames
summary(DataOrdF)

polychor(DataOrdF$age, DataOrdF$h1, std.err=T)

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
varnames <- colnames(dataR12)
summary(dataR12)





# ------------------------------------------------------------------------------
# Model 1: CCC Saturated for Factor Score @ Time 1 + 1 Cig SX 
# ------------------------------------------------------------------------------

dataR12$age_T1 <- ifelse(is.na(dataR12$age_T1), dataR12$age_T2, dataR12$age_T1)
dataR12$age_T1[is.na(dataR12$age_T1)] <- mean(dataR12$age_T1, na.rm = TRUE)

selVars_1 <- c("h63_T1", "marsx_T1", "h63_T2","marsx_T2")
dataR12[,c("age_T1")] <- dataR12[,c("age_T1")]	# Select definition variables
mzData_1 <- subset(dataR12, zyg==3, c(selVars_1,"age_T1"))
dzData_1 <- subset(dataR12, zyg==4, c(selVars_1,"age_T1"))
summary(mzData_1)

## Note to self, need to go back and sum Sx counts in SPSS and redo data file, create new variable "sx" and you ought to be set. JC 10.26.2012


vars	<- c( "h63","marsx", "age")
nv 		<- 2		# number of variables per twin
ntv 	<- nv*2		# number of variables per pair

nth 	<- 2		# Number of thresholds
lth     <- -0.438    	# start values for first threshold
ith     <-  0.5    	# start values for threshold increments

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
var1th     <- paste("var1","_th",1:nth, sep="") 
var2th     <- paste("var2","_th",1:nth, sep="")
var3th		<- paste("var3","_th",1:nth, sep="")
var4th		<- paste("var4","_th",1:nth, sep="")
dz1th     <- paste("DZ1","_th",1:nth, sep="")
dz2th     <- paste("DZ2","_th",1:nth, sep="")
dz3th		<- paste("DZ23","_th",1:nth, sep="")
aLabs     <- paste("a",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")
cLabs     <- paste("c",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")
eLabs     <- paste("e",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")

# thFree    <- c(rep(T,nth1),rep(F,nth-nth1),rep(T,nth2),rep(F,nth-nth2),rep(T,nth3),rep(F,nth-nth3),rep(T,nth4),rep(F,nth-nth4))

# ACE Model with one overall Threshold
# Matrices to store a, c, and e Path Coefficients
a    <-mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F, F,T), labels=aLabs, lbound=-0.99, ubound=0.99, values=c(0.5, 0, 0, 0.5), name="a" ) 
c    <-mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F, F, T), labels=cLabs, lbound=-0.99, ubound=0.99, values=c(0.5, 0, 0, 0.5), name="c" ) 
e    <-mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F, F, T), labels=eLabs, lbound=-0.99, ubound=0.99, values=c(0.5, 0, 0, 0.5), name="e" ) 
	
# Algebra to generate Matrices to hold A, C, and E computed Variance Components
A    <-mxAlgebra( expression=a %*% t(a), name="A" )
C    <-mxAlgebra( expression=c %*% t(c), name="C" ) 
E    <-mxAlgebra( expression=e %*% t(e), name="E" )

# Algebra to compute Total Variance
V    <-mxAlgebra( expression=A+C+E, name="V" )

# Matrices for expected Means & Thresholds (on liabilities) 
expMean <-mxMatrix( type="Zero", nrow=1, ncol=ntv, name="expMean" )

B4age	<- mxMatrix( type="Full", nrow=1, ncol=2, free=T, labels=c("Bagev1", "Bagev2"), values=0.3, lbound=-3, ubound=3, name="B4age" ) 
age   	<- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age_T1"), name="age")
i	    <- mxMatrix( type="Unit", nrow=2, ncol=1, name="i"  )

Inc     <- mxMatrix( type="Lower", nrow=nth, ncol=nth, free=F, values=1, name="Inc" )
Thre    <-mxMatrix( type="Full", nrow=2, ncol=nv, free=c(T,F, T,T), values=c(-1.1,0,  1.1,0.2), labels=cbind(var1th,var2th), lbound=thLB, ubound=4, name="Thre")
ExpThre   <- mxAlgebra( expression= cbind( ( Inc %*% Thre  -  (age %x% B4age) %x%i ), 
										   ( Inc %*% Thre  -  (age %x% B4age) %x%i ) ), name="ExpThre" )

# Note that I've taken out the age adjustment, but left B4age in, thus getting mx Red status code. 


# Algebra to generate Matrices to hold Parameter Estimates and Derived Variance Components
rowVars  <-rep(c('vars', nv))
colVars  <-rep(c('A','C','E','SA','SC','SE'),each=nv)
estVars  <-mxAlgebra( expression=cbind(A,C,E,A/V,C/V,E/V), dimnames=list(rowVars,colVars), name="Vars" )
                                          
## Beta matrix

b    <-mxMatrix( type="Full", nrow=nv, ncol=nv, free=c(F,T, F,F), labels=c("b11","b21", "b12","b22"), values=c(0,0.8, 0.0,0.0), lbound=-0.99, ubound=1, name="b" )
j    <-mxMatrix( type="Full", nrow=2, ncol=2, free=F, values=1, name="j" ) 
k    <-mxMatrix( type="Full", nrow=nv, ncol=nv,free = F, values=1, name="k" ) 
DoC     <-mxAlgebra( expression=j%x%(k-b), name="DoC" )


# Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
expCovMZ    <-mxAlgebra( expression= DoC%&%(rbind( cbind(A+C+E , A+C),
                                          cbind(A+C   , A+C+E))), name="expCovMZ" )
expCovDZ    <-mxAlgebra( expression= DoC%&%(rbind( cbind(A+C+E     , 0.5%x%A+C),
                                          cbind(0.5%x%A+C , A+C+E))), name="expCovDZ" )
                                          
# Constraint on variance of the liability of Binary variables (assumed to have a SND) 
Unv1	<-mxMatrix( type="Unit", nrow=nv, ncol=1, name="Unv1" )
Var1	<-mxConstraint( expression=diag2vec(V)==Unv1, name="Var1" )

# Data objects for Multiple Groups
dataMZ   <-mxData( observed=mzData_1, type="raw" )
dataDZ   <-mxData( observed=dzData_1, type="raw" )

# Objective objects for Multiple Groups
objMZ	<-mxFIMLObjective( covariance="expCovMZ", means="expMean", dimnames=selVars_1, thresholds="ExpThre" )
objDZ	<-mxFIMLObjective( covariance="expCovDZ", means="expMean", dimnames=selVars_1, thresholds="ExpThre" )

# Combine Groups
pars	<-list( a, c, e, A, C, E, V, expMean,B4age, i, Inc, Thre,b, j, k, DoC, Unv1 )
#B4age, i,
defs	<- list( age,ExpThre)
modelMZ	<-mxModel( pars, defs, expCovMZ, dataMZ, objMZ, name="MZ" )
#defs
modelDZ	<-mxModel(pars, defs, expCovDZ, dataDZ, objDZ, name="DZ" )
minus2ll<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj	<-mxAlgebraObjective( "m2LL" )
ciVC	<-mxCI(c('b'), interval = 0.95)
#fitfunc<-mxFitFunctionML()
Model_1_Chol <-mxModel( "Model_1_Chol", pars, modelMZ, modelDZ, Var1, minus2ll, obj, ciVC )

# -------------------------------------------------
# RUN GENETIC MODEL
# -------------------------------------------------

# Run ACE Model
More_precise <- mxOption(Model_1_Chol,"Function precision", 1e-100) # set the precision
Model_1_Chol_Fit  <- mxRun( More_precise, intervals=F, unsafe=T )
Model_1_Chol_Summ <- summary( Model_1_Chol_Fit )
Model_1_Chol_Summ

parameterSpecifications(Model_1_Chol_Fit)

tableFitStatistics(Model_1_Chol_Fit)


# ------------------------------------------------------------------------------
# Model 2: No correlated liability from PGD to sx count, only Beta 3 pathway remains 
# ------------------------------------------------------------------------------


Model_1_b3   <- mxModel( Model_1_Chol, name="Model_1_b3" )
Model_1_b3   <- omxSetParameters( Model_1_b3, labels=c("a41","c41","e41", "a42","c42", "e42"), free=FALSE, values=0 )
more_precise22<-mxOption(Model_1_b3, "Function precision", 1e-100)
Model_1_b3_Fit <-mxRun(more_precise22)
Model_1_b3_Summ<-summary(Model_1_b3_Fit)
Model_1_b3_Summ

tableFitStatistics(Model_1_Chol_Fit, Model_1_b3_Fit)

# ------------------------------------------------------------------------------
# Model 2a: No correlated liability from either PGD to initiation, only Beta 3 pathway remains 
# ------------------------------------------------------------------------------


Model_1a_b3   <- mxModel( Model_1_b3, name="Model_1a_b3" )
Model_1a_b3   <- omxSetParameters( Model_1a_b3, labels=c("a31","c31","e31","a32","c32", "e32"), free=FALSE, values=0 )
more_precise3<-mxOption(Model_1a_b3, "Function precision", 1e-250)
Model_1a_b3_Fit <-mxRun(more_precise3)
Model_1a_b3_Summ<-summary(Model_1a_b3_Fit)
Model_1a_b3_Summ


SatNested <- list(Model_1_b3_Fit, Model_1a_b3_Fit)
tableFitStatistics(Model_1_Chol_Fit, SatNested)

# ------------------------------------------------------------------------------
# Model 2b: Constraining C to be equal from both PGD timpoints to initiaiton
# ------------------------------------------------------------------------------

Model_1b_b3   <- mxModel( Model_1_b3, name="Model_2b" )
Model_1b_b3   <- omxSetParameters( Model_1b_b3, labels=c("c31","e31","c32","e32"), free=FALSE, values=0 )
more_precise4<-mxOption(Model_1b_b3, "Function precision", 1e-250)
Model_1b_b3_Fit <-mxRun(more_precise4)
Model_1b_b3_Summ<-summary(Model_1b_b3_Fit)
Model_1b_b3_Summ


SatNested <- list(Model_1_b3_Fit, Model_1a_b3_Fit, Model_1b_b3_Fit)
tableFitStatistics(Model_1_Chol_Fit, SatNested)


# Model 2b.1: Constraining C at timepoint 1

Model_1b1   <- mxModel( Model_1_b3, name="Model_2b.1" )
Model_1b1  <- omxSetParameters( Model_1b1, labels=c("c31","e31","e32"), free=FALSE, values=0 )
more_precise4<-mxOption(Model_1b1, "Function precision", 1e-250)
Model_1b1_Fit <-mxRun(more_precise4)
Model_1b1_Summ<-summary(Model_1b1_Fit)
Model_1b1_Summ


SatNested <- list(Model_1_b3_Fit, Model_1a_b3_Fit, Model_1b_b3_Fit, Model_1b1_Fit)
tableFitStatistics(Model_1_Chol_Fit, SatNested)

# Model 2b.2: Constraining C at timepoint 2

Model_1b2   <- mxModel( Model_1_b3, name="Model_2b.2" )
Model_1b2   <- omxSetParameters( Model_1b2, labels=c("e31","c32","e32"), free=FALSE, values=0 )
more_precise4<-mxOption(Model_1b2, "Function precision", 1e-250)
Model_1b2_Fit <-mxRun(more_precise4)
Model_1b2_Summ<-summary(Model_1b2_Fit)
Model_1b2_Summ


SatNested <- list(Model_1_b3_Fit, Model_1a_b3_Fit, Model_1b_b3_Fit, Model_1b1_Fit, Model_1b2_Fit)
tableFitStatistics(Model_1_Chol_Fit, SatNested)


# ------------------------------------------------------------------------------
# Model 2c: Constraining A to be equal from both PGD timpoints to initiaiton
# ------------------------------------------------------------------------------

Model_1c_b3   <- mxModel( Model_1_b3, name="Model_2c" )
Model_1c_b3   <- omxSetParameters( Model_1c_b3, labels=c("a31","e31","a32","e32"), free=FALSE, values=0 )
more_precise5<-mxOption(Model_1c_b3, "Function precision", 1e-250)
Model_1c_b3_Fit <-mxRun(more_precise5)
Model_1c_b3_Summ<-summary(Model_1c_b3_Fit)
Model_1c_b3_Summ


SatNested <- list(Model_1_b3_Fit, Model_1a_b3_Fit, Model_1b_b3_Fit, Model_1c_b3_Fit)
#tableFitStatistics(Model_1_Chol_Fit, SatNested)

# Model 2c.1: Constraining A from Timepoint 1

Model_1c1   <- mxModel( Model_1_b3, name="Model_2c.1" )
Model_1c1   <- omxSetParameters( Model_1c1, labels=c("a31","e31","e32"), free=FALSE, values=0 )
more_precise5<-mxOption(Model_1c1, "Function precision", 1e-250)
Model_1c1_Fit <-mxRun(more_precise5)
Model_1c1_Summ<-summary(Model_1c1_Fit)
Model_1c1_Summ


SatNested <- list(Model_1_b3_Fit, Model_1a_b3_Fit, Model_1b_b3_Fit, Model_1c_b3_Fit, Model_1c1_Fit)
#tableFitStatistics(Model_1_Chol_Fit, SatNested)


# Model 2c.2: Constraining A from timepoint 2

Model_1c2   <- mxModel( Model_1_b3, name="Model_2c.2" )
Model_1c2   <- omxSetParameters( Model_1c2, labels=c("e31","a32","e32"), free=FALSE, values=0 )
more_precise5<-mxOption(Model_1c2, "Function precision", 1e-250)
Model_1c2_Fit <-mxRun(more_precise5)
Model_1c2_Summ<-summary(Model_1c2_Fit)
Model_1c2_Summ


SatNested2 <- list(Model_1_b3_Fit, Model_1a_b3_Fit, Model_1b_b3_Fit,Model_1b1_Fit, Model_1b2_Fit, Model_1c_b3_Fit, Model_1c1_Fit, Model_1c2_Fit)
tableFitStatistics(Model_1_Chol_Fit, SatNested2)


sink("CCC_Mar.out", append=F, split=T)
sink()
# ------------------------------------------------------------------------------
# Model 3: No correlated liability 
# ------------------------------------------------------------------------------


## Set everything so there's only a11, a22, a33 and such
# ACE Model with one overall Threshold
# Matrices to store a, c, and e Path Coefficients
a    <-mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,F,F, F,T,F,F, F,F,T,F, F,F,F,T), labels=aLabs, values=c(0.1,0,0,0, 0,0.1,0,0, 0,0,0.8,0,  0,0,0,0.6), lbound=-0.99, ubound=0.99, name="a" ) 
c    <-mxMatrix( type="Lower", nrow=nv,
 ncol=nv, free=c(T,F,F,F, F,T,F,F, F,F,T,F, F,F,F,T), labels=cLabs, values=c(0.5,0,0,0, 0,0.1,0,0, 0,0,0.1,0,  0,0,0,0.1), lbound=-0.99, ubound=0.99, name="c" ) 
e    <-mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,F,F, F,T,F,F, F,F,T,F, F,F,F,T), labels=eLabs, values=c(0.8,0,0,0, 0,0.6,0,0, 0,0,0.5,0,  0,0,0,0.7), lbound=-0.99, ubound=0.99, name="e" ) 
	
A    <-mxAlgebra( expression=a %*% t(a), name="A" )
C     <-mxAlgebra( expression=c %*% t(c), name="C" ) 
E     <-mxAlgebra( expression=e %*% t(e), name="E" )

# Algebra to compute Total Variance
V     <-mxAlgebra( expression=A+C+E, name="V" )

# Matrices for expected Means & Thresholds (on liabilities) 
expMean    <-mxMatrix( type="Zero", nrow=1, ncol=ntv, name="expMean" )

B4age	<- mxMatrix( type="Full", nrow=1, ncol=2, free=T, labels=c("Bagev1","Bagev2"), lbound=-3, ubound=3, name="B4age" ) 
age   	<- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("age_T1"), name="age")
i	    <- mxMatrix( type="Unit", nrow=2, ncol=1, name="i"  )

Inc     <- mxMatrix( type="Lower", nrow=nth, ncol=nth, free=F, values=1, name="Inc" )
Thre    <-mxMatrix( type="Full", nrow=2, ncol=nv, free=c(T,T, T,T, T,T, T,T), values=c(-1.1,0.1,  1.1,0.5, 0.5,0.5, 0.5,0.5), labels=cbind(var1th,var2th,var3th,var4th), lbound=thLB, ubound=4, name="Thre")
ExpThre   <- mxAlgebra( expression= cbind( ( Inc %*% Thre ), 
										   ( Inc %*% Thre  ) ), name="ExpThre" )
										   
#rowVars  <-rep('vars',nv)
#colVars  <-rep(c('A','C','E','SA','SC','SE'),each=nv)
#estVars  <-mxAlgebra( expression=cbind(A,C,E,A/V,C/V,E/V), name="Vars", dimnames=list(rowVars,colVars))

b    <-mxMatrix( type="Full", nrow=nv, ncol=nv, free=c(F,T,T,T, F,F,T,T, F,F,F,T, F,F,F,F), labels=c("b11","b21","b31","b41", "b12", "b22", "b32", "b42","b13", "b23", "b33","b43", "b14", "b24", "b34", "b44"), values=c(0,0.8,0.5,0.5, 0,0,0.2,0.3, 0,0,0,0.5, 0,0,0,0), lbound=-0.99, ubound=1, name="b" )
j    <-mxMatrix( type="Iden", nrow=2, ncol=2, name="j" ) 
k    <-mxMatrix( type="Iden", nrow=nv, ncol=nv, name="k" ) 
DoC     <-mxAlgebra( expression=j%x%(k-b), name="DoC" )

# Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
expCovMZ    <-mxAlgebra( expression= DoC%&%(rbind( cbind(A+C+E , A+C),
                                          cbind(A+C   , A+C+E))), name="expCovMZ" )
expCovDZ    <-mxAlgebra( expression= DoC%&%(rbind( cbind(A+C+E     , 0.5%x%A+C),
                                          cbind(0.5%x%A+C , A+C+E))), name="expCovDZ" )

# Constraint on variance of the liability of Binary variables (assumed to have a SND) 
Unv1	<-mxMatrix( type="Unit", nrow=nv, ncol=1, name="Unv1" )
Var1	<-mxConstraint( expression=diag2vec(V)==Unv1, name="Var1" )

# Data objects for Multiple Groups
dataMZ   <-mxData( observed=mzData_1, type="raw" )
dataDZ   <-mxData( observed=dzData_1, type="raw" )

# Objective objects for Multiple Groups
objMZ	<-mxFIMLObjective( covariance="expCovMZ", means="expMean", dimnames=selVars_1, thresholds="ExpThre" )
objDZ	<-mxFIMLObjective( covariance="expCovDZ", means="expMean", dimnames=selVars_1, thresholds="ExpThre" )

# Combine Groups
pars	<-list( a, c, e, A, C, E, V, Inc, Thre, i, B4age, Unv1, b, j, k, DoC )
defs	<- list(age,ExpThre)
modelMZ	<-mxModel(  pars, defs, expMean, expCovMZ, dataMZ, objMZ, name="MZ" )
modelDZ	<-mxModel(  pars, defs, expMean, expCovDZ, dataDZ, objDZ, name="DZ" )
minus2ll<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj	<-mxAlgebraObjective( "m2LL" )
ciVC	<-mxCI(c('b'), interval = 0.95)
Model_3_nocor <-mxModel( "Model_3_nocor", pars, Var1, modelMZ, modelDZ, minus2ll, obj, ciVC )

# Run ACE model
More_precise2 <- mxOption(Model_3_nocor,"Function precision", 1e-100) # set the precision
Model_3_nocor_Fit  <- mxRun( More_precise2, intervals=F )
Model_3_nocor_Summ <- summary( Model_3_nocor_Fit )
Model_3_nocor_Summ

SatNested<-list(Model_1_b3_Fit, Model_1a_b3_Fit, Model_3_nocor_Fit)
tableFitStatistics(Model_1_Chol_Fit,SatNested)

# Model 3.1: No PGDt1 to sx count

Model_2_nob41     <- mxModel( Model_3_nocor_Fit, name="Model_2_nob41" )
Model_2_nob41_Fit    <- omxSetParameters( Model_2_nob41, free=c(F,T,T,F, F,F,T,T, F,F,F,T, F,F,F,F), labels=c("b11","b21","b31","b41", "b12", "b22", "b32", "b42","b13", "b23", "b33","b43", "b14", "b24", "b34", "b44"), values=c(0,0.5,0.5,0, 0,0,0.5,0.5, 0,0,0,0.5, 0,0,0,0) )
Model_2_nob41_Fit<- mxOption(Model_2_nob41_Fit,"Function precision", 1e-100) # set the precision
Model_2_nob41_Fit   <- mxRun(Model_2_nob41_Fit)
Model_2_nob41_Fit_Summ   <- summary(Model_2_nob41_Fit)
Model_2_nob41_Fit_Summ

SatNested<-list(Model_1_b3_Fit, Model_1a_b3_Fit, Model_3_nocor_Fit, Model_2_nob41_Fit)
tableFitStatistics(Model_1_Chol_Fit, SatNested)

# Model 3.2: No PGDt2 to Sx count

Model_2_nob42     <- mxModel( Model_2_nob41_Fit, name="Model_2_nob42" )
Model_2_nob42_Fit    <- omxSetParameters( Model_2_nob42, free=c(F,T,T,F, F,F,T,F, F,F,F,T, F,F,F,F), labels=c("b11","b21","b31","b41", "b12", "b22", "b32", "b42","b13", "b23", "b33","b43", "b14", "b24", "b34", "b44"), values=c(0,0.5,0.5,0, 0,0,0.5,0, 0,0,0,0.5, 0,0,0,0) )
Model_2_nob42_Fit<- mxOption(Model_2_nob42_Fit,"Function precision", 1e-100) # set the precision
Model_2_nob42_Fit   <- mxRun(Model_2_nob42_Fit)
Model_2_nob42_Fit_Summ   <- summary(Model_2_nob42_Fit)
Model_2_nob42_Fit_Summ

SatNested<-list(Model_1_b3_Fit, Model_1a_b3_Fit, Model_3_nocor_Fit, Model_2_nob41_Fit, Model_2_nob42_Fit)
tableFitStatistics(Model_1_Chol_Fit, SatNested)



# Model 3.3: No PGD to PGD

Model_2_nob21     <- mxModel( Model_2_nob42_Fit, name="Model_2_nob21" )
Model_2_nob21_Fit    <- omxSetParameters( Model_2_nob21, free=c(F,F,T,F, F,F,T,F, F,F,F,T, F,F,F,F), labels=c("b11","b21","b31","b41", "b12", "b22", "b32", "b42","b13", "b23", "b33","b43", "b14", "b24", "b34", "b44"), values=c(0,0,0.5,0, 0,0,0.5,0, 0,0,0,0.5, 0,0,0,0) )
Model_2_nob21_Fit<- mxOption(Model_2_nob21_Fit,"Function precision", 1e-100) # set the precision
Model_2_nob21_Fit   <- mxRun(Model_2_nob21_Fit)
Model_2_nob21_Fit_Summ   <- summary(Model_2_nob21_Fit)
Model_2_nob21_Fit_Summ

SatNested<-list(Model_1_b3_Fit, Model_1a_b3_Fit, Model_3_nocor_Fit, Model_2_nob41_Fit, Model_2_nob42_Fit, Model_2_nob21_Fit)
tableFitStatistics(Model_1_Chol_Fit, SatNested)

# Model 3.4 No PGDt2 to initiation

Model_2_nob32     <- mxModel( Model_2_nob21_Fit, name="Model_2_nob32" )
Model_2_nob32_Fit    <- omxSetParameters( Model_2_nob32, free=c(F,F,T,F, F,F,F,F, F,F,F,T, F,F,F,F), labels=c("b11","b21","b31","b41", "b12", "b22", "b32", "b42","b13", "b23", "b33","b43", "b14", "b24", "b34", "b44"), values=c(0,0,0.5,0, 0,0,0,0, 0,0,0,0.5, 0,0,0,0) )
Model_2_nob32_Fit<- mxOption(Model_2_nob32_Fit,"Function precision", 1e-100) # set the precision
Model_2_nob32_Fit   <- mxRun(Model_2_nob32_Fit)
Model_2_nob32_Fit_Summ   <- summary(Model_2_nob32_Fit)
Model_2_nob32_Fit_Summ

SatNested<-list(Model_1_b3_Fit, Model_1a_b3_Fit, Model_3_nocor_Fit, Model_2_nob41_Fit, Model_2_nob42_Fit, Model_2_nob21_Fit, Model_2_nob32_Fit)
tableFitStatistics(Model_1_Chol_Fit, SatNested)

# Model 3.4a No PGDt1 to initiation

Model_2_nob31     <- mxModel( Model_2_nob32_Fit, name="Model_2_nob31" )
Model_2_nob31_Fit    <- omxSetParameters( Model_2_nob31, free=c(F,F,F,F, F,F,F,F, F,F,F,T, F,F,F,F), labels=c("b11","b21","b31","b41", "b12", "b22", "b32", "b42","b13", "b23", "b33","b43", "b14", "b24", "b34", "b44"), values=c(0,0,0,0, 0,0,0,0, 0,0,0,0.5, 0,0,0,0) )
Model_2_nob31_Fit<- mxOption(Model_2_nob31_Fit,"Function precision", 1e-100) # set the precision
Model_2_nob31_Fit   <- mxRun(Model_2_nob31_Fit)
Model_2_nob31_Fit_Summ   <- summary(Model_2_nob31_Fit)
Model_2_nob31_Fit_Summ

SatNested<-list(Model_1_b3_Fit, Model_1a_b3_Fit, Model_3_nocor_Fit, Model_2_nob41_Fit, Model_2_nob42_Fit, Model_2_nob21_Fit, Model_2_nob32_Fit, Model_2_nob31_Fit)
tableFitStatistics(Model_1_Chol_Fit, SatNested)

# Model 3.5 No initiation to sx count


Model_2_nob43     <- mxModel( Model_2_nob32_Fit, name="Model_2_nob43" )
Model_2_nob43_Fit    <- omxSetParameters( Model_2_nob43, free=c(F,F,F,F, F,F,F,F, F,F,F,F, F,F,F,F), labels=c("b11","b21","b31","b41", "b12", "b22", "b32", "b42","b13", "b23", "b33","b43", "b14", "b24", "b34", "b44"), values=c(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0) )
Model_2_nob43_Fit<- mxOption(Model_2_nob43_Fit,"Function precision", 1e-100) # set the precision
Model_2_nob43_Fit   <- mxRun(Model_2_nob43_Fit)
Model_2_nob43_Fit_Summ   <- summary(Model_2_nob43_Fit)
Model_2_nob43_Fit_Summ

SatNested<-list(Model_1_b3_Fit, Model_1a_b3_Fit, Model_3_nocor_Fit, Model_2_nob41_Fit, Model_2_nob42_Fit, Model_2_nob21_Fit, Model_2_nob32_Fit, Model_2_nob31_Fit, Model_2_nob43_Fit)
tableFitStatistics(Model_1_Chol_Fit, SatNested)