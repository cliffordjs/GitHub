
# ------------------------------------------------------------------------------
# Polychoric correlations for MM3 PGD, SI, and SUDs
## ------------------------------------------------------------------------------						
# Revision History
# ------------------------------------------------------------------------------		
#  
# James Clifford -- 01 17 2014
#
# ------------------------------------------------------------------------------

# Require the packages

rm(list = ls(all = TRUE))
require(MASS)
require(OpenMx)
require(psych)
require(polycor)

# Source the GenEpiHelperFunctions and PolychoricCorrelations

source('~/Desktop/Super wicked important R files/GenEpiHelperFunctions.R', chdir = TRUE)
source('~/Desktop/Super wicked important R files/PolychoricCorrelations_v1.R', chdir = TRUE)

# Name variables

allvars<-c("FamNo", "ID", "Zyg", "age", "pgdt1", "pgdt2", "marin", "marsx")

# Call data, NAs = 99

data<-read.table("~/Desktop/ResearchProjects/MM3/Peers_substance_DOC/1. Initial Data/Data files/mm3_all_data_1.17.2014.dat", header=F,na.strings="99", col.names=allvars)
head(data)
describe(data)

# MX Factor the data, PGD = 4 levels, SUD = 3 levels, SI = 2 levels

 Part1 <- mxFactor( x=data[,c(5:7)], levels=c(0:1)) 			
 Part2 <- mxFactor( x=data[,c(9:10)], levels=c(0:3),ordered=T)  
 Part3<-mxFactor(x=data[,11:12],levels=c(0:2),ordered=T)  
 Part4<-mxFactor(x=data[,8], levels=c(0:2), ordered=T)
 # Repack the data post mxfactor
 ordat<-cbind(Part1, Part2,Part3,Part4)
 head(ordat)
 summary(ordat)
 # use Mike's polypairwise script
 results<-polypairwise(ordat)
 results
 
 # Add back in FamNo, ID, Zyg, and Age
 data2<- cbind(data[,c(1:4)], ordat)
 summary(data2)
 
 ## Prepare twin data by splitting into twin 1 (T1) and twin 2 (T2)
 
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

dataR12 <- twindat(dat=data2, famid= "FamNo", twinid= "ID", zygosity= "Zyg")
summary(dataR12)





# ------------------------------------------------------------------------------
# Model 1: CCC Saturated for Factor Score @ Time 1 + 1 Cig SX 
# ------------------------------------------------------------------------------

dataR12$age_T1 <- ifelse(is.na(dataR12$age_T1), dataR12$age_T2, dataR12$age_T1)
dataR12$age_T1[is.na(dataR12$age_T1)] <- mean(dataR12$age_T1, na.rm = TRUE)
dataR12$age_T1<-dataR12$age_T1/100
selVars_1 <- c("pgdt1_T1", "pgdt2_T1", "MarIn_T1", "marsx_T1", "pgdt1_T2", "pgdt2_T2", "MarIn_T2","marsx_T2")
dataR12[,c("age_T1")] <- dataR12[,c("age_T1")]	# Select definition variables
mzData_1 <- subset(dataR12, Zyg==3, c(selVars_1,"age_T1"))
dzData_1 <- subset(dataR12, Zyg==4, c(selVars_1,"age_T1"))
summary(mzData_1)

## Note to self, need to go back and sum Sx counts in SPSS and redo data file, create new variable "sx" and you ought to be set. JC 10.26.2012


vars	<- c("pgdt1", "pgdt2", "MarIn","marsx", "age")
nv 		<- 4		# number of variables per twin
ntv 	<- nv*2		# number of variables per pair

nth 	<- 3		# Number of thresholds
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
a    <-mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T,T, F,T,T,T,  F,F,T,F, F,F,F,T), labels=aLabs, lbound=-0.99, ubound=0.99, values=c(0.6,0.3,0.5,0.5,  0.0,0.8,0.7,0.5, 0.0,0.0,0.7,0.0, 0.0,0.0,0.0,0.7), name="a" ) 
c    <-mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T,T, F,T,T,T,  F,F,T,F, F,F,F,T), labels=cLabs, lbound=-0.99, ubound=0.99, values=c(0.1,0.1,0.1,0.1,  0.0,0.1,0.1,0.1, 0.0,0.0,0.1,0.0, 0.0,0.0,0.0,0.1), name="c" ) 
e    <-mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T,T, F,T,T,T,  F,F,T,F, F,F,F,T), labels=eLabs, lbound=-0.99, ubound=0.99, values=c(0.5,0.3,0.5,0.5,  0.0,0.4,0.5,0.5, 0.0,0.0,0.5,0.0, 0.0,0.0,0.0,0.5), name="e" ) 
	
# Algebra to generate Matrices to hold A, C, and E computed Variance Components
A    <-mxAlgebra( expression=a %*% t(a), name="A" )
C    <-mxAlgebra( expression=c %*% t(c), name="C" ) 
E    <-mxAlgebra( expression=e %*% t(e), name="E" )

# Algebra to compute Total Variance
V    <-mxAlgebra( expression=A+C+E, name="V" )

# Matrices for expected Means & Thresholds (on liabilities) 
expMean <-mxMatrix( type="Zero", nrow=1, ncol=ntv, name="expMean" )

B4age	<- mxMatrix( type="Full", nrow=1, ncol=2, free=T, labels=c("Bagev1", "Bagev2"), values=0.4, lbound=-3, ubound=3, name="B4age" ) 
age   	<- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age_T1"), name="age")
i	    <- mxMatrix( type="Unit", nrow=3, ncol=2, name="i"  )

Inc     <- mxMatrix( type="Lower", nrow=nth, ncol=nth, free=F, values=1, name="Inc" )
Thre    <-mxMatrix( type="Full", nrow=3, ncol=nv, free=c(T,T,T, T,T,T, T,F,F, T,T,F), values=c(-1.1,0.1,0.5,  1.1,0.2,0.5, 0.5,0,0, 0.5,0.5,0), labels=cbind(var1th,var2th,var3th,var4th), lbound=thLB, ubound=1, name="Thre")
ExpThre   <- mxAlgebra( expression= cbind( ( Inc %*% Thre  -  (age %x% B4age) %x%i ), 
										   ( Inc %*% Thre  -  (age %x% B4age) %x%i ) ), name="ExpThre" )

# Algebra to generate Matrices to hold Parameter Estimates and Derived Variance Components
rowVars  <-rep(c('vars', nv))
colVars  <-rep(c('A','C','E','SA','SC','SE'),each=nv)
estVars  <-mxAlgebra( expression=cbind(A,C,E,A/V,C/V,E/V), dimnames=list(rowVars,colVars), name="Vars" )
                                          
## Beta matrix

b    <-mxMatrix( type="Full", nrow=nv, ncol=nv, free=c(F,F,F,F, F,F,F,F, F,F,F,T, F,F,F,F), labels=c("b11","b21","b31","b41", "b12", "b22", "b32", "b42","b13", "b23", "b33","b43", "b14", "b24", "b34", "b44"), values=c(0,0,0,0, 0,0,0,0, 0,0,0,0.2, 0,0,0,0), lbound=-0.99, ubound=1, name="b" )
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
pars	<-list( a, c, e, A, C, E, V, expMean, Inc, B4age, i,Thre, Unv1,b, j, k, DoC )
#B4age, i,
defs	<- list( age,ExpThre)
modelMZ	<-mxModel(pars, defs, expMean, expCovMZ, dataMZ, objMZ, name="MZ")
#defs
modelDZ	<-mxModel(pars, defs, expMean, expCovDZ, dataDZ, objDZ, name="DZ")
minus2ll<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL")
obj	<-mxAlgebraObjective( "m2LL")
ciVC	<-mxCI(c('b','a','c','e'), interval = 0.95)
Model_1_Chol <-mxModel("Model_1_Chol", pars,Var1, modelMZ, modelDZ, minus2ll, obj, ciVC)

# -------------------------------------------------
# RUN GENETIC MODEL
# -------------------------------------------------

# Run ACE Model
More_precise <- mxOption(Model_1_Chol,"Function precision", 1e-100) # set the precision
Model_1_Chol_Fit  <- mxRun( More_precise, intervals=F )
Model_1_Chol_Summ <- summary( Model_1_Chol_Fit )
Model_1_Chol_Summ

parameterSpecifications(Model_1_Chol)
tableFitStatistics(Model_1_Chol_Fit)
# ------------------------------------------------------------------------------
# Model 2: No correlated liability, only Beta 3 pathway remains 
# ------------------------------------------------------------------------------


Model_1a_b3   <- mxModel( Model_1_Chol, name="Model 2" )
Model_1a_b3   <- omxSetParameters( Model_1a_b3, labels=c("a31","c31","e31","a41","c41", "e41"), free=FALSE, values=0 )
more_precise3<-mxOption(Model_1a_b3, "Function precision", 1e-100)
Model_1a_b3_Fit <-mxRun(more_precise3, intervals=F)
Model_1a_b3_Summ<-summary(Model_1a_b3_Fit)
Model_1a_b3_Summ

SatNested<-list(Model_1a_b3_Fit)
tableFitStatistics(Model_1_Chol_Fit, SatNested)

# ------------------------------------------------------------------------------
# Model 3: No correlated liability from either PGD to initiation, only Beta 3 pathway remains 
# ------------------------------------------------------------------------------


Model_1b_b3   <- mxModel( Model_1a_b3, name="Model 3" )
Model_1b_b3   <- omxSetParameters( Model_1b_b3, labels=c("a31","c31","e31","a32","c32", "e32"), free=FALSE, values=0 )
Model_1b_b3<-omxSetParameters(Model_1b_b3, labels=c("b32"), free =T, values=0.5)
more_precise3<-mxOption(Model_1b_b3, "Function precision", 1e-250)
Model_1b_b3_Fit <-mxRun(more_precise3)
Model_1b_b3_Summ<-summary(Model_1b_b3_Fit)
Model_1b_b3_Summ


SatNested <- list(Model_1a_b3_Fit, Model_1b_b3_Fit)
tableFitStatistics(Model_1_Chol_Fit, SatNested)

# ------------------------------------------------------------------------------
# Model 4: Constraining C to be equal from both PGD timpoints to initiaiton
# ------------------------------------------------------------------------------

Model_1c_b3   <- mxModel( Model_1b_b3, name="Model 4" )
Model_1c_b3   <- omxSetParameters( Model_1c_b3, labels=c("a21", "c21", "e21"), free=FALSE, values=0 )
Model_1c_b3<- omxSetParameters(Model_1c_b3, labels=c("b21"), free =T, values=0.5)
more_precise4<-mxOption(Model_1c_b3, "Function precision", 1e-250)
Model_1c_b3_Fit <-mxRun(more_precise4)
Model_1c_b3_Summ<-summary(Model_1c_b3_Fit)
Model_1c_b3_Summ


SatNested <- list(Model_1a_b3_Fit, Model_1b_b3_Fit, Model_1c_b3_Fit)
tableFitStatistics(Model_1_Chol_Fit, SatNested)


 
 