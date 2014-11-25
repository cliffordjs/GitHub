# ------------------------------------------------------------------------------
# Function for computing polychoric/polyserial/pearson correlations
# ------------------------------------------------------------------------------

polychoricMatrix <- function(data, useDeviations=TRUE) {
nvar <- dim(data)[[2]]
ncontinuous <- 0
nordinal <- 0
nthresh <- vector(mode="integer",nvar)
isord <- vector(mode="logical",nvar)
nameList <- names(data)
ordnameList <- vector(mode="character",nvar)
contnameList <- vector(mode="character",nvar)

# Make variables into factors and label correlation parameters
correlationLabels <- matrix(NA,nrow=nvar,ncol=nvar)
for (i in 1:nvar) 
{
    if (is.factor(data[,i]))
    {
        nordinal <- nordinal + 1
        nthresh[nordinal] <- length(table(data[,i]))-1
#       I think we can avoid this
#        data[,i] <- mxFactor(data[,i], c(0:nthresh[i]))
        ordnameList[nordinal] <- nameList[i]
        isord[i] <- TRUE
    }
    else 
    {
        ncontinuous <- ncontinuous + 1
        nthresh[i] <- 0
        contnameList[ncontinuous] <- nameList[i]
        isord[i] <- FALSE
    }
# Label correlation parameters
    for (k in 1:nvar) 
        {
        if (i > k) {
                    correlationLabels[i,k] <- paste("r",i,k)
                    correlationLabels[k,i] <- paste("r",i,k)
                    }
        }
}
if (nordinal>0) {ordnameList<-ordnameList[1:nordinal]} else {ordnameList <- NULL}
if (ncontinuous>0) {contnameList<-contnameList[1:ncontinuous]} else {contnameList <- NULL}
maxnthresh <- max(nthresh)

# Populate matrix with threshold deviations, starting threshold 1 at -1 and putting maximum at +1
# for efficiency, we could take a better guess to start with
minthresh <- -.5
maxthresh <- .5

# Construct either threshold deviation matrix or threshold direct estimate matrix - as long as there's at least one ordinal variable
if (nordinal > 0)
{
    if (useDeviations)
        {
            thresholdDeviationValues <- matrix(0,nrow=maxnthresh, ncol=nordinal)
            thresholdDeviationValues[1,] <- minthresh
            thresholdDeviationLbounds <- matrix(nrow=maxnthresh, ncol=nordinal)
            thresholdDeviationLabels <- matrix(nrow=maxnthresh, ncol=nordinal)
            thresholdDeviationLabels[1,] <- paste("ThresholdDeviation ", 1, 1:nordinal)
            thresholdDeviationFree <- matrix(F,nrow=maxnthresh, ncol=nordinal)
            thresholdDeviationFree[1,] <- TRUE
            iordvar <- 0
        for (i in 1:nvar) 
            { 
            if (isord[i]) 
                {
                    iordvar <- iordvar + 1
                    if(nthresh[iordvar]>1)
                    {
                        for (j in 2:nthresh[iordvar]) 
                            {
                                thresholdDeviationValues[j,iordvar] <- (maxthresh - minthresh) / nthresh[iordvar]
                                thresholdDeviationLbounds[j,iordvar] <- .001
                                thresholdDeviationLabels[j,iordvar] <- paste("ThresholdDeviation ", j, iordvar)
                                thresholdDeviationFree[j,iordvar] <- TRUE
                            }
                    }
                }
            }
        }
    else
        {
            thresholdDirectEstimatesValues <- matrix(0,nrow=maxnthresh, ncol=nordinal)
            thresholdDirectEstimatesLbounds <- matrix(-Inf,nrow=maxnthresh, ncol=nordinal)
            thresholdDirectEstimatesLabels <- matrix(nrow=maxnthresh, ncol=nordinal)
            thresholdDirectEstimatesFree <- matrix(F,nrow=maxnthresh, ncol=nordinal)
            thresholdDirectEstimatesValues[1,] <- minthresh
            thresholdDirectEstimatesLabels[1,] <- paste("ThresholdDirectEstimates ", 1, 1:nordinal)
            thresholdDirectEstimatesFree[1,] <- TRUE
            iordvar <- 0
            for (i in 1:nvar) 
                { 
                    if (isord[i]) 
                    {
                        iordvar <- iordvar + 1
                        if(nthresh[iordvar]>1)
                            {
                                for (j in 2:nthresh[iordvar]) 
                                    {
                                        thresholdDirectEstimatesValues[j,iordvar] <- minthresh + (j-1) * ((maxthresh - minthresh) / nthresh[iordvar])
                                        thresholdDirectEstimatesLabels[j,iordvar] <- paste("ThresholdDirectEstimate ", j, iordvar)
                                        thresholdDirectEstimatesFree[j,iordvar] <- TRUE
                                    }
                            }
                    }
                }
        }
}
nameList <- names(data)
tnames <- paste("Threshold",1:maxnthresh,sep='')

# Define the model
model <- mxModel('model')
model <- mxModel(model, mxMatrix("Stand", name = "R", nrow = nvar, ncol = nvar, free=TRUE, labels=correlationLabels, dimnames=list(nameList, nameList)))
model <- mxModel(model, mxMatrix("Full", name = "M", nrow = 1, ncol = nvar, free=!isord, dimnames = list('Mean', nameList)))
model <- mxModel(model, mxMatrix("Diag", name = "StdDev", nrow = nvar, ncol = nvar, free=!isord, values=1, lbound=.01, dimnames=list(nameList, nameList)))
model$expCov <- mxAlgebra(StdDev %&% R, dimnames=list(nameList,nameList))

# Algebra to compute Threshold matrix
if (nordinal > 0)
{
    if (useDeviations) 
    {
        # For Multiplication
        model <- mxModel(model, mxMatrix("Lower", name="UnitLower", nrow = maxnthresh, ncol = maxnthresh, free=F, values=1))
        # Threshold differences:
        model <- mxModel(model, mxMatrix("Full", name="thresholdDeviations", nrow = maxnthresh, ncol = nordinal, free=thresholdDeviationFree,   values=thresholdDeviationValues, lbound=thresholdDeviationLbounds, labels = thresholdDeviationLabels))
        model <- mxModel(model, mxAlgebra(UnitLower %*% thresholdDeviations, dimnames=list(tnames,ordnameList), name="thresholds"))
        }
    else 
    {
        model <- mxModel(model, mxMatrix("Full", name="thresholds", ncol = nordinal, nrow = maxnthresh, free=thresholdDirectEstimatesFree, values=thresholdDirectEstimatesValues, lbound=thresholdDirectEstimatesLbounds, labels = thresholdDirectEstimatesLabels))
        dimnames(model$thresholds)=list(tnames,ordnameList)
    }
}

# Define the objective function
if (nordinal > 0)
{
    objective <- mxFIMLObjective(covariance="expCov", means="M", thresholds="thresholds", threshnames=ordnameList)
}
else
{
    objective <- mxFIMLObjective(covariance="expCov", means="M")
}

# Define the observed covariance matrix
dataMatrix <- mxData(data, type='raw')

# Add the objective function and the data to the model
model <- mxModel(model, objective, dataMatrix)

# Run the job
model <- mxRun(model)

# Populate seMatrix for return
seMatrix <- matrix(NA,nvar,nvar)
k<-0
for (i in 1:nvar){
    for (j in i:nvar){
        if(i != j) {
            k <- k+1
            seMatrix[i,j] <- model@output$standardErrors[k]
            seMatrix[j,i] <- model@output$standardErrors[k]
        }
    }
}
# Add dimnames to thresholds, which oddly are not in model$thresholds' output
if(nordinal > 0) 
{
    if(useDeviations)
    {
        thresholds <- matrix(model@output$algebras$model.thresholds, nrow=maxnthresh, ncol=nordinal, dimnames=list(tnames,ordnameList))     
    }
    else
    {
        thresholds <- matrix(model@output$matrices$model.thresholds, nrow=maxnthresh, ncol=nordinal, dimnames=list(tnames,ordnameList))     
    }
}
else
{
    thresholds <- NULL
}
# Return results      
return(list(polychorics=model$expCov@result, thresholds=thresholds, polychoricStandardErrors=seMatrix, Minus2LogLikelihood=model@output$Minus2LogLikelihood, Hessian=model@output$calculatedHessian, estHessian=model@output$estimatedHessian,estimatedModel=model))
}


# ------------------------------------------------------------------------------
# Pairwise wrapper
# ------------------------------------------------------------------------------

require(OpenMx)

 # Pairwise wrapper
 polypairwise <- function (data, useDeviations=TRUE) {
     nvar <- dim(data)[[2]]
     ncor <- nvar*(nvar-1)/2
     pairCorrelationMatrix <- matrix(diag(,nvar),nvar,nvar,dimnames=list(names(data),names(data)))
     pairErrors <- matrix(0,ncor,1)
     pairCount <- 0
     namelist <- NULL
     for (var1 in 1:(nvar-1)) {
         for (var2 in (var1+1):(nvar)) {
             pairCount <- pairCount + 1
             print(pairCount)
             tempResult <- polychoricMatrix(data[,c(var1,var2)], useDeviations)
             pairCorrelationMatrix[var1,var2] <- tempResult$polychorics[2,1]
             pairCorrelationMatrix[var2,var1] <- pairCorrelationMatrix[var1,var2]
             pairErrors[pairCount] <- tempResult$polychoricStandardErrors[2,1]
             namelist <- c(namelist,paste(names(data[var1]),names(data[var2]),sep="-"))
         }
     }
     dimnames(pairErrors) <- list(namelist,"est(SE)")
     return(list(R=pairCorrelationMatrix,SE=pairErrors))
 }