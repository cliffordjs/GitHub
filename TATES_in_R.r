rm(list=ls(all=TRUE))

###############################################################################################
#
#    #######      #      #######  ######    ####
#       #        # #        #     #        #
#       #       #   #       #     ###       ###
#       #      # # # #      #     #            #
#       #     #       #     #     ######   ####
#
#   This script runs the TATES procedure, described in
#
#   Van der Sluis, S., Posthuma, D., & Dolan, C.V.
#       Dealing with genetically heterogeneous phenotypes: 
#       TATES, a rapid and powerful trait-based test based on the GATES procedure 
#
#  The TATES procedure is inspired by the Gates procedure described in: 
#     Li, M-X, Gui, H-S., Kwan, J.S.H., & Sham, P.C. (2011). 
#     GATES: A Rapid and Powerful Gene-Based Association Test 
#     Using Extended Simes Procedure, American Journal of Human Genetics, 88,283-293.
#
###############################################################################################
#
#  nvar = number of variables
#  nsnp = number of SNPs
#
# Required:
# 1. A file with the full, symmetrical correlation matrix between the nvar variables
#      This file should NOT have a header.
#      This file should have dimensions nvar*nvar.
#      See file "Example_cor" for an example file for nvar=12.
# 2. A file with the p-values
#      This file should NOT have a header.
#      This file should have the following sctructure:
#      Column 1: Chromosome number
#      Column 2: SNP name / rs number
#      Column 3: the p-values of the relations of all SNPs with variable 1
#      Column 4: the p-values of the relations of all SNPs with variable 2
#      etc
#      This file should thus have dimensions nsnp*(nvar+2)
#      Note that the order of the variables in this p-value file has to match 
#         the order of the variables in the correlation matrix file!
#      See file "Example_pvals" for an example file for nvar=12 and nsnp=100.
#
###############################################################################################

# set working directory 
# i.e., that contains the file with the correlation matrix 
# and the file with the p-values
setwd("/Users/ngillespie/Documents/ngillespie/papers/report_CA_FTND_gwas/2b.brain_volume/TATES/")  # work

# chose name for the results file
namefile=as.character("GWAS_Ramyg_d48_ftnd_tates_RESULTS.txt")            # name simulated data file

               
# ---------------------------------------------------------------------------------------------
# ----------------------------- READ IN p-values and correlations -----------------------------
# ---------------------------------------------------------------------------------------------
nvar=3  # nr of variables 
nsnp=2541501  # nr of tested SNPs

getwd()
cormat<-read.table("GWAS_Ramyg_d48_ftnd_tates_short_0.30_correlations.txt",header=F,sep = "\t") # nvar*nvar
pval<-read.table("GWAS_Ramyg_d48_ftnd_tates_short_0.30.txt",header=F)    # nsnp*(nvar+2)

r=as.matrix(cormat)
pval2=as.matrix(pval[,3:dim(pval)[2]]) # exclude CHR number and rs-number columns

# ---------------------------------------------------------------------------------------------
# ------------------------------------ run TATES procedure ------------------------------------
# ---------------------------------------------------------------------------------------------

# read in the beta weights from the 6th order polynomial
betaa=c(-0.0007609278,-0.0022965148,0.6226249243,0.0148755138,
         0.1095155903,-0.0218930325,0.2178970393)

# define the results matrix containing PT,
# i.e., the TATES trait-based p-value
synp=matrix(0,nsnp,1)

# ------------------- start TATES loop

for (isnp in 1:nsnp) {

# 1. ---------------- get the p-values 

ps=pval2[isnp,]
tmp=sort(ps,index.return=T)  # gives location of the sorted p-values
pj=tmp$x # sorted p-values
iorder=tmp$ix # index

# 2.----------------- get the correlation matrix between the variables 

r2=matrix(0,nvar,nvar)
# order correlation matrix according to the rankorder of the p-values
r2=r[iorder,iorder]

# 3. ---------------- weight symptom correlations [r2] with regression weights 
# matrix ro contains the predicted correlations between p-values,
# i.e., predicted from the correlations between variables
ro=diag(nvar)
for (i1 in 1:nvar)  {  
for (i2 in 1:i1) {
if (i1>i2) {
er=r2[i1,i2]
ro[i1,i2]=ro[i2,i1]= betaa[7]*er^6+betaa[6]*er^5+betaa[5]*er^4+betaa[4]*er^3+betaa[3]*er^2+betaa[2]*er+betaa[1]
}}}

# 4. ---------------- determine eigen values based on entire p-value matrix 
# get Mall (Me based on all p-values)
alllam=eigen(ro[1:nvar,1:nvar])$values #eigenvalues of the ro matrix
mepj=nvar
for (i1 in 1:nvar) {
mepj=mepj-(alllam[i1]>1)*(alllam[i1]-1) }

# 5. ---------------- determine eigen values top x p-values 
# sellam is eigenvalues based on varying nr of top SNPs

mej=matrix(c(seq(1,nvar,1)),nvar,1,byrow=T)

for (j in 1:nvar) { 
sellam=eigen(ro[1:j,1:j])$values #eigenvalues of the ro matrix
id=j
# subtract 1-eigenvalue for those eigenvalues >1    # page 284 Li et al., 2011
for (i1 in 1:id) {
mej[j,1]=mej[j,1]-(sellam[i1]>1)*(sellam[i1]-1)
}
}

# 6. ---------------- weight sorted p-values with eigenvalues ratio 

pg=matrix(0,nvar,1)
for (i in 1:nvar) {
pg[i,1]=(mepj/mej[i,1])*pj[i]
}

pg=pg[iorder]  # p-values back in original order so that pval[1] corresponds to item[1] etc !!


# 7. ----------------------  write out results 
synp[isnp,]=min(pg)
                                       
}

# ---------------------------------------------------------------------------------------------
# ------------------------------------ create results file ------------------------------------
# ---------------------------------------------------------------------------------------------
# The results file is similar to the p-value file that was read in, 
# but it has a header and 1 column is added at the end of this file, 
# which contains the TATES trait-based p-values.
# This way, the variable-specific p-values are readily at hand 
# so that one can examine "the origin/cause" of low TATES values 

res=cbind(pval,synp)

x=matrix(0,nvar,1)
for (i in 1:nvar) {
x[i,1]=c(paste('VAR',as.character(i),sep=""))
}



colname=c("CHR","SNP",as.vector(x),"TATES")
write(colname,file=namefile,nc=(nvar+3))
write(t(res),file="tmp",nc=(nvar+3))
file.append(namefile,"tmp")
file.remove("tmp")

