# Networks-for-renal-tubule-segments

A Matlab Tool for computing adjacency weighted matrices for networks by using the ORA-Quantile procedure proposed in Pagliarini et al. 
and other threshold approaches described in the same paper. 

function [Matrices] = ORAShrinkagePartialCorr_V5(data,ShrCov,VarInModel)

Input:

data: a n*p matrix, where n is the number of samples and p the number of
variables.
ShrCov: name of the method for the shrinkage estimator of the covariance
matrix. If ShrCov = 'Sha' the the method is the one proposed in J. Schaefer and K. Strimmer.  2005,
otherwise the method is the one proposed in Chen et al. 2010.
VarInModel: name of the variables/nodes. This input can be empty.

Output: a struct array that contains these fields:

GGMatrixCorrThre: the matrix associated with the output of the ORA-Quantile Algorithm;
CorrMatrixThreLocalFDR: the correlation matrix obtained by using Local-FDR threshold;
GGMatrixThreLocalFDR: the matrix associated with the Graphical Gaussian Model obtained by using Local-FDR threshold;
GGMatrix: the matrix associated with the Graphical Gaussian Model;
CorrMatrix: the correlation matrix;
GGMatrixThre: 
%CorrMatrixThre: threshold correlation matrix based on pval t-test;
GGMatrixThre: graphical Gaussian model matrix obtained by using the Quantile thresold;
GGMatrix_zscore: Graphical Gaussian Model Matrix z-score threshold. Values with absolute value of z-score < 1.96 are set to 0;
GGMatrix_2sigma: Graphical Gaussian Model Matrix 2-sigma threshold. Values with absolute value < mu+2*sigma are set to 0, 
where mu is the mean and sigm is the standard deviation;
PvalCorrMatrix: matrix of p-values associated with correlation matrix computed by applying t-test;
VarInModel: name of the variables/nodes, if this input is not empty.
