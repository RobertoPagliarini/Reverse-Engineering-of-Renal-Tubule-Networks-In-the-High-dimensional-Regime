function [Matrices] = ORASchafferShrinkageNetworkInference(data,OraSch,LocalFdr,VarInModel)

Computing Adjacency Weighted Matrices for networks by using the ORA-Quantile procedure proposed in Pagliarini et al., 
shrinkage covariace matrix computation plus Local-FDR (Pagliarini et al.), 
and quantile correlation and partical correlazion plus network deconvolution. 

Input:

data: a n*p matrix, where n is the number of samples and p the number of
variables.

OraSch: name of the method for the shrinkage estimator of the covariance
matrix. If OraSch = 'Sha' the the method is the one proposed in J. Schaefer and K. Strimmer.  2005,
OraSch = 'Ora' the method is the one proposed in Chen et al. 2010.

LocalFdr: name of the method for network pruning. If LocalFdr = 'Q'the
approach is correlation/partial correlation quantile test as described in
Pagliarini et al. If LocalFdr = 'L'the approach is the one described in J.
Schäfer, K. Strimmer 2005. If LocalFdr = 'N'the approach is correlation/partial correlation quantile test as described in
Pagliarini et al. plus network deconvolution developed in Feizi et al.

VarInModel: name of the variables/nodes. This input can be empty.

Output: a struct array that contains these fields:

GGMatrix: the matrix associated with the Graphical Gaussian Model;
GGMatrixThre: the pruned GGMatrix;
GGMatrixCorrThre: the matrix associated with the output of the ORA-Quantile Algorithm, or with the select pruning method;
CorrMatrix: the correlation matrix;
OraSch: the input OraSch parameter;
LocalFdr: the input LocalFdr parameter;
VarInModel: name of the variables/nodes, if this input is not empty.

