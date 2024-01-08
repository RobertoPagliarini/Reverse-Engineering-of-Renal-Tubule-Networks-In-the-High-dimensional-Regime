function [Matrices] = ORAShrinkagePartialCorr_V5(data,ShrCov,VarInModel)
%Computing Adjacency Weighted Matrices for networks by using the ORA-Quantile procedure proposed in Pagliarini et al. 
%and other threshold approached described in the same paper. 
%
%Input:
%
%data: a n*p matrix, where n is the number of samples and p the number of
%variables.
%ShrCov: name of the method for the shrinkage estimator of the covariance
%matrix. If ShrCov = 'Sha' the the method is the one proposed in J. Schaefer and K. Strimmer.  2005,
%otherwise the method is the one proposed in Chen et al. 2010.
%VarInModel: name of the variables/nodes. This input can be empty.
%
%Output: a struct array that contains these fields:
%
%GGMatrixCorrThre: the matrix associated with the output of the ORA-Quantile Algorithm;
%CorrMatrixThreLocalFDR: the correlation matrix obtained by using Local-FDR threshold;
%GGMatrixThreLocalFDR: the matrix associated with the Graphical Gaussian Model obtained by using Local-FDR threshold;
%GGMatrix: the matrix associated with the Graphical Gaussian Model;
%CorrMatrix: the correlation matrix;
%GGMatrixThre: 
%CorrMatrixThre: threshold correlation matrix based on pval t-test;
%GGMatrixThre: graphical Gaussian model matrix obtained by using the Quantile thresold;
%GGMatrix_zscore: Graphical Gaussian Model Matrix z-score threshold. Values with absolute value of z-score < 1.96 are set to 0;
%GGMatrix_2sigma: Graphical Gaussian Model Matrix 2-sigma threshold. Values with absolute value < mu+2*sigma are set to 0, 
%where mu is the mean and sigm is the standard deviation;
%PvalCorrMatrix: matrix of p-values associated with correlation matrix computed by applying t-test;
%VarInModel: name of the variables/nodes, if this input is not empty.

if nargin < 3 
    
    VarInModel = [];
    
else
    
    VarInModel = {VarInModel};
    
end

%Compute covariance matrix

if strcmp(ShrCov,'Sha')
    
    %Shrinkage estimate of a covariance matrix, using optimal shrinkage
    %coefficient proposed in J. Schaefer and K. Strimmer.  2005
    [CovMat,lambda] = covshrinkKPM(data,1);
    
else
    
    %Shrinkage estimate of a covariance matrix, using optimal shrinkage
    %coefficient proposed in Chen et al., IEEE Trans. on Sign. Proc., Volume 58, Issue 10, October 2010.
    [CovMat,lambda] = shrinkage_cov(data);
    
end

[mx,nx] = size(CovMat);

%Compute pseudo inverse of covariance matrix
PrecMat = pinv(CovMat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Computing the correlation matrix
CorrMatrix =  diag(1./sqrt(diag(CovMat)))*CovMat*diag(1./sqrt(diag(CovMat)));

%Pvalues of Corr Matrix based on t-test
[nd, md] = size(data);

t = CorrMatrix.*sqrt((nd-2)./(1-CorrMatrix.^2));

PvalCorrMatrix = 2*tcdf(-abs(t),nd-2);

%Computing the Graphical Gaussian model matrix
GGMatrix =  diag(1./sqrt(diag(PrecMat)))*PrecMat*diag(1./sqrt(diag(PrecMat)));

%Creating local FDR threshold matrices for correlaitons and graphical model
GGMatrixThreLocalFDR = sparse(zeros(nx,nx));

CorrMatrixThreLocalFDR = sparse(zeros(nx,nx));

%Creating vecotrs associated with local FDR matrices
VecGGMatrix = zeros((nx^2-nx)/2,1);
    
IndexVecGGMatrix = zeros((nx^2-nx)/2,2);

VecCorrMatrix = zeros((nx^2-nx)/2,1);
    
IndexVecCorrMatrix = zeros((nx^2-nx)/2,2);
     
tmp = 0;
 
for i = 1:nx

    for j = i+1:nx

        tmp = tmp+1;
            
        VecGGMatrix(tmp,1) = GGMatrix(i,j); 
        
        VecCorrMatrix(tmp,1) = CorrMatrix(i,j); 
        
        IndexVecGGMatrix(tmp,1) = i;
            
        IndexVecGGMatrix(tmp,2) = j;
               
    end
    
end

%%%%%Quantile threshold
quantGG05 = quantile(VecGGMatrix,[0.05 0.95]);

quantCor05 = quantile(VecCorrMatrix,[0.05 0.95]);

A = GGMatrix <= quantGG05(1);

B = GGMatrix >= quantGG05(2);

C = CorrMatrix <= quantCor05(1);

D = CorrMatrix >= quantCor05(2);

E = (A|B)&(C|D);
 
%Computing the Graphical Gaussian model matrix quantile thresold
GGMatrixThre =sparse(E.*GGMatrix);

%Computing the Graphical Gaussian model matrix quantile thresold plus
%correlation weight, the output of our ORA-quantile algorithm
GGMatrixCorrThre = sparse(E.*CorrMatrix);
    
%Threshold correlation matrix based on pval t-test
CorrMatrixThre = sparse(CorrMatrix.*(PvalCorrMatrix <= 0.05));
     
%%%%Local Fdr
[fdr] = locfdrEd(VecGGMatrix);

%Computile Graphical Gaussian outputs based on local-FDR by using 0.2 as threshold 
for i = 1:length(fdr)
        
    if fdr(i) < 0.2
            
        GGMatrixThreLocalFDR(IndexVecGGMatrix(i,1),IndexVecGGMatrix(i,2)) = VecGGMatrix(i);
            
        GGMatrixThreLocalFDR(IndexVecGGMatrix(i,2),IndexVecGGMatrix(i,1)) = VecGGMatrix(i);
        
        CorrGGMatrixThreLocalFDR(IndexVecGGMatrix(i,1),IndexVecGGMatrix(i,2)) = VecCorrMatrix(i);
            
        CorrGGMatrixThreLocalFDR(IndexVecGGMatrix(i,2),IndexVecGGMatrix(i,1)) = VecCorrMatrix(i);
            
       
    end
    
end

%Graphical Gaussian Model Matrix z-score threshold
[nx,mx] = size(GGMatrix);

VecGGMatrix = [];

VecGGMatrix = reshape(GGMatrix,nx*nx,1);

[z,mu,sigma] = zscore(VecGGMatrix);

GGMatrix_z = zeros(nx,nx);

GGMatrix_2sigma = zeros(nx,nx);

%Graphical Gaussian Model Matrix 2-sigma threshold
AM = GGMatrix;

AM(abs(AM) < 1.96) = 0;

GGMatrix_z = AM; 
 
AM = GGMatrix;

AM(abs(AM) < mu+2*sigma) = 0;

GGMatrix_2sigma = AM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Matrices = struct('GGMatrixCorrThre',GGMatrixCorrThre,'CorrMatrixThreLocalFDR',CorrMatrixThreLocalFDR,'GGMatrixThreLocalFDR',GGMatrixThreLocalFDR,...
                  'GGMatrix',GGMatrix,'CorrMatrix',CorrMatrix,'GGMatrixThre',GGMatrixThre,'CorrMatrixThre',CorrMatrixThre,...
                  'GGMatrix_zscore',GGMatrix_z,'GGMatrix_2sigma',GGMatrix_2sigma,'PvalCorrMatrix',PvalCorrMatrix,'VarInModel',VarInModel);


