function [Matrices] = ORASchafferShrinkageNetworkInference(data,OraSch,LocalFdr,VarInModel)
%Computing Adjacency Weighted Matrices for networks by using the ORA-Quantile procedure proposed in Pagliarini et al., 
%shrinkage covariace matrix computation plus Local-FDR (Pagliarini et al.), 
% and quantile correlation and partical correlazion plus network deconvolution. 
%
%Input:
%
%data: a n*p matrix, where n is the number of samples and p the number of
%variables.
%OraSch: name of the method for the shrinkage estimator of the covariance
%matrix. If OraSch = 'Sha' the the method is the one proposed in J. Schaefer and K. Strimmer.  2005,
%OraSch = 'Ora' the method is the one proposed in Chen et al. 2010.
%LocalFdr: name of the method for network pruning. If LocalFdr = 'Q'the
%approach is correlation/partial correlation quantile test as described in
%Pagliarini et al. If LocalFdr = 'L'the approach is the one described in J.
%Schäfer, K. Strimmer 2005. If LocalFdr = 'N'the approach is correlation/partial correlation quantile test as described in
%Pagliarini et al. plus network deconvolution developed in Feizi et al.
%VarInModel: name of the variables/nodes. This input can be empty.
%
%Output: a struct array that contains these fields:
%
%GGMatrix: the matrix associated with the Graphical Gaussian Model;
%GGMatrixThre: the pruned GGMatrix;
%GGMatrixCorrThre: the matrix associated with the output of the ORA-Quantile Algorithm, or with the select pruning method;
%CorrMatrix: the correlation matrix;
%OraSch: the input OraSch parameter;
%LocalFdr: the input LocalFdr parameter;
%VarInModel: name of the variables/nodes, if this input is not empty.

if nargin < 4 
    
    VarInModel = [];
    
else
    
    VarInModel = {VarInModel};
    
end

%Compute covariance matrix
if strcmp(OraSch,'Ora')

    %Shrinkage estimate of a covariance matrix, using optimal shrinkage
    %coefficient proposed in Chen et al., IEEE Trans. on Sign. Proc., Volume 58, Issue 10, October 2010.
    [CovMat] = shrinkage_cov(data);
    
else
    
    if strcmp(OraSch,'Sch')
        
        %Shrinkage estimate of a covariance matrix, using optimal shrinkage
        %coefficient proposed in J. Schaefer and K. Strimmer.  2005
        [CovMat] = covshrinkKPM(data,1);
        
    end
    
end

%Compute dimension of covariance matrix
[mx,nx] = size(CovMat);

%Compute pseudo inverse of covariance matrix
PrecMat = pinv(CovMat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Computing the correlation matrix
CorrMatrix =  sparse(diag(1./sqrt(diag(CovMat)))*CovMat*diag(1./sqrt(diag(CovMat))));

%Pvalues of Corr Matrix
%[nd, md] = size(data);

%Computing the GG matrix
GGMatrix =  sparse(diag(1./sqrt(diag(PrecMat)))*PrecMat*diag(1./sqrt(diag(PrecMat))));

% CorrMatrixThre = sparse(zeros(nx,nx));

%%%%%%%Removing zero from diagonal
V = zeros(nx,1);

CorrMatrix = CorrMatrix - diag(diag(CorrMatrix)) + diag(V);

GGMatrix = GGMatrix - diag(diag(GGMatrix)) + diag(V);

%Costruction of vectors for filtering by applying quantile test
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

%Correlation and partial correlation quantile test as described in
%Pagliarini et al
if strcmp(LocalFdr,'Q')

    quantGG05 = quantile(VecGGMatrix,[0.05 0.95]);

    quantCor05 = quantile(VecCorrMatrix,[0.05 0.95]);

    A = GGMatrix <= quantGG05(1);

    B = GGMatrix >= quantGG05(2);

    C = CorrMatrix <= quantCor05(1);

    D = CorrMatrix >= quantCor05(2);

    E = (A|B)&(C|D);
    
    GGMatrixThre =sparse(E.*GGMatrix);
    
    GGMatrixCorrThre = sparse(E.*CorrMatrix);
   
end

if strcmp(LocalFdr,'L')
    
    %Creating threshold matrices
    GGMatrixThre = sparse(zeros(nx,nx));

    GGMatrixCorrThre = sparse(zeros(nx,nx));
     
    %%%%Local Fdr
    [fdr] = locfdrEd(VecGGMatrix);

    for i = 1:length(fdr)
        
        if fdr(i) < 0.2
            
            GGMatrixThre(IndexVecGGMatrix(i,1),IndexVecGGMatrix(i,2)) = VecGGMatrix(i);
            
            GGMatrixThre(IndexVecGGMatrix(i,2),IndexVecGGMatrix(i,1)) = VecGGMatrix(i);
        
            GGMatrixCorrThre(IndexVecGGMatrix(i,1),IndexVecGGMatrix(i,2)) = VecCorrMatrix(i);
            
            GGMatrixCorrThre(IndexVecGGMatrix(i,2),IndexVecGGMatrix(i,1)) = VecCorrMatrix(i);
            
        end
    
    end
    
end

%Correlation and partial correlation quantile test as described in
%Pagliarini et al plus network deconvolution described in paper by paper by
%S. Feizi et al.
if strcmp(LocalFdr,'N')

    quantGG05 = quantile(VecGGMatrix,[0.05 0.95]);

    quantCor05 = quantile(VecCorrMatrix,[0.05 0.95]);

    A = GGMatrix <= quantGG05(1);

    B = GGMatrix >= quantGG05(2);

    C = CorrMatrix <= quantCor05(1);

    D = CorrMatrix >= quantCor05(2);

    E = (A|B)&(C|D);
    
    GGMatrixThre =sparse(E.*GGMatrix);
    
    GGMatrixCorrThre = sparse(E.*CorrMatrix);
    
    mat_nd=ND(abs(full(GGMatrixCorrThre)));
    
    VecND = (reshape(mat_nd,nx*nx,1));
    
    quantND = quantile(VecND,0.95);
    
    E = mat_nd >= quantND;
    
    GGMatrixCorrThre = sparse(E.*CorrMatrix);
        
end
    
            
    
Matrices = struct('GGMatrix',GGMatrix,'GGMatrixThre',GGMatrixThre,'GGMatrixCorrThre',GGMatrixCorrThre,'CorrMatrix',CorrMatrix,...
                  'OraSch',OraSch,'LocalFdr',LocalFdr,'VarInModel',VarInModel);


