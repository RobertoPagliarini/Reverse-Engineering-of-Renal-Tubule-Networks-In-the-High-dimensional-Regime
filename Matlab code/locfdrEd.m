function [ffdr f0 f x]=locfdrEd(zz,bre,pct,pct0,disfit)
%zz: A vector of summary statistics, one for each case under
%          simultaneous consideration. In a microarray experiment there
%          would be one component of zz for each gene, perhaps a
%          t-statistic comparing gene expression levels under two
%          different conditions. The calculations assume a large number
%          of cases, say at least length(zz) > 100.



%bre: Number of breaks in the discretization of the z-score axis,
%          set to 120 by default. This can also be a vector of
%          breakpoints fully describing the discretization.



%pct:  Excluded tail proportions of zz's when fitting f(z);
%          pct=1/1000 by default; pct=0 includes full range of zz's; pct
%          can also be a 2-vector, describing the fitting range.

%pct0: 2-length vector. Included proportion of range (from the maximum of the distribution)
%      used in fitting null density f0(z); [percentile to the left, percentile to the right]
%      default value is [2/3 2/3].

%disfit: distribution to use while fitting the histogram; default is
%poisson. Other can be 'gamma'

if ~exist('bre') || isempty(bre)
    bre=120;
end
if ~exist('pct0') || isempty(pct0)
    pct0=[2/3 2/3];
end
if ~exist('pct') || isempty(pct)
    pct=1/1000;
end
if ~exist('disfit') || isempty(disfit)
    disfit='poisson';
end
%begin excluding proportion of data
v=prctile(zz,[100*pct 100*(1-pct)]);
lo=v(1);
hi=v(2);
zzz=zz;
zzz(find(zzz < lo))=lo;
zzz(find(zzz > hi))=hi;
%end excluding proportion of data;


%begin the histogram
breaks=linspace(lo,hi,bre);
[y,x]=hist(zzz,breaks);
%y = y./sum(y);
%x=(x(2:end)+x(1:end-1))./2; %promediando las particiones para hacer glmfit
%end the histogram

%begin suavidad en los extremos
if pct>0
    y(1)=min(y(1),1); y(end)=min(y(end),1);
end
%end suavidad en los extremos

%begin the mixture density and the generalized linear model
%y=y/sum(y); %in units of probability
%the knots over the suitable percentiles of the data
nknots=10; 
%knots=linspace(0,100,nknots);  knots=prctile(zzz,knots);
knots=linspace(lo,hi,nknots);
%the basis spline
x_sp=spcol(knots,3,x);
%f=glmfit(x_sp,y,'poisson','log');
if strcmp(disfit,'gamma');y(y<1)=1;end % for gamma fitting y cannot be 0
f=glmfit(x_sp,y,disfit,'log');
f=glmval(f,x_sp,'log');

l=log(f);

%end the mixture density and the generalized linear model
%xmax=x(l==max(l));
imax=find(l==max(l)); xmax=x(imax);
%%%%Eduardo;
if length(xmax)>1; xmax=xmax(1);end
%%%%

%begin the zero neigborhood
if xmax<lo; %taking into account one-sided distributions
   lo0=lo;
else
%     lo0=prctile(zz(find(zz<xmax)),100*(1-pct0)); % from the data
     lo0=prctile(zz(find(zz<xmax)),100*(1-pct0(1))); % from the data,
%      kk=cumsum(f(1:imax));kk=kk./kk(end); %for defining percentil from f instead from the data
%      [pp,ipp]=min(abs(kk-(1-pct0))); lo0=x(ipp);    % from the modeling f of the data
end
%hi0=prctile(zz(find(zz>xmax)),100*(pct0)); %from the data
hi0=prctile(zz(find(zz>xmax)),100*(pct0(2))); %from the data
% kk=cumsum(f(imax:end));kk=kk./kk(end); %for defining percentil from f instead from the data
% [pp,ipp]=min(abs(kk-pct0)); hi0=x(imax+ipp-1);    % from the modeling f of the data
nx=length(x);
%los indices de los bins que estan dentro de la condicion de vecindad zero
i0=1:nx;
i0=i0(find(x > lo0 & x < hi0));
%los bins que estan dentro de la vecindad
x0=x(i0);
y0=l(i0);
%end the zero neighborhood

%begin linear estimatiation of log(zz) dentro de la vecindad usando un poly
%de 2nd order
X0=[(x0-xmax)' (x0-xmax).^2',repmat(1,size(x0'))];
coeff=regress(y0,X0);
%coeff1=robustfit(X0(:,1:2),y0); %is this better?
%end linear estimation of log(zz) dentro de la vecindad


%begin the empirical fdr evaluation in all the range
X00=[(x-xmax)', (x-xmax).^2',repmat(1,size(x'))]; %base cuadrada
%estimation of the log of fo(0)
l0=X00*coeff;
%evaluation of the exponential
f0=exp(l0);
%the empirical fdr
fdr=f0./f;
fdr(find(fdr>1))=1;
%end the empirical fdr evaluation in all the range

%begin interpolation
% este tenia nan
% ffdr= interp1(x,fdr,zz);
% este lo hace mejor
ffdr= spline(x,fdr,zz);

%disp('End')