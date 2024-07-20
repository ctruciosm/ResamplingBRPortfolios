function [diff, delta, sigi, pvalpw, pvalboot] = robustVariance_onesided(data,boot,M,bl,kernel,extsim,bootMat)
% [diff, delta, sigi, pvalpw, pvalboot] = ...
% robustVariance(data,boot,M,bl,kernel,extsim,bootMat)
%  minimal input:    robustVariance(data) <- takes some minutes to compute
% Tests whether the difference of two variances is statistically significantly 
% different from zero based on Ledoit & Wolf (2010). Tested as two-sided 
% hypothesis by simulating M datasets with the circular block bootstrap, 
% taking critical values as the empirical quantiles of the simulated datasets.
% The Quadratic Spectral Kernel is used as default in the prewhitened kernel variance 
% estimator. The R implementation uses the Gallant kernel.
%
% Inputs:
%   data        [Tx2] matrix of excess returns (not in log)
%   boot        1 if you want to carry out bootstrap inference (default). 
%               0 if you merely want the prewhitened HAC results (faster).
%               boot inference is more flexible but computationally costly.
%   M           number of bootstrap iterations to compute critical values
%   bl          block size in Circular Block Bootstrap.
%               Use routine optimalblrobustVariance.m to determine optimal
%               block size.
%   kernel      'G' for Gallant/Parzen, 'QS' for Quadratic spectral (default)
%   extsim      1 if the indices matrix bootMat in the circluar block bootstrap is fed
%               in rather than simulated in robustVariance itself, 0 else
%               useful to achieve comparability of results based on other implementations
%   bootMat     if extims==1: exogenous indices matrix in circular block bootstrap 
%               of size [MxT] where M is number of CBB iterations
%
% Outputs:
%   diff        Difference of sample variances
%   delta       Difference of log sample variances
%   sigi        Robust prewhitened standard errors of the two sample variances
%   pvalpw      prewhitened kernel p-value of the difference of log variances
%   pvalboot    bootstrap p-value of the difference of log variances
%
% c2010 Dan Wunderli, Institute for Empirical Research in Economics, U Zurich

%Defaults
if not(ismember('stud',who)), stud=1; end
if not(ismember('boot',who)), boot=1; end
if not(ismember('M',who)), M=5000; end
if not(ismember('bl',who)), fprintf('%s \n','Computation of optimal block size:'), bl=optimalblrobustVariance(data,0,1000,199); fprintf('%s \n','Testing of H0:'), end
if not(ismember('extsim',who)), extsim=0; end
if not(ismember('bootMat',who)), bootMat=0; end
if not(ismember('kernel',who)), kernel='QS'; end
H0=0;

%Start
t=size(data,1);

%Computation of studentized test statistic and Generation of 
%Circular Block Bootstrap Index Matrix
X=data;

diff = var(X(:,1)) - var(X(:,2));
delta = log(var(X(:,1))) - log(var(X(:,2)));

%%%needed input arguments for andmon implementation of Mike Cliff, V1.1
gmmopt.prt             = 0;
gmmopt.aminfo.p        = 1;
gmmopt.aminfo.q        = 0;
gmmopt.aminfo.vardum   = 0;
gmmopt.aminfo.kernel   = kernel;
gmmopt.aminfo.nowhite  = 0;
gmmopt.aminfo.diagdum  = 0;
gmmopt.plot            = 0;

nu = [mean(X(:,1)); mean(X(:,2)); mean(X(:,1).^2); mean(X(:,2).^2)];
nux = [X(:,1)-nu(1) X(:,2)-nu(2) X(:,1).^2-nu(3) X(:,2).^2-nu(4)];
nabla = [-(2*nu(1))/(nu(3)-nu(1)^2) (2*nu(2))/(nu(4)-nu(2)^2) 1/(nu(3)-nu(1)^2) -1/(nu(4)-nu(2)^2)]';

%%%Prewhiten data with VAR(1) model, estimate HAC kernel estimator
%%%using AR(1) models as univariate approximating parametric models
dataVAR = vare2(nux,1,0);          %alternatively by vare(nux,1)
PSI = horzcat(dataVAR(1:4).beta)';
[U,S,V]=svd(PSI);
for i=1:size(PSI,2)
    if S(i,i)>0.97
        S(i,i)=0.97;
    elseif S(i,i)<-0.97
        S(i,i)=-0.97;
    end
end
PSI = U*S*V';

u = (nux(1+1:end,:)' - PSI*nux(1:end-1,:)')'; %equals V.star in R implementation

[~,Z] = andmon6cvm2(gmmopt,u,PSI);    %covariance matrix by kernel-based HAC estimator of Andrews-Monahan (1992). andmon6 is an altered version of Mike Cliff's andmon function Version 1.1
varde = (nabla'*((eye(size(u,2)) - PSI) \Z/ (eye(size(u,2)) - PSI)')*nabla);

z_T = (delta - H0)./(sqrt(varde/t));     %studentization of data
sigi = sqrt(varde)./sqrt(t);                %HAC std estimate
r = mod(t,bl); L = floor(t/bl);
if boot==1
    Xbootind = zeros(M,t);
    
    if extsim==0
        for m = 1:M                         %generation of M cbb matrices X_T*m, 1<=m<=M
            [ind,L,r] = cbb_seq(t,bl);
            Xbootind(m,:) = ind;
        end
        csvwrite('CBBMatMatlab.csv',Xbootind);
    elseif extsim==1
        Xbootind=bootMat;
    else error('only extsim=1 or extsim=0 accepted')
    end
    
    %bsstat is a matrix with corresponding studentized test statistics for each bootstrap iteration (row)
    %muboot as optional output contains CBB simulated excess returns of two assets,
    %sigboot as optional output is bootstrap std estimate of difference of two Sharpe ratios
    [bsstat] = bsstats11_onesided(X,Xbootind,delta,bl,L,r,t);
    
    pvalboot=(sum(bsstat>= z_T)+1)/(M+1);
end

pvalpw = 2*normcdf(-(delta-H0)/sigi,0,1);

fprintf('%s \n',horzcat('Difference in variances = ',num2str(diff)))
fprintf('%s \n',horzcat('Difference in log variances = ',num2str(delta)))
fprintf('%s \n',horzcat('HAC pw standard errors of log difference = ',num2str(sigi)))
fprintf('%s \n',horzcat('HAC pw p-value of log difference = ',num2str(pvalpw)))
if boot==1
    fprintf('%s \n',horzcat('Bootstrap p-value = ',num2str(pvalboot)))
end

