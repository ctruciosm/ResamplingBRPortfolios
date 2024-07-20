function [optbl blcand] = optimalblrobustVariance(data,silent,M,MBM,kernel,blcand,avbl,alph,Tstart,extsimSB,extSB,extsimMB,extMB)
% [optbl blcand] = optimalblrobustVariance(data,silent,M,MBM,kernel,blcand,avbl,alph,Tstart,extsimSB,extSB,extsimMB,extMB)
% minimal input is optimalblrobustVariance(data)
% In robustVariance, confidence regions are simulated. To find the optimal block size,
% one first simulates datasets by a parametric VAR(1) bootstrap with stationary bootstrapping
% the residuals. Then one constructs confidence regions by circular block bootstrap
% constructed critical values and observes how often rejections of the null hypothesis occur.
% The optimal block size is the one that minimizes the deviation of the simulated
% coverage from the desired coverage 1-alpha i.e. minimizes the deviation of the
% desired error rate alpha from the simulated error rate errr.
%
% Inputs:
%   data        [Tx2] matrix of excess returns (NOT in log already)
%   silent      1 for yes, 0 for no. Default is not silent, so that
%               estimated time left in seconds and current error rates are printed
%   M           number of first order bootstrap iterations (stationary bootstrap)
%   MBM         number of second order bootstrap iterations (circular blocks bootstrap)
%   kernel      'QS' for quadratic spectral kernel (default), 'G' for Gallant kernel
%   blcand      candidate block sizes, e.g. [1 3 6 10]
%   avbl        average blocksize in first order bootstrap (stationary bootstrap)
%   alph        upper bound for type I error, used in testing the hypothesis
%   Tstart      in the parametric bootstrap, e.g. Tstart=1000 to make starting values irrelevant
%   extsimSB    1 if the indices matrix extSB in the stationary bootstrap is fed
%               in rather than simulated in optimalblrobustVariance itself, 0 else
%   extSB       exogenous indices matrix in stationary bootstrap of size
%               [MxT] if extsimSB==1
%   extsimMB    1 if the indices matrix extSB in the circular blocks bootstrap is fed
%               in rather than simulated in optimalrobustVariance itself, 0 else
%   extMB       exogenous indices matrix in circular block bootstrap of
%               size [MBMxT] if extsimMB==1
%
% Outputs:
% (during execution, estimated time left in seconds and current error rates are printed)
%   optbl:  the optimal blocksize so that the simulated coverage is closest to the desired one
%   blcand: blcand vector that was put into the algorithm
%
% ?2010 Dan Wunderli, Institute for Empirical Research in Economics, U Zurich

if not(ismember('silent',who)), silent=0; end
if not(ismember('MBM',who)), MBM=199; end
if not(ismember('M',who)), M=1000; end
if not(ismember('alph',who)), alph=0.05; end
if not(ismember('avbl',who)), avbl=5; end
if not(ismember('blcand',who)), blcand=[1 3 6 10]; end
if not(ismember('kernel',who)), kernel='QS'; end
if not(ismember('Tstart',who)), Tstart=50; end
if not(ismember('extsimSB',who)), extsimSB=0; end
if not(ismember('extSB',who)), extSB=0; end
if not(ismember('extsimMB',who)), extsimMB=0; end
if not(ismember('extMB',who)), extMB=0; end

varlag=1;
t = size(data,1)-varlag; T = size(data,1); S = size(data,2);
% datac = data - repmat(mean(data,1),[T 1]);

theta_0star = log(var(data(:,1))) - log(var(data(:,2)));

Pt = vare2(data,1,1);           %varsvdfit numerically more advisable here
                                %vare2 used for compatability reasons
yhat=horzcat(Pt(1:S).yhat);
A=horzcat(Pt(1:S).beta)';
% A = A(:,1:end-1);             %cuts off estimate of constant

resid = data(1+varlag:end,:) - yhat;
residc = resid;
% residc = resid - repmat(mean(resid,1),[t 1]);   %!!!better numerical
% centering?

%% needed input arguments for andmon implementation of Mike Cliff, V1.1
gmmopt.prt             = 0;
gmmopt.aminfo.p        = 1;
gmmopt.aminfo.q        = 0;
gmmopt.aminfo.vardum   = 0;
gmmopt.aminfo.kernel   = kernel;
gmmopt.aminfo.nowhite  = 0;
gmmopt.aminfo.diagdum  = 0;
gmmopt.plot            = 0;

%%%%%%%%%%% Simulation of data sets by Stationary Bootstrap
Out = zeros(1,length(blcand));
tic
for m=1:M       %%% first level bootstrap
    timeb=toc;
    if extsimSB==0
        sbb_seq = sbbseq(t,Tstart,avbl);        % stat. BB sequence of length T+Tstart
        % indmatSB(m,:) = sbb_seq;
    else sbb_seq=extSB(m,:);
    end
    residstar = residc(sbb_seq,:);
    
    xold = zeros(S,1); xboot = zeros(t+Tstart,S);
    
    for j=1:(t+Tstart)
        xnew = A(:,3)+ A(:,[1 2])*xold + residstar(j,:)';
        xboot(j,:) = xnew';
        xold = xnew;
    end
    xboot = xboot(Tstart:t+Tstart,:);           % throw away first Tstart-1 valued
        
    X = xboot;
    mubi = log(var(X(:,1))) - log(var(X(:,2)));
    
    nu = [mean(X(:,1)); mean(X(:,2)); mean(X(:,1).^2); mean(X(:,2).^2)];
    nux = [X(:,1)-nu(1) X(:,2)-nu(2) X(:,1).^2-nu(3) X(:,2).^2-nu(4)];
    nabla = [-(2*nu(1))/(nu(3)-nu(1)^2) (2*nu(2))/(nu(4)-nu(2)^2) 1/(nu(3)-nu(1)^2) -1/(nu(4)-nu(2)^2)]';
    
    %%% Prewhiten data with VAR(1) model, estimate HAC kernel estimator
    %%% using AR(1) models as univariate approximating parametric models
    dataVAR = vare2(nux,1,0);          % alternatively by vare(nux,1)
    THETA = horzcat(dataVAR(1:4).beta)';
    [U,Q,V]=svd(THETA);
    for i=1:size(THETA,2)
        if Q(i,i)>0.97
            Q(i,i)=0.97;
        elseif Q(i,i)<-0.97
            Q(i,i)=-0.97;
        end
    end
    THETA = U*Q*V';
    
    
    u = (nux(1+1:end,:)' - THETA*nux(1:end-1,:)')'; % equals V.star in R implementation
    
    [~,Z] = andmon6cvm2(gmmopt,u,THETA);    % covariance matrix by kernel-based HAC estimator of Andrews-Monahan (1992). 
                                            % andmon6cvm2 is an altered version of Mike Cliff's andmon function
                                            % Version 1.1
    % recoloring
    vard = (nabla'*((eye(size(u,2))-THETA)\Z/(eye(size(u,2))-THETA)')*nabla);
    
    sigbl = sqrt(vard./T);
    
    %loop over candidate block sizes
    for b=1:length(blcand)
        blc=blcand(b);
        
        % second-level bootstrap by CB Bootstrap
        Xbootind=zeros(MBM,T);
        for i=1:MBM
            if extsimMB==0
                [ind,L,r] = cbb_seq(T,blc);
                Xbootind(i,:) = ind;
            elseif extsimMB==1
                Xbootind(i,:) = extMB(i,:);
            end
        end
        
        % computation of simulated coverage and ghat(b)
        [bsstatistics] = bsstats11(X,Xbootind,mubi,blc,floor(T/blc),mod(T,blc),T);
        
        [~,rej] = critvalues1(bsstatistics,mubi,sigbl,theta_0star,alph);
        
        Out(b) = Out(b) + rej;
    end
    timee=toc;
    FWE = sum(Out,1)./m;
    if silent==0
        time=['Estimated time remaining: ' num2str((timee-timeb)*(M-(m-1))) 's, error rates: ' num2str(FWE)];
        fprintf('%s \n',time)
    end
end

errorrates = sum(Out,1)./M

optbl = max(blcand(abs(FWE-alph)==min(abs(FWE-alph))))

csvwrite('optimalblrobustVariance.csv',[blcand errorrates optbl]);