function results = vare2(y,nlag,const)
%Vector Autoregression Estimation
%y is [number of obs x variables] matrix
%nlag is number of lags
%const=1 is include intercept

[nobs neqs] = size(y);

results.meth = 'vare';

% adjust nobs to feed the lags
nobse = nobs - nlag;

% nvar adjusted for constant term 
 k = neqs*nlag + 1;
 nvar = k;

results.nvar = nvar;

xlag = mlag(y,nlag);

results.nobs = nobse;
results.neqs = neqs;
results.nlag = nlag;


% form x-matrix
if const==1
xmat = [xlag(nlag+1:nobs,:) ones(nobs-nlag,1)];
else
xmat = [xlag(nlag+1:nobs,:)];
end;


% pull out each y-vector and run regressions
for j=1:neqs;

 yvec = y(nlag+1:nobs,j);
 res = ols(yvec,xmat);
 results(j).beta  = res.beta;      % bhats
 results(j).resid = res.resid;     % resids 
 sigu = res.resid'*res.resid;
 results(j).yhat = xmat*results(j).beta;       % yhats
 results(j).y    = yvec;           % actual y
end; 

 



