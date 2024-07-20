function [band,Z,alpha,num,den,bandoverall,thet,sigma2] = andmon6cvm2(gmmopt,u,THETA)
%function [band,D,vari,bandoverall] = andmon6(gmmopt,u,THETA)
%
% ANDREWS-MONAHAN HAC ESTIMATOR
%
% Modified version by:
% Dan Wunderli, Institute for Empirical Research in Economics, Zurich
% University, dan.wunderli@mac.com
% 11/1/07
%
% Original version:
% WRITTEN BY:  Mike Cliff, Purdue Finance,  mcliff@mgmt.purdue.edu
% VERSION 1.1  
% CREATED: 11/6/99

%========================================================================
%   INITIALIZATIONS
%========================================================================

eflag = 0;
p = gmmopt.aminfo.p;
q = gmmopt.aminfo.q;
vardum = gmmopt.aminfo.vardum;
diagdum = gmmopt.aminfo.diagdum;
nobs = size(u,1);
amdum=1;

if p > 0
  if p > 1, error('ANDMON written for p<=1 in ARMA(p,q)/VAR(p)');  end
  if q > 1, error('ANDMON written for ARMA(p,q), q <=1 when p = 1'); end
end

switch gmmopt.aminfo.kernel
    case 'H',  const = 0.6611; ktype = 2; kname = 'Truncated (Hansen)';
    case 'NW', const = 1.1447; ktype = 1; kname = 'Bartlett (Newey-West)';
    case 'G',  const = 2.6614; ktype = 2; kname = 'Parzen (Gallant)';
    case 'QS', const = 1.3221; ktype = 2; kname = 'Quadratic-Spectral';
    case 'TH', const = 1.7462; ktype = 2; kname = 'Tukey-Hanning';
    otherwise error('Need kernel in ANDMON');
end

north = size(u,2);  
In = eye(north);
alpha = zeros(north,1);
band = zeros(north,1);
theta = zeros(north,1);
sigma = zeros(north,1);
num = zeros(north,1);
den = zeros(north,1);

interc=1;
  for i = 1:north
        tempout = ar2(u(:,i),p,interc);         %Estimate the AR model

        thet(i) = tempout.beta(interc+1);       %rho_a in equation (3.6)

    	sigma2(i) = tempout.sige;               %sigma_a in equation (3.6)
        
        % --- AR(1) Case --------------------------------------------
        num(i) = 4*sigma2(i)^2*thet(i)^2/((1-thet(i))^8);
        den(i) = sigma2(i)^2/((1-thet(i))^4);
        alpha(i) = num(i)/den(i);
        band(i) = const*(alpha(i)*(nobs))^0.2;
  end
 bandoverall = const*(sum(num)/sum(den)*nobs)^0.2;

%========================================================================
%   FINAL CALCULATIONS
%========================================================================

% --- The Matrix D for 'Recoloring' --------------------------------------
D = zeros(1,north);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following is taken partly from the gmm toolbox, function gmmS.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxlags = nobs-1;
Z = zeros(north,north);

for lag = 0:maxlags

    Gamma = u(1:nobs-lag,:)'*u(1+lag:nobs,:);	% Vectorized (fast!!) as in (2.3) of AndMon paper

    Gamma = Gamma/nobs;
    if lag >= 1, Gamma = Gamma + Gamma'; end

    x = lag/bandoverall;
    switch gmmopt.aminfo.kernel
        case 'G'
            % --- Parzen (Gallant) --------------------------
            if x < 0.5
                wt = 1 - 6*x^2 + 6*x^3;
            elseif x < 1
                wt = 2*(1-x)^3;
            else
                wt = 0;
            end
        case 'QS'
            % --- Quadratic Spectral -------------------------------------------
            term = 6*pi*x/5;
            if lag == 0
                wt = 1;
            else
                wt = 25*(sin(term)/term - cos(term))/(12*pi^2*x^2);
            end
    end
    Z = Z + wt*Gamma*nobs/(nobs-size(u,2));
end