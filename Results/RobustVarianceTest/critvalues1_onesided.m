function [d_K,rej] = critvalues1_onesided(bsstat,mu,sigi,theta_0,alph)
%calculates the critical value d_K from the bootstrap matrix of studentized
%test statistics

zma = bsstat;

zmat = zma;              
d_K = quantileR(zmat,1-alph,1);   

%invert generalized confidence regions, i.e. reject hypothesis
rej = zeros(1);
    if  d_K <= (mu-theta_0)/sigi
        rej = 1;
    end

rej;