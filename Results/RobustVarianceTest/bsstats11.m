function [bsstatistics,muboot,sigboot] = bsstats11(X,Xbootind,theta,bl,L,r,t) 
%Input:
%Xbootind is an index matrix, rows are bootstrap iterations, columns
%variables. bl is block size, l is floor(t/bl), R is mod(t/bl), t as usual.
%Output:
%bsstatistics is a matrix with corresponding studentized test statistics for each bootstrap
%iteration (row) and each variable (column)
%muboot is the test statistic w*_(T,r_s), sigboot the estimated sigma* for
%each bs-run and iteration

M=size(Xbootind,1);
muboot = zeros(M,1); sigboot = zeros(M,1); zma = zeros(M,1);
zeta = zeros(L,4);
for m=1:M
    bsinds=Xbootind(m,:);
    muboot(m,:) = log(var(X(bsinds,1))) - log(var(X(bsinds,2)));
    nu = [mean(X(bsinds,1)); mean(X(bsinds,2)); mean(X(bsinds,1).^2); mean(X(bsinds,2).^2)];
    nabla = [-(2*nu(1))/(nu(3)-nu(1)^2) (2*nu(2))/(nu(4)-nu(2)^2) 1/(nu(3)-nu(1)^2) -1/(nu(4)-nu(2)^2)]';
    Vt = [X(bsinds,1)-nu(1) X(bsinds,2)-nu(2) X(bsinds,1).^2-nu(3) X(bsinds,2).^2-nu(4)];
    Vt = Vt(1:t-r,:);
    for J=1:L
        %zeta(J,:) = mean(Vt(((J-1)*bl+1):((J-1)*bl+bl),:),1);
        zeta(J,:) = sum(Vt(((J-1)*bl+1):((J-1)*bl+bl),:),1)./sqrt(bl);
    end
    sigb = (t/(t-4))*zeta'*zeta./L;
    %for g=1:S
    %sigboot(m,g) = sqrt(1/L*1/bl*zeta(:,g)'*zeta(:,g))/sqrt(t-r);
    %end
    sigboot(m) = sqrt(nabla'*sigb*nabla/t);
    zma(m) = abs(muboot(m) - theta)/sigboot(m);

end

bsstatistics = zma;

