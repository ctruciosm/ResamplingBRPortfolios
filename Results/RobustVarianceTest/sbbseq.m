function [sequence,bls] = sbbseq(T,Tstar,avbl)
%Generates a Stationary BB index sequence of blocks of length bl from time series of
%length T. blave is average block length used as mean of the geometric
%distribution generating the random block sizes.
%If the SBB index sequence 'overfills',
%sbb_seq cuts off so that a sequence of length T results.
bls = geornd(1/avbl,1,T+Tstar);
i = length(bls);
while sum(bls)<T+Tstar
    bls(1+i) = geornd(1/avbl,1);
    i=i+1;
end

while sum(bls(1:end-1))>T+Tstar
    bls = bls(1:end-1);
end
bls(end) = T+Tstar-sum(bls(1:end-1));

startpoints = unidrnd(T,[length(bls) 1]);
indexsequence = [1:1:T,1:1:max(bls)];
sequence = zeros(1,T+Tstar);
where=1;

for i = 1:length(bls)
    start = startpoints(i);
    if where <= sum(bls(1:i))
        sequence(where:where+bls(i)-1) = indexsequence(start:(start+bls(i)-1));
    end
    where = 1+sum(bls(1:i));
end
