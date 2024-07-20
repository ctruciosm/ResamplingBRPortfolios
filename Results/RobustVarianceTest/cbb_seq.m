function [sequence,L,r] = cbb_seq(T,bl)
%Generates a CBB index sequence of blocks of length bl from time series of
%length T. If T/bl is not an integer, the cbb index sequence 'overfills'
%and cuts off so that a sequence of length T results.
L = floor(T/bl);
r=mod(T,bl);             %check whether T/bl is an integer
if r==0
startpoints = unidrnd(T,[L,1]);
indexsequence = [1:1:T,1:1:bl];
sequence = zeros(1,T);
    for i = 1:L
        start = startpoints(i);
        sequence(((i-1)*bl+1):(i*bl)) = indexsequence(start:(start+bl-1));
    end
else startpoints = unidrnd(T,[L+1,1]);
    indexsequence = [1:1:T,1:1:bl];
    sequence = zeros(1,T+(bl-r));
    for i = 1:(L+1)
        start = startpoints(i);
        sequence(((i-1)*bl+1):(i*bl)) = indexsequence(start:(start+bl-1));
    end
    sequence = sequence(1:T);
end