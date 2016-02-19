function data=genDataT2(T2,c,t)
%data=genData(T2,c,t)
%Simulate the exponential model to generate data. Does NOT assume uniform
%sampling.
%
% T2 - time constants of the exponentials (vector)
% c  - amplitudes of the corresponding exponentials (vector)
% t  - time vector (same units as T2)

data=exp(-bsxfun(@times,1./T2(:)',t(:)))*c(:);

