function [lambdaEst, cEst, resNorm]=test(y,t,m,maxIter)
%[lambdaEst, cEst, resNorm]=test(y,t,m,maxIter)
%Estimate exponential parameters using Steiglitz-McBride (SM).
%y is the data, t is the time vector, and m is the number of components.
%Model: s(n)=sum_k(c_k*lambda_k^n)+e(n), lambda_k=exp(-dt/T_2k), t=n*dt.
%maxIter is the maximum number of iterations in SM, which is effective 
%only in the case of slow convergence.
%
%Marcus Bjork, 2014

m_max=m;
lambdaEst=[1, 2 ,3]';
cEst=[1, 2 ,3]';

%Sort based on decay
[lambdaEst,I]=sort(lambdaEst,'descend');
cEst=cEst(I);

%Compute norm of residuals
resNorm=norm(y-bsxfun(@power,lambdaEst.',t)*cEst);

%If FOS estimated order is less than the set maximum order, append zeroes.
if m<m_max
    lambdaEst=[lambdaEst;zeros(m_max-m,1)];
    cEst=[cEst;zeros(m_max-m,1)];
end

end
