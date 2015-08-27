function [lambdaEst, cEst, resNorm]=EASI_SM_BIC(y,t,m,maxIter)
%[lambdaEst, cEst, resNorm]=EASI_SM_BIC(y,t,m,maxIter)
%Estimate exponential parameters using Steiglitz-McBride (SM).
%y is the data, t is the time vector, and m is the maximum number of components.
%Model: s(n)=sum_k(c_k*lambda_k^n)+e(n), lambda_k=exp(-dt/T_2k), t=n*dt.
%maxIter is the maximum number of iterations in SM, which is effective 
%only in the case of slow convergence.
%
%Marcus Bjärk, 2014

m_max=m;
N=length(y);
BIC=zeros(m_max,1);
lambdaEst=zeros(m_max);
cEst=zeros(m_max);

for r=1:m_max
  [b,a]=stmcbFlex(y,[1;zeros(length(y)-1,1)],r,maxIter);
  [cEst(1:r,r),lambdaEst(1:r,r)] = residue(b,a);

  %Compute BIC criterion
  BIC(r)=N*log(norm(y-bsxfun(@power,lambdaEst(1:r,r).',t)*cEst(1:r,r)).^2)+2*r*log(N);
end

%Choose the BIC solution
[~,ind]=min(BIC);
cEst=cEst(:,ind);  
lambdaEst=lambdaEst(:,ind);

%Compute norm of residuals
resNorm=norm(y-bsxfun(@power,lambdaEst.',t)*cEst);

%Sort based on decay
[lambdaEst,I]=sort(lambdaEst,'descend');
cEst=cEst(I);
end

function [b,a]=stmcbFlex(y,u,m,iter)
%A more flexible implementation of Stieglitz-McBride, compared to the MATLAB built-in stmcb() function,
%that does not assume unit impulse input.

%Initial A estimate
a=1;
ab_old=zeros(2*m,1);
for k=1:iter
    yf=filter(1,a,y);
    uf=filter(1,a,u);
    
    %Least squares (confirmed equal to stmcb for white noise)
    ab=[-toeplitz([0;yf(1:end-1)],zeros(m,1)) toeplitz(uf,[uf(1);zeros(m-1,1)])]\yf;
    a=[1; ab(1:m)];
    %Stabilization needed??
    rootsA=roots(a);
    indic=abs(rootsA)>=1;
    if any(indic) %Mirror
        rootsA(indic)=1./rootsA(indic);
        a=poly(rootsA);
    elseif norm(ab-ab_old)/norm(ab_old)<1e-8
        break;
    end
    ab_old=ab;
end

b=ab(m+1:end);

end
