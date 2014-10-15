function [lambdaEst, cEst, resNorm]=EASI-SM_FOS(y,t,m,maxIter)
%[lambdaEst, cEst, resNorm]=EASI-SM(y,t,m,maxIter)
%Estimate exponential parameters using Steiglitz-McBride (SM).
%y is the data, t is the time vector, and m is the number of components.
%Model: s(n)=sum_k(c_k*lambda_k^n)+e(n), lambda_k=exp(-dt/T_2k), t=n*dt.
%maxIter is the maximum number of iterations in SM, which is effective 
%only in the case of slow convergence.
%
%Marcus Bjork, 2014

N=length(y);
finished=false;
m_max=m;
while ~finished

  [b,a]=stmcbFlex(y,[1;zeros(N-1,1)],m,maxIter);
  [cEst,lambdaEst] = residue(b,a);

  %Feasibility-based order selection (FOS)
  if any(lambdaEst<0 | lambdaEst>1 | ~isreal(lambdaEst)) 
      m=m-1;
  else
      %Compute vandermonde matrix
      vMonde=bsxfun(@power,lambdaEst.',t);
      if cond(vMonde)<1e3  %Only solutions with OK conditioning are passed
          if any(cEst<=0) && m>1 %If any cEst is negative or zero, also reduce order
              m=m-1;
          else
              finished=true;
          end
      else %If poor conditioning, reduce order (will never happen for r=1 since it is perfectly conditioned)
          m=m-1;
      end
  end
end

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

function [b,a]=stmcbFlex(y,u,m,iter)
%A more flexible implementation of Stieglitz-McBride, compared to the MATLAB built-in stmcb() function,
that does not assume unit impulse input.

%Initial A estimate
a=1;
ab_old=zeros(2*m,1);
for k=1:iter
    yf=filter(1,a,y);
    uf=filter(1,a,u);
   
    
    %Least squares (confirmed equal to stmcb for white noise)
    %ab=[-toeplitz([zeros(m-1,1);yf(1:end-1)],zeros(m,1)) toeplitz([uf(1:end-1); zeros(m-1,1)],[uf(1);zeros(m-1,1)])]\yf;
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
