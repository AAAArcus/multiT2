function CRBvec=CRB(c_true,T2_true,noiseVar,t)
%CRBvec=CRB(c_true,T2_true,noiseVar,TE)
%Computes the numerical Cramér-Rao lower bound (variance) of [T2 c] under 
%the Gaussian assumption, given true parameter values, noise variances, and 
%sampling times t.
%Model: s(t)=sum_k(c_k*exp(-t/T_2k)+e(t)
%Works for vector noiseVar and returns a matrix with each parameter CRB
%along a column (row vector for scalar noiseVar).
%
%Marcus Björk, 2014

%Force column vectors
t=t(:);
c_true=c_true(:);
T2_true=T2_true(:);

%Jacobian
lambda=exp(bsxfun(@times,-t,1./T2_true'));
J=[bsxfun(@times,t,(c_true./T2_true.^2)').*lambda lambda];

%Return diagonal of CRB matrix (inverse of FIM)
CRBvec=noiseVar(:)*diag(inv((J'*J)))';
