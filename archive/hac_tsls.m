function [beta,vbeta,se_beta] = hac_tsls(y,x,w,nma,ikern);
%{
     Modified by MWW, 3-8-14

     Procedure for estimating the regression y = xbeta+ u using 2SLS
     The procedure produces the TSLS estimate of b
     and a hetero/autocorrelation consistent estimate of
     the autocorrelation matrix.

Input:
     y = tx1
     x = txk (regressors)
     w = tx1 (instruments)
     nma=truncation parameter (nma=0, White SEs) 
     ikern = kernel indicator
             1 => triangular
             2 => rectangular

Output:
     Beta = OLS estimate of beta (kx1)
    VBeta = Robust estimate of covariance matrix of beta (kxk)
            (Note this is computed used PINV, if X'X is singular

%}
xhat = w*(inv(w'*w)*(w'*x));
xx=xhat'*xhat;
xxi = inv(xx);
beta=xxi*(xhat'*y);
u=y-x*beta;

z = xhat.*repmat(u,1,size(x,2));
v=zeros(size(x,2),size(x,2));


% Form Kernel 
kern=zeros(nma+1,1);
for ii = 0:nma;
 kern(ii+1,1)=1;
 if nma > 0;
  if ikern == 1; 
    kern(ii+1,1)=(1-(ii/(nma+1))); 
  end;
 end;
end;

% Form Hetero-Serial Correlation Robust Covariance Matrix 
for ii = -nma:nma;
 if ii <= 0; 
    r1=1; 
    r2=size(z,1)+ii; 
 else; 
    r1=1+ii; 
    r2=size(z,1); 
 end;
 v=v + kern(abs(ii)+1,1)*(z(r1:r2,:)'*z(r1-ii:r2-ii,:));
end;

vbeta=xxi*v*xxi;
se_beta = sqrt(diag(vbeta));


end