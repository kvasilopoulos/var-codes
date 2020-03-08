function [c,der_c] = upperchol_id_minusj(AL,Sigma,Gamma)
%  -Zero restriction for the two shock-two instrument SVARIV model (-jth shock)
%  -Syntax:
%    [ c,der_c,e_j ] = zerorestriction(A,Sigma,Gamma)
%  -Inputs:     
%       AL: matrix of autoregressive coefficients        (n x np)
%    Sigma: covariance matrix of residuals               (n x n)
%    Gamma: covariance matrix between inst and shocks    (n x 2)
%  -Output:
%        c: vector defining the zero restriction         (n x 1)
%    der_c: matrix of derivatives of c                   (n^2p + (n*(n+1)/2) + 2n x n)
%
% This version: August 4th, 2017
% Last edited by José Luis Montiel-Olea

%% 1) Define the vector of zero restrictions implied by the lower Chol decomposition

n                    = size(Sigma,1);

aux                  = eye(n);

Sigmainv             = (Sigma\aux);

[IRFSVARIV, dIRFdmu] = IRFSVARIV_2inst_j(AL,Sigma,Gamma,...
                      1,-1,2,@upperchol_id);

c                    = Sigmainv*IRFSVARIV(:,1);

%% 2) Define the vector of derivatives

p           = size(AL,2)/n;

dcdvecA     = dIRFdmu(1:(n^2)*p,:,1)*Sigmainv;

M           = vechtovec(n);

dcdvecSigma = dIRFdmu(((n^2)*p)+1:((n^2)*p)+(n*(n+1)/2),:,1)*Sigmainv... 
             -M'*kron(Sigmainv,Sigmainv)*kron(IRFSVARIV(:,1),aux);

dcdvecGamma = dIRFdmu(((n^2)*p)+(n*(n+1)/2)+1:end,:,1)*Sigmainv;

der_c        = [ dcdvecA; dcdvecSigma; dcdvecGamma ];

end



function [IRFSVARIV,dIRFdmu] = IRFSVARIV_2inst_j(AL,Sigma,Gamma,hori,x,nvar,c)
%  -Computes IRFs identified using two external instruments possibly
%  correlated with the two target shocks
%  -Syntax:
%    [ IRFSVARIV, dIRFdmu ] = IRFSVARIV(AL,Gamma,hori,x,nvar)
%  -Inputs:     
%       AL: matrix of autoregressive coefficients        (n x np)
%    Gamma: covariances between zt and etat              (n x 2)
%     hori: number of horizons to evaluate IRFs          (1 x 1)  
%        x: scale                                        (1 x 1)
%     nvar: normalizing variable                         (1 x 1)
%        c: zero restriction                             c(AL,Sigma,Gamma)
%        j: shock over which the restriction is imposed  (1 or 2)
%  -Output:
%IRFSVARIV: vector of IRFs                               (n x hori+1)    
%dIRFdmu  : vector of derivatives of the IRFs            (n^2p + n(n+1)/2 + 2n x n x hori+1)   
% 
% This version: July 17th, 2017
% Last edited by José Luis Montiel-Olea

%% 1) Main definitions and IRFs

n            = size(Sigma,1);

j            = nvar;

p            = size(AL,2)/n;

[crest,der_c]= c(AL,Sigma,Gamma);

Bcircj       = (Gamma*[0,-1;1,0]*Gamma')*crest;

B1           = x*Bcircj./Bcircj(j,1);

Cauxsim      = [eye(n),MARep(AL,p,hori)]; 
      
Csim         = reshape(Cauxsim,[n,n,hori+1]);

IRFSVARIV    = reshape(sum(bsxfun(@times,Csim,B1'),2),[n,hori+1]);

if nargout > 1

%% Auxiliary Section 1 
%  Derivative of Bcircj w.r.t. mu = [vec(A), vech(Sigma), vec(Gamma)];

dBcircjdvecmu    = der_c*(Gamma*[0,-1;1,0]*Gamma')';

T1               = zeros(2*n,n);

T1(1:2:(2*n)-1,:)= eye(n);

T2               = zeros(2*n,n);

T2(2:2:(2*n),:)  = eye(n);

T                = [T1,T2];

dBcircjdvecGamma = kron( [0,-1;1,0]*Gamma'*crest , eye(n)) ...
                 + T'*kron(crest, [0,-1;1,0]'*Gamma');
             
dBcircdmu        = dBcircjdvecmu + ...
                   [zeros(((n^2)*p)+(n*(n+1)/2),n);dBcircjdvecGamma];  

%% Auxiliary Section 2
%  Derivative of B_j with respect to mu

aux              = eye(n);

dB_jdmu          = dBcircdmu*[(x*eye(n))-(aux(:,j)*B1')]./(Bcircj(j,1));

%% 2) Derivative of the (k,i,j)-th IRF coefficient

[G,~]    = Gmatrices(AL,MARep(AL,p,hori),p,hori,n); 
% 3D array of n^2 times (n^2 p) matrices

dIRFdmu  = zeros(((n^2)*p)+(n*(n+1)/2)+2*n,n,hori+1);

for i_hori = 1: hori+1
    
   for i_var = 1: n
       
      dIRFdmu(:,i_var,i_hori) = dB_jdmu*Csim(:,:,i_hori)'*aux(:,i_var) ...
                              + [ kron(B1',aux(:,i_var)')*G(:,:,i_hori), ...
                                 zeros(1,(n*(n+1)/2)+(2*n))]';       
   end
    
end

end
end

