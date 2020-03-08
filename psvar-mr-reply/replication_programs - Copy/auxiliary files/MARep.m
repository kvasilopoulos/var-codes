function [C] = MARep(AL,p,hori)
% Transforms the A(L) parameters of a reduced-form VAR into MA coeffs.
% Syntax:
%     [C] = MARep(AL,p,hori)
% Inputs:
%     AL: VAR model coefficients
%      p: number of lags in the VAR model
%   hori: forecast horizon
% Outputs:
%      C: MA representation coefficients
% 
% This version: June 6th, 2017

%% Reshape AL into a 3-D array

n     = size(AL,1);

vecAL = reshape(AL,[n,n,p]); 


%% Initialize the value of the auxiliary array vecALrevT

vecALrevT = zeros(n,n,hori);

for i = 1:hori
    
    if i<(hori-p)+1
        
        vecALrevT(:,:,i) = zeros(n,n);
        
    else
        
        vecALrevT(:,:,i) = vecAL(:,:,(hori-i)+1)';
        
    end
    
end

vecALrevT = reshape(vecALrevT,[n,n*hori]);


%% MA coefficients

C = repmat(vecAL(:,:,1),[1,hori]);

for i = 1:hori-1
    
    C(:,(n*i)+1:(n*(i+1))) = [eye(n),C(:,1:n*i)] ...
                           * vecALrevT(:,(hori*n-(n*(i+1)))+1:end)';  

end

end