function VARci = doDeltaMethod_2shock(VAR,clevel,NWlags)

scale = -1;
[Whatall,~,~] = CovAhat_Sigmahat_Gamma(VAR.p,VAR.X(:,[end end-size(VAR.DET,2)+1:end-1 1:VAR.n*VAR.p]),VAR.m,VAR.res',NWlags);                

AL          = VAR.bet(1:end-size(VAR.DET,2),:)';
Gamma       = VAR.res'*VAR.m./VAR.T; 
Caux        = [eye(VAR.n),MARep(AL,VAR.p,VAR.irhor-1)]; 
C        = reshape(Caux,[VAR.n,VAR.n,VAR.irhor]); 

[G,~]    = Gmatrices(AL,Caux,VAR.p,VAR.irhor,VAR.n); 

for shock = 1:4

    if (shock==1) || (shock==2)
    nvar = 1;
    elseif (shock==3) || (shock==4)
    nvar = 2;
    end

    if shock==1
    [crest,der_c] = upperchol_id_minusj(AL,VAR.Sigma,Gamma);
    elseif (shock==2)||(shock==3)
    [crest,der_c] = lowerchol_id(AL,VAR.Sigma,Gamma);    
    elseif shock==4
    [crest,der_c] = lowerchol_id_minusj(AL,VAR.Sigma,Gamma); 
    end

    Bcircj         = (Gamma*[0,-1;1,0]*Gamma')*crest;
    B1             =  scale*Bcircj./Bcircj(nvar,1);
    IRFSVARIV(:,:,shock)    = reshape(sum(bsxfun(@times,C,B1'),2),[VAR.n,VAR.irhor]);

   
    %% Auxiliary Section 1 
    %  Derivative of Bcircj w.r.t. mu = [vec(A), vech(Sigma), vec(Gamma)];

        dBcircjdvecmu    = der_c*(Gamma*[0,-1;1,0]*Gamma')';
        T1               = zeros(2*VAR.n,VAR.n);
        T1(1:2:(2*VAR.n)-1,:)= eye(VAR.n);
        T2               = zeros(2*VAR.n,VAR.n);
        T2(2:2:(2*VAR.n),:)  = eye(VAR.n);
        T                = [T1,T2];
        dBcircjdvecGamma = kron( [0,-1;1,0]*Gamma'*crest , eye(VAR.n)) ...
                         + T'*kron(crest, [0,-1;1,0]'*Gamma');
        dBcircdmu        = dBcircjdvecmu + ...
                           [zeros(((VAR.n^2)*VAR.p)+(VAR.n*(VAR.n+1)/2),VAR.n);dBcircjdvecGamma];  


        aux              = eye(VAR.n);
        dB_jdmu          = dBcircdmu*[(scale*eye(VAR.n))-(aux(:,nvar)*B1')]./(Bcircj(nvar,1));

        for i_hori = 1: VAR.irhor
           for i_var = 1: VAR.n
              V = dB_jdmu*C(:,:,i_hori)'*aux(:,i_var) ...
                                      + [ kron(B1',aux(:,i_var)')*G(:,:,i_hori), ...
                                      zeros(1,(VAR.n*(VAR.n+1)/2)+(2*VAR.n))]'; 
              S = (V'*Whatall*V./VAR.T).^.5;     
             for jj = 1:length(clevel)
                critval         = norminv(1-((1-clevel(jj)/100)/2),0,1);
                VARci.irsH(i_hori,i_var,shock,jj) =  IRFSVARIV(i_var,i_hori,shock) + critval*S; 
                VARci.irsL(i_hori,i_var,shock,jj) =  IRFSVARIV(i_var,i_hori,shock) - critval*S;        
             end    
           end
        end
end

