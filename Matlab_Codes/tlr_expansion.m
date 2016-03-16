function tlr_exp_net = tlr_expansion(N,p,cond)
% TLR_EXPANSION finds the Taylor Series expansion for URC
% given N and p. It should be kept in mind which Taylor
% Series expansion you want--conditional or unconditional.

syms z;
prod = 1;
Gsum = 0;
M_max = N; 
n_Tx = 150; 

for M=1:M_max
    K = N-M;
    PNNm = 1-(2^K-1)/(2^N-1);
    prod = prod * ( (PNNm)*z*(1-p) ) / ( 1-z+PNNm*z*(1-p)  );
    GNM = prod;
    Gsum = Gsum + nchoosek(N,M)*p^M*(1-p)^(N-M)*GNM;
end

G = Gsum/(1-(1-p)^N);
if ~cond
    tlr = taylor(Gsum,n_Tx);
    alphas = sym2poly(tlr);
    tlr_exp_net = alphas(end:-1:1);
    tlr_exp_net(1) = (1-p)^N;
else
    tlr = taylor(G,n_Tx);
    alphas = sym2poly(tlr);
    tlr_exp_net = alphas(end:-1:1);
end

end