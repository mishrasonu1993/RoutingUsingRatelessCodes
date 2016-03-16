function [ tlr_exp_net ] = taylorexp_newTx( N, p, cond )
% TLR_EXPANSION finds the Taylor Series expansion for EURC
% given N and p. It should be kept in mind which Taylor
% Series expansion you want--conditional or unconditional.

syms z;
pdt = 1;
Gsum = nchoosek(N,1)*p*(1-p)^N;
M_max = N;
n_Tx = 150;

for M = 1 : M_max
    K = N-M;
    PNm = 1-(2^K-1)/(2^N-1);
    pdt = pdt * ( (PNm)*z*(1-p) ) / ( 1-z+PNm*z*(1-p)  );
    GNM = pdt;
    Gsum = Gsum + nchoosek(N+1,M+1)*p^(M+1)*(1-p)^(N-M)*GNM;
end

G = z*Gsum/(1-(1-p)^N);

if ~cond
    tlr = taylor(z*Gsum,n_Tx);
    alphas = sym2poly(tlr);
    tlr_exp_net = alphas(end:-1:1);
    tlr_exp_net(1) = (1-p)^N;
else
    tlr = taylor(G, n_Tx);
    alphas = sym2poly(tlr);
    tlr_exp_net = alphas(end:-1:1);
end
end