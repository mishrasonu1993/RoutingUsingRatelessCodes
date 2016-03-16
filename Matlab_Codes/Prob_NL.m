function Prob = Prob_NL(K_thresh,index)
% PROB_NL finds the probability that the concerned
% node has <= K_thresh packets

p = 0.50; % erasure probability of source to node
N = 40; % total symbols per message
n_Tx = 120; % eq to infinity
M_max = N; % max packs lost in N Tx
M_min = N-K_thresh; % min packs lost in N Tx
syms z;
Prob = zeros(1,numel(index));

for i=1:numel(index)
    L = index(i);
    prod = 1;
    low_lim = N-K_thresh;
    up_lim = low_lim-1;
    
    for m = M_min:M_max
        while up_lim~=m
            up_lim = up_lim+1;
            n = up_lim;
            K = N-n;
            PNNn = 1-(2^K-1)/(2^N-1);
            Hn = ( (PNNn)*z*(1-p) ) / ( 1-z+PNNn*z*(1-p)  );
            prod = prod * Hn;
        end
        tlr = taylor(prod,n_Tx);
        alphas = sym2poly(tlr);
        gamma(m+1,:) = alphas(end:-1:1);
        sum_gam = sum(gamma(m+1,L+2:n_Tx));
        Prob(i) = Prob(i) + nchoosek(N,m)*p^m*(1-p)^(N-m)*sum_gam;
    end
end
end