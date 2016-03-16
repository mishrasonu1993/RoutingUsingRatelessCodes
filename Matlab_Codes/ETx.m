function EX = ETx(N,p,K_thresh,bts,k)
% ETx finds the expected number of transmissions required
% given that the receiver has more than K_thesh symbols but
% has not ACKed yet, i.e K_thresh<K_2<N
% p: erasure probability of link to best forwarder
% K_thresh: after this E(N-K,N,p) becomes less
% bts: uncondtional Taylor expansion pertaining to p
% k: k additional transmissions is over

nTx = 120;
P = zeros(N-K_thresh-1);

for j=K_thresh+1:N-1
    i = j-K_thresh;
    P(i) = 1 - ( Prob_NL(j,k) + sum(bts(1:k+1)));
end

for i=1:numel(K_thresh+1:N-2)
    P(i) = P(i)-P(i+1);
end

E = ExtraTrans(N,p);
EX = E(N-K_thresh-1:-1:1).*P / sum(P);
    
end
