function E = ExtraTrans (N,ps)
E = zeros(1,N);

E(1) = 1/((1-ps)*(1-(2^(N-1)-1)/(2^N-1)));
%E(1) = 2/(1-ps);
for M = 2:N   % can this create problema in later indexing
    K = N-M;
    %PNM = 1-2^(-N+K);   % earlier analysis
    PNM = 1-(2^K-1)/(2^N-1);      % new analysis
    E(M) = 1/((1-ps)*PNM) + E(M-1);        
end
end

