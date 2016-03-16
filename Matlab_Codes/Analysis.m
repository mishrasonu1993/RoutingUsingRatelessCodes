function ntrans = Analysis(p,N,version)
%N = 40;  %no of packets in one message
K = 1; %no of bits in one packet
ps = 1 - (1-p)^K;
%ps = p; %29 Nov
Sum = 0;

E = ExtraTrans(N,ps);

if strcmp(version,'old')
    for M = 1:N
        Sum = Sum + E(M)*nchoosek(N,M)*(ps)^M*(1-ps)^(N-M);
    end
    ntrans = N + Sum;

elseif strcmp(version,'new')
    E_new(1)= (1-p) + p*( 1+ E(1));
    for i=2:N
        E_new(i) = 1+ (1-p)*E(i-1) + p*E(i);
    end
    
    for M = 1:N
        Sum = Sum + E_new(M)*nchoosek(N,M)*(ps)^M*(1-ps)^(N-M);
    end
    ntrans = N + Sum;
end

end
