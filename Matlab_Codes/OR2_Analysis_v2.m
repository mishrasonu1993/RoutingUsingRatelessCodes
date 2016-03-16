% It finds the cost of OR-2 Routing Algorithm
% als: Taylor coeffiencts of source to forwarder 1
% bts: Taylor coefficients of source to forwarder 2
% forwarder 2 has lesser CiD

N = 40; % no of symbols per message
C1D = 147.1682; % cost from forwarder 1 to destination
C2D = 123.7450; % cost from forwarder 2 to destination
n_Tx = 120; % max transmissions eq to infinity
p = 0.55; % erasure probability from source to forwarder having less CiD

Cost = 0;
P1 = zeros(1,n_Tx);
P2 = zeros(1,n_Tx);
P21 = zeros(1,n_Tx);
P22 = zeros(1,n_Tx);
margin = C1D-C2D;
% if node 2 has more than K_thresh packets, it'll be selected as the
% forwarder
E = ExtraTrans (N,p);
K_thresh = N - ( find(E<margin,1,'last')+1);
%EX = ETx(N,p,K_thresh,bts);

for i=0:n_Tx-1
    j=i
    P21(j+1) = bts(j+1)*( sum(als(j+1:end)));
    P22(j+1) = als(j+1)*( 1-( Prob_NL(K_thresh, j) + sum(bts(1:j+1))));
    % node 1 ACKs and node 2 has <= K_thresh packets
    P1(i+1) = als(i+1)*Prob_NL(K_thresh, i);
    
    EX = ETx(N,p,K_thresh,bts,i);
    Cost = Cost + (i+N+C1D)*P1(i+1) + (j+N+C2D)*P21(j+1) + (j+N+EX+C2D)*P22(j+1);
end

add1=sum(P1);
add2=sum(P2);