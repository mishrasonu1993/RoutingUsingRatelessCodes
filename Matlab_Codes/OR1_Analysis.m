% It finds the cost of OR-1 Routing Algorithm
% als: Taylor coeffiencts of source to forwarder 1
% bts: Taylor coefficients of source to forwarder 2
% forwarder 2 has lesser CiD

N = 40; % no of symbols per message
C1D = 43.2937; % cost from forwarder 1 to destination
C2D = 43.2937; % cost from forwarder 2 to destination
n_Tx = 300; % max transmissions eq to infinity
p = 0;

Cost = 0;
P1 = zeros(1,n_Tx);
P2 = zeros(1,n_Tx);
margin = C1D-C2D; %83.2134; %39.9197; %C1D-C2D;
L_thresh = Set_Thresh(N,p,margin,bts);
if isempty(L_thresh)
    L_thresh = Inf;
end

for i=0:n_Tx-1
    j=i;
    if i<L_thresh
        P2(j+1) = bts(j+1)*sum(als(j+1:end)); 
    else P2(j+1) = bts(j+1)*sum(als(L_thresh+1:end));
    end    
    if i<L_thresh
        P1(i+1) = als(i+1)*sum(bts(i+2:end));
    end    
    Cost = Cost + (i+N+C1D)*P1(i+1) + (j+N+C2D)*P2(j+1);
end

add1=sum(P1);
add2=sum(P2);
Cost=Cost/(add1+add2);
