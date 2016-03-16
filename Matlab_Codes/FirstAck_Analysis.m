% It finds the cost of First-ACK Routing Algorithm
% als1 or als: Taylor coeffiencts of source to forwarder 1 
% als2 or bts: Taylor coefficients of source to forwarder 2
% als3: Taylor coefficients of source to forwarder 3
% forwarder 3 (2 case of J=2) has the least CiD

N = 40; % no of symbols per message
C1D = 166.4268; % cost from forwarder 1 to destination
C2D = 92.4593; % cost from forwarder 2 to destination
% C3D = 43.2937; % cost from forwarder 3 to destination
n_Tx = 300; % max transmissions eq to infinity
Cost = 0;

% % J=2 Breaking ties in favor of the node having better C_{iD}
% for i=0:n_Tx-1
%     j=i;
%     P1(i+1) = als(i+1)*sum(bts(i+2:end));
%     P2(j+1) = bts(j+1)*sum(als(j+1:end));
%     
%     Cost = Cost + (i+N+C1D)*P1(i+1) + (j+N+C2D)*P2(j+1);
% end

% % J=2 Breaking ties randomly
for i=0:n_Tx-1
    j=i;
    P1(i+1) = als(i+1)*sum(bts(i+2:end)) + (1/2)*als(i+1)*bts(i+1);
    P2(j+1) = bts(j+1)*sum(als(j+2:end)) + (1/2)*als(i+1)*bts(i+1);
    
    Cost = Cost + (i+N+C1D)*P1(i+1) + (j+N+C2D)*P2(j+1);
end

% J=3 Breaking ties randomly
% for i=0:n_Tx-1
%     j=i;
%     k=j;
%     P1(i+1) = als1(i+1)*sum(als2(j+2:end))*sum(als3(k+2:end)) + (1/2)*als1(i+1)*als2(j+1)*sum(als3(k+2:end)) + (1/2)*als1(i+1)*als3(k+1)*sum(als2(j+2:end)) + (1/3)*als1(i+1)*als2(j+1)*als3(k+1);
%     P2(j+1) = als2(j+1)*sum(als1(i+2:end))*sum(als3(k+2:end)) + (1/2)*als2(j+1)*als1(i+1)*sum(als3(k+2:end)) + (1/2)*als2(j+1)*als3(k+1)*sum(als1(i+2:end)) + (1/3)*als1(i+1)*als2(j+1)*als3(k+1);
%     P3(j+1) = als3(k+1)*sum(als1(i+2:end))*sum(als2(j+2:end)) + (1/2)*als3(k+1)*als1(i+1)*sum(als2(j+2:end)) + (1/2)*als2(j+1)*als3(k+1)*sum(als1(i+2:end)) + (1/3)*als1(i+1)*als2(j+1)*als3(k+1);
%     
%     Cost = Cost + (i+N+C1D)*P1(i+1) + (j+N+C2D)*P2(j+1) + (k+N+C3D)*P3(j+1);
% end

add1=sum(P1);
add2=sum(P2);
Cost=Cost/(add1+add2);
% add3=sum(P3);
% Cost=Cost/(add1+add2+add3);
