function L_thresh = Set_Thresh(N,p,margin, bts)
% SET_THERSH sets the threshold after (including) which 
% forwarder 1 is never selected, i.e C*S2|L < margin
% L_thresh = argmin_j(C*S2|j<margin)

% tlr_exp_net = tlr_expansion(N,p,1);
% alphas = tlr_exp_net(2:70);

alphas = bts(2:end)/sum(bts(2:end));

cost = -1*ones(1,11);

for L=0:300
    betas = alphas(L+1:end)./sum(alphas(L+1:end));
    n=numel(betas);
    cost(L+1) = sum([1:n].*betas);
end
L_thresh = find(cost<margin,1)-1; 
% minus 1 just because in MATLAB indexing starts from 1,
% not 0. 
end