function cost = C_star_SDL(bts)

alphas = bts(2:end)/sum(bts(2:end));

cost = -1*ones(1,11);

for L=0:300
    betas = alphas(L+1:end)./sum(alphas(L+1:end));
    n=numel(betas);
    cost(L+1) = sum([1:n].*betas);
end
end


