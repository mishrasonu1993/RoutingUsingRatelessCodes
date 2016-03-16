% It finds the Taylor series expansion coefficients for each link 
% INPUT: a 2-D matrix containing the erausre rates of the links
% Since it was a time consuming to run this code every time we changed the
% network, I obtained the taylor expansions for all probability values from
% 0 to 0.85 and stored them. So this code became redundant

file_name = 'ip_mesh_v1.txt';
fp=fopen(file_name,'r');
if(~fp)
    printf('Error:File does not exist');
	exit(0);
else
    pe = fscanf(fp,'%f');
end
fclose(fp);
pe = reshape(pe, sqrt(size(pe,1)), sqrt(size(pe,1)))';

no_nodes = size(pe,1);

N = 100;
n_Tx = 40;   % max 40 transmissions will be required
M_max = 20;  % max 20 packs lost in N transmissions
% nlinks = 3*no_nodes;  

l=0;
for i = 1:no_nodes
    for j = i+1:no_nodes
        p = pe(i,j);
        if p<0.95
            l = l+1;
            prob_map(l)=p;
        end
    end
end

nlinks = l; % total number of links in the network
form='';
for i=1:nlinks, form = strcat(form,'%f\t'); end;
fp = fopen('prob_map.txt','w');
fprintf(fp,form,prob_map);
fclose(fp);     

tlr_expansion_mat = zeros(nlinks,n_Tx);
for i=1:nlinks
    p = prob_map(i);
    syms z;
    prod = 1;
    Gsum=0;
    for M=1:M_max;
        K = N-M;
        PNNm = 1-(2^K-1)/(2^N-1);
        prod = prod * ( (PNNm)*z*(1-p) ) / ( 1-z+PNNm*z*(1-p)  ); 
        GNM = prod;
        Gsum = Gsum + nchoosek(N,M)*p^M*(1-p)^(N-M)*GNM;
    end
    G = Gsum/(1-(1-p)^N);
    tlr = taylor(G,n_Tx);
    alphas = sym2poly(tlr);
    tlr_exp_net = alphas(end:-1:1);
    tlr_expansion_mat(i,:) = tlr_exp_net;
end
       
form = '';
for i=1:nlinks
    for j=1:n_Tx, form = strcat(form,'%f\t'); end       
    form = strcat(form,'\n');
end
fp = fopen('tlr_exp_mat_mesh_v1.txt','w');
fprintf(fp,form,tlr_expansion_mat');
fclose(fp);   