% It finds the mean link costs of all the links in the network.
% INPUT: a 2-D matrix containing the erausre rates of the links
% OUTPUT: a 2-D matrix containing the mean link costs of the links

file_name = 'ip_2h_journal2.txt';
fp=fopen(file_name,'r');
if(~fp)
    printf('Error:File does not exist');
	exit(0);
else
    pe = fscanf(fp,'%f');
end
fclose(fp);
pe = reshape(pe, sqrt(size(pe,1)), sqrt(size(pe,1)))';
N = 40;
version = 'old';
no_nodes = size(pe,1);
mean_lnk_cst = 9999*ones(no_nodes, no_nodes);
for i=1:no_nodes
    for j=1:no_nodes
        if i==j, continue;end
        mean_lnk_cst(i,j) = Analysis(pe(i,j),N,version);        
    end
end

form='';
for i=1:no_nodes
    for j=1:no_nodes
        form = strcat(form,'%f\t');
    end       
    form = strcat(form,'\n');
end

fp = fopen('lcost_ip_2h_journal2.txt','w');
fprintf(fp,form,mean_lnk_cst');
fclose(fp);
