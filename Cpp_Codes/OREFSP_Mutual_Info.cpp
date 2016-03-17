#define _CRT_SECURE_NO_DEPRECATE
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <stdlib.h>
#include <queue>
#include <deque>
#include <algorithm>
#include <vector>
#include <time.h>                      // define time()
#include "randomc.h"                   // define classes for random number generators
#include "stocc.h"                     // define random library classes


#ifndef MULTIFILE_PROJECT
// If compiled as a single file then include these cpp files, 
// If compiled as a project then compile and link in these cpp files.
   #include "mersenne.cpp"             // code for random number generator
   #include "stoc1.cpp"                // random library source code
   #include "userintf.cpp"             // define system specific user interface
#endif

# define nmessgs 1000                 // total number of messages
# define npacks 40                    // no of packets per message 
# define no_nodes 11                  // no of nodes in network 
# define M_max 5                      // max no of lost packs
# define n_Tx 120                     // max no of additional retrans
# define nlinks 23                    // no of links in network 


using namespace std;
int no_dest;                           //no of destinations
float pe[no_nodes][no_nodes], mean_cost[no_nodes], mean_link_cost[no_nodes][no_nodes], 
	  tlr_coeffs[M_max][n_Tx][nlinks], prob_map[nlinks];
/* pe:erasure probability, mean_cost:cost to destinaion, mean_link_cost:cost of links,
   tlr_coeffs:taylor coefficients, prob_map: probability map*/

int linr_indep(bool[npacks][npacks][no_nodes], bool[npacks], int, int);  
void cost2D_est(int, int); 
float addnl_Tx_est(int, float, int, float);
void mean_link_cost_calc();
//void tlr_coeffs_read();
double nchoosek(int,int);
void reset_all(int[no_nodes], bool[no_nodes], int[no_nodes], bool[no_nodes], bool[npacks][npacks][no_nodes], int);

int main()
{
	char file_input, file_name[50], link, modify, ext,  input, fname[50];
	FILE *fp, *fptr;
	int source, node1, node2, destination;
	int messgID, ntrans, fwd, nxtfwdset[no_nodes], partfilledfwds[no_nodes], nfwds, pfwds, packID, trans;
	bool rec_matrix[npacks][npacks][no_nodes], vect[npacks], ACK[no_nodes], ackfwd[no_nodes], fwd_select=false;
	int rec_vect[no_nodes], bestfwd, sfwd=-1, freq1=0, freq2=0, cost1=0, cost2=0;
	float ackmincost, nackmincost;
	int seed = (int)time(0);            // random seed
	StochasticLib1 sto(seed);           // make instance of random library

	printf("Do you want to save the results in a file?-y/n: ");
	cin >> input;

	switch(input)
	{
		case 'y':
			printf("Enter the file name in which you want to store the results\n(format:filename.txt): ");
			cin >> fname;
			fptr = fopen(fname,"w");
			break;
		case 'n':
			break;
		default:
			break;
	}

	printf("Do you want to input ETX values from a file?-y/n: ");
	cin >> file_input;

	switch(file_input)
	{
		case 'n':
		{
			printf("Enter the file name in which you want to store network's pe value\n(format:filename.txt): ");
			cin >> file_name;
			for(int i=0;i<no_nodes;i++)
			{
				for(int j=0;j<no_nodes;j++)
				{
					if(i==j) pe[i][j]=0;
					else 
					{
						printf("Is there a link between %d and %d?-y/n: ",i,j);
						cin >> link;
						if(link=='y')
						{
							printf("Enter the symbol erasure probability between %d and %d: ", i,j);
							cin >> pe[i][j];
						}
						else pe[i][j]=0.95;
					}
				}
			}
			break;
		}
		case 'y':
		{
			printf("Enter the name of file from which we input the network data\n(format:filename.txt): ");
			cin >> file_name;
			fp=fopen(file_name,"r");
			if(!fp)
			{
				printf("Error:File doesn't exist");
				exit(0);
			}
			else
			{
				for(int i=0;i<no_nodes;i++)
				{
					for(int j=0;j<no_nodes;j++)
					{
						fscanf(fp,"%f",&pe[i][j]);
					}
				}
				fclose(fp);
				printf("Do you want to modify any link's pe value in \"%s\"-y/n: ",file_name);
				cin >> modify;

				while(modify=='y')
				{
					printf("Enter the nodes which you wish to change the links' pe\n(format: node1 node2): ");
					cin >> node1 >> node2;
					printf("Previous value=%f, Enter new value",pe[node1][node2]);
					cin >> pe[node1][node2];
					printf("Do you want to change another link's pe?-y/n: ");
					cin >> modify;
				}
			}
			break;
		}
		default: printf("Invalid Response -- system will now exit"); exit(0);
	}
    
    fp = fopen(file_name,"w");

    for(int i=0;i<no_nodes;i++)
    {
    	for(int j=0;j<no_nodes;j++)
    		fprintf(fp,"%f\t",pe[i][j]);
    	fprintf(fp,"\n");
    }
    fclose(fp);
    printf("Enter Source node Id (Id range-<0-%d>): ",no_nodes-1);
    cin >> source;
    printf("Enter no.of destinations: ");
    cin >> no_dest;
    
    for(int i=0;i<no_dest;i++)
    {
    	printf("Enter the destination Id %d/%d (Range (0-%d)): ",i+1,no_dest,no_nodes-1);
    	cin >> destination;
    }

    mean_link_cost_calc();              // Calculate link costs
    cost2D_est(source, destination);    // Estimate costs from nodes to destination
    // tlr_coeffs_read();                  // read taylor coefficients from a file

    printf("\n--------------------------------------------------------\n");
    printf("\t\tsimulation begins\t\t");
    printf("\n--------------------------------------------------------\n");

    messgID = 0;

    while(messgID < nmessgs)
    {
    	trans = 0;
    	sfwd = -1;
    	fwd = source;                   // initially, source has the message

    	for (int i=0; i<no_nodes; ++i)
    	{
    		rec_vect[i] = 0;         // all nodes have zero packets
    		ackfwd[i] = false;       // none of the nodes have sent any ACK
    		nxtfwdset[i] = 0;        // downstream forwarders aren't chosen yet
    		partfilledfwds[i] = 0;   // no fwd is partially filled
    		ACK[i] = false;          // none of the nodes have received any ACK
    	}

    	for (int i=0; i<npacks; ++i)
    	{
    		for (int j=0; j<npacks; ++j)
    			for (int k=0; k<no_nodes; ++k)
    				rec_matrix[i][j][k] = false;  // all nodes have empty buckets
		}

		while (fwd!=destination)        // until the destination has the message
		{
			int j = 0, k = 0;
			for (int i=0; i<no_nodes; ++i)  // selecting the potential forwarders
			{
				if (i==fwd) continue;
				if (pe[fwd][i] < 0.94 && mean_cost[i] < mean_cost[fwd])
				{
					nxtfwdset[j]=i;
					++j;
				}
				if (rec_vect[i])
				{
					partfilledfwds[k]=i;
					++k;
				}
			}

			nfwds = j;
			printf("Forwarder set of node %d : node ID",fwd);

			for (int i=0; i<nfwds; ++i)
			{
				printf(" %d,", nxtfwdset[i]);
			}
			
			printf("\n");
			pfwds = k;
			if (pfwds)
			{
				printf("Partially Filled:");
				for (int i=0; i<pfwds; ++i)
				{
					printf(" %d,", partfilledfwds[i]);
				}
				printf("\n");
			}

			ntrans = 0;                 // no of transmissions is zero
			ACK[fwd] = false;			// fwd hasn't received any ACK
			fwd_select = false;			// actual forwarder hasn't been selected
			bestfwd = -1;

			for (packID = 0; packID < npacks; ++packID)  // initial transmissions of N packets
			{
				ntrans++;
				for (int i=0; i<nfwds; ++i)
				{
					int f = nxtfwdset[i];
					int r = sto.Binomial(1, pe[fwd][f]);   // sampling binomial random variable
					if (!r)
					{
						int found = 0;
						for (int j=0; j<pfwds; ++j)
						{
							if (f==partfilledfwds[j])
							{
								found = 1;
								break;
							}
						}
						if (!found)
						{
							rec_vect[f] = rec_vect[f] + 1;
							int j = rec_vect[f];
							if (j==npacks)
							{
								ACK[fwd]=true;
								ackfwd[f] = true;
								printf("node : %d sends ACK\n",f);
							}
							rec_matrix[packID][j-1][f] = true;
						}
						else
						{
							for (int p=0; p<npacks; ++p)
							{
								vect[p] = false;
							}
							vect[packID] = true;
							int li = linr_indep(rec_matrix, vect, f, rec_vect[f]);  // check the independence of the vectors
							if (!li) continue;
							rec_vect[f] = rec_vect[f] + 1;
							int j = rec_vect[f];
							for (int k=0; k<npacks; ++k)
							{
								rec_matrix[k][j-1][f] = vect[k];  // insert the vector into the bucket
							}

							if (j==npacks)    // if received all the npacks, send the ACK
							{
								ACK[fwd]=true;
								ackfwd[f] = true;
								printf("node : %d sends ACK\n",f);
							}
						}
					}
				}
				if (ACK[fwd])
					break;
			}

			if (ACK[fwd])               // is it wise to continue transmissions
			{
				ackmincost = 9999;
				nackmincost = 9999;
				for (int i=0; i<nfwds; ++i)
				{
					int f = nxtfwdset[i];
					if (ackfwd[f])
					{
						if (mean_cost[f] < ackmincost)
						{
							ackmincost = mean_cost[f];
							bestfwd = f;
						}
					}
					else
					{
						float linccost = mean_link_cost[fwd][f];
						float addcost = addnl_Tx_est(rec_vect[f],pe[fwd][f],ntrans, linccost);    // addnl retrans others need to receive the message
						if (mean_cost[f] + addcost < nackmincost)
							nackmincost = mean_cost[f] + addcost;
					}
				}
				if (ackmincost < nackmincost)  // if it holds, use the bestfwd as the actual forwarder
				{
					fwd_select = true;
					fwd = bestfwd;
					if (fwd!=3)
					{
						sfwd = fwd;
					}
					if (fwd==1) 
                    {
                        freq1++; cost1+=ntrans;
                    }
					if (fwd==2) 
                    {
                        freq2++; cost2+=ntrans;
                    }
					reset_all(rec_vect, ackfwd, nxtfwdset, ACK, rec_matrix, fwd);
				}
			}

			while(!fwd_select)			// if no forwarder has been selected yet
			{
				ntrans = ntrans + 1;
				int sum = 0;
				while(sum==0)
				{
					for (int i=0; i<npacks; ++i)
					{
						vect[i] = sto.Binomial(1, 0.5);
						sum = sum + vect[i];
					}
				}

				for (int i=0; i<nfwds; ++i)		// send to all the potential forwarders
				{
					int f = nxtfwdset[i];
					if(ackfwd[f] == 1) continue;  // if the forwarder has already received all npacks
					int r = sto.Binomial(1, pe[fwd][f]);  // sample the binomial random variable
					if (!r)
					{
						int li = linr_indep(rec_matrix, vect, f, rec_vect[f]);  // check the independence of the vectors  
						if (!li) continue;
						rec_vect[f] = rec_vect[f] + 1;
						j = rec_vect[f];
						for (int k=0; k<npacks; ++k)
						{
							rec_matrix[k][j-1][f] = vect[k];  // insert the vector into the bucket
						}
						if (j==npacks)    // if received all the npacks, send the ACK
						{
							ACK[fwd]=true;
							ackfwd[f] = true;
							printf("node : %d sends ACK\n",f);
						}
					}
				}

				if (ACK[fwd])			// check if it wise to continue transmission
				{
					ackmincost = 9999;
					nackmincost = 9999;
					for (int i=0; i<nfwds; ++i)
					{
						int f = nxtfwdset[i];
						if (ackfwd[f])
						{
							if (mean_cost[f] < ackmincost)
							{
								ackmincost = mean_cost[f];
								bestfwd = f;
							}
						}
						else
						{
							float linccost = mean_link_cost[fwd][f];
							float addcost = addnl_Tx_est(rec_vect[f],pe[fwd][f],ntrans, linccost);
							if (mean_cost[f] + addcost < nackmincost)
								nackmincost = mean_cost[f] + addcost;
						}
					}
					if (ackmincost < nackmincost)
					{
						fwd_select = true;
						fwd = bestfwd;
						
						if (fwd!=3)
						{
							sfwd = fwd;
						}
						if (fwd==1) 
                        {
                            freq1++; cost1+=ntrans;
                        }
						if (fwd==2) 
                        {
                            freq2++; cost2+=ntrans;
                        }
						reset_all(rec_vect, ackfwd, nxtfwdset, ACK, rec_matrix, fwd);
					}
				}
			}
			trans = trans + ntrans;
		}
		messgID ++;
		printf("Send: message ID: %d\n",messgID);
		fprintf(fptr,"%d\t%d\t%d\n", messgID, sfwd, trans);
	}
	printf("\nDone: total messages sent : %d\n",messgID);
	printf("freq1: %d, cost1: %d; freq2: %d, cost2: %d\n",freq1,cost1,freq2,cost2);
	printf("\nPress any key and Enter for exit: ");
	cin>>ext;
}

/*
 * this function calculates the mean cost from nodes (including source) to the destination
 */
void cost2D_est(int source, int destination) 
{/*
	std::deque<int> upnodes;
	int node;

	for (int i=0; i<no_nodes; ++i)
		mean_cost[i]=9999;

	mean_cost[destination] = 0;
	node = destination;
	upnodes.push_back(node);

	while (1)
	{
		float min_cost = 9999;
		for (int i=0; i<no_nodes; ++i)
		{
			if (pe[node][i]<0.94 && mean_cost[i]<mean_cost[node])
			{
				int dwnode = i;
				float cost = mean_link_cost[node][dwnode] + mean_cost[dwnode];
				if (cost<min_cost)
					min_cost = cost;
				mean_cost[node] = min_cost;
			}
		}

		for (int i=0; i<no_nodes; ++i)
		{
			if (pe[i][node]<0.85 && mean_cost[i]>mean_cost[node] && (std::find(upnodes.begin(), upnodes.end(), i) == upnodes.end() || node==destination))
				upnodes.push_back(i);
		}
		upnodes.pop_front();
		if (upnodes.empty())
			break;
		node = upnodes.front();
	}
 */	
	char file_name[50];
	FILE *fp;
	printf("Enter the name of file where costs to destinations are stored\n(format:filename.txt): ");
	cin >> file_name;
	fp=fopen(file_name,"r");
	if(!fp)
	{
		printf("Error:File doesn't exist");
		exit(0);
	}
	else
	{
		for(int i=0;i<no_nodes;i++)
			fscanf(fp,"%f",&mean_cost[i]);
		fclose(fp);
	}
}

/*
 * this class is used for getting the rank of boolean matrices
 */

class BooleanMatrix
{
	vector< vector<bool> > mat; //boolean matrix
	int n, m;           //size of matrix nxm
	int rank;           //rank of the matrix
	
	public:

    /*Constructor
     * Required Parameters:
     * M ==> boolean matrix
     * n ==> number of rows
     * m ==> number of columns
     */
    template <size_t size_m>
    BooleanMatrix(bool M[][size_m], int n, int m)
    {
    	this -> n = n;
    	this -> m = m;
    	for (int i = 0; i < n; i++)
    	{
    		vector<bool> row(m);
    		for (int j = 0; j < m; j++)
    			row[j] = M[i][j];
    		mat.push_back(row);
    	}
    	gaussElimination();
    }

    /* Does Gauss Elimination with partial pivoting on the matrix */
    void gaussElimination()
    {
    	rank = n;
    	for (int i = 0; i < n; i++)
    	{
    		if (!mat[i][i])
    		{
    			int j;
    			for (j = i+1; j < n && !mat[j][i]; j++);
    			if (j == n)
    			{
    				rank--;
                    continue;
                }
                else
                	for (int k = i; k < m; k++)
                	{
                		bool t = mat[i][k];
                		mat[i][k] = mat[j][k];
                		mat[j][k] = t;
                	}
            }
            for (int j = i+1; j < n; j++)
            {
              	if (mat[j][i])
               	{
               		for (int k = i; k < m; k++)
               			mat[j][k] = mat[j][k] - mat[i][k];
               	}
            }
        }
    }

    /* Get the row rank of the boolean matrix
     * If you require the rank of the matrix, make sure that n > m.
     * i.e. if n < m, call the constructor over the transpose.
     */
    int getRank()
    {
        return rank;
    }
};

/*
 * this functions checks if vect is lindearly independent of the cols of matrix
 */
int linr_indep(bool matrix[npacks][npacks][no_nodes], bool vect[npacks], int f, int cols)
{
	int rank, li=0;
	const int m=cols+1;
	printf("node %d : %d\n",f,cols);
	if (cols>=npacks)
		printf("here is the error");
	bool M1[npacks][npacks]={false};
	for (int i=0; i<npacks; ++i)
	{
		int j=0;
		while (j<cols)
		{
			M1[i][j] = matrix[i][j][f];
			++j;
		}
		M1[i][j] = vect[i];
	}

	BooleanMatrix booleanMatrix1(M1, npacks, npacks);
	rank = booleanMatrix1.getRank();
	if (rank==(cols+1))
		li = 1;
	return(li);
}

/*
 * this function estimates the additional number of retransmissions required on the link
 * having erasure probability p to successfully receive all the npacks packets in the message
 */
float addnl_Tx_est(int K, float p, int ntrans, float linccost)
{
	if (ntrans>npacks)
	{
		int M;                      //K: no of packs rec already has; M: no of errors
		M = npacks-K;
		if (M==0)
		{
			printf("message has already been been received");
			exit(0);
	    }
		float E[100]={0};
		E[0]= 1/((1-p)*(1-(pow(2.0,npacks-1)-1)/(pow(2.0,npacks)-1)));  //one erasure
		for (int m=2; m<=M; ++m)	
		{
			float PNm = 1-((pow(2.0,K)-1)/(pow(2.0,npacks)-1));
			E[m-1] =  1/((1-p)*PNm) + E[m-2];
		}
		return E[M-1];
	}
	else 
	{
		float addnl_Tx = linccost-ntrans;
		return addnl_Tx;
	}
}

/*
 * this function calculates the mean costs of the links
 */
void mean_link_cost_calc()
{
	char file_name[50];
	FILE *fp;
	printf("Enter the name of file where mean link costs are stored\n(format:filename.txt): ");
	cin >> file_name;
	fp=fopen(file_name,"r");
	if(!fp)
	{
		printf("Error:File doesn't exist");
		exit(0);
	}
	else
	{
		for(int i=0;i<no_nodes;i++)
		{
			for(int j=0;j<no_nodes;j++)
			{
				fscanf(fp,"%f",&mean_link_cost[i][j]);
			}
		}
		fclose(fp);
	}
} 

/*
 * this function calculates nCk
 */
double nchoosek(int n, int k)
{
	if (k > n) return 0;
	if (k * 2 > n) k = n-k;
	if (k == 0) return 1;

	double result = n;
	for( int i = 2; i <= k; ++i ) 
	{
		result *= (n-i+1);
		result /= i;
	}
	return result;
}


void reset_all(int rec_vect[no_nodes], bool ackfwd[no_nodes], int nxtfwdset[no_nodes], bool ACK[no_nodes], bool rec_matrix[npacks][npacks][no_nodes], int fwd)
{
	for (int i=0; i<no_nodes; ++i)
	{
		ackfwd[i] = false;       // none of the nodes have sent any ACK
		nxtfwdset[i] = 0;        // downstream forwarders aren't chosen yet
		ACK[i] = false;          // none of the nodes have received any ACK
		if (!(pe[fwd][i] < 0.94 && mean_cost[i] < mean_cost[fwd]))
			rec_vect[i] = 0;         // all nodes have zero packets
	}

	for (int k=0; k<no_nodes; ++k)
	{
		if (!(pe[fwd][k] < 0.94 && mean_cost[k] < mean_cost[fwd]))
		{
			for (int i=0; i<npacks; ++i)
			{
				for (int j=0; j<npacks; ++j)
					rec_matrix[i][j][k] = false;  // all nodes have empty buckets
			}
		}
	}
}