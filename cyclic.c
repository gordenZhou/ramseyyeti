#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>

#include "fifo.h"	/* for taboo list */


#define MAXSIZE (512)

#define TABOOSIZE (1000)
#define BIGCOUNT (9999999)

/***
 *** example of very simple search for R(6,6) counter examples
 ***
 *** starts with a small randomized graph and works its way up to successively
 *** larger graphs one at a time
 ***
 *** uses a taboo list of size #TABOOSIZE# to hold and encoding of and edge
 *** (i,j)+clique_count
 ***/

/*
 * PrintGraph
 *
 * prints in the right format for the read routine
 */
void PrintGraph(int *g, int gsize,FILE *ofp)
{
	int i;
	int j;

	fprintf(ofp,"%d\n",gsize);

	for(i=0; i < gsize; i++)
	{
		for(j=0; j < gsize; j++)
		{
			if (j>=i) fprintf(ofp,"%d",g[i*gsize+j]);
            else fprintf(ofp,"%d",g[i*gsize+j]);
		}
		fprintf(ofp,"\n");
	}

	return;
}

/*
 * CopyGraph 
 *
 * copys the contents of old_g to corresponding locations in new_g
 * leaving other locations in new_g alone
 * that is
 * 	new_g[i,j] = old_g[i,j]
 */
void CopyGraph(int *old_g, int o_gsize, int *new_g, int n_gsize)
{
	int i;
	int j;

	/*
	 * new g must be bigger
	 */
	if(n_gsize < o_gsize)
		return;

	for(i=0; i < o_gsize; i++)
	{
		for(j=0; j < o_gsize; j++)
		{
			new_g[i*n_gsize+j] = old_g[i*o_gsize+j];
		}
	}

	return;
}

int minlin(int a,int b){
	if (a<b)
		return a;
	return b;
}

int maxlin(int a,int b){
	if (a>b)
		return a;
	return b;
}

void outputtofile(char* ostr){
    FILE* ofp = fopen("../log/log.txt","a");
    fprintf("%s",ostr);
    fclose(ofp);
}
/*
 *** Read in a current solution from file
 ***
 ***
 */
int* readinsolution(int gsize){
    int* g = (int *)malloc(gsize*gsize*sizeof(int));
    memset(g,0,gsize*gsize*sizeof(int));
    
    FILE* ifp = fopen("./98.txt","r");
    int val = 0;
    int i;
    for (i=0;i<gsize*gsize;i++){
        fscanf(ifp,"%d",&val);
        g[i] = val;
    }
    fclose(ifp);
    return g;
}
/*
 ***
 *** returns the number of monochromatic cliques in the graph presented to
 *** it
 ***
 *** graph is stored in row-major order
 *** only checks values above diagonal
 */

int CliqueCount(int *g,
	     int gsize)
{
    int i;
    int j;
    int k;
    int l;
    int m;
    int n;
    int count=0;
    int sgsize = 6;
    
    for(i=0;i < gsize-sgsize+1; i++)
    {
	for(j=i+1;j < gsize-sgsize+2; j++)
        {
	    for(k=j+1;k < gsize-sgsize+3; k++) 
            { 
		if((g[i*gsize+j] == g[i*gsize+k]) && 
		   (g[i*gsize+j] == g[j*gsize+k]))
		{
		    for(l=k+1;l < gsize-sgsize+4; l++) 
		    { 
			if((g[i*gsize+j] == g[i*gsize+l]) && 
			   (g[i*gsize+j] == g[j*gsize+l]) && 
			   (g[i*gsize+j] == g[k*gsize+l]))
			{
			    for(m=l+1;m < gsize-sgsize+5; m++) 
			    {
				if((g[i*gsize+j] == g[i*gsize+m]) && 
				   (g[i*gsize+j] == g[j*gsize+m]) &&
				   (g[i*gsize+j] == g[k*gsize+m]) && 
				   (g[i*gsize+j] == g[l*gsize+m])) {
					for(n=m+1; n < gsize-sgsize+6; n++)
					{
						if((g[i*gsize+j]
							== g[i*gsize+n]) &&
						   (g[i*gsize+j] 
							== g[j*gsize+n]) &&
						   (g[i*gsize+j] 
							== g[k*gsize+n]) &&
						   (g[i*gsize+j] 
							== g[l*gsize+n]) &&
						   (g[i*gsize+j] 
							== g[m*gsize+n])) {
			      					count++;
			      				return count;
						}
					}
				}
			    }
			}
		    }
		}
	    }
         }
     }
    return(count);
}
void randomgenrow(int gsize){
    int *r = (int *)malloc(gsize*sizeof(int));
    int ones = gsize/2;
    int zeros = gsize - ones;
    int j = 0;
    int k = 0;
    int tmp;
    int same;
    int a = 1;
    int step = 0;
    int tt = 0;
    FILE* logfp = fopen("../log/cycliclog.txt","a");
    while (a == 1){

        r[0] = 0;
        ones = gsize/2;
        zeros = gsize - ones;
        j = 1;
        
        while (j < gsize/2){
            if (zeros == 0){
                r[j] = 1;
                r[gsize-j] = 1;
                ones--;
                if (j != gsize - j) ones--;
                j++;
                continue;
            }else if (ones == 0){
                r[j] = 0;
                r[gsize-j] = 0;
                if (j != gsize -j) zeros--;
                j++;
                continue;
            }
            
            tmp = rand() % 2;
            if (tmp == 0){
                r[j] = 0;
                r[gsize-j] = 0;
                zeros--;
            }else{
                r[j] = 1;
                r[gsize-j] = 1;
                ones--;
            }
            j++;
        }
        for (step = 1;step < gsize/5;step++){
            same = 1;
            for (k=1;k<5;k++){
                if (r[k*step] == r[(k+1)*step])
                    same++;
            }
            if (same >= 5)
                break;
        }
        if (same >= 5 && tt %100 == 0){
            fprintf(logfp,"same: %d, step: %d tt: %d\n",same,step,tt);
        }
        tt++;
    }
    for (j=0;j<gsize;j++)
        fprintf(logfp,"%d ",r[j]);
    fclose(logfp);
    free(r);
}
int* iterategenrowmid(int gsize, FILE* ofp,int* mid,int startd){
	//fprintf(ofp,"Start Mid %d\n",startd);
    int n = gsize/2;
    int k = n/2;

    int s=0,news=0;
    int f=0;

    int end;

    int i,j,d,ii;
    int same;

    int *ret = (int *)malloc(startd*sizeof(int));
    int *row = (int *)malloc((n+1)*sizeof(int));
    int *wholerow = (int *)malloc(gsize*sizeof(int));
    int *newret;
    s = 1;
    for(i=0;i<startd;i++)
        ret[i] = mid[i];

    int x,step,kk;
    for (d=startd;d<k;d++){
        news = 0;

        for (i=0;i<s;i++){
        	f = ret[i*d+d-1];
            
            for (x=0;x<n+1;x++)
                row[x] = 0;
        	for (x=0;x<d;x++)
        		row[ret[i*d+x]] = 1;

        	if (d == k-1){
	        	for(kk=0;kk<=n;kk++)
					wholerow[kk] = row[kk];
				for(kk=gsize-1;kk>n;kk--)
					wholerow[kk] = row[gsize-kk];
			}

            same = 0;
            for(j=f+1;j<n-k+2+d;j++){
				row[j] = 1;
				if (d == k-1){
					wholerow[j] = row[j];
					wholerow[gsize-j] = row[j];

					end = gsize/5;
					for(step=1;step<=end;step++){
						same = 1;
						for(kk=1;kk<5;kk++)
							if (wholerow[kk*step] == wholerow[(kk+1) * step])
								same++;
						if (same >= 5)
							break;

						same = 1;
	                    for (kk=1;kk<4;kk++)
	                    	if (wholerow[kk*step] == wholerow[(kk+1)*step])
	                            same++;
	                	if (same >= 4)
	                		for (ii=1;ii<=gsize-1-4*step;ii++)
	                			if (wholerow[ii] == wholerow[step] && wholerow[ii+step] == wholerow[step] && wholerow[ii+2*step] == wholerow[step]
	                				&& wholerow[ii+3*step] == wholerow[step] && wholerow[ii+4*step] == wholerow[step]){
	                				same = 5;
	                				break;
	                			}
	                	if (same >= 5)
	                    	break;
					}
					
					row[j] = 0;
					wholerow[j] = row[j];
					wholerow[gsize-j] = row[j];
					if (same>=5)
						continue;
				}
				else{
					end = j/5;

	                for(step=1;step <= end;step++){
	                    same = 1;
	                    for(kk=1;kk<5;kk++)
	                        if (row[kk*step] == row[(kk+1)*step])
	                            same++;
	                    if (same >= 5)
	                    	break;

	                    same = 1;
	                    for (kk=1;kk<4;kk++)
	                    	if (row[kk*step] == row[(kk+1)*step])
	                            same++;
	                	if (same >= 4)
	                		for (ii=1;ii<=j-4*step;ii++)
	                			if (row[ii] == row[step] && row[ii+step] == row[step] && row[ii+2*step] == row[step]
	                				&& row[ii+3*step] == row[step] && row[ii+4*step] == row[step]){
	                				same = 5;
	                				break;
	                			}
	                	if (same >= 5)
	                    	break;
	                }
	                step = 0;
					for (kk=0;kk<j;kk++){
						if (row[kk] == row[kk+1]){
							step++;
							if (step > 8)
								break;
						}
						else
							step = 1;			
					}
	                row[j] = 0;
	                if (same >= 5 || step > 7)
	                    continue;
	            }
                news++;
            }
            //if (news > 10000000)
            //	break;
        }
        newret = (int *)malloc((news)*(d+1)*sizeof(int));
        news = 0;

        for (i=0;i<s;i++){
        	f = ret[i*d+d-1];
            
            for (x=0;x<n+1;x++)
                row[x] = 0;
        	for (x=0;x<d;x++)
        		row[ret[i*d+x]] = 1;

        	if (d == k-1){
	        	for(kk=0;kk<=n;kk++)
					wholerow[kk] = row[kk];
				for(kk=gsize-1;kk>n;kk--)
					wholerow[kk] = row[gsize-kk];
			}
            same = 0;
            for(j=f+1;j<n-k+2+d;j++){
            	
				row[j] = 1;
				if (d == k-1){
					wholerow[j] = row[j];
					wholerow[gsize-j] = row[j];

					end = gsize/5;
					for(step=1;step<=end;step++){
						same = 1;
						for(kk=1;kk<5;kk++)
							if (wholerow[kk*step] == wholerow[(kk+1) * step])
								same++;
						if (same >= 5)
							break;


						same = 1;
	                    for (kk=1;kk<4;kk++)
	                    	if (wholerow[kk*step] == wholerow[(kk+1)*step])
	                            same++;
	                	if (same >= 4)
	                		for (ii=1;ii<=gsize-1-4*step;ii++)
	                			if (wholerow[ii] == wholerow[step] && wholerow[ii+step] == wholerow[step] && wholerow[ii+2*step] == wholerow[step]
	                				&& wholerow[ii+3*step] == wholerow[step] && wholerow[ii+4*step] == wholerow[step]){
	                				same = 5;
	                				break;
	                			}
	                	if (same >= 5)
	                    	break;
					}
					row[j] = 0;
					wholerow[j] = row[j];
					wholerow[gsize-j] = row[j];
					if (same>=5)
						continue;
				}
				else{
					end = j/5;
	                for(step=1;step <= end;step++){

	                    same = 1;
	                    for(kk=1;kk<5;kk++)
	                        if (row[kk*step] == row[(kk+1)*step])
	                            same++;
	                    if (same >= 5)
	                    	break;

	                    same = 1;
	                    for (kk=1;kk<4;kk++)
	                    	if (row[kk*step] == row[(kk+1)*step])
	                            same++;
	                	if (same >= 4)
	                		for (ii=1;ii<=j-4*step;ii++)
	                			if (row[ii] == row[step] && row[ii+step] == row[step] && row[ii+2*step] == row[step]
	                				&& row[ii+3*step] == row[step] && row[ii+4*step] == row[step]){
	                				same = 5;
	                				break;
	                			}
	                	if (same >= 5)
	                    	break;
	                }
	                step = 0;
					for (kk=0;kk<j;kk++){
						if (row[kk] == row[kk+1]){
							step++;
							if (step > 8)
								break;
						}
						else
							step = 1;			
					}
	                row[j] = 0;
	                if (same >= 5 || step > 7)
	                    continue;
				}
                for (x=0;x<d;x++)
                	newret[news*(d+1)+x] = ret[i*d+x];
                newret[news*(d+1)+d] = j;
                news++;
        	}
        	//if (news > 10000000)
            //	break;
        }
        free(ret);
        ret = newret;
        s = news;
        //printf("d:%d s:%d\n",d,s);
        //fflush(stdout);
    }
    /*fprintf(ofp,"(%d,%d) = %d\n",n,k,s);
    printf("(%d,%d) = %d\n",n,k,s);
    for(i=0;i<3;i++){
    	for(j=0;j<d;j++){
			fprintf(ofp,"%d ",ret[d*i+j]);
			printf("%d ",ret[d*i+j]);
    	}
    	fprintf(ofp,"\n");
    	printf("\n");
    }*/
    ret[0] = s;
    return ret;
}

int* iterategenrow(int gsize,FILE* ofp){
	// set a cap for the news, then we can have different range of combinations.
	fprintf(ofp,"Start Iterating %d\n",gsize);
    int n = gsize/2;
    int k = n/2;

    int s=0,news=0;
    int f=0;

    int end;

    int i,j,d,ii;
    int same;

    
    int *ret = (int *)malloc(1*sizeof(int));
    int *row = (int *)malloc((n+1)*sizeof(int));
    int *newret;
    s = 0;
    for(i=5;i<6;i++){
        ret[s] = i;
        s++;
    }

    int x,step,kk;
    for (d=1;d<k;d++){
        //newret = []
        news = 0;

        for (i=0;i<s;i++){
        	f = ret[i*d+d-1];
            
            for (x=0;x<n+1;x++)
                row[x] = 0;
        	for (x=0;x<d;x++)
        		row[ret[i*d+x]] = 1;

            same = 0;
            for(j=f+1;j<n-k+2+d;j++){
				row[j] = 1;
				end = j/5;
                for(step=1;step <= end;step++){
                    same = 1;
                    for(kk=1;kk<5;kk++)
                        if (row[kk*step] == row[(kk+1)*step])
                            same++;
                    if (same >= 5)
                    	break;

                    same = 1;
                    for (kk=1;kk<4;kk++)
                    	if (row[kk*step] == row[(kk+1)*step])
                            same++;
                	if (same >= 4)
                		for (ii=1;ii<=j-4*step;ii++)
                			if (row[ii] == row[step] && row[ii+step] == row[step] && row[ii+2*step] == row[step]
                				&& row[ii+3*step] == row[step] && row[ii+4*step] == row[step]){
                				same = 5;
                				break;
                			}
                	if (same >= 5)
                    	break;
                }
                step = 0;
				for (kk=0;kk<j;kk++){
					if (row[kk] == row[kk+1]){
						step++;
						if (step > 8)
							break;
					}
					else
						step = 1;			
				}
                row[j] = 0;
                if (same >= 5 || step > 7)
                    continue;
            
                news++;
            }

            //if (news > 10000000)
            //	break;
        }
        newret = (int *)malloc((news)*(d+1)*sizeof(int));
        news = 0;

        for (i=0;i<s;i++){
        	f = ret[i*d+d-1];
            
            for (x=0;x<n+1;x++)
                row[x] = 0;
        	for (x=0;x<d;x++)
        		row[ret[i*d+x]] = 1;
            same = 0;

            for(j=f+1;j<n-k+2+d;j++){
            	
				row[j] = 1;
				
				end = j/5;
                for(step=1;step <= end;step++){

                    same = 1;
                    for(kk=1;kk<5;kk++)
                        if (row[kk*step] == row[(kk+1)*step])
                            same++;
                    if (same >= 5)
                    	break;

                    same = 1;
                    for (kk=1;kk<4;kk++)
                    	if (row[kk*step] == row[(kk+1)*step])
                            same++;
                	if (same >= 4)
                		for (ii=1;ii<=j-4*step;ii++)
                			if (row[ii] == row[step] && row[ii+step] == row[step] && row[ii+2*step] == row[step]
                				&& row[ii+3*step] == row[step] && row[ii+4*step] == row[step]){
                				same = 5;
                				break;
                			}
                	if (same >= 5)
                    	break;
                }
                step = 0;
				for (kk=0;kk<j;kk++){
					if (row[kk] == row[kk+1]){
						step++;
						if (step > 8)
							break;
					}
					else
						step = 1;			
				}
                row[j] = 0;
                if (same >= 5 || step > 7)
                    continue;

                for (x=0;x<d;x++)
                	newret[news*(d+1)+x] = ret[i*d+x];
                newret[news*(d+1)+d] = j;
                news++;
        	}
        	//if (news > 10000000)
            //	break;
        }
        free(ret);
        ret = newret;
        s = news;
        fprintf(ofp,"d:%d s:%d\n",d,s);
        if (s>10000000)
        	break;
        //fflush(stdout);
    }
    d++;
    fprintf(ofp,"(%d,%d) = %d\n",n,k,s);
    //printf("(%d,%d) = %d\n",n,k,s);
    for(i=0;i<3;i++){
    	for(j=0;j<d;j++){
			fprintf(ofp,"%d ",ret[d*i+j]);
			//printf("%d ",ret[d*i+j]);
    	}
    	fprintf(ofp,"\n");
    	//printf("\n");
    }
    fflush(ofp);
    ret[0] = s;
    ret[1] = d;
    return ret;
}



int main(int argc,char *argv[])
{
	int *g;
	int gsize;
	int count=0;
	int i,j;
	
    int k;
    
    gsize = 99;
    int n = gsize/2;
    k = n/2;

    //randomgenrow(99);
    FILE *logfp = fopen("../log/cycliclog.txt","a");
	FILE *ofp = fopen("../data/cyclicr66.txt","a");
    int* ret = iterategenrow(gsize,logfp);

    int* finalret;
    int s = ret[0];
    int d = ret[1];
    ret[0] = 5;
    ret[1] = 6;
    if (1 == 2){
       return 1;
    }
     
	/*
	 * start with graph of size 8
	 */
	g = (int *)malloc(gsize*gsize*sizeof(int));
	if(g == NULL) {
		exit(1);
	}
	int si,try,finals;
	int* mid = (int *)malloc(d*sizeof(int));

	for (si=0;si<s;si++){
		for (i=0;i<d;i++)
			mid[i] = ret[si*d+i];

		finalret = iterategenrowmid(gsize,logfp,mid,d);
		finals = finalret[0];
		finalret[0] = 5;

		if (si % 1000 == 0){
			fprintf(logfp,"si/s: %d/%d\n",si,s);
			fflush(logfp);
		}	
			
		for (try=0;try<finals;try++){
			memset(g,0,gsize*gsize*sizeof(int));
			for (i=0;i<gsize;i++)
				g[i] = 0;
			
			for (j=0;j<k;j++)
				g[finalret[try*k+j]] = 1;
			for (i=gsize-1;i>n;i--)
				g[i] = g[gsize-i];

			for (i=1;i<gsize;i++){
				for (j=i+1;j<gsize;j++)
					g[i*gsize+j] = g[j-i];
			}

			count = CliqueCount(g,gsize);
			if (count == 0){
				fprintf(logfp,"Eureka!  Counter-example found! si/s: %d/%d try/finals:%d/%d\n",si,s,try,finals);
				fflush(logfp);
				PrintGraph(g,gsize,ofp);
			}
		}
		free(finalret);
	}
	free(ret);
	fclose(ofp);
	fclose(logfp);
	return(0);

}
