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

/*
 ***
 *** returns the number of monochromatic cliques in the graph presented to
 *** it
 ***
 *** graph is stored in row-major order
 *** only checks values above diagonal
 */

int CliqueCountFourD(int *g,
	     int gsize,int i,int j,int precount)
{
    int k;
    int l;
    int m;
    int n;
    int count=0;
    int sgsize = 4;
    
    int newvalue = g[i*gsize + j];
    int prevalue = 1-newvalue;

    int oldcount = 0;
    for(k=0;k < gsize-sgsize+1; k++) 
    {
	if(k != i && k != j && (prevalue == g[minlin(i,k)*gsize+maxlin(i,k)]) && 
	   	(prevalue == g[minlin(j,k)*gsize+maxlin(j,k)]))
	{
	    for(l=k+1;l < gsize-sgsize+2; l++) 
	    { 
		if(l!=i && l != j && (prevalue == g[minlin(i,l)*gsize+maxlin(i,l)]) && 
		   (prevalue == g[minlin(j,l)*gsize+maxlin(j,l)]) && 
		   (prevalue == g[k*gsize+l]))
		{
		    for(m=l+1;m < gsize-sgsize+3; m++) 
		    {
			if(m != i && m!=j && (prevalue == g[minlin(i,m)*gsize+maxlin(i,m)]) && 
			   (prevalue == g[minlin(j,m)*gsize+maxlin(j,m)]) &&
			   (prevalue == g[k*gsize+m]) && 
			   (prevalue == g[l*gsize+m])) {
				for(n=m+1; n < gsize-sgsize+4; n++)
				{
					if(n != i && n!=j && (prevalue
						== g[minlin(i,n)*gsize+maxlin(i,n)]) &&
					   (prevalue 
						== g[minlin(j,n)*gsize+maxlin(j,n)]) &&
					   (prevalue 
						== g[k*gsize+n]) &&
					   (prevalue 
						== g[l*gsize+n]) &&
					   (prevalue 
						== g[m*gsize+n])) {
		      					oldcount++;
					}
				}
			}
		    }
		}
	    }
	}
    }


    int newcount = 0;
    for(k=0;k < gsize-sgsize+1; k++) 
    {
	if(k != i && k != j && (newvalue == g[minlin(i,k)*gsize+maxlin(i,k)]) && 
	   	(newvalue == g[minlin(j,k)*gsize+maxlin(j,k)]))
	{
	    for(l=k+1;l < gsize-sgsize+2; l++) 
	    { 
		if(l!=i && l != j && (newvalue == g[minlin(i,l)*gsize+maxlin(i,l)]) && 
		   (newvalue == g[minlin(j,l)*gsize+maxlin(j,l)]) && 
		   (newvalue == g[k*gsize+l]))
		{
		    for(m=l+1;m < gsize-sgsize+3; m++) 
		    {
			if(m != i && m!=j && (newvalue == g[minlin(i,m)*gsize+maxlin(i,m)]) && 
			   (newvalue == g[minlin(j,m)*gsize+maxlin(j,m)]) &&
			   (newvalue == g[k*gsize+m]) && 
			   (newvalue == g[l*gsize+m])) {
				for(n=m+1; n < gsize-sgsize+4; n++)
				{
					if(n != i && n!=j && (newvalue
						== g[minlin(i,n)*gsize+maxlin(i,n)]) &&
					   (newvalue 
						== g[minlin(j,n)*gsize+maxlin(j,n)]) &&
					   (newvalue 
						== g[k*gsize+n]) &&
					   (newvalue 
						== g[l*gsize+n]) &&
					   (newvalue 
						== g[m*gsize+n])) {
		      					newcount++;
					}
				}
			}
		    }
		}
	    }
	}
    }

    count = precount - oldcount + newcount;
    return(count);
}

/*
 ***
 *** returns the number of monochromatic cliques in the graph presented to
 *** it
 ***
 *** graph is stored in row-major order
 *** only checks values above diagonal
 */

void* CliqueCountLin(int *g,
	     int gsize,void *edge_list)
{
    edge_list = FIFOResetEdge(edge_list);

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
			      				if (!FIFOFindEdge(edge_list,i,j)) FIFOInsertEdge(edge_list,i,j);
			      				if (!FIFOFindEdge(edge_list,i,k)) FIFOInsertEdge(edge_list,i,k);
			      				if (!FIFOFindEdge(edge_list,i,l)) FIFOInsertEdge(edge_list,i,l);
			      				if (!FIFOFindEdge(edge_list,i,m)) FIFOInsertEdge(edge_list,i,m);
			      				if (!FIFOFindEdge(edge_list,i,n)) FIFOInsertEdge(edge_list,i,n);
			      				if (!FIFOFindEdge(edge_list,j,k)) FIFOInsertEdge(edge_list,j,k);
			      				if (!FIFOFindEdge(edge_list,j,l)) FIFOInsertEdge(edge_list,j,l);
			      				if (!FIFOFindEdge(edge_list,j,m)) FIFOInsertEdge(edge_list,j,m);
			      				if (!FIFOFindEdge(edge_list,j,n)) FIFOInsertEdge(edge_list,j,n);
			      				if (!FIFOFindEdge(edge_list,k,l)) FIFOInsertEdge(edge_list,k,l);
			      				if (!FIFOFindEdge(edge_list,k,m)) FIFOInsertEdge(edge_list,k,m);
			      				if (!FIFOFindEdge(edge_list,k,n)) FIFOInsertEdge(edge_list,k,n);
			      				if (!FIFOFindEdge(edge_list,l,m)) FIFOInsertEdge(edge_list,l,m);
			      				if (!FIFOFindEdge(edge_list,l,n)) FIFOInsertEdge(edge_list,l,n);
			      				if (!FIFOFindEdge(edge_list,m,n)) FIFOInsertEdge(edge_list,m,n);
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
    return(edge_list);
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

    int *ret = (int *)malloc((n-k)*sizeof(int));
    int *row = (int *)malloc((n+1)*sizeof(int));
    int *wholerow = (int *)malloc(gsize*sizeof(int));
    int *newret;
    for(i=2;i<n-k+2;i++){
        ret[i-2] = i;
    }
    s = n-k;

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
	                row[j] = 0;

	                if (same >= 5)
	                    continue;
	            }
                news++;
            }
            if (news > 1000000)
            	break;
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
	                row[j] = 0;
	                if (same >= 5)
	                    continue;
				}
                for (x=0;x<d;x++)
                	newret[news*(d+1)+x] = ret[i*d+x];
                newret[news*(d+1)+d] = j;
                news++;
        	}
        	if (news > 1000000)
            	break;
        }
        free(ret);
        ret = newret;
        s = news;
    }
    fprintf(ofp,"(%d,%d) = %d\n",n,k,s);
    printf("(%d,%d) = %d\n",n,k,s);
    for(i=0;i<3;i++){
    	for(j=0;j<k;j++){
			fprintf(ofp,"%d ",ret[k*i+j]);
			printf("%d ",ret[k*i+j]);
    	}
    	fprintf(ofp,"\n");
    	printf("\n");
    }
    ret[0] = s;
    return ret;
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
    int s = ret[0];
    ret[0] = 1; 
    if (1 == 1){
       return 1;
    }
     
	/*
	 * start with graph of size 8
	 */
	
	g = (int *)malloc(gsize*gsize*sizeof(int));
	if(g == NULL) {
		exit(1);
	}
    /*
	 * start out with random
	 */
	
	int try;
	for (try=0;try<s;try++){
		memset(g,0,gsize*gsize*sizeof(int));
		for (i=0;i<gsize;i++)
			g[i] = 0;
		
		for (j=0;j<k;j++)
			g[ret[try*k+j]] = 1;
		for (i=gsize-1;i>n;i--)
			g[i] = g[gsize-i];

		for (i=1;i<gsize;i++){
			for (j=i+1;j<gsize;j++)
				g[i*gsize+j] = g[j-i];
		}

		count = CliqueCount(g,gsize);
		if (try % 1000 == 0)
			fprintf(logfp,"%d %d\n",try,count);
		if (count == 0){
			fprintf(logfp,"Eureka!  Counter-example found! %d\n",try);
			PrintGraph(g,gsize,ofp);
		}
		if (s-try<10){
			for (j=0;j<k;j++){
				fprintf(logfp,"%d ",ret[try*k+j]);
			}
			fprintf(logfp,"\n");
		}
	}
	fclose(ofp);
	fclose(logfp);
	return(0);

}
