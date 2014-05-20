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
            else fprintf(ofp,"%d",g[j*gsize+i]);
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
    FILE* ofp = fopen("./log.txt","a");
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

void iterategenrow(int gsize){
    int n = gsize/2;
    int k = n/2;
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
                zeros--;
            }else{
                r[j] = 1;
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
        if (same >= 5){
            printf("same: %d, step: %d\n",same,step);
        }
    }
    for (j=0;j<gsize;j++)
        printf("%d ",r[j]);
    free(r);
}

int main(int argc,char *argv[])
{
    srand(time(NULL));
	int *g;
	int *new_g;
	int gsize;
	int count;
	int i;
	int j;
	int best_count;
	int best_i;
	int best_j;
	void *taboo_list;
    void *edge_list;
	int psen = 0;
	int precount = 0;
	int checkcount = 0;
    int ffct = 0;
    int rbk = 0;
    int ini;
    int lottery;
    int k;
    
    int* carray = (int *)malloc(50*sizeof(int));
    int* iarray = (int *)malloc(50*sizeof(int));
    int* jarray = (int *)malloc(50*sizeof(int));
    
    randomgenrow(99);
    if (1 == 1)
        return 1;
	/*
	 * start with graph of size 8
	 */
	gsize = 99;
	g = (int *)malloc(gsize*gsize*sizeof(int));
	if(g == NULL) {
		exit(1);
	}
    /*
	 * start out with random
	 */
	memset(g,0,gsize*gsize*sizeof(int));
    for (i=0;i<gsize;i++){
        for(j=i+1;j<gsize;j++)
            g[i*gsize+j] = rand()%2;
    }

    /*
     *  start from a file
     */
    //gsize = 98;
    //g = readinsolution(98);

	/*
	 * make a fifo to use as the taboo list
	 */
    
	taboo_list = FIFOInitEdge(TABOOSIZE);
	if(taboo_list == NULL) {
		exit(1);
	}
    edge_list = FIFOInitEdge(gsize*gsize/2);
    if(edge_list == NULL) {
		exit(1);
	}
	
    FILE *ofp = fopen("./r66.txt","w");

	/*
	 * while we do not have a publishable result
	 */
	while(gsize < 101)
	{
		/*
		 * find out how we are doing
		 */
		count = CliqueCount(g,gsize);
		precount = count;

		/*
		 * if we have a counter example
		 */
		if(count == 0)
		{
			printf("Eureka!  Counter-example found!\n");
			PrintGraph(g,gsize,ofp);
			/*
			 * make a new graph one size bigger
			 */
			new_g = (int *)malloc((gsize+1)*(gsize+1)*sizeof(int));
			if(new_g == NULL){
                fclose(ofp);
				exit(1);
            }
			/*
			 * copy the old graph into the new graph leaving the
			 * last row and last column alone
			 */
			CopyGraph(g,gsize,new_g,gsize+1);

			/*
			 * zero out the last column and last row
			 */
			for(i=0; i < (gsize+1); i++)
			{
				new_g[i*(gsize+1) + gsize] = rand() % 2; // last column
				new_g[gsize*(gsize+1) + i] = 0; // last row
			}

			/*
			 * throw away the old graph and make new one the
			 * graph
			 */
			free(g);
			g = new_g;
			gsize = gsize+1;
            rbk = 0;

			/*
			 * reset the taboo list for the new graph
			 */
			taboo_list = FIFOResetEdge(taboo_list);

			/*
			 * keep going
			 */
            break;
			//continue;
		}

		/*
		 * otherwise, we need to consider flipping an edge
		 *
		 * let's speculative flip each edge, record the new count,
		 * and unflip the edge.  We'll then remember the best flip and
		 * keep it next time around
		 *
		 * only need to work with upper triangle of matrix =>
		 * notice the indices
		 */
		best_count = precount;//BIGCOUNT;

		edge_list = CliqueCountLin(g,gsize,edge_list);
		psen = 0;
        for(i=0;i<50;i++){
            carray[i] = BIGCOUNT;
            iarray[i] = -1;
            jarray[i] = -1;
        }
		for(i=0; i < gsize; i++)
		{
			for(j=i+1; j < gsize; j++)
			{
				if (!FIFOFindEdge(edge_list,i,j) || FIFOFindEdge(taboo_list,i,j))
					continue;

				//FIFODeleteEdge(edge_list,i,j);
				/*
				 * flip it
				 */
				g[i*gsize+j] = 1 - g[i*gsize+j];
				//count = CliqueCount(g,gsize);
				count = CliqueCountFourD(g,gsize,i,j,precount);
				//checkcount = CliqueCount(g,gsize);
				//if (count != checkcount){
				//	printf("count: %d checkcount: %d\n", count,checkcount);
				//	exit(1);
				//}
				psen++;
				/*
				 * is it better and the i,j,count not taboo?
				 */
                ini = 49;
                for(;ini>=0;ini--){
                    if (count < carray[ini]){
                        if (carray[ini] < BIGCOUNT && ini<49){
                            carray[ini+1] = carray[ini];
                            iarray[ini+1] = iarray[ini];
                            jarray[ini+1] = jarray[ini];
                        }
                    }else{
                        break;
                    }
                }
                if (ini<49){
                    carray[ini+1] = count;
                    iarray[ini+1] = i;
                    jarray[ini+1] = j;
                }
 				/*if((count < best_count) &&
					!FIFOFindEdge(taboo_list,i,j))
//					!FIFOFindEdgeCount(taboo_list,i,j,count))
				{
					best_count = count;
					best_i = i;
					best_j = j;
				}*/
                

				/*
				 * flip it back
				 */
				g[i*gsize+j] = 1 - g[i*gsize+j];
			}
		}
        if (carray[0] < best_count){
            ini = 0;
            best_count = carray[0];
            best_i = iarray[0];
            best_j = jarray[0];
        }else{
            lottery = 0;
            k = 50;
            for (ini = 0;ini<50;ini++)
                if (carray[ini] < BIGCOUNT)
                    lottery = lottery + 10;
            lottery = rand() % lottery;
            for (ini = 0;ini<50;ini++){
                lottery = lottery-10;
                k--;
                if (lottery <=0){
                    best_count = carray[ini];
                    best_i = iarray[ini];
                    best_j = jarray[ini];
                    break;
                }
            }
        }

		if(best_count == BIGCOUNT) {
			printf("no best edge found, terminating\n");
            fclose(ofp);
			exit(1);
		}
		
		/*
		 * keep the best flip we saw
		 */
		g[best_i*gsize+best_j] = 1 - g[best_i*gsize+best_j];

		/*
		 * taboo this graph configuration so that we don't visit
		 * it again
		 */
		//count = CliqueCount(g,gsize);
		FIFOInsertEdge(taboo_list,best_i,best_j);
//		FIFOInsertEdgeCount(taboo_list,best_i,best_j,count);
        ffct = FIFOCount(taboo_list);
        
		if (ffct % 10 == 1 || (ffct >= gsize * maxlin(6,gsize/10) && psen >= gsize*3)){
            FILE* logfp = fopen("./log.txt","a");
            fprintf(logfp,"ce sz: %d, b_ct: %d, b_eg: (%d,%d), new c: %d psen: %d q sz:%d ini:%d\n",
			gsize,
			best_count,
			best_i,
			best_j,
			g[best_i*gsize+best_j],
			psen,
			ffct,
            ini);
            fclose(logfp);
        }
        if (ini > 0){
            outputtofile("Reset 100 Edges\n");
            fflush(stdout);
            for (ini=0;ini<100;ini++){
                i = rand()%gsize;
                j = rand()%gsize;
                if (i<=j)
                    g[i*gsize+j] = 1 - g[i*gsize+j];
                else
                    g[j*gsize+i] = 1 - g[j*gsize+i];
            }
            taboo_list = FIFOResetEdge(taboo_list);
            PrintGraph(g,gsize,ofp);
        }
        
		if ((1==0) && (ffct >= gsize * maxlin(6,gsize/10) && psen >= gsize*3)){
			// Roll Back all the changes
			for(i=0; i < gsize; i++)
			{
				for(j=i+1; j < gsize; j++)
				{
					if (FIFOFindEdge(taboo_list,i,j))
						g[i*gsize+j] = 1-g[i*gsize+j];

				}
			}
			taboo_list = FIFOResetEdge(taboo_list);
			// Regenerate from n to n + 1
			for(i=0; i < (gsize); i++)
			{
				new_g[i*(gsize) + gsize-1] = rand() % 2; // last column
				new_g[(gsize-1)*gsize + i] = 0; // last row
			}
            rbk++;
            if (rbk % 10 == 0)
                printf("--Regenerate: %d\n",rbk);
		}
		/*
		 * rinse and repeat
		 */
	}

	FIFODeleteGraph(taboo_list);
    fclose(ofp);

	return(0);

}
