/*
 #  File            : condat_fast_tv.c 
 #
 #  Version History : 1.0, Feb. 2012
 #
 #  Author          : Laurent Condat, PhD, CNRS research fellow in France.
 #
 #  Description     : This file contains an implementation in the C language
 #                    of algorithms described in the research paper:
 #	
 #                    L. Condat, "A Direct Algorithm for 1D Total Variation
 #                    Denoising", preprint hal-00675043, 2012.
 #
 #                    This implementation comes with no warranty: due to the
 #                    limited number of tests performed, there may remain
 #                    bugs. In case the functions would not do what they are
 #                    supposed to do, please email the author (contact info
 #                    to be found on the web).
 #
 #                    If you use this code or parts of it for any purpose,
 #                    the author asks you to cite the paper above or, in 
 #                    that event, its published version. Please email him if 
 #                    the proposed algorithms were useful for one of your 
 #                    projects, or for any comment or suggestion.
 #
 #  Usage rights    : Copyright Laurent Condat.
 #                    This file is distributed under the terms of the CeCILL
 #                    licence (compatible with the GNU GPL), which can be
 #                    found at the URL "http://www.cecill.info".
 #
 #  This software is governed by the CeCILL license under French law and
 #  abiding by the rules of distribution of free software. You can  use,
 #  modify and or redistribute the software under the terms of the CeCILL
 #  license as circulated by CEA, CNRS and INRIA at the following URL :
 #  "http://www.cecill.info".
 #
 #  As a counterpart to the access to the source code and rights to copy,
 #  modify and redistribute granted by the license, users are provided only
 #  with a limited warranty  and the software's author,  the holder of the
 #  economic rights,  and the successive licensors  have only  limited
 #  liability.
 #
 #  In this respect, the user's attention is drawn to the risks associated
 #  with loading,  using,  modifying and/or developing or reproducing the
 #  software by the user in light of its specific status of free software,
 #  that may mean  that it is complicated to manipulate,  and  that  also
 #  therefore means  that it is reserved for developers  and  experienced
 #  professionals having in-depth computer knowledge. Users are therefore
 #  encouraged to load and test the software's suitability as regards their
 #  requirements in conditions enabling the security of their systems and/or
 #  data to be ensured and,  more generally, to use and operate it in the
 #  same conditions as regards security.
 #
 #  The fact that you are presently reading this means that you have had
 #  knowledge of the CeCILL license and that you accept its terms.
*/


#include <math.h>
#include <time.h>			
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


/* 
This function implements the 1D total variation denoising 
algorithm described in the paper referenced above. 
If output=input, the process is performed in place. Else, 
the values of input are left unchanged. 
lambda must be nonnegative. lambda=0 is admissible and 
yields output[k]=input[k] for all k. 
If width<=0, nothing is done. 
*/
void TV1D_denoise(float* input, float* output, const int width, const float lambda) {
	if (width>0) {				/*to avoid invalid memory access to input[0]*/
		int k=0, k0=0;			/*k: current sample location, k0: beginning of current segment*/
		float umin=lambda, umax=-lambda;	/*u is the dual variable*/
		float vmin=input[0]-lambda, vmax=input[0]+lambda;	/*bounds for the segment's value*/
		int kplus=0, kminus=0; 	/*last positions where umax=-lambda, umin=lambda, respectively*/
		const float twolambda=2.0*lambda;	/*auxiliary variable*/
		const float minlambda=-lambda;		/*auxiliary variable*/
		for (;;) {				/*simple loop, the exit test is inside*/
			while (k==width-1) {	/*we use the right boundary condition*/
				if (umin<0.0) {			/*vmin is too high -> negative jump necessary*/
					do output[k0++]=vmin; while (k0<=kminus);
					umax=(vmin=input[kminus=k=k0])+(umin=lambda)-vmax;
				} else if (umax>0.0) {	/*vmax is too low -> positive jump necessary*/
					do output[k0++]=vmax; while (k0<=kplus);
					umin=(vmax=input[kplus=k=k0])+(umax=minlambda)-vmin;
				} else {
					vmin+=umin/(k-k0+1); 
					do output[k0++]=vmin; while(k0<=k); 
					return;
				}
			}
			if ((umin+=input[k+1]-vmin)<minlambda) {		/*negative jump necessary*/
				do output[k0++]=vmin; while (k0<=kminus);
				vmax=(vmin=input[kplus=kminus=k=k0])+twolambda;
				umin=lambda; umax=minlambda;
			} else if ((umax+=input[k+1]-vmax)>lambda) {	/*positive jump necessary*/
				do output[k0++]=vmax; while (k0<=kplus);
				vmin=(vmax=input[kplus=kminus=k=k0])-twolambda;
				umin=lambda; umax=minlambda;
			} else { 	/*no jump necessary, we continue*/
				k++;
				if (umin>=lambda) {		/*update of vmin*/
					vmin+=(umin-lambda)/((kminus=k)-k0+1);
					umin=lambda;
				} 
				if (umax<=minlambda) {	/*update of vmax*/
					vmax+=(umax+lambda)/((kplus=k)-k0+1);
					umax=minlambda;
				} 	
			}
		}
	}
}


/* 
This function implements the fused lasso signal approximator,
described in the paper referenced above. 
If output=input, the process is performed in place. Else, the 
values of input are left unchanged. 
lambda and mu must be nonnegative. If mu=0, the function does 
the same as TV1D_denoise.
If width<=0, nothing is done. 
*/
void fused_lasso(float* input, float* output, const int width, const float lambda, const float mu) {
	if (width>0) {				/*to avoid invalid memory access to input[0]*/
		int k=0, k0=0;			/*k: current sample location, k0: beginning of current segment*/
		float umin=lambda, umax=-lambda;	/*u is the dual variable*/
		float vmin=input[0]-lambda, vmax=input[0]+lambda;	/*bounds for the segment's value*/
		int kplus=0, kminus=0; 	/*last positions where umax=-lambda, umin=lambda, respectively*/
		const float twolambda=2.0*lambda;	/*auxiliary variable*/
		const float minlambda=-lambda;		/*auxiliary variable*/
		for (;;) {				/*simple loop, the exit test is inside*/
			while (k==width-1) {	/*we use the right boundary condition*/
				if (umin<0.0) {			/*vmin is too high -> negative jump necessary*/
					vmin = vmin>mu ? vmin-mu : vmin<-mu ? vmin+mu : 0.0;  
					do output[k0++]=vmin; while (k0<=kminus);
					umax=(vmin=input[kminus=k=k0])+(umin=lambda)-vmax;
				} else if (umax>0.0) {	/*vmax is too low -> positive jump necessary*/
					vmax = vmax>mu ? vmax-mu : vmax<-mu ? vmax+mu : 0.0;
					do output[k0++]=vmax; while (k0<=kplus);
					umin=(vmax=input[kplus=k=k0])+(umax=minlambda)-vmin;
				} else {
					vmin+=umin/(k-k0+1);
					vmin = vmin>mu ? vmin-mu : vmin<-mu ? vmin+mu : 0.0;  
					do output[k0++]=vmin; while(k0<=k); 
					return;
				}
			}
			if ((umin+=input[k+1]-vmin)<minlambda) {		/*negative jump necessary*/
				vmin = vmin>mu ? vmin-mu : vmin<-mu ? vmin+mu : 0.0; 
				do output[k0++]=vmin; while (k0<=kminus);
				vmax=(vmin=input[kplus=kminus=k=k0])+twolambda;
				umin=lambda; umax=minlambda;
			} else if ((umax+=input[k+1]-vmax)>lambda) {	/*positive jump necessary*/
				vmax = vmax>mu ? vmax-mu : vmax<-mu ? vmax+mu : 0.0;
				do output[k0++]=vmax; while (k0<=kplus);
				vmin=(vmax=input[kplus=kminus=k=k0])-twolambda;
				umin=lambda; umax=minlambda;
			} else { 	/*no jump necessary, we continue*/
				k++;
				if (umin>=lambda) {		/*update of vmin*/
					vmin+=(umin-lambda)/((kminus=k)-k0+1);
					umin=lambda;
				} 
				if (umax<=minlambda) {	/*update of vmax*/
					vmax+=(umax+lambda)/((kplus=k)-k0+1);
					umax=minlambda;
				} 	
			}
		}
	}
}


#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define floattype double
/* 
This function implements 1D total variation denoising by
the taut string algorithm. It was adapted from the Matlab 
code written by Lutz Duembgen, which can be found at 
http://www.imsv.unibe.ch/content/staff/personalhomepages/
duembgen/software/multiscale_densities/index_eng.html
Note: numerical blow-up if floattype=float for N>=10^6
because the algorithm is based on the running sum of the 
signal values.
*/
void TV1D_denoise_tautstring(float* input, float* output, int width, const float lambda) {
	width++;
	int* index_low=calloc(width,sizeof(int));
	floattype* slope_low=calloc(width,sizeof(floattype));
	int* index_up=calloc(width,sizeof(int));
	floattype* slope_up=calloc(width,sizeof(floattype));
	int* index=calloc(width,sizeof(int));
	floattype* z=calloc(width,sizeof(floattype));
	floattype* y_low=malloc(width*sizeof(floattype));
	floattype* y_up=malloc(width*sizeof(floattype));
	int s_low=0;
	int c_low=0;
	int s_up=0;
	int c_up=0;
	int c=0;
	int i=2;
	y_low[0]=y_up[0]=0;
	y_low[1]=input[0]-lambda;
	y_up[1]=input[0]+lambda;
	for (;i<width;i++) {
		y_low[i]=y_low[i-1]+input[i-1];
		y_up[i]=y_up[i-1]+input[i-1];
	}
	y_low[width-1]+=lambda;
	y_up[width-1]-=lambda;
	slope_low[0] = INFINITY;
	slope_up[0] = -INFINITY;
	z[0]=y_low[0];
	for (i=1;i<width;i++) {
		index_low[++c_low] = index_up[++c_up] = i;
		slope_low[c_low] = y_low[i]-y_low[i-1];
		while ((c_low > s_low+1) && (slope_low[MAX(s_low,c_low-1)]<=slope_low[c_low])) {
			index_low[--c_low] = i;
			if (c_low > s_low+1) 
				slope_low[c_low] = (y_low[i]-y_low[index_low[c_low-1]]) /
					(i-index_low[c_low-1]);
			else
				slope_low[c_low] = (y_low[i]-z[c]) / (i-index[c]);
		}
		slope_up[c_up] = y_up[i]-y_up[i-1];
		while ((c_up > s_up+1) && (slope_up[MAX(c_up-1,s_up)]>=slope_up[c_up])) {
			index_up[--c_up] = i;
			if (c_up > s_up+1)
				slope_up[c_up] = (y_up[i]-y_up[index_up[c_up-1]]) /
					(i-index_up[c_up-1]);
			else
				slope_up[c_up] = (y_up[i]-z[c]) / (i-index[c]);
		}
		while ((c_low==s_low+1) && (c_up>s_up+1) && (slope_low[c_low]>=slope_up[s_up+1])) {
			index[++c] = index_up[++s_up];
			z[c] = y_up[index[c]];
			index_low[s_low] = index[c];
			slope_low[c_low] = (y_low[i]-z[c]) / (i-index[c]);
		}
		while ((c_up==s_up+1) && (c_low>s_low+1) && (slope_up[c_up]<=slope_low[s_low+1])) {
			index[++c] = index_low[++s_low];
			z[c] = y_low[index[c]];
			index_up[s_up] = index[c];
			slope_up[c_up] = (y_up[i]-z[c]) / (i-index[c]);
		}
	}
	for (i=1;i<=c_low-s_low;i++) 
		z[c+i]=y_low[index[c+i]=index_low[s_low+i]];
	c = c + c_low-s_low;
	int j=0;
	float a;
	i=1;
	while (i<=c) {
		a = (z[i]-z[i-1]) / (index[i]-index[i-1]);
		while (j<index[i]) {
			output[j] = a;
			j++;
		}
		i++;
	}
	free(index_low); free(slope_low);
	free(index_up);  free(slope_up);
	free(index);     free(z);
	free(y_low);     free(y_up);
}


/*
Function used for saving the signals in the .cimg format with floating point values.
See http://cimg.sourceforge.net
*/
void save_cimg(const char *filename,int width, float* buffer1, float* buffer2, float* buffer3) {
	FILE *file = fopen(filename,"wb");
	fprintf(file,"%u float little_endian\n",1);    
	fprintf(file,"%u %u %u %u\n",width,1,1,3);
	fwrite(buffer1,sizeof(float),width,file);
	fwrite(buffer2,sizeof(float),width,file);
	fwrite(buffer3,sizeof(float),width,file);
	fclose(file);
}


/*
The main function generates the Figure 2 of the paper.
If not already done, you should consider installing G'MIC: see http://gmic.sourceforge.net
Depending on the settings of your machine, you may just type on the command line:
gcc condat_fast_tv.c -O1 -lm -lgsl -o condat_fast_tv -I/usr/local/include/; ./condat_fast_tv; gmic signal.cimg -plot 1,1
*/
int main(int argc, char* argv[]) {
  	clock_t time1, time2;
	int i, width;
	const gsl_rng_type *T;
  	gsl_rng *r;
  	double aux;
  	float* x0;
  	float* y;
  	float* x;	
	width=1000;
	gsl_rng_env_setup();
	gsl_rng_default_seed=40;	
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	x0=malloc(width*sizeof(float));
	y=malloc(width*sizeof(float));
	x=malloc(width*sizeof(float));
	x0[0]=10;	
	for (i=1; i<width; i++) {
		aux=gsl_rng_uniform(r);
		if (aux<0.95) aux=0.0;	
		else aux=gsl_ran_gaussian(r,4.0); 
		x0[i]=x0[i-1]+aux;
	}
	for (i=0; i<width; i++) 
		x[i]=y[i]=x0[i]+gsl_ran_gaussian(r,1.0);	
	time1=clock();
	TV1D_denoise(x,x,width,2.0);
	/*TV1D_denoise_tautstring(x,x,width,2.0);*/
	time2=clock();
	printf("computation time: %f\n",(double)(time2-time1)/CLOCKS_PER_SEC);
	save_cimg("signal.cimg",width,y,x0,x);
	free(x0);
	free(y);
	free(x);
	return 0;
}
