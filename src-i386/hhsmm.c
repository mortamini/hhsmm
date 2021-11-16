#include "stdio.h"
#include <R.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <Rinternals.h>
#include <stddef.h>

#define int8 unsigned char

void checkmem(void *x) {
  if(x==NULL) error("Out of memory.");
}

int min(int a,int b) {
  if(a<b) return(a);
  else return(b);
}

void simulate_markov(double *start, double *a,int *nstates,int *state,int *T,int *nseq) {
  int t,i,n,*s=NULL;
  int K = *nstates;
  int N = *nseq;
  double tmp;
  GetRNGstate();
  for(n=0;n<N;n++) {
    if(n==0) s = state;
    else s = s+T[n-1];
    tmp=unif_rand();
    i=0;
    while(tmp>start[i])	i++;
    s[0]=i+1;
    for(t=1;t<T[n];t++) {
      tmp=unif_rand();
      i=0;
      while(tmp>a[i*K+s[t-1]-1])	i++;
      s[t]=i+1;
    }
  }
  PutRNGstate();
}

void **alloc_matrix(int nrow,int ncol,int size) {
  int i;
  void **x = calloc(nrow,sizeof(void *));
  checkmem(x);
  for(i=0;i<nrow;i++)	{
    x[i]=calloc(size,ncol);
    checkmem(x[i]);
  }
  return(x);
}

void free_matrix(int nrow,int ncol,void **x) {
  int i;
  for(i=0;i<nrow;i++)	free(x[i]);
  free(x);
}


void forward(double *a,double *start,double *p0,double *d,double *D,int *timelength,int *nstates,int *M,
	     double **F0,double *N0,double **si0,int *nsequences,int *totallength,double *semi) {
  double obs;
  int u,i,j,t,n;
  int T;
  int J = *nstates;
  double **p,**F,**si,*N;
  int offset;
  int ln = *totallength;
  int nseq = *nsequences;
  double prod;
//
  p = (double **)calloc(J,sizeof(double *));
  F = (double **)calloc(J,sizeof(double *));
  si = (double **)calloc(J,sizeof(double *));

/*
    printf("J = %d     T = %d\n",J,T);
    print_matrix(J,J,a);
    print_matrix(1,J,start);

    for(j=0;j<J;j++)  printf("M[%d] = %d\t",j,M[j]);
    printf("\nd = \n");
    print_matrix(J,100,d);
    printf("p = \n");
    print_matrix(J,T,p);
    printf("D = \n");
    print_matrix(J,100,D);
*/
  for(j=0;j<J;j++) {
    p[j]  = p0 + j*ln;
    F[j]  = F0[j];
    si[j] = si0[j];
  }
  N = N0;

  for(n=0;n<nseq;n++) {
    T = timelength[n];
    if(n>0){
      offset = timelength[n-1];
      for(j=0;j<J;j++) {
	p[j]  += offset;
	F[j]  += offset;
	si[j] += offset;
      }
      N += offset;
    }
    for(t=0;t<T;t++) {
      N[t]=0;
      for(j=0;j<J;j++) {
		if(semi[j]==1){
			F[j][t]=0;
			obs=p[j][t];
			if(t<T-1) {
	  			for(u=1;u<=min(t+1,M[j]);u++) {
	    				if(u<t+1) {
	      				F[j][t] += obs*d[j*M[j]+u-1]*si[j][t-u+1];
	      				N[t] += obs*D[j*M[j]+u-1]*si[j][t-u+1];
		 				prod = p[j][t-u]/N[t-u];
		 				if(isnan(prod)) prod = 1;
		 				if(isinf(prod)) prod = 1e10;
	      				obs*= prod;
		 				if(isnan(obs)) obs = 1;
		 				if(isinf(obs)) obs = 1e300;
	    				}
	    				else {
	      				F[j][t] += obs*d[j*M[j]+t]*start[j];
	      				N[t] += obs*D[j*M[j]+t]*start[j];
	    				}
	  			}
			}
			else{
	  			for(u=1;u<=min(t+1,M[j]);u++) {
	    				if(u<T) {
	      				F[j][T-1] += obs*D[j*M[j]+u-1]*si[j][T-u];
		 				prod = p[j][T-1-u]/N[T-1-u];
		 				if(isnan(prod)) prod = 1;
		 				if(isinf(prod)) prod = 1e10;
	      				obs*= prod;
		 				if(isnan(obs)) obs = 1;
		 				if(isinf(obs)) obs = 1e300;
	    				}
	    				else F[j][T-1] += obs*D[j*M[j]+T-1]*start[j];
	  			}
	  			N[T-1] += F[j][T-1];
			}
		  }//  semi states
		  else{
		  	if(t == 0){
				F[j][t] +=  p[j][t]*start[j];
         		}
		     else{
			 	F[j][t] +=  p[j][t]*si[j][t];
			}
       		N[t] += F[j][t];
      	}//   markov states
	  }
      for(j=0;j<J;j++){
		if(N[t]==0) N[t]=1e-300;
		F[j][t] /= N[t];
		F[j][t]+=1e-300;
      }
      if(t<T-1) {
		for(j=0;j<J;j++){
	  		si[j][t+1]=0;
	  		for(i=0;i<J;i++)
	    			si[j][t+1]+=a[j*J+i]*F[i][t];
		}
     }
   }
 }
  free(si);
  free(p);
  free(F);
  /*
    printf("\nF=\n");
    print_matrix2(J,T,F);
    print_matrix(1,T,N);
  */
}

void forward_backward(double *a,double *start,double *p0,double *d,double *D,int *timelength,int *nstates,int *M,double *L10,
	      double *N0,double *eta,double *F1,double *statein,double *gamma,int *nsequences,int *totallength,double *Gret
		  ,double *semi) {

  double **L0,**G0,**F0,**si0,**L,**G,obs,**si,**F,**num,*den,*N,**p,**L1;
  int u,i,j,k,t,T,n;
  int J = *nstates;
  int offset;
  int ln = *totallength;
  int nseq = *nsequences;
  double prod;
  double frac;
  double term;

  F0 =   (double **)alloc_matrix(J,ln,sizeof(double));
  si0  = (double **)alloc_matrix(J,ln,sizeof(double));
  den  = (double *)calloc(J,sizeof(double));
  num  = (double **)alloc_matrix(J,J,sizeof(double));

  forward(a,start,p0,d,D,timelength,nstates,M,F0,N0,si0,nsequences,totallength,semi);
  /*
    for(j=0;j<J;j++) memcpy(F1+j*T,F[j],T*sizeof(double));
    for(j=0;j<J;j++) memcpy(statein+j*T,si[j],T*sizeof(double));
    printf("J = %d     T = %d\n",J,T);
    for(j=0;j<J;j++)  printf("M[%d] = %d\t",j,M[j]);
    printf("\n");
    print_matrix(1,T,N);
    print_matrix2(J,50,F);
    printf("\n");
    print_matrix2(J,50,si);
  */

  //	L1 =  (double **)alloc_matrix(J,T,sizeof(double));
  L0  =  (double **)alloc_matrix(J,ln,sizeof(double));
  G0  =  (double **)alloc_matrix(J,ln,sizeof(double));

  p = (double **)calloc(J,sizeof(double *));
  F = (double **)calloc(J,sizeof(double *));
  G = (double **)calloc(J,sizeof(double *));
  L = (double **)calloc(J,sizeof(double *));
  si = (double **)calloc(J,sizeof(double *));
  L1 = (double **)calloc(J,sizeof(double *));

  /*
    printf("J = %d     T = %d\n",J,T);
    print_matrix(J,J,a);
    print_matrix(1,J,start);

    for(j=0;j<J;j++)  printf("M[%d] = %d\t",j,M[j]);
    printf("\nd = \n");
    print_matrix(J,100,d);
    printf("p = \n");
    print_matrix(J,T,p);
    printf("D = \n");
    print_matrix(J,100,D);
  */
  for(j=0;j<J;j++) {
    p[j]  = p0 + j*ln;
    F[j]  = F0[j];
    si[j] = si0[j];
    G[j] = G0[j];
    L[j] = L0[j];
    L1[j] = L10+j*ln;
  }
  N = N0;

  for(n=0;n<nseq;n++) {
    T = timelength[n];
    if(n>0) {
      offset = timelength[n-1];
      for(j=0;j<J;j++) {
	p[j]  += offset;
	F[j]  += offset;
	si[j] += offset;
	G[j] += offset;
	L[j] += offset;
	L1[j] += offset;
      }
      N += offset;
    }

    for(j=0;j<J;j++) L[j][T-1] = F[j][T-1];

    if(n==0) {
      for(j=0;j<J;j++)
		for(u=1;u<M[j];u++)
	  		eta[j*M[j]+u-1] = 0;
    }
    for(t=T-2;t>=0;t--) {
      for(j=0;j<J;j++) {
		  if(semi[j]==1){
			G[j][t+1] = 0;
			obs = 1;
			for(u=1;u<=min(T-1-t,M[j]);u++) {
		 		prod = p[j][t+u]/N[t+u];
		 		if(isnan(prod)) prod = 1;
		 		if(isinf(prod)) prod = 1e10;
	      obs*= prod;
		 if(isnan(obs)) obs = 1;
		 if(isinf(obs)) obs = 1e300;
		 if(u<T-1-t) {
		 	frac = L1[j][t+u]/F[j][t+u];
		 	if(isnan(frac)) frac = 1;
		 	if(isinf(frac)) frac = 1e10;
		    G[j][t+1] += frac*obs*d[j*M[j]+u-1];
			if(isinf(G[j][t+1])) G[j][t+1]=1e300;
			term = frac*obs*d[j*M[j]+u-1]*si[j][t+1];
			if(isnan(term)) term = 1e-300;
	    		eta[j*M[j]+u-1] += term;
			if(isinf(eta[j*M[j]+u-1])) eta[j*M[j]+u-1]=1e300;
	  	}
	  	else {
	    		G[j][t+1] += obs*D[j*M[j]+T-1-t-1];
			if(isinf(G[j][t+1])) G[j][t+1]=1e300;
			term = obs*D[j*M[j]+T-1-t-1]*si[j][t+1];
			if(isnan(term)) term = 1e-300;
	    		eta[j*M[j]+u-1] += term;
			if(isinf(eta[j*M[j]+u-1])) eta[j*M[j]+u-1]=1e300;
	  	}
	  if(t==0){
	    if(u>T-1){
		 	frac = L1[j][t+u]/F[j][t+u];
		 	if(isnan(frac)) frac = 1;
		 	if(isinf(frac)) frac = 1e10;
			term = frac*obs*d[j*M[j]+u-1]*start[j];
			if(isnan(term)) term = 1e-300;
			eta[j*M[j]+u-1] += term;
			if(isinf(eta[j*M[j]+u-1])) eta[j*M[j]+u-1]=1e300;
		}
	    else{
			term = obs*d[j*M[j]+u-1]*start[j];
			if(isnan(term)) term = 1e-300;
			eta[j*M[j]+u-1] += term;
			if(isinf(eta[j*M[j]+u-1])) eta[j*M[j]+u-1]=1e300;
		}
	  }
	}
	} // semi states
	else{
		G[j][t+1] = 0;
		//obs = 1;
	    //obs *= p[j][t+1]/N[t+1];
		if(si[j][t+1]==0) si[j][t+1]=1e-300;
		G[j][t+1] = L[j][t+1]/si[j][t+1];
		if(isinf(G[j][t+1])) G[j][t+1]=1e300;
		if(isnan(G[j][t+1])) G[j][t+1]=1e-300;
     }// markove state
   }
      for(j=0;j<J;j++){
		L1[j][t] = 0;
		for(k=0;k<J;k++) L1[j][t] += G[k][t+1]*a[k*J+j];
		L1[j][t] *= F[j][t];
		if(semi[j]==1){
			L[j][t] = L1[j][t] + L[j][t+1] - G[j][t+1]*si[j][t+1];
		}
		else{
			L[j][t] = L1[j][t];
		}
		if(((isnan(L[j][t])) || (L[j][t]==0))) L[j][t]=1e-300;
		if(((isnan(L1[j][t])) || (L1[j][t]==0))) L1[j][t]=1e-300;
      }
    }
  }
  //reset pi
  for(i=0;i<J;i++) {
    start[i]=0;
    den[i] = 0;
    for(j=0;j<J;j++) num[i][j] = 0;
  }
  //new estimates for a and pi
  for(i=0;i<J;i++){
    for(n=0;n<nseq;n++) {
      T = timelength[n];
      if(n==0) {
		for(j=0;j<J;j++) {
	  		F[j]  = F0[j];
	  		G[j] = G0[j];
	  		L[j] = L0[j];
	  		L1[j] = L10+j*ln;
		}
     }
     else {
		offset = timelength[n-1];
		for(j=0;j<J;j++) {
	  		F[j]  += offset;
	  		G[j] += offset;
	  		L[j] += offset;
	  		L1[j] += offset;
		}
    }
      start[i]+=L[i][0];
      for(t=0;t<=(T-2);t++){
		if(semi[i]==1) den[i] += L1[i][t];
		else den[i] += L[i][t];
	 }
      	for(j=0;j<J;j++){
		for(t=0;t<=(T-2);t++) num[i][j] += G[j][t+1]*a[j*J+i]*F[i][t];
      }
    }
  }
  for(i=0;i<J;i++) {
    	start[i]/=nseq;
    	for(j=0;j<J;j++){
		if(((num[i][j]==0) && (den[i]==0))){
			a[j*J+i]=0;
		}
		else a[j*J+i]=num[i][j]/den[i];
	}
  }
  /*
    printf("\n\nG=\n");
    print_matrix2(J,T,G);
    printf("\nL=\n");
    print_matrix2(J,T,L);
    printf("\nL1=\n");
    print_matrix(J,T,L1);
    printf("\n");
  */
  for(j=0;j<J;j++) {
    memcpy(gamma+j*ln,L0[j],ln*sizeof(double));
    memcpy(F1+j*ln,F0[j],ln*sizeof(double));
    memcpy(Gret+j*ln,G0[j],ln*sizeof(double));
  }
  //	free_matrix(J,T,(void **)L1);
  free_matrix(J,ln,(void **)F0);
  free_matrix(J,ln,(void **)G0);
  free_matrix(J,ln,(void **)si0);
  free_matrix(J,ln,(void **)L0);
  free_matrix(J,J,(void **)num);
  free(den);
  free(p);
  free(F);
  free(G);
  free(L);
  free(si);
  free(L1);
}


void viterbi(double *a,double *start,double *p0,double *d0,double *D0,int *timelength,int *nstates,int *M,
	     double *alpha0,int *statehat,int *psi_state0,int *psi_time0, double *semi) {

  double obs;
  int u,i,j,t,T;
  int J = *nstates;

  double **p,**d,**D,**alpha,**si,**si0,tmp_max=-1e300,tmp1;
  int *q;
  int **psi_time;
  int **psi_state;
  T = *timelength;

  si0 = (double **)alloc_matrix(J,T,sizeof(double));
  psi_time = (int **)calloc(J,sizeof(int *));
  psi_state = (int **)calloc(J,sizeof(int *));
  p = (double **)calloc(J,sizeof(double *));
  d = (double **)calloc(J,sizeof(double *));
  D = (double **)calloc(J,sizeof(double *));
  alpha = (double **)calloc(J,sizeof(double *));
  si = (double **)calloc(J,sizeof(double *));

  checkmem(si0);
  checkmem(psi_time);
  checkmem(psi_state);
  checkmem(p);
  checkmem(d);
  checkmem(D);
  checkmem(alpha);
  checkmem(si);


  for(j=0;j<J;j++) {
    d[j] = d0 + j*M[j];
    D[j] = D0 + j*M[j];
    p[j]  = p0 + j*T;
    alpha[j]  = alpha0 + j*T;
    si[j] = si0[j];
    psi_time[j] = psi_time0 + j*T;
    psi_state[j] = psi_state0 + j*T;
  }

  for(t=0;t<T;t++) {
    for(j=0;j<J;j++) {
		if(semi[j]==1){// semi Markov
      		obs=0;
      		if(t<T-1) {
				for(u=1;u<=min(t+1,M[j]);u++) {
	  				if(u<t+1) {
	    					tmp1 = obs+d[j][u-1]+si[j][t-u+1];
	    					if(u==1 || tmp_max < tmp1) {
	      					tmp_max = tmp1;
	      					psi_time[j][t]=u;
	    					}//if
	    					obs += p[j][t-u];
	  				}
	  				else {
	    					tmp1 = obs+d[j][t]+start[j];
	    					if(u==1 || tmp_max < tmp1) {
	      					tmp_max = tmp1;
	      					psi_time[j][t]=u;
	    					}//if
	  				}//if else
				}// for u
				alpha[j][t] = tmp_max+p[j][t];
      		}
      		else{
				for(u=1;u<=min(t+1,M[j]);u++) {
	  				if(u<T) {
	    					tmp1 = obs+D[j][u-1]+si[j][t-u+1];
	    					if(u < 2000)
	      					if(u==1 || tmp_max < tmp1) {
								tmp_max = tmp1;
								psi_time[j][t]=u;
	      					}//if
	    						obs += p[j][T-1-u];
	  				}
	  				else {
	    					tmp1 = obs+D[j][T-1]+start[j];
	    					if(u==1 || tmp_max < tmp1) {
	      					tmp_max = tmp1;
	      					psi_time[j][t]=u;
	    					}//if
	  				}//if else
				}//for u
			alpha[j][t] = tmp_max+p[j][t];
      	}//if else
	}	// semi Markov States
	else{
		obs = 0;
		if(t == 0){
			alpha[j][t] = p[j][t] + start[j];
		}
		else{
			alpha[j][t] = p[j][t] + si[j][t];
		}
	}// Markov states
  }//for j
  if(t<T-1) {
      for(j=0;j<J;j++){
		i=0;
		si[j][t+1]=a[j*J+i]+alpha[i][t];
		psi_state[j][t+1]=0;
		for(i=1;i<J;i++)
	  		if(i!=j || semi[j]==0) {
	    			tmp1 = a[j*J+i]+alpha[i][t];
	    			if(si[j][t+1] <= tmp1) {
	      			si[j][t+1] = tmp1;
	      			psi_state[j][t+1]=i;
	    			}//if
	  		}//if
    }//for j
    }//if t
  }// for t
  //and now we backtrack!
  	q=statehat;
	q[T-1] = 1;
  //q[T-1] = J-1;
  //for(j=0;j<J;j++) if(semi[j]==0) q[T-1] = j;
  for(j=0;j<J;j++) if(alpha[q[T-1]][T-1] <= alpha[j][T-1]) q[T-1] = j;
  u=1;
  for(t=T-2;t>=0;t--) {
    if(u < psi_time[q[t+u]][t+u]) {
      q[t] = q[t+u];
      u++;
    }
    else {
      q[t] = psi_state[q[t+u]][t+u];
      u=1;
    }
  }

  free(si);
  free_matrix(J,T,(void **)si0);
  free(p);
  free(alpha);
  free(psi_time);
  free(psi_state);
}


