// to build: /usr/bin/gcc main.c finish.c output.c advance.c rhs.c halo_update.c SWE_testcase5.c vmath.c Input.c profiling.c -O3 -I. -o pdex

#ifdef CACHEQ
#include "cq.h"   // use cq_malloc()
#define getTime() cq_nstime()
#endif
#ifndef CQ_POOL
#define CQ_POOL(x)
#endif
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <floating_types.h>
#include <SWE.h>
#include <NS.h>
#include <rcm.h>
#include <profiling.h>
#include <vmath.h>
#include <memory.h>
// Only need _VERBOSE_ to debug initialization
#undef _VERBOSE_
#define NCheckAnswers 100

//void RK3_advance(fType dt, fType *t, int NPts, fType *S, void fswe( const fType t, const fType* S, fType* K));
#ifdef CACHEQ
void RK4_advance(fType dt, timing_struct CQ_POOL(4)* t, int NPts, fType CQ_POOL(2)*S_cq,fType CQ_POOL(3)*K_cq,fType CQ_POOL(2)*dS,fType CQ_POOL(3)*Stmp, SWE_struct* SWE_cq);
fType CQ_POOL(2)* reorder_2D_fp_arr_fpga(fType CQ_POOL(2)* var, int* mapping, int dim1, int dim1_stride, int dim2);
#else
void RK4_advance(fType dt, fType *t, int NPts, fType *S, void fswe( const fType t, const fType* S, fType* K));
#endif
void output(const fType t, const int NState, const fType *S);
void history(const char* pdecase, fType *S, fType *S_obs, int NPts);

fType* reorder_2D_fp_arr(fType* var, int* mapping, int dim1, int dim1_stride, int dim2);

void finish(fType *S, fType *K);

//supported RHSs

void fsimple( const fType t, const fType *S, fType* K );
void   faero( const fType t, const fType *S, fType* K );
void    fswe( const fType t, const fType* S, fType* K );
void fstraka( const fType t, const fType* S, fType* K );

// SWE operators...

void SWE_init_bin(char* inputFile, char* testcase, SWE_struct *SWE, fType *S);
void SWE_get_init_filename(char *inputFile, int NNodes, char *testcase);
void SWE_profiler(timing_struct local_timer, SWE_struct SWE, int RK_Order, int Nsteps);
void SWE_correctness(SWE_struct SWE, char *inputFile, fType* S, fType* S_ref);
void SWE_final_conditions(char* inputFile, fType* S_ref);

// Global variables

//SWE_struct CQ_POOL(2)* SWE;
SWE_struct*  SWE;
NS_struct* NS;

timing_struct local_timer;

int main(int argc, char *argv[])
{

  int NState;
  int NNodes;
  int NPts;
  int Nsteps;

  int RK_Order;
  fType dt;
  fType t;
  //fType* S=NULL;
  fType* S_obs=NULL;
  fType* S_ref=NULL;
  //fType* K=NULL;
  int *mapping=NULL;
  char inputFile[30];
#ifdef CACHEQ
        fType CQ_POOL(2)* S=NULL;
	fType CQ_POOL(3)* K=NULL;
	fType CQ_POOL(2)* dS=NULL;
        fType CQ_POOL(3)* Stmp=NULL;
#else
	fType* S=NULL;
	fType* K=NULL;
#endif

#ifndef CACHEQ
  void (*fptr)( const fType t, const fType* S, fType* K) = NULL; //function pointer to assign RHS for specified problem
#endif
  if(argc==2)
    {
      if(! strcmp(argv[1], "simple")){
	NNodes = 1;    // number of nodes
	NState = 1;    // number of state variables in SW PDE
	NPts   = NState*NNodes;
	S = (fType*) malloc(NPts * sizeof(fType));
	Nsteps = 100;  // number of timesteps
	RK_Order = 4;  // Runge Kutta Order
	dt = 0.1;      // timestep in seconds
	t = 0.0;
	S[0] = 1.0;
#ifndef CACHEQ
	fptr = &fsimple;
#endif
      }
      else if(! strcmp(argv[1], "aero")){
	NNodes = 1;    // number of nodes
	NState = 2;    // number of state variables in SW PDE
	NPts   = NState*NNodes;
	S = (fType*) malloc(NPts * sizeof(fType));
	Nsteps = 100;  // number of timesteps
	RK_Order = 4;  // Runge Kutta Order
	dt = 0.1;      // timestep in seconds
	t = 0.0;
	S[0] = 2000.0;
	S[1] = 0.0;

#ifndef CACHEQ
	fptr = &faero;
#endif
      }
      else if(! strcmp(argv[1], "swe")){

	int RCM = 1;            // Reverse Cuthill McKee reorder 1 = Yes; 2 = No
	NNodes = 2562;         // Resolution in points
        t = 0;        // Start time 
#ifndef CACHEQ
	fptr = &fswe;
#endif
	Nsteps = 100;           // number of timesteps
	RK_Order = 4;           // Runge Kutta Order

	// Build init file name

	char testcase[4];
	strcpy(testcase,"tc5"); // Select testcase
	SWE_get_init_filename(inputFile, NNodes, testcase);
#ifdef _VERBOSE_
	printf("inputFile = %s\n", inputFile);
#endif
	NState = 4;
	NPts = NNodes*NState;
        SWE   = (SWE_struct*) malloc(sizeof(SWE_struct));
#ifdef CACHEQ
	S = (fType CQ_POOL(2) *) cq_malloc(2,NPts * sizeof(fType));
        K = (fType CQ_POOL(3) *) cq_malloc(3,NPts * sizeof(fType));
	dS = (fType CQ_POOL(2)*) cq_malloc(2,NPts * sizeof(fType));
        Stmp = (fType CQ_POOL(3)*) cq_malloc(3,NPts * sizeof(fType));
#else
	S     = (fType*) malloc(NPts * sizeof(fType));
	K     = (fType*) malloc(NPts * sizeof(fType));
#endif
	S_obs = (fType*) calloc(NPts, sizeof(fType));
	S_ref = (fType*) calloc(NPts, sizeof(fType));

	SWE_init_bin(inputFile, testcase, SWE, S);
	dt = SWE->dt;           // timestep in seconds

	if (RCM == 1){
	  printf("\n==============================================================\
                  \n\nReordering RBF-FD stencil matrix using Reverse Cuthill McKee\
                  \n\nINITIAL stencil matrix:\n");
	  rcm_print_max_bandwidth(SWE->NNodes, SWE->NNbrs_Padded, SWE->NNbrs, SWE->idx);
	  mapping = rcm_mapping(SWE->NNodes, SWE->NNbrs_Padded, SWE->NNbrs, SWE->idx);
	  rcm_check_mapping(mapping, SWE->NNodes);
	  SWE_reorder_nodes(SWE, mapping);
	  printf("\nREORDERED (RCM) stencil matrix:\n");
	  rcm_print_max_bandwidth(SWE->NNodes, SWE->NNbrs_Padded, SWE->NNbrs, SWE->idx);
#ifdef CACHEQ
	  S = reorder_2D_fp_arr_fpga(S, mapping, SWE->NNodes, NState, NState);
#else
	  S = reorder_2D_fp_arr(S, mapping, SWE->NNodes, NState, NState);
#endif
	}
	else{
	  for(int i=0; i<NNodes; i++){
	    mapping[i]=i;
	  }
	}

	// Initialize verification array for SWE case

	SWE_final_conditions(inputFile, S_ref);
	printf("\nReordering verification State Vector");
	S_ref = reorder_2D_fp_arr(S_ref, mapping, SWE->NNodes, NState, NState);
	printf(" ...done\n");
#ifdef _VERBOSE_
	output(t,NState,S);
#endif

      }
      else if(! strcmp(argv[1], "straka")){
        NS= (NS_struct*) malloc(sizeof(NS_struct));
	NNodes = 3000; // number of RBF nodes
	NState = 4;    // number of state variables in SW PDE
	NPts   = NState*NNodes;
	S = (fType*) malloc(NPts * sizeof(fType));
	Nsteps = 100; // number of timesteps
	RK_Order = 4; // Runge Kutta Order
	dt = 2;       // timestep in seconds
	t = 0.0;
#ifndef CACHEQ
	fptr = &fstraka;
#endif
      }
      else{
	printf("ERROR IN MAIN: Unknown PDE case %s\n",argv[1]);
	exit(-1);
      }
    }
  else{
    printf("ERROR IN MAIN: Invalid num of command line arguments\n");
    printf("Correct usage is: ./pdex {simple,aero,swe,straka}\n");
    exit(-1);
  }

   init_timer(&local_timer);

#ifdef CACHEQ
        CQ_SWE_struct CQ_POOL(2)* SWE_cq = (CQ_SWE_struct CQ_POOL(2)*) cq_malloc(2,sizeof(CQ_SWE_struct));	
        memcpy((void*) SWE_cq, SWE, sizeof(*SWE));
	
        SWE_cq->x = cq_canonicalize_pointer((void*)SWE->x);
        SWE_cq->y = cq_canonicalize_pointer((void*)SWE->y);
        SWE_cq->z = cq_canonicalize_pointer((void*)SWE->z);
        SWE_cq->f = cq_canonicalize_pointer((void*)SWE->f);

        SWE_cq->ghm = cq_canonicalize_pointer((void*)SWE->ghm);
        SWE_cq->gradghm = cq_canonicalize_pointer((void*)SWE->gradghm);

        SWE_cq->p_u = cq_canonicalize_pointer((void*)SWE->p_u);
        SWE_cq->p_v = cq_canonicalize_pointer((void*)SWE->p_v);
        SWE_cq->p_w = cq_canonicalize_pointer((void*)SWE->p_w);

        SWE_cq->idx = cq_canonicalize_pointer((void*)SWE->idx);

        SWE_cq->Dx = cq_canonicalize_pointer((void*)SWE->Dx);
        SWE_cq->Dy = cq_canonicalize_pointer((void*)SWE->Dy);
        SWE_cq->Dz = cq_canonicalize_pointer((void*)SWE->Dz);
        SWE_cq->L = cq_canonicalize_pointer((void*)SWE->L);
        

        CQ_timing_struct CQ_POOL(4)* cq_timing =
                            (CQ_timing_struct CQ_POOL(4)*) cq_malloc(4,sizeof(CQ_timing_struct));
	cq_timing->t_rhs = 0.0;
	cq_timing->t_update = 0.0;
	uint64_t t_start;
	uint64_t t_main = 0.0;
	uint64_t init_rhs;
	uint64_t init_update;
#endif


#ifdef _VERBOSE_
  output(t,NState,S);
#endif

  // Integrate the equations...
  printf("\n==============================================================\n\n");
  printf("Advancing PDE equations %d timesteps with Runge Kutta Order = %d\n",Nsteps,RK_Order);

  if (RK_Order == 3){
    for(int it=0; it<Nsteps; it++){
#ifndef CACHEQ
      //RK3_advance(dt,&t,NPts,S,fptr);
#endif
      if (it == NCheckAnswers-1){
	printf("\nSaving state variable After %d iterations ",it+1);
	history(argv[1],S,S_obs,NPts);
	printf("...done\n");
      }
#ifdef _VERBOSE_
      output(t,NState,S);
#endif
    }
  }
  else if (RK_Order == 4){
#ifdef CACHEQ
    init_rhs = cq_timing->t_rhs;
    init_update = cq_timing->t_update;
    //printf("\nbcq_timing->t_rhs: %lu ",(cq_timing->t_rhs));
    //printf("\nf bcq_timing->t_rhs: %f ",(cq_timing->t_rhs)*1e-9);
#endif
    for(int it=0; it<Nsteps; it++){
#ifdef CACHEQ
      _cq_ignorelcd(); 
      //init_rhs = cq_timing->t_rhs;
      //init_update = cq_timing->t_update;
      t_start = cq_nstime();
      //printf("\nbcq_timing->t_rhs: %f ",(cq_timing->t_rhs)*1e-9);
      //printf("\nbcq_timing->t_update: %f ",(cq_timing->t_update)*1e-9);   
      //RK4_advance(dt,cq_timing,NPts,S,fswe); //S-Size of NPts
      RK4_advance(dt,(timing_struct*)cq_timing,NPts,S,K,dS,Stmp,(SWE_struct*)SWE_cq); //S-Size of NPts
      t_main += cq_nstime() - t_start;
      //printf("\niter %d time: %f ",it+1,(cq_nstime() - t_start)*1e-9); //%ld
      //printf("\ntotal time: %f ",(t_main)*1e-9); //%ld
      //printf("\ncq_timescale: %d ",cq_timescale);
      //printf("\ncq_timestamp: %ld ",cq_timestamp());
      //printf("\niter %d acq_timing->t_rhs: %f ",it,((cq_timing->t_rhs)-init_rhs)*1e-9);
      //printf("\nacq_timing->t_update: %f ",((cq_timing->t_update)-init_update)*1e-9);
      if (it==NCheckAnswers-1){
        printf("\nSaving state variable After %d iterations ",it+1);
        history(argv[1],S,S_obs,NPts);
        printf("...done\n");
      }
      //output(t,NState,S);
#else
      RK4_advance(dt,&t,NPts,S,fswe); //S-Size of NPts

      if (it==NCheckAnswers-1){
	printf("\nSaving state variable After %d iterations ",it+1);
	history(argv[1],S,S_obs,NPts);
	printf("...done\n");
      }
      //output(t,NState,S);
#endif
#ifdef _VERBOSE_
      output(t,NState,S);
#endif
    }
  }
  else{
    printf("ERROR IN MAIN: RK_Order = %i not supported\n",RK_Order);
    exit(-1);
  }

#ifndef CACHEQ
  if(! strcmp(argv[1], "swe")){
    SWE_profiler(local_timer,*SWE,RK_Order,Nsteps);
    SWE_correctness(*SWE, inputFile, S_obs, S_ref);
    }


  finish(S,S_ref);
#else
  double t_rhs_d,t_update_d;
  printf("Nsteps= %d cq_timing->t_rhs = %lu cq_timing->t_update = %lu\n",Nsteps,cq_timing->t_rhs,cq_timing->t_update);
  t_rhs_d = (double)cq_timing->t_rhs;
  t_update_d = (double)cq_timing->t_update;
  //memcpy(cq_timing->t_rhs,t_rhs_d,sizeof(cq_timing->t_rhs));
  //memcpy(cq_timing->t_update,t_update_d,sizeof(cq_timing->t_update));  
 SWE_correctness(*SWE, inputFile, S_obs, S_ref);
 if(! strcmp(argv[1], "swe")){
  double rk_flop_per_step = RK_Order*((2*SWE->NNbrs-1)*4*(3+1) + 17*3 + 31)*SWE->NNodes;
  //double t_rhs_per_step = 1e-9*t_rhs_d/Nsteps;
  double t_rhs_per_step = 1e-9*t_main/Nsteps;
  double rhs_Gflops = 1e-9*(rk_flop_per_step/t_rhs_per_step);
  double update_flop_per_step = (4+2*RK_Order)*(RK_Order*SWE->NNodes);
  double t_update_per_step = 1e-9*t_update_d/Nsteps;
  double update_Gflops = 1.e-9*(update_flop_per_step/t_update_per_step);
  printf("\n==============================================================\n\n");
  printf("Benchmark results:\n\n");
  //printf("Nsteps= %d cq_timing->t_rhs = %f cq_timing->t_update = %f\n",Nsteps,t_rhs_d,t_update_d);
  printf("time in rhs = %f Gflops = %f rk_flops/step = %f\n",t_rhs_per_step,rhs_Gflops,rk_flop_per_step);
  //printf("time in rhs = %f Gflops = %f rk_flops/step = %f\n",t_rhs_per_step,rhs_Gflops,rk_flop_per_step);
  //printf("time in update = %f Gflops = %f\n",t_update_per_step,update_Gflops);
 }
#endif
/*
//#ifdef CACHEQ 
  free(SWE_cq->x);
  free(SWE_cq->y);
  free(SWE_cq->z);
  free(SWE_cq->f);
  free(SWE_cq->ghm);
  free(SWE_cq->p_u);
  free(SWE_cq->p_v);
  free(SWE_cq->p_w);
  free(SWE_cq->gradghm);
  free(SWE_cq->Dx);
  free(SWE_cq->Dy);
  free(SWE_cq->Dz);
  free(SWE_cq->L);
  free(SWE_cq->idx);
  free(SWE_cq);
  cq_free(K_cq);
  cq_free(S_cq);
  cq_free(K1_cq);
  cq_free(Stmp_cq);
  cq_free(dS_cq);
  
//#endif
*/

  if(! strcmp(argv[1], "straka")){
     free(NS);
  }
  //free(S);
  free(S_obs);
  free(SWE->x);
  free(SWE->y);
  free(SWE->z);
  free(SWE->f);
  free(SWE->ghm);
  free(SWE->p_u);
  free(SWE->p_v);
  free(SWE->p_w);
  free(SWE->gradghm);
  free(SWE->Dx);
  free(SWE->Dy);
  free(SWE->Dz);
  free(SWE->L);
  free(SWE->idx);
  free(SWE);
  free(K);
  //free(dS);
  //free(Stmp);
  free(mapping);

}
