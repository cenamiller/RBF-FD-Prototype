#include <stdlib.h>
#include <floating_types.h>
#include <profiling.h>
//#include <advance.h>

#ifdef CACHEQ
#include <cq.h>
#define getTime()         (cq_nstime() / 1e18)
#else
#define CQ_POOL(x)
#define CLANG_NOINLINE
#endif

extern timing_struct local_timer;
#include <SWE.h>
extern SWE_struct* SWE;

void halo_update(fType* S);
/*
void RK3_advance(const fType dt, fType *t, const int NPts, fType *S, void f ( const fType t, const fType* S, fType* K)){

  const fType sixth = 1.0/6.0;

  fType *K1;
  fType *K2;
  fType *K3;

  fType *K = (fType*) malloc(NPts * sizeof(fType));
  fType *dS = (fType*) malloc(NPts * sizeof(fType));
  fType *Stmp1 = (fType*) malloc(NPts * sizeof(fType));
  fType *Stmp2 = (fType*) malloc(NPts * sizeof(fType));

  double t_start = getTime();
  K1 = K;
  f((*t), S, K1);
  local_timer.t_rhs +=  (getTime() - t_start);

  t_start = getTime();
  for(int i=0; i<NPts; i++){
    Stmp1[i] =  S[i] + 0.5*dt*K1[i]; //2 flop/site
    Stmp2[i] =  S[i] - dt*K1[i]; //2 flop/site
    dS[i] = K1[i];
  }
  local_timer.t_update +=  (getTime() - t_start);

  t_start = getTime();
  K2 = K;
  f((*t)+0.5*dt, Stmp1, K2);
  local_timer.t_rhs +=  (getTime() - t_start);

  t_start = getTime();
  for(int i=0; i<NPts; i++){
    Stmp2[i] += 2.0*dt*K2[i]; //2 flop/site
    dS[i] += 4.0*K1[i];  //2 flop/site
  }
  local_timer.t_update +=  (getTime() - t_start);

  t_start = getTime();
  K3 = K;
  f((*t)+dt, Stmp2, K3);
  local_timer.t_rhs +=  (getTime() - t_start);

  t_start = getTime();
  for(int i=0; i<NPts; i++){
    S[i] += (sixth*dt)*(dS[i]+K3[i]); //2 flop/site
  }
  local_timer.t_update +=  (getTime() - t_start);

  (*t) += dt;

  free(K);
  free(dS);
  free(Stmp1);
  free(Stmp2);

}
*/

#if 1

//void RK4_advance(fType dt, CQ_timing_struct *t, int NPts, fType *S, void fswe( const CQ_timing_struct t, const fType* S, fType* K)){
CLANG_NOINLINE
void RK4_advance(fType dt, timing_struct CQ_POOL(4) * restrict cq_timing, int NPts, fType CQ_POOL(2) * restrict S,fType CQ_POOL(3) * restrict K, fType CQ_POOL(2) * restrict dS,fType CQ_POOL(3) * restrict Stmp, SWE_struct CQ_POOL(2)* SWE){

  const int X_COMP_ID = 0;
  const int Y_COMP_ID = 1;
  const int Z_COMP_ID = 2;

  const fType sixth = 1.0/6.0;
  fType CQ_POOL(3)* K1;
  fType CQ_POOL(3)* K2;
  fType CQ_POOL(3)* K3;
  fType CQ_POOL(3)* K4;
  //fType CQ_POOL(3)* dS = (fType CQ_POOL(3)*) cq_malloc(3,NPts * sizeof(fType));
  //fType CQ_POOL(3)* Stmp = (fType CQ_POOL(3)*) cq_malloc(3,NPts * sizeof(fType));

  //uint64_t t_start = getTime();
  //uint64_t t_start = cq_nstime();
  K1 = K;
  //fswe((*t),S,K1);
  // Dimensions
  
  const int NNodes = 2562; //SWE->NNodes;
  const int NState = 4; //SWE->NState;
  const int NNbrs  = 31; //SWE->NNbrs;
  const int NNbrs_Padded = 32; //SWE->NNbrs_Padded;

  // 3-D node coordinates (x,y,z)

  const fType* x = SWE->x;
  const fType* y = SWE->y;
  const fType* z = SWE->z;

  // spherical projetion operators

  const fType* p_u = SWE->p_u;
  const fType* p_v = SWE->p_v;
  const fType* p_w = SWE->p_w;

  // planetary variables
  const fType  gh0 = SWE->gh0;          // reference geopotential (gh0)
  const fType* f = SWE->f;              // coreolis term
  const fType* ghm = SWE->ghm;          // mountain geopotential height
  const fType* gradghm = SWE->gradghm;  // gradient of geopotential height

  // stencil and derivative operators

  const int*  idx = SWE->idx;
  const fType* Dx = SWE->Dx;
  const fType* Dy = SWE->Dy;
  const fType* Dz = SWE->Dz;
  const fType* L  = SWE->L;

  fType u;
  fType v;
  fType w;
  fType h;

  int iu;
  int iv;
  int iw;
  int ih;
  uint64_t t_start = cq_nstime();
  uint64_t new_rhs = 0;
  uint64_t new_up = 0;
  //t->t_rhs +=  (cq_nstime() - t_start-t_start1);
  //t->t_update += (cq_nstime() - t_start-t_start1);


  for(int i=0; i<NNodes; i++){
    _cq_ignorelcd();
    _cq_ignorelcd_pool(3);
    // initialize temporary variables for spatial differentiation at respective node
    fType du_dx = 0.0;
    fType dv_dx = 0.0;
    fType dw_dx = 0.0;
    fType dh_dx = 0.0;

    fType du_dy = 0.0;
    fType dv_dy = 0.0;
    fType dw_dy = 0.0;
    fType dh_dy = 0.0;

    fType du_dz = 0.0;
    fType dv_dz = 0.0;
    fType dw_dz = 0.0;
    fType dh_dz = 0.0;

    fType Lu = 0.0;
    fType Lv = 0.0;
    fType Lw = 0.0;
    fType Lh = 0.0;

    for (int j = 0; j < NNbrs; j++) {
      //_cq_unroll(NNbrs);
      const int dm_lind = i*NNbrs_Padded + j;
      const int inbr = idx[dm_lind];
      iu = U_SV_ID + NState*inbr;
      iv = V_SV_ID + NState*inbr;
      iw = W_SV_ID + NState*inbr;
      ih = H_SV_ID + NState*inbr;
    // temporary variables to hold state variable data for current node
      fType u = S[iu];
      fType v = S[iv];
      fType w = S[iw];
      fType h = S[ih];
      // get neighbor node's state vars
      // update sums
      du_dx += Dx[dm_lind] * u;
      dv_dx += Dx[dm_lind] * v;
      dw_dx += Dx[dm_lind] * w;
      dh_dx += Dx[dm_lind] * h;

      du_dy += Dy[dm_lind] * u;
      dv_dy += Dy[dm_lind] * v;
      dw_dy += Dy[dm_lind] * w;
      dh_dy += Dy[dm_lind] * h;

      du_dz += Dz[dm_lind] * u;
      dv_dz += Dz[dm_lind] * v;
      dw_dz += Dz[dm_lind] * w;
      dh_dz += Dz[dm_lind] * h;

      Lu += L[dm_lind] * u;
      Lv += L[dm_lind] * v;
      Lw += L[dm_lind] * w;
      Lh += L[dm_lind] * h;
    }

    iu = U_SV_ID + NState*i;
    iv = V_SV_ID + NState*i;
    iw = W_SV_ID + NState*i;
    ih = H_SV_ID + NState*i;

    u = S[iu];
    v = S[iv];
    w = S[iw];
    h = S[ih];
    fType rhs_u = - ((u * du_dx) + (v * du_dy) + (w * du_dz) + (f[i] * ((y[i] * w) - (z[i] * v))) + dh_dx);
    fType rhs_v = - ((u * dv_dx) + (v * dv_dy) + (w * dv_dz) + (f[i] * ((z[i] * u) - (x[i] * w))) + dh_dy);
    fType rhs_w = - ((u * dw_dx) + (v * dw_dy) + (w * dw_dz) + (f[i] * ((x[i] * v) - (y[i] * u))) + dh_dz);

    // vector linear ids
    int vect_idx = 3*i + X_COMP_ID;
    int vect_idy = 3*i + Y_COMP_ID;
    int vect_idz = 3*i + Z_COMP_ID;

    // evaluate projections and apply hyperviscosity
    K1[iu] = (p_u[vect_idx] * rhs_u) + (p_u[vect_idy] * rhs_v) + (p_u[vect_idz] * rhs_w) + Lu;
    K1[iv] = (p_v[vect_idx] * rhs_u) + (p_v[vect_idy] * rhs_v) + (p_v[vect_idz] * rhs_w) + Lv;
    K1[iw] = (p_w[vect_idx] * rhs_u) + (p_w[vect_idy] * rhs_v) + (p_w[vect_idz] * rhs_w) + Lw;
    // ----------------- Use RBF-FD Approximations to Calculate the RHS of the Geopotential Equation ---------------------//

    K1[ih] = - ((u * (dh_dx - gradghm[vect_idx])) + (v * (dh_dy - gradghm[vect_idy])) + (w * (dh_dz - gradghm[vect_idz]))
            + ((h + gh0 - ghm[i]) * (du_dx + dv_dy + dw_dz))) + Lh;
  }


  //t->t_rhs +=  (getTime()() - t_start);
  new_rhs +=  (cq_nstime() - t_start);
  //cq_timing->t_rhs +=  (cq_nstime() - t_start);
  
  t_start = cq_nstime();
  //t_start = getTime();
  for(int i=0; i<NPts; i++){
    Stmp[i] =  S[i] + 0.5*dt*K1[i]; // 2 flop/site
    dS[i] = K1[i];
  }
  //halo_update(Stmp);
  new_up +=  (cq_nstime() - t_start);
  //cq_timing->t_update +=  (cq_nstime() - t_start);
  //t->t_update +=  (getTime() - t_start);
  
  t_start = cq_nstime();
  //t_start = getTime();
  K2 = K;
  //fswe((*t)+0.5*dt, Stmp, K2);
    
  for(int i=0; i<NNodes; i++){
    _cq_ignorelcd();
    _cq_ignorelcd_pool(3);
    // initialize temporary variables for spatial differentiation at respective node
    fType du_dx = 0.0;
    fType dv_dx = 0.0;
    fType dw_dx = 0.0;
    fType dh_dx = 0.0;

    fType du_dy = 0.0;
    fType dv_dy = 0.0;
    fType dw_dy = 0.0;
    fType dh_dy = 0.0;

    fType du_dz = 0.0;
    fType dv_dz = 0.0;
    fType dw_dz = 0.0;
    fType dh_dz = 0.0;

    fType Lu = 0.0;
    fType Lv = 0.0;
    fType Lw = 0.0;
    fType Lh = 0.0;

    for (int j = 0; j < NNbrs; j++) {
      //_cq_unroll(31);
      const int dm_lind = i*NNbrs_Padded + j;
      const int inbr = idx[dm_lind];
      iu = U_SV_ID + NState*inbr;
      iv = V_SV_ID + NState*inbr;
      iw = W_SV_ID + NState*inbr;
      ih = H_SV_ID + NState*inbr;
    // temporary variables to hold state variable data for current node
      fType u = Stmp[iu];
      fType v = Stmp[iv];
      fType w = Stmp[iw];
      fType h = Stmp[ih];
      // get neighbor node's state vars
      // update sums
      du_dx += Dx[dm_lind] * u;
      dv_dx += Dx[dm_lind] * v;
      dw_dx += Dx[dm_lind] * w;
      dh_dx += Dx[dm_lind] * h;

      du_dy += Dy[dm_lind] * u;
      dv_dy += Dy[dm_lind] * v;
      dw_dy += Dy[dm_lind] * w;
      dh_dy += Dy[dm_lind] * h;

      du_dz += Dz[dm_lind] * u;
      dv_dz += Dz[dm_lind] * v;
      dw_dz += Dz[dm_lind] * w;
      dh_dz += Dz[dm_lind] * h;

      Lu += L[dm_lind] * u;
      Lv += L[dm_lind] * v;
      Lw += L[dm_lind] * w;
      Lh += L[dm_lind] * h;
    }

    iu = U_SV_ID + NState*i;
    iv = V_SV_ID + NState*i;
    iw = W_SV_ID + NState*i;
    ih = H_SV_ID + NState*i;

    u = Stmp[iu];
    v = Stmp[iv];
    w = Stmp[iw];
    h = Stmp[ih];
    fType rhs_u = - ((u * du_dx) + (v * du_dy) + (w * du_dz) + (f[i] * ((y[i] * w) - (z[i] * v))) + dh_dx);
    fType rhs_v = - ((u * dv_dx) + (v * dv_dy) + (w * dv_dz) + (f[i] * ((z[i] * u) - (x[i] * w))) + dh_dy);
    fType rhs_w = - ((u * dw_dx) + (v * dw_dy) + (w * dw_dz) + (f[i] * ((x[i] * v) - (y[i] * u))) + dh_dz);

    // vector linear ids
    int vect_idx = 3*i + X_COMP_ID;
    int vect_idy = 3*i + Y_COMP_ID;
    int vect_idz = 3*i + Z_COMP_ID;

    // evaluate projections and apply hyperviscosity
    K2[iu] = (p_u[vect_idx] * rhs_u) + (p_u[vect_idy] * rhs_v) + (p_u[vect_idz] * rhs_w) + Lu;
    K2[iv] = (p_v[vect_idx] * rhs_u) + (p_v[vect_idy] * rhs_v) + (p_v[vect_idz] * rhs_w) + Lv;
    K2[iw] = (p_w[vect_idx] * rhs_u) + (p_w[vect_idy] * rhs_v) + (p_w[vect_idz] * rhs_w) + Lw;
    // ----------------- Use RBF-FD Approximations to Calculate the RHS of the Geopotential Equation ---------------------//

    K2[ih] = - ((u * (dh_dx - gradghm[vect_idx])) + (v * (dh_dy - gradghm[vect_idy])) + (w * (dh_dz - gradghm[vect_idz]))
            + ((h + gh0 - ghm[i]) * (du_dx + dv_dy + dw_dz))) + Lh;
  }

  new_rhs +=  (cq_nstime() - t_start);
  //cq_timing->t_rhs +=  (cq_nstime() - t_start);
  //t->t_rhs +=  (getTime() - t_start);
 
  t_start = cq_nstime();
  //t_start = getTime();
  for(int i=0; i<NPts; i++){
    Stmp[i] = S[i] + 0.5*dt*K2[i];  // 2 flop/site
    dS[i] += 2.0*K2[i]; // 2 flop/site
  }
  //halo_update(Stmp);
  new_up +=  (cq_nstime() - t_start);
  //cq_timing->t_update +=  (cq_nstime() - t_start);
  //t->t_update +=  (getTime() - t_start);
  
  t_start = cq_nstime();
  //t_start = getTime();
  K3 = K;
  //fswe((*t)+0.5*dt, Stmp, K3);
  
  
  for(int i=0; i<NNodes; i++){
    _cq_ignorelcd();
    _cq_ignorelcd_pool(3);
    // initialize temporary variables for spatial differentiation at respective node
    fType du_dx = 0.0;
    fType dv_dx = 0.0;
    fType dw_dx = 0.0;
    fType dh_dx = 0.0;

    fType du_dy = 0.0;
    fType dv_dy = 0.0;
    fType dw_dy = 0.0;
    fType dh_dy = 0.0;

    fType du_dz = 0.0;
    fType dv_dz = 0.0;
    fType dw_dz = 0.0;
    fType dh_dz = 0.0;

    fType Lu = 0.0;
    fType Lv = 0.0;
    fType Lw = 0.0;
    fType Lh = 0.0;

    for (int j = 0; j < NNbrs; j++) {
      //_cq_unroll(31);
      const int dm_lind = i*NNbrs_Padded + j;
      const int inbr = idx[dm_lind];
      iu = U_SV_ID + NState*inbr;
      iv = V_SV_ID + NState*inbr;
      iw = W_SV_ID + NState*inbr;
      ih = H_SV_ID + NState*inbr;
    // temporary variables to hold state variable data for current node
      fType u = Stmp[iu];
      fType v = Stmp[iv];
      fType w = Stmp[iw];
      fType h = Stmp[ih];
      // get neighbor node's state vars
      // update sums
      du_dx += Dx[dm_lind] * u;
      dv_dx += Dx[dm_lind] * v;
      dw_dx += Dx[dm_lind] * w;
      dh_dx += Dx[dm_lind] * h;

      du_dy += Dy[dm_lind] * u;
      dv_dy += Dy[dm_lind] * v;
      dw_dy += Dy[dm_lind] * w;
      dh_dy += Dy[dm_lind] * h;

      du_dz += Dz[dm_lind] * u;
      dv_dz += Dz[dm_lind] * v;
      dw_dz += Dz[dm_lind] * w;
      dh_dz += Dz[dm_lind] * h;

      Lu += L[dm_lind] * u;
      Lv += L[dm_lind] * v;
      Lw += L[dm_lind] * w;
      Lh += L[dm_lind] * h;
    }

    iu = U_SV_ID + NState*i;
    iv = V_SV_ID + NState*i;
    iw = W_SV_ID + NState*i;
    ih = H_SV_ID + NState*i;

    u = Stmp[iu];
    v = Stmp[iv];
    w = Stmp[iw];
    h = Stmp[ih];
    fType rhs_u = - ((u * du_dx) + (v * du_dy) + (w * du_dz) + (f[i] * ((y[i] * w) - (z[i] * v))) + dh_dx);
    fType rhs_v = - ((u * dv_dx) + (v * dv_dy) + (w * dv_dz) + (f[i] * ((z[i] * u) - (x[i] * w))) + dh_dy);
    fType rhs_w = - ((u * dw_dx) + (v * dw_dy) + (w * dw_dz) + (f[i] * ((x[i] * v) - (y[i] * u))) + dh_dz);

    // vector linear ids
    int vect_idx = 3*i + X_COMP_ID;
    int vect_idy = 3*i + Y_COMP_ID;
    int vect_idz = 3*i + Z_COMP_ID;

    // evaluate projections and apply hyperviscosity
    K3[iu] = (p_u[vect_idx] * rhs_u) + (p_u[vect_idy] * rhs_v) + (p_u[vect_idz] * rhs_w) + Lu;
    K3[iv] = (p_v[vect_idx] * rhs_u) + (p_v[vect_idy] * rhs_v) + (p_v[vect_idz] * rhs_w) + Lv;
    K3[iw] = (p_w[vect_idx] * rhs_u) + (p_w[vect_idy] * rhs_v) + (p_w[vect_idz] * rhs_w) + Lw;
    // ----------------- Use RBF-FD Approximations to Calculate the RHS of the Geopotential Equation ---------------------//

    K3[ih] = - ((u * (dh_dx - gradghm[vect_idx])) + (v * (dh_dy - gradghm[vect_idy])) + (w * (dh_dz - gradghm[vect_idz]))
            + ((h + gh0 - ghm[i]) * (du_dx + dv_dy + dw_dz))) + Lh;
  }
  
  // t->t_rhs +=  (getTime() - t_start);
  new_rhs +=  (cq_nstime() - t_start);
  //cq_timing->t_rhs +=  (cq_nstime() - t_start);
  
  t_start = cq_nstime();
  //t_start = getTime();
  for(int i=0; i<NPts; i++){
    _cq_ignorelcd();
    Stmp[i] = S[i] + dt*K3[i];  // 2 flop/site
    dS[i] += 2.0*K3[i];         // 2 flop/site
  }
  //halo_update(Stmp);
  new_up +=  (cq_nstime() - t_start);
  //cq_timing->t_update +=  (cq_nstime() - t_start);
  //t->t_update +=  (getTime() - t_start);

  t_start = cq_nstime();
  //t_start = getTime();
  K4 = K;
  //fswe((*t)+dt, Stmp, K4);
  
  for(int i=0; i<NNodes; i++){
    _cq_ignorelcd();
    _cq_ignorelcd_pool(3);
    // initialize temporary variables for spatial differentiation at respective node
    fType du_dx = 0.0;
    fType dv_dx = 0.0;
    fType dw_dx = 0.0;
    fType dh_dx = 0.0;

    fType du_dy = 0.0;
    fType dv_dy = 0.0;
    fType dw_dy = 0.0;
    fType dh_dy = 0.0;

    fType du_dz = 0.0;
    fType dv_dz = 0.0;
    fType dw_dz = 0.0;
    fType dh_dz = 0.0;

    fType Lu = 0.0;
    fType Lv = 0.0;
    fType Lw = 0.0;
    fType Lh = 0.0;

    for (int j = 0; j < NNbrs; j++) {
      //_cq_unroll(31);
      const int dm_lind = i*NNbrs_Padded + j;
      const int inbr = idx[dm_lind];
      iu = U_SV_ID + NState*inbr;
      iv = V_SV_ID + NState*inbr;
      iw = W_SV_ID + NState*inbr;
      ih = H_SV_ID + NState*inbr;
    // temporary variables to hold state variable data for current node
      fType u = Stmp[iu];
      fType v = Stmp[iv];
      fType w = Stmp[iw];
      fType h = Stmp[ih];
      // get neighbor node's state vars
      // update sums
      du_dx += Dx[dm_lind] * u;
      dv_dx += Dx[dm_lind] * v;
      dw_dx += Dx[dm_lind] * w;
      dh_dx += Dx[dm_lind] * h;

      du_dy += Dy[dm_lind] * u;
      dv_dy += Dy[dm_lind] * v;
      dw_dy += Dy[dm_lind] * w;
      dh_dy += Dy[dm_lind] * h;

      du_dz += Dz[dm_lind] * u;
      dv_dz += Dz[dm_lind] * v;
      dw_dz += Dz[dm_lind] * w;
      dh_dz += Dz[dm_lind] * h;

      Lu += L[dm_lind] * u;
      Lv += L[dm_lind] * v;
      Lw += L[dm_lind] * w;
      Lh += L[dm_lind] * h;
    }

    iu = U_SV_ID + NState*i;
    iv = V_SV_ID + NState*i;
    iw = W_SV_ID + NState*i;
    ih = H_SV_ID + NState*i;

    u = Stmp[iu];
    v = Stmp[iv];
    w = Stmp[iw];
    h = Stmp[ih];
    fType rhs_u = - ((u * du_dx) + (v * du_dy) + (w * du_dz) + (f[i] * ((y[i] * w) - (z[i] * v))) + dh_dx);
    fType rhs_v = - ((u * dv_dx) + (v * dv_dy) + (w * dv_dz) + (f[i] * ((z[i] * u) - (x[i] * w))) + dh_dy);
    fType rhs_w = - ((u * dw_dx) + (v * dw_dy) + (w * dw_dz) + (f[i] * ((x[i] * v) - (y[i] * u))) + dh_dz);

    // vector linear ids
    int vect_idx = 3*i + X_COMP_ID;
    int vect_idy = 3*i + Y_COMP_ID;
    int vect_idz = 3*i + Z_COMP_ID;

    // evaluate projections and apply hyperviscosity
    K4[iu] = (p_u[vect_idx] * rhs_u) + (p_u[vect_idy] * rhs_v) + (p_u[vect_idz] * rhs_w) + Lu;
    K4[iv] = (p_v[vect_idx] * rhs_u) + (p_v[vect_idy] * rhs_v) + (p_v[vect_idz] * rhs_w) + Lv;
    K4[iw] = (p_w[vect_idx] * rhs_u) + (p_w[vect_idy] * rhs_v) + (p_w[vect_idz] * rhs_w) + Lw;
    // ----------------- Use RBF-FD Approximations to Calculate the RHS of the Geopotential Equation ---------------------//

    K4[ih] = - ((u * (dh_dx - gradghm[vect_idx])) + (v * (dh_dy - gradghm[vect_idy])) + (w * (dh_dz - gradghm[vect_idz]))
            + ((h + gh0 - ghm[i]) * (du_dx + dv_dy + dw_dz))) + Lh;
  }
  
  //t->t_rhs +=  (getTime() - t_start);
  new_rhs += cq_nstime() - t_start;
  //cq_timing->t_rhs += cq_nstime() - t_start;
  //t->t_rhs -= t_start1;
  t_start = cq_nstime();
  //t_start = getTime();
  for(int i=0; i<NPts; i++) {
    S[i] += (sixth*dt)*(dS[i]+K4[i]); //2 flop/site
  }
  //halo_update(S);
  new_up +=  (cq_nstime() - t_start);
  //cq_timing->t_update +=  (cq_nstime() - t_start);
  //t->t_update -= t_start1;
  //t->t_update +=  (getTime() - t_start);
  cq_timing->t_rhs = new_rhs;
  cq_timing->t_update = new_up;
  //free((CQ_POOL(3)*)dS);
  //free((CQ_POOL(3)*)Stmp);
  //free((void*)dS);
  //free((void*)Stmp);
  //(*t) += dt;
  //t->t_rhs += dt;
}

#else

void RK4_advance(const fType dt, fType *t, int NPts, fType *S, void f ( const fType t, const fType* S, fType* K)){

  const fType sixth = 1.0/6.0;

  fType *K1;
  fType *K2;
  fType *K3;
  fType *K4;

  fType *dS = (fType*) malloc(NPts * sizeof(fType));
  fType *Stmp = (fType*) malloc(NPts * sizeof(fType));
  fType *K = (fType*) malloc(NPts * sizeof(fType));

  double t_start = getTime();
  K1 = K;
  f((*t),S,K1);
  local_timer.t_rhs +=  (getTime() - t_start);

  t_start = getTime();
  for(int i=0; i<NPts; i++){
    Stmp[i] =  S[i] + 0.5*dt*K1[i]; // 2 flop/site
    dS[i] = K1[i];
  }
  halo_update(Stmp);
  local_timer.t_update +=  (getTime() - t_start);

  t_start = getTime();
  K2 = K;
  f((*t)+0.5*dt, Stmp, K2);
  local_timer.t_rhs +=  (getTime() - t_start);

  t_start = getTime();
  for(int i=0; i<NPts; i++){
    Stmp[i] = S[i] + 0.5*dt*K2[i];  // 2 flop/site
    dS[i] += 2.0*K2[i]; // 2 flop/site
  }
  halo_update(Stmp);
  local_timer.t_update +=  (getTime() - t_start);

  t_start = getTime();
  K3 = K;
  f((*t)+0.5*dt, Stmp, K3);
  local_timer.t_rhs +=  (getTime() - t_start);

  t_start = getTime();
  for(int i=0; i<NPts; i++){
    Stmp[i] = S[i] + dt*K3[i];  // 2 flop/site
    dS[i] += 2.0*K3[i];         // 2 flop/site
  }
  halo_update(Stmp);
  local_timer.t_update +=  (getTime() - t_start);

  t_start = getTime();
  K4 = K;
  f((*t)+dt, Stmp, K4);
  local_timer.t_rhs +=  (getTime() - t_start);

  t_start = getTime();
  for(int i=0; i<NPts; i++)
    S[i] += (sixth*dt)*(dS[i]+K4[i]); //2 flop/site
  halo_update(S);
  local_timer.t_update +=  (getTime() - t_start);

  free(K);
  free(dS);
  free(Stmp);

  (*t) += dt;
}
#endif
