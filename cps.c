/* 
 * cps.c --- written by Yasuko Matsubara
 *       --- updated by Takato Honda
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"
#include "hmm.h"
#include "tool.h"
#include "cubemarker.h"

#define DBG  1
#define S0  1
#define S1 -1

MarkerWS *ws;
double CPSearch();
double _search_aux();
void AllocCPS();
void FreeCPS();

void _printP(double *P, int k){
    int i;
    for(i=0;i<k;i++){
        fprintf(stderr, "%f ", P[i]);
    }
    fprintf(stderr, "\n");
}

void AllocCPS(CPS *cps, int maxk, int maxlen){
    fprintf(stderr, "alloc segment assignment...(k:%d,len:%d)\n", maxk, maxlen);
    cps->Pu = dvector(0, maxk);
    cps->Pv = dvector(0, maxk);
    cps->Pi = dvector(0, maxk);
    cps->Pj = dvector(0, maxk);
    cps->Su = imatrix(0, maxk, 0, maxlen);
    cps->Sv = imatrix(0, maxk, 0, maxlen);
    cps->Si = imatrix(0, maxk, 0, maxlen);
    cps->Sj = imatrix(0, maxk, 0, maxlen);
    cps->nSu = ivector(0, maxk);
    cps->nSv = ivector(0, maxk);
    cps->nSi = ivector(0, maxk);
    cps->nSj = ivector(0, maxk);
}

void FreeCPS(CPS *cps, int maxk, int maxlen){
    free_dvector(cps->Pu, 0, maxk);
    free_dvector(cps->Pv, 0, maxk);
    free_dvector(cps->Pi, 0, maxk);
    free_dvector(cps->Pj, 0, maxk);
    free_imatrix(cps->Su, 0, maxk, 0, maxlen);
    free_imatrix(cps->Sv, 0, maxk, 0, maxlen);
    free_imatrix(cps->Si, 0, maxk, 0, maxlen);
    free_imatrix(cps->Sj, 0, maxk, 0, maxlen);
    free_ivector(cps->nSu, 0, maxk);
    free_ivector(cps->nSv, 0, maxk);
    free_ivector(cps->nSi, 0, maxk);
    free_ivector(cps->nSj, 0, maxk);
}

double CPSearch(SegBox *Sx, SegBox *s0, SegBox *s1, MarkerWS *wsd){
    ws = wsd;
    ResetStEd(s0, ws->T.ntrip);
    ResetStEd(s1, ws->T.ntrip);
    double lh = 0.0;
    for(int i=0;i<Sx->totalnSb;i++){
        lh += _search_aux(Sx->Sb[i].st, Sx->Sb[i].len, Sx->Sb[i].trip, s0, s1);
    }
    return lh;
} 

int _findMax(double *P, int k){
    int i;
    int loc = -1;
    double max = -INF;
    for(i=0;i<k;i++){
        if(max < P[i]){ max = P[i]; loc = i; }
    }
    return loc; 
}

int _copy_path(int *from, int nfrom, int *to){
    int i;
    for(i=0;i<nfrom;i++){
        to[i]=from[i];
    }
    return nfrom;
}

void _reset_npaths(int *nS, int k){
    int i;
    for(i=0;i<k;i++){
        nS[i] = 0;
    }
}

void _print_path(FILE *fp, int *S, int len){
    int i;
    fprintf(fp, "---path---\n");
    for(i=0;i<len;i++){
        fprintf(fp, "%d\n", S[i]);
    }
    fprintf(fp, "----------\n");
}

double _estimate_b(HMM *m, int i, int t, int trip){
    double b=0;
    b = pdfL(m,i,ws->T.x[trip].O[t]);
    return b;
}

double _search_aux(int st, int len, int trip, SegBox *s0, SegBox *s1){
    int i,j,u,v,t;
    int k0 = s0->model.k;
    int k1 = s1->model.k;
    ///--- setting ---///
    // float **O = ws->T.x[trip].O;
    HMM *m0 = &s0->model;
    HMM *m1 = &s1->model;
    double d0 = s0->delta;
    double d1 = s1->delta; 
    ///--- setting ---///
    double *Pu = ws->cps.Pu;
    double *Pv = ws->cps.Pv;
    double *Pi = ws->cps.Pi;
    double *Pj = ws->cps.Pj;
    int **Su = ws->cps.Su;
    int **Sv = ws->cps.Sv;
    int **Si = ws->cps.Si;
    int **Sj = ws->cps.Sj;  
    int *nSu = ws->cps.nSu;
    int *nSv = ws->cps.nSv;
    int *nSi = ws->cps.nSi;
    int *nSj = ws->cps.nSj;
    _reset_npaths(nSu, k0);
    _reset_npaths(nSv, k0);
    _reset_npaths(nSi, k1);
    _reset_npaths(nSj, k1);
    // for swap
    // HMM *hmm;
    int **imat;
    int *ivec;
    double *dvec;
    ///--- setting ---///

    if(d0<=0 || d1<=0){
        error("delta == 0","cpsearch");
    }

    // (t == 0);
    t=st;
    for(v=0;v<k0;v++){
        Pv[v] = log(d1) + log(m0->pi[v]+ZERO) + _estimate_b(m0,v,t, trip);
    }
    for(j=0;j<k1;j++){
        Pj[j] = log(d0) + log(m1->pi[j]+ZERO) + _estimate_b(m1,j,t, trip);
    }
    
    // (t >= 1);
    for(t=st+1;t<st+len;t++){
        ///---  Pu(t)  ---///
        // regime-switch
        int maxj = _findMax(Pj, k1);

        for(u=0;u<k0;u++){
            // if switch
            double maxPj = Pj[maxj] + log(d1) + log(m0->pi[u]+ZERO) + _estimate_b(m0,u,t, trip);

            // else
            double maxPv = -INF; int maxv;

            for(v=0;v<k0;v++){
                double val = log(1.0-d0) + Pv[v] + log(m0->A[v][u]+ZERO) + _estimate_b(m0,u,t, trip);
                if(val > maxPv){ maxPv = val; maxv = v; }
            }//v

            // if switch
            if(maxPj > maxPv){
                Pu[u] = maxPj;
                nSu[u] = _copy_path(Sj[maxj], nSj[maxj], Su[u]);
                Su[u][nSu[u]] = t; nSu[u]++;
            }else{
                Pu[u] = maxPv;
                nSu[u] = _copy_path(Sv[maxv], nSv[maxv], Su[u]);
            }
        }//u

        ///---  Pi(t)  ---///
        // regime-switch
        int maxv = _findMax(Pv, k0);

        for(i=0;i<k1;i++){
            // if switch
            double maxPv = Pv[maxv] + log(d0) + log(m1->pi[i]+ZERO) + _estimate_b(m1,i,t,trip);
            //else
            double maxPj = -INF; int maxj;
            for(j=0;j<k1;j++){
                double val = log(1.0-d1) + Pj[j] + log(m1->A[j][i]+ZERO) + _estimate_b(m1,i,t, trip);
                //update
                if(val > maxPj){ maxPj = val; maxj = j; }
            }//j

            //if switch
            if(maxPv > maxPj){
                Pi[i] = maxPv;
                nSi[i] = _copy_path(Sv[maxv], nSv[maxv], Si[i]);
                Si[i][nSi[i]] = t; nSi[i]++;
            }else{
                Pi[i] = maxPj;
                nSi[i] = _copy_path(Sj[maxj], nSj[maxj], Si[i]);
            }
        }//i

        // swap Pu Pv; Pi Pj;
        dvec = Pu; Pu = Pv; Pv = dvec;
        dvec = Pi; Pi = Pj; Pj = dvec;
        // swap Su Sv; Si Sj;
        imat = Su; Su = Sv; Sv = imat;
        imat = Si; Si = Sj; Sj = imat;
        // swap nSu nSv; nSi nSj;
        ivec = nSu; nSu = nSv; nSv = ivec;
        ivec = nSi; nSi = nSj; nSj = ivec;
        // each trip
    }//t

    // check the best path
    int maxv = _findMax(Pv, k0);
    int maxj = _findMax(Pj, k1);

    int *path;
    int npath;
    int firstID;
    double lh = 0;
    if(Pv[maxv] > Pj[maxj]){
        // fixPath(maxv);
        path = Sv[maxv];
        npath = nSv[maxv];
        firstID=pow(-1,npath)*S0;
        lh = Pv[maxv];
    }else{
        path = Sj[maxj];
        npath = nSj[maxj];
        firstID=pow(-1,npath)*S1;
        lh = Pj[maxj];
    }
    // add paths
    int curSt=st; int nxtSt;
    for(i=0;i<npath;i++){
        nxtSt = path[i];
        if(firstID*pow(-1,i)==S0){
            AddStEd(s0, curSt, nxtSt-curSt, trip, ws->lmax);
        }else{
            AddStEd(s1, curSt, nxtSt-curSt, trip, ws->lmax);
        }
        curSt = nxtSt;
    }//i
    if(firstID*pow(-1,npath)==S0){
        AddStEd(s0, curSt, st+len-curSt, trip, ws->lmax);
    }else{
        AddStEd(s1, curSt, st+len-curSt, trip, ws->lmax);
    }
    double costC = -lh/log(2.0);

    return costC;
}



