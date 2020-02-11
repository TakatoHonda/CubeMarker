#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include "nrutil.h"
#include "hmm.h"
#include "tool.h"
#include "cubemarker.h"

//! debug
#define DBG_VSPLIT FALSE
#define DBG_HSPLIT FALSE
#define DBG_RM FALSE
#define DBG_RGMREF FALSE
#define SNAP FALSE

//! for experiments
#define RN TRUE  // RemoveSegmentNoise
#define REFINE FALSE  // RegimeRefinement
#define HSP TRUE  // H-Split
#define INI_SP_REST 2
#define INFER_ITER_MIN 2
#define INFER_ITER_MAX 10

//! minimum/maximum value for optimization
#define MAXR 16  // maximum number of regimes
#define MINK 1  // minimum number of K
#define LM 0.1  // DEFAULT: 0.1 for sampling

#define SPLIT_R  3.0
#define REGIME_R  0.03  // for V/H-Split
#define SEGMENT_R  0.03  // for RemoveSegmentNoise
#define INTG_R  0.03     // for RegimeRefinement

//! output option
#define PRINTHMM 1

/// main functions ///
void _marker();
void _multiSplit();
void _vsplit();
void _vassignment();
void _hsplit();
void _hassignment();

/// find regimes & model estimation ///
void _vsplit_aux();
void _hsplit_aux();
void _estimateHMM_k(SegBox *s, int k);
void _estimateHMM();
double _findCentroid();
double _computeDiffVR();
int _findOptSeedLen();
void _removeSegmentNoise();
void _regimeRefinement();
void _selectLargest();

/// compute Vit & MDL ///
double _viterbi();
double _MDLtotal();
void _computeLhMDL();
double _computeLhCost();
double _MDL();

/// output ///
void _printSplit();

/// segment set --- get & release ///
void _setSegment();
SegBox *_getS();
void _releaseS();

/// memory allocation ///
void AllocMarker();

//! workspace
MarkerWS *ws;

/// main functions ///
void Marker(MarkerWS *pws){
    ws = pws;
    fprintf(stdout, "---------- \n");
    fprintf(stdout, "i|r|m|Cost \n");
    fprintf(stdout, "---------- \n");
    _marker();
    printf("\n");
}

void _marker(){
    // initialize regime
    SegBox *Sx = _getS("", "");
    for(int i=0; i<ws->T.ntrip; i++){
        AddStEd(Sx, 0, ws->T.x[i].m, i, ws->lmax);
    }
    _estimateHMM(Sx);
    Push(Sx, &ws->C);

    for(int iter=1; TRUE; iter++){
        // save snapshot for each iteration 
        #if(SNAP)
            char dir[100];
            char buf[32];
            IntToString(buf, iter);
            strcat(strcat(strcpy(dir, ws->fout), buf), "/");
            SaveMarker(ws, dir);
            CopyFile(ws->fin, strcat(dir, "input")); 
        #endif

        // compute total cost
        ws->costT = _MDLtotal(&ws->Opt, &ws->C, ws->T.ntrip, iter);
        
        // get target regime in stack C
        SegBox *Sx = Pop(&ws->C);
        
        // if converged: break
        if(Sx == NULL) break;

        // create new regimes s0 and s1
        SegBox *s0 = _getS(Sx->label, "0");
        SegBox *s1 = _getS(Sx->label, "1"); 

        // try to split
        _multiSplit(Sx, s0, s1, iter);

        // splitting result evaluation
        bool splitFlag = FALSE;
        int split = 0; int keep = 0;
        for(int i=0; i<Sx->ntrip; i++){
            // if # of regimes > maximum # of regimes
            if(ws->Opt.idx + ws->C.idx + 1 >= MAXR){
                keep = INF;
                break;
            }
            // if cannot split
            if(s0->costT == INF || s1->costT == INF){
                keep = INF;
                break;
            }
            // if cannot split
            if(Sx->len[Sx->trip[i]] == ws->lmax){
                split = INF;
                break;
            }

            // compare cost for each regime
            double costT_s01 = 0.0;
            double costT_sx = Sx->costT_s[Sx->trip[i]];
            if (s0->thash[Sx->trip[i]]) costT_s01 += s0->costT_s[Sx->trip[i]];
            if (s1->thash[Sx->trip[i]]) costT_s01 += s1->costT_s[Sx->trip[i]];
            if(costT_s01 != 0.0 && costT_s01 + costT_sx * REGIME_R < costT_sx){
                split++;
            }else{
                keep++;
            }
        }
        // compare total cost
        if (s0->costT + s1->costT < Sx->costT) split++;

        // split or not
        if (split > keep) splitFlag = TRUE;

        // split regime
        if(splitFlag){
            Push(s0, &ws->C);
            Push(s1, &ws->C);
            _releaseS(Sx);
        // do not split
        }else{
            Push(Sx, &ws->Opt);
            _releaseS(s0);
            _releaseS(s1);
        }
    }

    // regime refinement algorithm (unused)
    if(REFINE){
        _regimeRefinement();
        ws->costT = _MDLtotal(&ws->Opt, &ws->C, ws->T.ntrip, -1);
    }
}

void _multiSplit(SegBox *Sx, SegBox *s0, SegBox *s1, int iter){
    // initialize s0 and s1
    int seedlen = (int)ceil(ws->lmax*LM);
    _findCentroid(Sx, s0, s1, NSAMPLE, seedlen);

    // generate initial model
    SegBox *s0_v = _getS("", ""); SegBox *s1_v = _getS("", "");
    SegBox *s0_h = _getS("", ""); SegBox *s1_h = _getS("", "");
    CopyStEd(s0, s0_v); CopyStEd(s1, s1_v);
    CopyStEd(s0, s0_h); CopyStEd(s1, s1_h);
    _estimateHMM_k(s0_v, MINK); _estimateHMM_k(s1_v, MINK);
    _estimateHMM_k(s0_h, MINK); _estimateHMM_k(s1_h, MINK);

    // try to V-Split and H-Split
    _vassignment(Sx, s0_v, s1_v, RN, FALSE);
    _hassignment(Sx, s0_h, s1_h);

    // recompute costs 
    if(s0_v->totalnSb > 0 && s1_v->totalnSb > 0){
        _estimateHMM_k(s0_v, MINK);
        _estimateHMM_k(s1_v, MINK);
        _computeLhMDL(s0_v);
        _computeLhMDL(s1_v);
    }
    if(s0_h->totalnSb > 0 && s1_h->totalnSb > 0){
        _estimateHMM_k(s0_h, MINK);
        _estimateHMM_k(s1_h, MINK);
        _computeLhMDL(s0_h);
        _computeLhMDL(s1_h);
    }

    // compare V cost and H cost
    double cost_v = s0_v->costT + s1_v->costT;
    double cost_h = s0_h->costT + s1_h->costT;
    if(cost_v > cost_h && HSP && (INI_SP_REST > iter)){
        if(s0_h->totalnSb==0 || s1_h->totalnSb==0) return;
        // do H-Split
        CopyStEd(s0_h, s0); CopyStEd(s1_h, s1);
        _selectLargest(s0); _selectLargest(s1);
        _estimateHMM_k(s0, MINK); _estimateHMM_k(s1, MINK);
        _hsplit(Sx, s0, s1);
    }else{
        // do V-Split
        if(s0_v->totalnSb==0 || s1_v->totalnSb==0) return;
        CopyStEd(s0_v, s0); CopyStEd(s1_v, s1);
        _selectLargest(s0); _selectLargest(s1);
        _estimateHMM_k(s0, MINK); _estimateHMM_k(s1, MINK);
        _vsplit(Sx, s0, s1);
    }

    _releaseS(s0_v); _releaseS(s1_v);
    _releaseS(s0_h); _releaseS(s1_h);
}

void _vsplit(SegBox *Sx, SegBox *s0, SegBox *s1){
    // main
    _vsplit_aux(Sx, s0, s1);
    if(s0->totalnSb==0 || s1->totalnSb==0) return;
    
    // final model estimation
    _estimateHMM(s0);
    _estimateHMM(s1);

    #if(DBG_VSPLIT)
        _printSplit(stderr, Sx, s0, s1);
    #endif
}

void _vsplit_aux(SegBox *Sx, SegBox *s0, SegBox *s1){
    SegBox *opt0 = _getS("",""), *opt1 = _getS("","");
    for(int i=0; i<INFER_ITER_MAX; i++){
        
        // Phase 1: Find cut-points (SegmentAssignment)
        _vassignment(Sx, s0, s1, RN, TRUE);
        if(s0->totalnSb==0 || s1->totalnSb==0) break;
        
        // if improving, update the opt segment set
        double diff = (opt0->costT+opt1->costT) - (s0->costT+s1->costT);
        if(diff > 0){ CopyStEd(s0, opt0); CopyStEd(s1, opt1); }
        
        // if not improving, then break iteration (efficient convergent)
        else if(i>=INFER_ITER_MIN) break;

        // Phase 2: Model Update
        _selectLargest(s0);
        _selectLargest(s1);
        _estimateHMM_k(s0, MINK);
        _estimateHMM_k(s1, MINK);
    }

    CopyStEd(opt0, s0); CopyStEd(opt1, s1);
    _releaseS(opt0); _releaseS(opt1);
}

void _vassignment(SegBox *Sx, SegBox *s0, SegBox *s1, int RM, bool debug){
    CPSearch(Sx, s0, s1, ws);

    if(RM) _removeSegmentNoise(Sx, s0, s1, debug);

    _computeLhMDL(s0);
    _computeLhMDL(s1);
}

void _hsplit(SegBox *Sx, SegBox *s0, SegBox *s1){
    // main
    _hsplit_aux(Sx, s0, s1);
    if(s0->totalnSb==0 || s1->totalnSb==0) return;
    
    // final model estimation
    _estimateHMM(s0);
    _estimateHMM(s1);

    #if(DBG_HSPLIT)
        _printSplit(stderr, Sx, s0, s1);
    #endif
}

void _hsplit_aux(SegBox *Sx, SegBox *s0, SegBox *s1){
    SegBox *opt0 = _getS("",""), *opt1 = _getS("","");
    for(int i=0; i<INFER_ITER_MAX; i++){
        
        // Phase 1: RegimeAssignment
        _hassignment(Sx, s0, s1);
        if(s0->totalnSb==0 || s1->totalnSb==0) break;
        
        // if improving, update the opt segment set
        double diff = (opt0->costT+opt1->costT) - (s0->costT+s1->costT);
        if(diff > 0){ CopyStEd(s0, opt0); CopyStEd(s1, opt1); }
        // if not improving, then break iteration (efficient convergent)
        else if(i>=INFER_ITER_MIN) break;

        // Phase 2: Model Update
        _selectLargest(s0);
        _selectLargest(s1);
        _estimateHMM_k(s0, MINK);
        _estimateHMM_k(s1, MINK);
    }

    CopyStEd(opt0, s0); CopyStEd(opt1, s1);
    _releaseS(opt0); _releaseS(opt1);
}

void _hassignment(SegBox *Sx, SegBox *s0, SegBox *s1){
    // similarity score for each sequence
    int w = Sx->ntrip;
    double cost_s0[ws->fmax], cost_s1[ws->fmax];
    for(int i=0; i<w; i++){
        cost_s0[Sx->trip[i]] = cost_s1[Sx->trip[i]] = 0.0;
    }

    // compute similarity
    for(int i=0; i<Sx->totalnSb; i++){
        SubS segment = Sx->Sb[i];
        int sid = segment.trip;  // sequence ID

        cost_s0[sid] += _computeLhCost(segment, s0);
        cost_s1[sid] += _computeLhCost(segment, s1);
    }

    // assign regimes
    ResetStEd(s0, ws->T.ntrip);
    ResetStEd(s1, ws->T.ntrip);
    for(int i=0; i<Sx->totalnSb; i++){
        SubS segment = Sx->Sb[i];
        int sid = segment.trip;
        int st  = segment.st;
        int len = segment.len;

        if(cost_s0[sid] < cost_s1[sid]){
            AddStEd(s0, st, len, sid, ws->lmax);
        }else{
            AddStEd(s1, st, len, sid, ws->lmax);
        }
    }

    // compute regime costs
    _computeLhMDL(s0);
    _computeLhMDL(s1);
}

double _computeLhCost(SubS segment, SegBox *s){
    double cost = 0;
    int st = segment.st;
    int len = segment.len;
    int sid = segment.trip;
    cost = _viterbi(&s->model, s->delta, st, len, &ws->vit, sid) - s->costC;
    return cost;
}

/// sub functions ///
void _selectLargest(SegBox *s){
    int id = FindMaxStEd(s);
    int st = s->Sb[id].st;
    int len = s->Sb[id].len;
    int trip = s->Sb[id].trip;
    ResetStEd(s, ws->T.ntrip);
    AddStEd(s, st, len, trip, ws->lmax);
}

int _findMinDiff(SegBox *s0, SegBox *s1, double *diffp){
    double min = INF;
    int loc = -1;
    
    for(int i=0;i<s0->totalnSb;i++){
        int st  = s0->Sb[i].st;
        int len = s0->Sb[i].len;
        int trip = s0->Sb[i].trip;
        double costC0 = _viterbi(&s0->model, s0->delta, st, len, &ws->vit, trip);
        double costC1 = _viterbi(&s1->model, s1->delta, st, len, &ws->vit, trip);
        double diff = costC1-costC0;
        if(min > diff) {loc = i; min = diff; }
    }//i
    *diffp = min;
    return loc;
}

double _scanMinDiff(SegBox *Sx, SegBox *s0, SegBox *s1){
    double diff;

    int loc0 = _findMinDiff(s0, s1, &diff);
    int loc1 = _findMinDiff(s1, s0, &diff);
    if(loc0==-1 || loc1==-1) return INF;

    SegBox *tmp0 = _getS("",""); SegBox *tmp1 = _getS("","");
    AddStEd(tmp0, s0->Sb[loc0].st, s0->Sb[loc0].len, s0->Sb[loc0].trip, ws->lmax);
    AddStEd(tmp1, s1->Sb[loc1].st, s1->Sb[loc1].len, s1->Sb[loc1].trip, ws->lmax);
    _estimateHMM_k(tmp0, MINK); _estimateHMM_k(tmp1, MINK);
   
    double costC = CPSearch(Sx, tmp0, tmp1, ws);
    _releaseS(tmp0); _releaseS(tmp1);
   
    return costC;
}

void _removeSegmentNoise_aux(SegBox *Sx, SegBox *s0, SegBox *s1, double per){
    if(per==0) return;
    int mprev=INF;
    double th = (ws->costT/ws->T.ntrip)*per; // threshold
    // double th = (ws->costT/Sx->ntrip)*per;
    // double th = (ws->costT)*per;  // threshold

    while(mprev > s0->totalnSb+s1->totalnSb){
        mprev = s0->totalnSb+s1->totalnSb;
       
        // find minimum segment
        double diff0 = INF, diff1 = INF;
        int loc0 = _findMinDiff(s0, s1, &diff0);
        int loc1 = _findMinDiff(s1, s0, &diff1);
        if(loc0==-1 || loc1==-1) return;

        // check convergence
        bool rmflag = FALSE;
        for(int i=0; i<Sx->ntrip; i++)
            if(s0->nSb[Sx->trip[i]]>1 && s1->nSb[Sx->trip[i]]>1) rmflag = TRUE;
        if(!rmflag) return;

        // check s0->s1 or s1->s0
        double min; int id;
        if(diff0 < diff1){ 
            min = diff0; id = 0;
        }else{ 
            min = diff1; id = 1;
        }

        // check remove or not
        if(min<th){
            // from s0 to s1
            if(id==0){
                AddStEd(s1, s0->Sb[loc0].st, s0->Sb[loc0].len, s0->Sb[loc0].trip, ws->lmax); 
                RemoveStEd(s0,loc0, s0->Sb[loc0].trip);
            // from s1 to s0
            }else{
                AddStEd(s0, s1->Sb[loc1].st, s1->Sb[loc1].len, s1->Sb[loc1].trip, ws->lmax);
                RemoveStEd(s1,loc1, s1->Sb[loc1].trip);
            }
            // printf("removeNoise_aux: %0.f < %0.f\n", min, th);
            // _printRemoveSegmentNoise(s0, s1);
            // printf("end============\n");
       }
    }//while
}

void _printRemoveSegmentNoise(SegBox *s0, SegBox *s1){
    int cur = -1, pre = -1;
    printf("s0:\n");
    for(int i=0; i<s0->totalnSb; i++){
        cur = s0->Sb[i].trip;
        if(cur!=pre){
            if(i != 0) printf("]\n");
            printf("%d [", cur);
        }
        int st = s0->Sb[i].st;
        int ed = st + s0->Sb[i].len;
        printf(" (%d-%d) ", st, ed);
        pre = cur;
    }
    printf("]\n");
    cur = pre = -1;
    printf("s1:\n");
    for(int i=0; i<s1->totalnSb; i++){
        cur = s1->Sb[i].trip;
        if(cur!=pre){
            if(i != 0) printf("]\n");
            printf("%d [", cur);
        }
        int st = s1->Sb[i].st;
        int ed = st + s1->Sb[i].len;
        printf(" (%d-%d) ", st, ed);
        pre = cur;
    }
    printf("]\n");
}

void _removeSegmentNoise(SegBox *Sx, SegBox *s0, SegBox *s1, bool debug){
    if(DBG_RM && debug){
        printf("========== _removeSegmentNoise ==========\n");
      _printRemoveSegmentNoise(s0, s1);
      printf("-----\n");
    }

    // remove noise or not
    bool rmflag = FALSE;
    for(int i=0; i<Sx->ntrip; i++)
        if(s0->nSb[Sx->trip[i]]>1 && s1->nSb[Sx->trip[i]]>1)
            rmflag = TRUE;
    if(!rmflag) return;

    // default pruning
    double per = SEGMENT_R;
    _removeSegmentNoise_aux(Sx, s0, s1, per);
    double costC=_scanMinDiff(Sx, s0, s1);

    // opt: optimal segment set
    SegBox *opt0= _getS("",""); SegBox *opt1 = _getS("","");
    CopyStEd(s0, opt0); CopyStEd(s1, opt1);
    double prev = costC;
    
    // find optimal pruning point
    while(per<=SEGMENT_R*10){
        if(costC>=INF) break;
        per*=2.0;

        _removeSegmentNoise_aux(Sx, s0, s1, per);

        // check convergence
        rmflag = FALSE;
        for(int i=0; i<Sx->ntrip; i++)
            if(s0->nSb[Sx->trip[i]] > 1 && s1->nSb[Sx->trip[i]] > 1) rmflag = TRUE;
        if(!rmflag) break;

        // if improved the cost
        costC=_scanMinDiff(Sx, s0, s1);
        if(prev <= costC) break;
        
        CopyStEd(s0, opt0);
        CopyStEd(s1, opt1);
        prev = costC;
    }

    CopyStEd(opt0, s0); CopyStEd(opt1,s1);
    _releaseS(opt0);  _releaseS(opt1);

    if(DBG_RM && debug){
        _printRemoveSegmentNoise(s0, s1);
        printf("========== _removeSegmentNoise =======end\n");
    }
}

void _removeSegBox(int loc, Stack *S){
    _releaseS(S->s[loc]);
    for(int i=loc; i<S->idx; i++){
        S->s[i] = S->s[i+1];
    }
    S->idx--;
}

double _findMinDiffVR(Stack *S, int *loc0, int *loc1){
    double best = INF;
    double cost = 0.0;
    for(int i=0; i<S->idx-1; i++){
        for(int j=i+1; j<S->idx; j++){
            SegBox *s0 = S->s[i];
            SegBox *s1 = S->s[j];
            
            cost = _computeDiffVR(s0, s1);
            // printf("cost: %0.f [%d %d]\n", cost, i,j);
            if(cost < best){
                best = cost;
                *loc0 = i;
                *loc1 = j;
            }
        }
    }
    return best;
}

double _regimeIntegration(Stack *S, int loc0, int loc1, SegBox *tmp){
    // regime # loc0
    for(int i=0; i<S->s[loc0]->totalnSb; i++){
        int st = S->s[loc0]->Sb[i].st;
        int len = S->s[loc0]->Sb[i].len;
        int trip = S->s[loc0]->Sb[i].trip;
        AddStEd(tmp, st, len, trip, ws->lmax);
    }
    // regime # loc1
    for(int i=0; i<S->s[loc1]->totalnSb; i++){
        int st = S->s[loc1]->Sb[i].st;
        int len = S->s[loc1]->Sb[i].len;
        int trip = S->s[loc1]->Sb[i].trip;
        AddStEd(tmp, st, len, trip, ws->lmax);
    }
    _estimateHMM(tmp);
    return tmp->costT;
}

void _regimeRefinement(){
    #if(DBG_RGMREF)
    printf("===== Regime Refinement =====\n");
    printf("r:%d m:%d\n", ws->Opt.idx, mSegment(&ws->Opt));
    printf("-----\n");
    #endif

    while(ws->Opt.idx > 1){
        // find similar regime set
        int loc0, loc1; // loc0 < loc1
        loc0 = loc1 = -1;
        _findMinDiffVR(&ws->Opt, &loc0, &loc1);
        
        if(loc0 == -1 || loc1 == -1) break;
        printf("similar V-regimes:(%d, %d)\n", loc0, loc1);
        
        // cost computation
        SegBox *tmp = _getS("", ws->Opt.s[loc0]->label);
        int m = mSegment(&ws->Opt);
        int r  = ws->Opt.idx;
        
        // old cost
        double alphaOld = log_s(r) + m*log_2(r) + FB*r*r;
        double oldCost = (ws->Opt.s[loc0]->costT+ws->Opt.s[loc1]->costT);

        // new cost
        double curCost = _regimeIntegration(&ws->Opt, loc0, loc1, tmp);
        int mNew = tmp->totalnSb;
        double alphaNew = log_s(r-1) + mNew*log_2(r-1) + FB*(r-1)*(r-1);
        printf("%0.f(%0.f+%0.f)  vs. %0.f(%0.f+%0.f+%0.f)\n",
            curCost+alphaNew, curCost, alphaNew, 
            oldCost+alphaOld*(1+INTG_R), oldCost, alphaOld, (oldCost+alphaOld)*INTG_R);

        // if decrease the cost
        // if(curCost + alphaNew < (oldCost+alphaOld) * (1+INTG_R)){
        if(curCost < oldCost * (1+INTG_R)){
        // if(ws->Opt.idx > MAXR){
            // printf("# INTEGRATION #\n");
            _removeSegBox(loc1, &ws->Opt);
            _removeSegBox(loc0, &ws->Opt);
            Push(tmp, &ws->Opt);
        }else{
            // printf("# GOOD ENOUGH #\n");
            break;
        }
        
    }

    #if(DBG_RGMREF)
    printf("r:%d m:%d\n", ws->Opt.idx, mSegment(&ws->Opt));
    printf("===== Regime Refinement =====\n");
    #endif
}

double _viterbi(HMM *phmm, double delta, int st, int len, VITERBI *vit, int trip){
    double Lh = ViterbiL(phmm, len, ws->T.x[trip].O+st, vit);
    if(delta<=0 || delta >= 1){
        error("not appropriate delta", "_computeLhMDL");
    }
    Lh += log(delta); // switch
    Lh += (len-1)*log(1.0-delta); // else (stay)
    double costC = -Lh/log(2.0);
    return costC;
}

void _estimateHMM_k(SegBox *s, int k){
    if(k<MINK) k=MINK;
    if(k>MAXK) k=MAXK;
    int wantKMEANS=1;
    int n = s->totalnSb;
    /// if nSb is too big, use small size
    if(n > MAXBAUMN){
        n=MAXBAUMN;
    }
    for(int i=0; i<n; i++){
        int st   = s->Sb[i].st;
        int len  = s->Sb[i].len;
        int trip = s->Sb[i].trip;
        Split(&ws->T.x[trip], st, len, &ws->x_tmp[i]);
    }

    // select longest segment
    // int n = 1;
    // int id = FindMaxStEd(s);
    // int st = s->Sb[id].st;
    // int len = s->Sb[id].len;
    // int trip = s->Sb[id].trip;
    // Split(&ws->T.x[trip], st, len, &ws->x_tmp[0]);

    s->model.k = k;
    s->model.n = 0;
    ResetHMM(&s->model, k, ws->d);
    BaumWelch(&s->model, n, ws->x_tmp, &ws->baum, wantKMEANS);
    s->delta = (double)s->totalnSb/(double)s->totallen; //ws->lmax;
    // s->delta = (double)s->totallen/(double)ws->lmax;
    if(s->delta <= 0) s->delta = ZERO;
    else if(s->delta >= 1) s->delta = ONE;
}

void _estimateHMM(SegBox *s){
    s->costT = INF;
    int optk = MAXK;
    if (s->totalnSb <= 0) return;
    for(int k=MINK; k<=MAXK; k++){
        double prev = s->costT;
        _estimateHMM_k(s, k);
        _computeLhMDL(s);
        if(s->costT > prev){ optk=k-1; break; }
    }
    if(optk<MINK){ optk=MINK; }
    if(optk>MAXK){ optk=MAXK; }
    _estimateHMM_k(s, optk);
    _computeLhMDL(s);
}

/**
 * @fn
 * select initial models s0 and s1 with sampling
 * @param Sx: target regime
 * @param s0: initial model
 * @param s1: initial model
 * @param nsamples: # of samples for initial model
 * @param seedlen: sample segment length
 */
double _findCentroid(SegBox *Sx, SegBox *s0, SegBox *s1, int nsamples, int seedlen){
    // double costMax = 0.0;
    double costMin = INF;

    /// keep best seeds
    int s0stB, s1stB, s0lenB, s1lenB, s0tripB, s1tripB; //best
    int s0stC, s1stC, s0lenC, s1lenC, s0tripC, s1tripC; //current

    /// make sample set
    UniformSet(Sx, seedlen, nsamples, &ws->U, ws->T.ntrip);
    /// start uniform sampling
    for(int iter1=0; iter1<ws->U.totalnSb; iter1++){
        for(int iter2=iter1+1; iter2<ws->U.totalnSb; iter2++){
            // get samples s0 and s1 from U
            UniformSampling(Sx, s0, s1, seedlen, iter1, iter2, &ws->U, ws->lmax, ws->T.ntrip);
            if(s0->totalnSb==0 || s1->totalnSb==0){ continue; } // not sufficient

            // copy positions
            s0stC=s0->Sb[0].st; s0lenC=s0->Sb[0].len; s0tripC=s0->Sb[0].trip;
            s1stC=s1->Sb[0].st; s1lenC=s1->Sb[0].len; s1tripC=s1->Sb[0].trip;

            // model initialization
            _estimateHMM_k(s0, MINK);
            _estimateHMM_k(s1, MINK);

            // segmentation and clustering
            // _vassignment(Sx, s0, s1, FALSE, FALSE);
            // _hassignment(Sx, s0, s1);
            // double curCost = s0->costT+s1->costT;
            double curCost = -_computeDiffVR(s0, s1);
            
            if(s0->totalnSb==0 || s1->totalnSb==0) continue;
            if(costMin > curCost){
            // if(costMax < curCost){
                // update best seeds
                costMin = curCost;
                s0stB=s0stC; s0lenB=s0lenC; s0tripB=s0tripC;
                s1stB=s1stC; s1lenB=s1lenC; s1tripB=s1tripC;
            }
        }
    }
    if(costMin==INF){
    // if(costMax==0.0){
        FixedSampling(Sx, s0, s1, seedlen, ws->lmax, ws->T.ntrip);
        return INF;
    }

    ResetStEd(s0, ws->T.ntrip); ResetStEd(s1, ws->T.ntrip);
    AddStEd(s0, s0stB, s0lenB, s0tripB, ws->lmax);
    AddStEd(s1, s1stB, s1lenB, s1tripB, ws->lmax);
    return costMin;
}

void _computeLhMDL(SegBox *s){
    int k = s->model.k;
    int d = s->model.d;
    // int *m = s->nSb;
    int *m = s->len;

    if(s->totalnSb==0){
        ResetRegimeCost(s, INF);
        return;
    }

    // compute costT, costM, costC and costLen
    ResetRegimeCost(s, 0.0);
    s->costM = costHMM(k, d);
    for(int i=0; i<s->totalnSb; i++){
        int st = s->Sb[i].st;
        int len = s->Sb[i].len;
        int trip = s->Sb[i].trip;
        double costC = _viterbi(&s->model, s->delta, st, len, &ws->vit, trip);
        double costLen = log_2(s->Sb[i].len);
        s->costC += costC;
        s->costC_s[trip] += costC;
        s->costLen += costLen;
        s->costLen_s[trip] += costLen;
    }
    for(int i=0; i<s->ntrip; i++){
        double costLen = m[s->trip[i]]*log_2(k);
        s->costLen += costLen;
        s->costLen_s[s->trip[i]] += costLen;
    }
    // s->costC /= s->ntrip;
    // s->costLen /= s->ntrip;
    for(int i=0; i<s->ntrip; i++){
        int sid = s->trip[i];
        s->costT_s[sid] = s->costC_s[sid] + s->costM + s->costLen_s[sid];
        s->costT += s->costC_s[sid] + s->costLen_s[sid];
    }
    // s->costT /= s->ntrip;
    s->costT += s->costM;
}

double _MDLtotal(Stack *Opt, Stack *C, int ntrip, int iter){
    int w = ntrip;
    int r = Opt->idx + C->idx;
    int m = mSegment(Opt) + mSegment(C);

    double cost = MDLSegment(Opt) + MDLSegment(C);
    double costT = cost + log_s(w) + log_s(r) + log_s(m)/w + m*log_2(r)/w + FB*r*r;
    fprintf(stdout, "%d %d %d %.0f\n", iter, r, m, costT); 

    return costT;
}

double _computeDiffVR(SegBox *s0, SegBox *s1){
    double cost01 = 0.0;
    double cost10 = 0.0;
    double cost11 = 0.0;
    double cost00 = 0.0;

    // for(int nsb_iter = 0; nsb_iter<s1->totalnSb; nsb_iter++){
    //     int st1 = s1->Sb[nsb_iter].st;
    //     int len1 = s1->Sb[nsb_iter].len;
    //     int trip1 = s1->Sb[nsb_iter].trip;
    //     cost10 += _viterbi(&s0->model, s0->delta, st1, len1, &ws->vit, trip1);
    //     cost11 += _viterbi(&s1->model, s1->delta, st1, len1, &ws->vit, trip1);
    // }

    // for(int nsb_iter = 0; nsb_iter<s1->totalnSb; nsb_iter++){
    //     int st0 = s0->Sb[nsb_iter].st;
    //     int len0 = s0->Sb[nsb_iter].len;
    //     int trip0 = s0->Sb[nsb_iter].trip;
    //     cost01 += _viterbi(&s1->model, s1->delta, st0, len0, &ws->vit, trip0);
    //     cost00 += _viterbi(&s0->model, s0->delta, st0, len0, &ws->vit, trip0);
    // }

    int id0 = FindMaxStEd(s0);
    int st0 = s0->Sb[id0].st;
    int len0 = s0->Sb[id0].len;
    int trip0 = s0->Sb[id0].trip;
    
    int id1 = FindMaxStEd(s1);
    int st1 = s1->Sb[id1].st;
    int len1 = s1->Sb[id1].len;
    int trip1 = s1->Sb[id1].trip;

    cost10 = _viterbi(&s0->model, s0->delta, st1, len1, &ws->vit, trip1);
    cost01 = _viterbi(&s1->model, s1->delta, st0, len0, &ws->vit, trip0);
    cost00 = _viterbi(&s0->model, s0->delta, st0, len0, &ws->vit, trip0);
    cost11 = _viterbi(&s1->model, s1->delta, st1, len1, &ws->vit, trip1);
    
    return cost10 + cost01 - cost11 - cost00;
}

void _printSplit(FILE *fp, SegBox *Sx, SegBox *s0, SegBox *s1){
    fprintf(fp, "[%s] %.0f(k:%d): {[%s] %.0f [%s] %.0f = %.0f(k:%d,%d)}\n",
        Sx->label, Sx->costT, Sx->model.k,
        s0->label, s0->costT,
        s1->label, s1->costT,
        s0->costT+s1->costT, s0->model.k, s1->model.k);

    /// segment
    fprintf(fp, "====\n");
    PrintStEd(fp, s0);
    fprintf(fp, "-------------\n");
    PrintStEd(fp, s1);
    fprintf(fp, "====\n");
}

void SaveMarker(MarkerWS *pws, char *fdir){
    ws = pws;
    FILE *fp;
    char filename[BUFSIZE];
    int total_r = ws->Opt.idx + ws->C.idx;

    mkdir(fdir, 0777);

    // print segment positions (segment.[d])
    for(int i=0; i<ws->Opt.idx; i++){
        sprintf(filename, "%ssegment.%d", fdir, i); 
        if(( fp = fopen (filename, "w") ) == NULL )
          error("can not open:",filename);
        PrintStEd(fp, ws->Opt.s[i]);
        fclose(fp);
    }
    for(int i=ws->Opt.idx; i<total_r; i++){
        sprintf(filename, "%ssegment.%d", fdir, i); 
        if(( fp = fopen (filename, "w") ) == NULL )
          error("can not open:",filename);
        PrintStEd(fp, ws->C.s[i-ws->Opt.idx]);
        fclose(fp);
    }

    // print regime labels (segment.labels)
    sprintf(filename, "%ssegment.labels",fdir); 
    if(( fp = fopen (filename, "w") ) == NULL )
        error("can not open:",filename);
    for(int i=0; i<ws->Opt.idx; i++){
        fprintf(fp, "%d\t\t%s\t\t%.0f\t\t%d \n",
        i, 
        ws->Opt.s[i]->label, 
        ws->Opt.s[i]->costT, 
        ws->Opt.s[i]->model.k);
    }
    for(int i=ws->Opt.idx; i<total_r; i++){
        fprintf(fp, "%d\t\t%s\t\t%.0f\t\t%d \n",
        i, 
        ws->C.s[i-ws->Opt.idx]->label, 
        ws->C.s[i-ws->Opt.idx]->costT, 
        ws->C.s[i-ws->Opt.idx]->model.k);
    }
    fclose(fp);
   
    #if(PRINTHMM)
        // print HMM params
        for(int i=0; i<ws->Opt.idx; i++){
            sprintf(filename, "%smodel.%d",fdir, i); 
            if(( fp = fopen (filename, "w") ) == NULL )
              error("can not open:",filename);
            PrintHMM(fp, &ws->Opt.s[i]->model);
            fclose(fp);
        }
        for(int i=ws->Opt.idx; i<total_r; i++){
            sprintf(filename, "%smodel.%d",fdir, i); 
            if(( fp = fopen (filename, "w") ) == NULL )
              error("can not open:",filename);
            PrintHMM(fp, &ws->C.s[i-ws->Opt.idx]->model);
            fclose(fp);
        }
    #endif
}

SegBox *_getS(char *parent, char *label){
    SegBox *s = Pop(&ws->S);
    if(s==NULL){
        error("too small ", "_getS()");
    }
    ResetStEd(s, ws->T.ntrip);
    ResetHMM(&s->model, MAXK, ws->d);
    sprintf(s->label, "%s%s", parent, label);
    s->delta = 1.0/(double)ws->lmax;
    s->opt = FALSE;
    return s;
}

void _releaseS(SegBox *s){ 
    Push(s, &ws->S);
}

void _allocSegBox(SegBox *s, int n){
    for(int i=0; i<ws->T.ntrip; i++){
        s->Sb = (SubS*)calloc((unsigned)n, sizeof(SubS));
        InitHMM(&s->model, MAXK, ws->d);
        s->nSb = (int *)calloc((unsigned)ws->T.ntrip, sizeof(int));
        s->trip  = (int*)calloc((unsigned)ws->T.ntrip, sizeof(int));
        s->thash  = (bool*)calloc((unsigned)ws->T.ntrip, sizeof(bool));
        s->len = (int *)calloc((unsigned)ws->T.ntrip, sizeof(int));
        s->costT_s = (double *)calloc((unsigned)ws->T.ntrip, sizeof(double));
        s->costC_s = (double *)calloc((unsigned)ws->T.ntrip, sizeof(double));
        s->costLen_s = (double *)calloc((unsigned)ws->T.ntrip, sizeof(double));
    }
}

/**
 * @fn
 * Allocate memories for our method
 * @param pws: work space
 */
void AllocMarker(MarkerWS *pws){
    fprintf(stderr, "memory allocation...\n");
    ws = pws;
    ws->maxss  = ws->T.ntrip + 120;
    ws->maxc   = (int) (log(ws->lmax) * ws->T.ntrip + 240);
    ws->maxseg = (int) ws->lmax*ws->T.ntrip*0.5;

    // Viterbi
    AllocViterbi(&ws->vit, ws->maxseg, ws->lmax, MAXK);
    AllocViterbi(&ws->vit2, ws->maxseg, ws->lmax, MAXK);
    ws->q = (int*) ivector(0, ws->lmax);
    
    // SegmentAssignment
    AllocCPS(&ws->cps, MAXK, ws->lmax);
    
    // BaumWelch
    int n = ws->maxseg;
    if(n > MAXBAUMN*ws->T.ntrip){ n = MAXBAUMN*ws->T.ntrip; }
    ws->x_tmp = (Input*)malloc(sizeof(Input)*n);
    AllocBaum(&ws->baum, n, ws->lmax, MAXK, ws->d);

    // Regimes
    ws->S.s = (SegBox **)calloc((unsigned)ws->maxc, sizeof(SegBox *));
    ws->S.idx = ws->maxc;
    for(int i=0;i<ws->maxc;i++){
        ws->S.s[i] = (SegBox *)calloc((unsigned)ws->maxc, sizeof(SegBox));
    }
    for(int i=0;i<ws->maxc;i++){
        ws->S.s[i]->Sb  = (SubS*)calloc((unsigned)ws->maxseg, sizeof(SubS));
        InitHMM(&ws->S.s[i]->model, MAXK, ws->d);
    }
    for(int i=0;i<ws->maxc;i++){
        ws->S.s[i]->trip  = (int*)calloc((unsigned)ws->T.ntrip, sizeof(int));
    }
    for(int i=0;i<ws->maxc;i++){
        ws->S.s[i]->thash  = (bool*)calloc((unsigned)ws->T.ntrip, sizeof(bool));
    }
    for(int i=0; i<ws->maxc; i++){
        ws->S.s[i]->nSb = (int *)calloc((unsigned)ws->T.ntrip, sizeof(int));
    }
    for(int i=0; i<ws->maxc; i++){
        ws->S.s[i]->len = (int *)calloc((unsigned)ws->T.ntrip, sizeof(int));
    }
    for(int i=0; i<ws->maxc; i++){
        ws->S.s[i]->costT_s = (double *)calloc((unsigned)ws->T.ntrip, sizeof(double));
    }
    for(int i=0; i<ws->maxc; i++){
        ws->S.s[i]->costC_s = (double *)calloc((unsigned)ws->T.ntrip, sizeof(double));
    }
    for(int i=0; i<ws->maxc; i++){
        ws->S.s[i]->costLen_s = (double *)calloc((unsigned)ws->T.ntrip, sizeof(double));
    }

    // candidate stack C and optimal stack Opt
    ws->C.s   = (SegBox**)calloc((unsigned)ws->maxc, sizeof(SegBox*));
    ws->Opt.s = (SegBox**)calloc((unsigned)ws->maxc, sizeof(SegBox*));
    ws->C.idx = 0;
    ws->Opt.idx = 0;

    /// for uniform sampling
    _allocSegBox(&ws->U, NSAMPLE*ws->maxss);
}


/// additional functions ///
int _findMinDiff_each(SegBox *s0, SegBox *s1, double *diffp, int tid){
    double min = INF;
    int loc = -1;
    
    // find trip [tid]
    int start = -1;
    for(int nsb_i=0; nsb_i<s0->nSb[tid]; nsb_i++){
        if(s0->Sb[nsb_i].trip == tid){
            start = nsb_i;
            break;
        }
    }
    if(start == -1) return -1;

    for(int i=start; i<s0->nSb[tid]; i++){
        int st  = s0->Sb[i].st;
        int len = s0->Sb[i].len;
        int trip = s0->Sb[i].trip;
        if(tid != trip) error("cannot find min diff\n", "_findMinDiff");
        double costC0 = _viterbi(&s0->model, s0->delta, st, len, &ws->vit, trip);
        double costC1 = _viterbi(&s1->model, s1->delta, st, len, &ws->vit, trip);
        double diff = costC1-costC0;
        if(min > diff) {loc = i; min = diff; }
    }
    *diffp = min;
    return loc;
}

double _scanMinDiff_each(SegBox *Sx, SegBox *s0, SegBox *s1, int tid){
    double diff;

    int loc0 = _findMinDiff_each(s0, s1, &diff, tid);
    int loc1 = _findMinDiff_each(s1, s0, &diff, tid);
    if(loc0==-1 || loc1==-1) return INF;

    SegBox *tmp0 = _getS("",""); SegBox *tmp1 = _getS("","");
    AddStEd(tmp0, s0->Sb[loc0].st, s0->Sb[loc0].len, s0->Sb[loc0].trip, ws->lmax);
    AddStEd(tmp1, s1->Sb[loc1].st, s1->Sb[loc1].len, s1->Sb[loc1].trip, ws->lmax);
    _estimateHMM_k(tmp0, MINK); _estimateHMM_k(tmp1, MINK);
   
    double costC = CPSearch(Sx, tmp0, tmp1, ws);
    _releaseS(tmp0); _releaseS(tmp1);
   
    return costC;
}