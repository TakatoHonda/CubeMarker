/* 
 * segbox.c --- written by Takato Honda
 * SegBox and Stack
 *
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

// output option
#define PRINTHMM 1
#define VITPATH 0
#define DBG 0
#define DBG_UNIEACH 0
#define DBG_RESET 0
#define DBG_UNISAMP 0
#define DBG_RMSTED 0


void PrintRegimes(char *fdir, Stack *Opt){
    FILE *fp;
    int i;
    char filename[BUFSIZE];
    /// print segment positions and trips
    for(i=0; i<Opt->idx; i++){
        sprintf(filename, "%s/segment.%d",fdir, i);
        if(( fp = fopen (filename, "w") ) == NULL ){
            error("can not open:",filename);
        }
        PrintStEd(fp, Opt->s[i]);
        fclose(fp);
    }

    /// print regime labels
    sprintf(filename, "%s/segment.labels",fdir);
    if(( fp = fopen (filename, "w") ) == NULL ){
        error("can not open:",filename);
    }
    for(i=0;i<Opt->idx;i++){
        fprintf(fp, "%d\t\t%s\t\t%.0f\t\t%d \n",
            i,
            Opt->s[i]->label,
            Opt->s[i]->costT,
            Opt->s[i]->model.k);
    }
    fclose(fp);

    #if(PRINTHMM)
    /// print HMM params
    for(i=0; i<Opt->idx; i++){
        sprintf(filename, "%s/model.%d",fdir, i);
        if(( fp = fopen (filename, "w") ) == NULL ){
            error("can not open:",filename);
        }
        PrintHMM(fp, &Opt->s[i]->model);
        fclose(fp);
    }//i
    #endif

    #if(VITPATH)
    /// print vit path
    _computeViterbiPath(&ws->Opt);
    #endif
}

void PrintStEd(FILE *fp, SegBox *s){
    int i;
    for(i=0;i<s->totalnSb;i++){
        fprintf(fp, "%d %d %d \n", 
            s->Sb[i].st, 
            s->Sb[i].st+s->Sb[i].len-1,
            s->Sb[i].trip);
    }
}

void _printResetStEd(SegBox *s, int w){
    // trips
    printf("s->ntrip: %d\n", s->ntrip);
    printf("[ ");
    for(int i=0; i<w; i++){
        int trip = s->trip[i];
        printf("(id:%d, n:%d, len:%d) ", trip, s->nSb[trip], s->len[trip]);
    }
    printf("]\n");
    
    // hash
    printf("hash[ ");
    for(int i=0; i<w; i++){
        printf("%d ", s->thash[i]);
    }
    printf("]\n");

    // segments
    printf("[ ");
    for(int i=0; i<s->totalnSb; i++){
        int st = s->Sb[i].st;
        int ed = st + s->Sb[i].len;
        int trip = s->Sb[i].trip;
        printf("(%d-%d, %d) ", st, ed, trip);
    }
    printf("]\n");

    //total
    printf("totalnSb: %d\n", s->totalnSb);
    printf("totallen: %d\n", s->totallen);
    printf("costC: %f\n", s->costC);
    printf("costT: %f\n", s->costT);
}

void ResetRegimeCost(SegBox *s, double cost){
    for(int i=0; i<s->ntrip; i++){
        s->costT_s[s->trip[i]] = cost; 
        s->costC_s[s->trip[i]] = cost; 
        s->costLen_s[s->trip[i]] = cost; 
    }
    s->costT = cost;
    s->costC = cost;
    s->costM = cost;
    s->costLen = cost;
}

void ResetStEd(SegBox *s, int w){ 
    #if(DBG_RESET)
    printf("=== Reset ===\n");
    _printResetStEd(s, w);
    printf("-----\n");
    #endif

    for(int i=0; i<w; i++){
        s->nSb[i] = 0;
        s->len[i] = 0;
        s->thash[i] = FALSE;
        s->costT_s[i] = INF;
        s->costC_s[i] = INF;
        s->costLen_s[i] = INF;
    }
    s->totalnSb = 0;
    s->totallen = 0;
    s->ntrip = 0;
    s->costC = INF;
    s->costT = INF;
    s->costM = INF;
    s->costLen = INF;

    #if(DBG_RESET)
    _printResetStEd(s, w);
    printf("=== Reset end\n");
    #endif
}

void _printRemoveStEd(SegBox *s, int id, int trip){
    printf("label: %s\n", s->label);
    printf("s->ntrip: %d\n", s->ntrip);
    printf("[ ");
    for(int i=0; i<s->ntrip; i++){
        printf("%d", s->trip[i]);
        printf("(hash:%d, nSb:%d) ", s->thash[s->trip[i]], s->nSb[s->trip[i]]);
    }
    printf("]\n");

    printf("[ ");
    for(int i=0; i<s->totalnSb; i++){
        int st = s->Sb[i].st;
        int ed = st + s->Sb[i].len;
        int trip = s->Sb[i].trip;
        printf("(%d-%d, %d)", st, ed, trip);
        if (i == id) printf("<-!! ");
        else printf(" ");
    }
    printf("]\n");
}

void RemoveStEd(SegBox *s, int id, int trip){
    #if(DBG_RMSTED)
        printf("========== RemoveStEd ==========\n");
        printf("RemoveStEd(s, id:%d, trip:%d)\n", id, trip);
        _printRemoveStEd(s, id, trip);
        printf("-----\n");
    #endif

    if(id >= s->totalnSb) {
        fprintf(stderr, "id: %d, s->totalnSb: %d\n", id, s->totalnSb);
        error("too large id","RemoveStEd");
    }

    int len = s->Sb[id].len;
    for(int i=id;i<s->totalnSb-1;i++){
        s->Sb[i].st=s->Sb[i+1].st;
        s->Sb[i].len=s->Sb[i+1].len;
        s->Sb[i].trip=s->Sb[i+1].trip;
    }//i
    s->totallen-=len;
    s->len[trip]-=len;
    s->totalnSb--;
    s->nSb[trip]--;

    // remove trip information
    if(s->nSb[trip] == 0){
        s->thash[trip] = FALSE;
        int loc = -1;
        for(loc=0; loc<s->ntrip; loc++){
            if(s->trip[loc] == trip) break;
        }
        for(int j=loc; j<s->ntrip-1; j++){
            s->trip[j] = s->trip[j+1];
        }
        s->ntrip--;
    }

    #if(DBG_RMSTED)
    _printRemoveStEd(s, id, trip);
    #endif
}

void RemoveVR(Stack *c, int id){
    if(id >= c->idx){
        error("too large id","RemoveRegime");
    }
    int i;
    for(i=id; i<c->idx-1; i++){
        c->s[i] = c->s[i+1];
    }
    c->idx--;
}

void _rm_overlap(SegBox *s){
    int i;
    int curr=INF;
    // check overlapped segments
    while(curr > s->totalnSb){
        curr = s->totalnSb;
        for(i=0;i<s->totalnSb-1;i++){
            int st0 = s->Sb[i].st;
            int ed0 = st0 + s->Sb[i].len;
            int trip0 = s->Sb[i].trip;
            int st1 = s->Sb[i+1].st;
            int ed1 = st1 + s->Sb[i+1].len;
            int trip1 = s->Sb[i+1].trip;
            int ed;
            if(ed0 > ed1) ed = ed0; else ed = ed1;
            // if overlapped
            if((ed0+1 >= st1) && (trip0 == trip1)){
                RemoveStEd(s, i+1, trip0);
                s->Sb[i].len = ed-st0;
            }//if
        }//i
    }//while
}

// generate subsequence in SegBox
void AddStEd(SegBox *s, int st, int len, int trip, int lmax){
    // if error ////////////////////////////////////
    if(len<=0){
        // fprintf(stderr, "### debug info ###\n");
        // fprintf(stderr, "st = %d\n", st);
        // fprintf(stderr, "len = %d\n", len);
        // fprintf(stderr, "trip = %d\n", trip);
        // fprintf(stderr, "lmax = %d\n", lmax);
        //error("len<=0", "addStEd"); 
        return;
    }
    if(st<0){  // st=0;
        error("st < 0", "addStEd"); //return;
    }
    if(st+len > lmax){
        // len = lmax - st;
        fprintf(stderr, "### debug info ###\n");
        fprintf(stderr, "st = %d\n", st);
        fprintf(stderr, "len = %d\n", len);
        fprintf(stderr, "trip = %d\n", trip);
        fprintf(stderr, "lmax = %d\n", lmax);
        error("st+len > lmax", "AddStEd");
    }
    // if error ////////////////////////////////////

    // find location to insert a segment ///////////
    int loc = -1; int i; int befseq = 0;
    for(i=0; i<s->totalnSb; i++){
        if( (s->Sb[i].st > st) && (s->Sb[i].trip == trip) ){
            break;
        // sequence changing point
        }else if( (befseq == trip) && (befseq != s->Sb[i].trip) ){
            break;
        }
        befseq = s->Sb[i].trip;
    }
    loc = i;
    // find location to insert a segment ///////////

    // add new segment
    for(int i=s->totalnSb-1;i>=loc;i--){
        s->Sb[i+1].st  = s->Sb[i].st;
        s->Sb[i+1].len = s->Sb[i].len;
        s->Sb[i+1].trip = s->Sb[i].trip;
    }
    s->Sb[loc].st = st;
    s->Sb[loc].len = len;
    s->Sb[loc].trip = trip;
    s->totalnSb++;
    s->nSb[trip]++;
    _rm_overlap(s);
    
    // length sum
    s->totallen = 0;
    for(int i=0; i<s->ntrip; i++){
        s->len[s->trip[i]] = 0;
    }
    for(int i=0; i<s->totalnSb; i++){
        s->len[s->Sb[i].trip] += s->Sb[i].len;
        s->totallen += s->Sb[i].len;
    }

    // if include the trip
    if(s->thash[trip] == TRUE) return;
    // if not include the trip
    // find location to insert a trip
    loc = -1; 
    for(i=0; i<s->ntrip; i++){
        if(s->trip[i] > trip){
            break;
        }
    }
    loc = i;
    // add new trip
    for(int i=s->ntrip-1;i>=loc;i--){
        s->trip[i+1] = s->trip[i];
    }
    s->trip[loc] = trip;
    s->ntrip++;
    s->thash[trip] = TRUE;  
}

/// allow overlapped segments
void AddStEd_ex(SegBox *s, int st, int len, int trip){
    s->totallen += len;
    s->len[trip] += len;
    int next = s->totalnSb;
    s->Sb[next].st = st; 
    s->Sb[next].len = len; 
    s->Sb[next].trip = trip;
    s->nSb[trip]++;
    s->totalnSb++;
    // including trips
    if(s->thash[trip] == FALSE){
        s->trip[s->ntrip] = trip;
        s->thash[trip] = TRUE;
        s->ntrip++;
    }
}

int Push(SegBox *s, Stack *C){
    C->s[C->idx] = s;
    C->idx++;
    return C->idx;
}

SegBox *Pop(Stack *C){
    if(C->idx==0) return NULL;
    C->idx--;
    return C->s[C->idx];
}

double MDLSegment(Stack *C){
    double cost = 0;
    for(int i=0; i<C->idx; i++){
        // cost += C->s[i]->costC;
        // cost += C->s[i]->costM;
        // cost += C->s[i]->costLen;
        cost += C->s[i]->costT;
    }
    return cost;
}

int mSegment(Stack *C){
    int i;
    int m = 0;
    for(i=0;i<C->idx;i++)
        m += C->s[i]->totalnSb;
    return m;
}

void CopyStEd(SegBox *from, SegBox *to){
    for(int i=0;i<from->totalnSb;i++){
        to->Sb[i].st = from->Sb[i].st;
        to->Sb[i].len = from->Sb[i].len;
        to->Sb[i].trip = from->Sb[i].trip;
    }
    for(int i=0; i<from->ntrip; i++){
        to->nSb[from->trip[i]] = from->nSb[from->trip[i]];
        to->trip[i] = from->trip[i];
        to->thash[from->trip[i]] = from->thash[from->trip[i]];
        to->len[from->trip[i]] = from->len[from->trip[i]];
        to->costT_s[from->trip[i]] = from->costT_s[from->trip[i]];
        to->costC_s[from->trip[i]] = from->costC_s[from->trip[i]];
        to->costLen_s[from->trip[i]] = from->costLen_s[from->trip[i]];
    }
    to->ntrip = from->ntrip;
    to->totalnSb   = from->totalnSb;
    to->totallen   = from->totallen;
    to->costT = from->costT;
    to->costC = from->costC;
    to->costLen = from->costLen;
    to->delta = from->delta;
    strcpy(to->label, from->label);
}

int IsSameStEd(SegBox *a, SegBox *b){
    int i;
    if(a->totalnSb != b->totalnSb)
        return FALSE;
    for(i=0;i<a->totalnSb;i++){
        if(a->Sb[i].st !=b->Sb[i].st){
            return FALSE;
        }
        if(a->Sb[i].len !=b->Sb[i].len){
            return FALSE;
        }
        if(a->Sb[i].trip !=b->Sb[i].trip){
            return FALSE;
        }
    }
    return TRUE;
}

int FindMinStEd(SegBox *s){
    int i;
    int id=-1; int min=INF;
    for(i=0;i<s->totalnSb;i++){
        if(s->Sb[i].len < min){ 
            id=i; min=s->Sb[i].len; 
        }
    }
    return id;
}

// find the longest segment
int FindMaxStEd(SegBox *s){
    int id=-1; int max=-INF;
    for(int i=0; i<s->totalnSb; i++){
        if(s->Sb[i].len > max){ 
            id=i; max=s->Sb[i].len; 
        }
    }
    return id;
}

void FixedSampling(SegBox *Sx, SegBox *s0, SegBox *s1, int len, int lmax, int w){
    int loc, r, trip;
    // init segments
    ResetStEd(s0, w);
    ResetStEd(s1, w);
    // segment s0
    loc = 0%Sx->totalnSb;
    r = Sx->Sb[loc].st;
    trip = Sx->Sb[loc].trip;
   if(r + len > lmax) len = lmax-r; // maybe OK
   AddStEd(s0, r, len, trip, lmax);
    // segment s1
   loc = 1%Sx->totalnSb; 
   r = Sx->Sb[loc].st+(int)(Sx->Sb[loc].len/2);
    if(r + len > lmax) len = lmax-r; // maybe OK
    AddStEd(s1, r, len, trip, lmax);
}

void RandomSampling(SegBox *Sx, SegBox *s0, SegBox *s1, int len, int lmax, int w){
    int loc, r, trip;
    // init segments
    ResetStEd(s0, w);
    ResetStEd(s1, w);
    // segment s0
    loc = (int) floor(getrand()*Sx->totalnSb);
    r = (int) floor(Sx->Sb[loc].st + getrand()*(Sx->Sb[loc].len-len));
    trip = Sx->Sb[loc].trip;
    AddStEd(s0, r, len, trip, lmax);
    // segment s1
    loc = (int) floor(getrand()*Sx->totalnSb);
    r = (int) floor(Sx->Sb[loc].st + getrand()*(Sx->Sb[loc].len-len));
    trip = Sx->Sb[loc].trip;
    AddStEd(s1, r, len, trip, lmax);
}

void PrintSegs(SegBox *s){
    printf("ntrip: %d\n", s->ntrip);
    printf("[ ");
    for(int i=0; i<s->ntrip; i++){
        printf("%d ",s->trip[i]);
    }
    printf("]\n\n");
    printf("totalnSb: %d\n", s->totalnSb);
    printf("[ ");
    for(int i=0; i<s->totalnSb; i++){
        int st = s->Sb[i].st;
        int ed = st + s->Sb[i].len;
        int trip = s->Sb[i].trip;
        printf("(%d-%d, %d) ", st, ed, trip);
    }
    printf("]\n");
    printf("===================================end\n");
}

void UniformSet(SegBox *Sx, int len, int trial, SegBox *U, int w){
    #if(DBG)
    printf("============ UniformSet ============\n");
    printf("Sx->len[%d]: %d\n", Sx->Sb[0].trip, Sx->len[Sx->Sb[0].trip]);
    #endif

    // create uniform blocks
    ResetStEd(U, w);
    int slideW = (int)ceil((Sx->totallen-len)/trial);
    for(int i=0; i<Sx->totalnSb; i++){
        if(U->totalnSb >= trial){
            #if(DBG)
                PrintSegs(U);
            #endif
            return;
        }
        int st = Sx->Sb[i].st;
        int ed = st+Sx->Sb[i].len;
        int trip=Sx->Sb[i].trip;
        // get trial (nsamples) samples
        for(int j=0; j<trial; j++){
            int next = st+j*slideW;
            if(next+len>ed){
                int st=ed-len; 
                if(st<0) st=0;
                AddStEd_ex(U, st, len, trip);
                break;
            }
            AddStEd_ex(U, next, len, trip);
        }//j
    }//i

    #if(DBG)
     PrintSegs(U);
    #endif
}

void _printUniformSampling(SegBox *s0, SegBox *s1){
    int st0 = s0->Sb[0].st;
    int st1 = s1->Sb[0].st;
    int ed0 = st0+s0->Sb[0].len;
    int ed1 = st1+s1->Sb[0].len;
    int trip0 = s0->Sb[0].trip;
    int trip1 = s1->Sb[0].trip;

    printf("s0: (%d-%d, %d)\n", st0, ed0, trip0);
    printf("s1: (%d-%d, %d)\n", st1, ed1, trip1);
}

void UniformSampling(SegBox *Sx, SegBox *s0, SegBox *s1, int len, int n1, int n2, SegBox *U, int lmax, int w){
    #if(DBG_UNISAMP)
    printf("========== UniformSampling ==========\n");
    #endif

    // init segments
    ResetStEd(s0, w);
    ResetStEd(s1, w);
    int i = (int) (n1 % U->totalnSb);
    int j = (int) (n2 % U->totalnSb);
    int st0 = U->Sb[i].st;
    int st1 = U->Sb[j].st;
    int trip0 = U->Sb[i].trip;
    int trip1 = U->Sb[j].trip;
    /// if overlapped, then, ignore
    if(abs(st0-st1)<len && trip0 == trip1) {
        #if(DBG_UNISAMP)
        printf("========== UniformSampling ========end\n");
        #endif
        return;
    }
    AddStEd(s0, st0, len, trip0, lmax);
    AddStEd(s1, st1, len, trip1, lmax);
    
    #if(DBG_UNISAMP)
    _printUniformSampling(s0, s1);
    printf("========== UniformSampling ========end\n");
    #endif
}

void GetSingleTrip(int tripid, SegBox *Sx, SegBox *singleTripSegBox, int lmax){
    // detect start position
    int start = 0;
    for(int i=0; i<Sx->totalnSb; i++) if(Sx->Sb[i].trip == tripid){start = i;break;}

    for(int i=start; i<start+Sx->nSb[tripid]; i++){
        int st = Sx->Sb[i].st;
        int len = Sx->Sb[i].len;
        int trip = Sx->Sb[i].trip;
        AddStEd(singleTripSegBox, st, len, trip, lmax);
    }
}