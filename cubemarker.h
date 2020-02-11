/* 
 * cubemarker.h --- written by Takato Honda 
 */

#define DEFAULT 1 
#if(DEFAULT)
    #define MAXK     8
    #define LSET     NO
    #define MAXBAUMN 3
    #define NSAMPLE  10
#else
    #define MAXK     4
    #define LSET NO  //YES
    #define MAXBAUMN 1
    #define NSAMPLE  5
#endif

typedef int bool;
#define TRUE  1
#define FALSE 0

typedef struct {
    int st;
    int len;
    int trip; // trip number  
}SubS; // SubSequence

typedef struct {
    SubS *Sb; // Sub sequence list
    int *nSb; // the number of SubS for each trip
    int totalnSb; // total number of SubS
    int *len; // length for each trip
    int totallen; // total length of regime
    int *trip; // including trips
    bool *thash; // hash table of trip
    int ntrip; // number of trips
    bool opt;

    double *costT_s;   // total cost
    double *costC_s;   // coding cost
    double *costLen_s; // length cost

    double costT;   // total cost
    double costC;   // coding cost
    double costM;   // model cost
    double costLen; // length cost

    int optimal;
    char label[BUFSIZE];
    HMM model;
    double delta;
}SegBox; // regime

typedef struct {
    SegBox **s; // regimes set
    int idx;    // number of regimes
}Stack; // regimes stack

typedef struct {
    Input *x;  // geographical complex sequences (tensor)
    int ntrip; // number of trips
}TripTensor;

typedef struct {
    int S[10000]; // cut points
    int regime[10000]; // regime numbers for each cut point
}CutPoint;

typedef struct {
    double *Pu, *Pv, *Pi, *Pj; // probabilities
    int **Su, **Sv, **Si, **Sj; // cut points
    int *nSu, *nSv, *nSi, *nSj; // number of cut points
}CPS;  //CutPointSearch

typedef struct {
    TripTensor T;
    char *fout;
    char *fin;
    int maxc; 
    int maxseg;
    int maxss; // maximum number of H-regimes
    int d;     // dimension
    int lmax; // maximum duration
    int fmax; // maximum number of sequences

    CPS cps; // CutPointSearch
    BAUM baum; // baum & viterbi
    int *q; //q[m] ... Viterbi path
    VITERBI vit;
    VITERBI vit2;
    Input *x_tmp;
    SegBox *s;

    // Stack tmp;
    Stack C;
    Stack Opt;
    Stack S;

    double costT;
    SegBox U; // uniform sampling (segment)
}MarkerWS; // cubemarker workspace

/// marker.c
void PrintParams();
void AllocMarker();
void Marker();
void SaveMarker();
void MergeSegment();

/// segbox.c
void PrintStEd(FILE *fp, SegBox *s);
void PrintRegimes(char *, Stack *);
void CopyStEd(SegBox *from, SegBox *to);
void ResetStEd(SegBox *s, int w);
void ResetRegimeCost(SegBox *s, double cost);
void AddStEd(SegBox *s, int st, int len, int trip, int lmax);
void AddStEd_ex(SegBox *s, int st, int len, int trip);
void RemoveStEd(SegBox *s, int id, int trip);
void RemoveVR(Stack *, int);
int Push(SegBox *s, Stack *C);
SegBox *Pop(Stack *C);
int mSegment(Stack *C);
double MDLSegment();
int IsSameStEd(SegBox *a, SegBox *b);
int FindMaxStEd();
int FindMinStEd();
void FixedSampling();
void RandomSampling();
void UniformSet(SegBox *Sx, int len, int trial, SegBox *U, int w);
void UniformSampling(SegBox *Sx, SegBox *s0, SegBox *s1, int len, int n1, int n2, SegBox *U, int lmax, int w);
void PrintSegs();
void GetSingleTrip();

/// cps.c
double CPSearch(SegBox *Sx, SegBox *s0, SegBox *s1, MarkerWS *wsd);
void AllocCPS();
void FreeCPS();

/// test.c 
void TestAllocSA();
