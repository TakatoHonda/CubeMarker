/*
 * tool.c --- written by Takato Honda
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "nrutil.h"
#include "hmm.h"
#include "tool.h"
#include "cubemarker.h"

#define DBG 0

SegBox *_getS();
void _estimateHMM();
void _estimateHMM_k();
void notfin();
void error();

double EuclideanDistance(Input *x, int trip0, int trip1){
    double distance = 0.0;
    int m = x->m; // length of a sequence 
    int d = x->d; // dimension
    // sequences
    float **T0 = x[trip0].O;
    float **T1 = x[trip1].O;
    for(int i=0; i<m; i++){
        for(int j=0; j<d; j++){
            distance += sqrt(pow(T0[i][j] - T1[i][j], 2.0));
        }
    }

    return distance;
}

//----------------------------------------//
// random
//---------------------------------------//
void setseed(void) {
    int seed = 1;
    // int seed = (int) getpid();
    srand((unsigned int)seed);
}

double getrand(void){
    return (double) rand()/RAND_MAX;
}

//----------------------------------------//
// MDL 
//---------------------------------------//
double log_2(int x){
    return log((double)x)/log(2.0);
}

double log_s(int x){
    return 2.0*log_2(x)+1.0;
}

double costHMM(int k, int d){
    double fB   = (double) FB; 
    double cost = (double) fB * (k + k*k + 2*k*d) + log_s(k);
    return cost;
}

void Split(Input *X, int st, int len, Input *x){
    if(st < 0) st = 0;
    if(X->m < st + len){
        len = X->m - st;
    }
    // pointer shift
    x->O = X->O + st;
    x->m = len;
    x->parent = X->id;
    x->st = st;
    sprintf(x->tag, "[%d] %s [%d-%d](%d)", x->id, X->tag, st, st+len, len);
}

//----------------------------------------//
// Input Output
//---------------------------------------//
void PrintIdx(FILE *fp, int *idx, int n){
    int i;
    for(i=0;i<n;i++){
        fprintf(fp, "[%d] \t : \t %d\n", i, idx[i]);
    }
}

void PrintInput(FILE *fp, Input *x){
    fprintf(fp, "id %d :  pid %d ", x->id, x->pid);
    fprintf(fp, "[%s] (%d)\n", x->tag, x->m);
}

void ReadSequence(FILE *fp, int m, int *d, float **O){
    int i,j;
    char c;
    for (i=0; i < m; i++){
        for (j=0; j < *d; j++){
            fscanf(fp,"%f", &O[i][j]);
            #if(DBG)
            fprintf(stderr, "%f", O[i][j]);
            #endif
            c = getc(fp);
            if(c == '\n'){
                *d = j+1;
                break;
            }
        }
        #if(DBG)
        fprintf(stderr, "\n");
        #endif
    }
}

int LoadSequence(char *fn, int *mx, int *d,float ***pO){

    // load length of the sequence
    int len = count_lines(fn);
    fprintf(stderr, "file: %s (len=%d)\n", fn, len);
    int m;
    if(len < *mx) m = len;
    else          m = *mx;

    *pO = (float**)matrix(0, m, 0, *d);
    FILE *fp = fopen(fn, "r");
    if (fp == NULL) error("cannot open", fn);
    
    ReadSequence(fp, m, d, *pO);
    fclose(fp);
    return m;
}

void LoadDB(char *fn, int n, int *mx, int *d, TripTensor *T){
    int i;
    // open filelist
    FILE *fp = fopen(fn, "r");
    char buf[BUFSIZE];
    if (fp == NULL) error("cannot open list", fn);

    // alloc memory for sequences
    Input *xlst = (Input*)malloc(sizeof(Input)*n);

    for(i=0;i<n;i++){
        // for scalability
        if(fscanf(fp, "%s", buf)==EOF){
            //error("error (filesize < n)", fn);
            break;
        }
        xlst[i].m = LoadSequence(buf, mx, d, &xlst[i].O);
        xlst[i].d = *d;
        strcpy(xlst[i].tag,buf);
        xlst[i].id=i;
    }//i
    fclose(fp);

    T->x = xlst;
    T->ntrip = i;
}

void PrintSequence(FILE *fp, int m, int d, float **O){
    int i,j;
    for(i=0;i<m;i++){
        for(j=0;j<d;j++){
            fprintf(fp,"%f ", O[i][j]);
        }
        fprintf(fp, "\n");
    }//i
}

void NormSequence(Input *xlst, int n, int dim){
    int i,j,d;
    fprintf(stderr, "Normalize sequences (0:%f) ... \n", MAXVAL);
    for(d=0;d<dim;d++){
        float min = INF, max=-INF;
        for(i=0;i<n;i++){
            Input *x = &xlst[i];
            /// find min & max values
            for(j=0;j<x->m;j++)
                if(x->O[j][d] > max)
                    max = x->O[j][d];
                else if(x->O[j][d] < min)
                    min = x->O[j][d];
        }//i
        for(i=0;i<n;i++){
            Input *x = &xlst[i];
            /// normalize -> [0:MAXVAL]
            for(j=0;j<x->m;j++){
                x->O[j][d] = MAXVAL*(x->O[j][d]-min) / (max-min);
            }
        }//i
    }//d
    //PrintSequence(stdout, x->m, dim, x->O);
}

void ZnormSequence(TripTensor *T, int dim){
    int i,j,d;
    fprintf(stderr, "Z-normalization... \n");
    for(d=0;d<dim;d++){

        float mean=0.0, std=0.0;
        int cnt=0;
        /// compute mean
        for(i=0;i< T->ntrip;i++){
            Input *x = &(T->x[i]);
            cnt+=x->m;
            for(j=0;j<x->m;j++){
                mean += x->O[j][d];
            }
        }//i
        mean /= cnt;
        
        /// compute std
        for(i=0;i<T->ntrip;i++){
            Input *x = &(T->x[i]);
            for(j=0;j<x->m;j++){
                std += pow(x->O[j][d]-mean,2.0);
            }
        }//i
        std = sqrt(std/cnt);

        //printf("std = %f\n",std);
        for(i=0;i<T->ntrip;i++){
            Input *x = &(T->x[i]);
            /// normalize 
            for(j=0;j<x->m;j++){
                x->O[j][d] = MAXVAL*(x->O[j][d]-mean) / std;
                //printf("%f\n", x->O[j][d]);
            }
    }//i

    }//
    //PrintSequence(stdout, T->x[0].m, dim, T->x[0].O);
}


void SaveModel(HMM *hmm, int g, char *fn){
    FILE *fp;
    int i;
    char filename[BUFSIZE];

    for(i=0;i<g;i++) {
        sprintf(filename, "%s.model.0.%d",fn, i);
        if( ( fp = fopen (filename, "w") ) == NULL ){
            error("can not open:",fn);
        }
        PrintHMM(fp, &(hmm[i]));
        fclose(fp);
    }// i
}

int Max_dvector(double *value, int len){
    int i;
    int max = -INF;
    for(i=0;i<len;i++){
        if(max < value[i]){
            max = value[i];
        }
    }
    return max;
}

int count_lines(char *filename){
    int j;
    char  read_buf[BUFSIZE];
    FILE  *fp;
    if ( ( fp = fopen ( filename , "r" ) ) == NULL ){
        fprintf(stderr, "  %s can not open\n" , filename);
        exit(1);
    }
    j = 0;
    while (fgets(read_buf,BUFSIZE,fp) != NULL ){
        j++;
    }
    fclose(fp);
    return j;
}

void CopyFile(char *fin, char *fout){
    FILE *fp_in;
    FILE *fp_out;
    char str[BUFSIZE];
    
    // make file 
    fp_out = fopen(fout, "w");

    // read copied file
    fp_in = fopen(fin, "r"); // .list file
    char *p = fgets(str, BUFSIZE, fp_in);
    while( p != NULL){
        fprintf(fp_out, "%s", str);
        p = fgets(str,BUFSIZE, fp_in);
    }

    fclose(fp_in);
    fclose(fp_out);    
}

void IntToString(char *str, int num){
    sprintf(str, "%d", num);
}

//----------------------------------------//
// system 
//---------------------------------------//
void error(char *msg, char *name){
    fprintf(stderr, "error: %s [%s] \n", msg, name);
    exit(1);
}

// getrusage_sec()
double getrusage_sec(){
    struct rusage t;
    struct timeval tv;
    getrusage(RUSAGE_SELF, &t);
    tv = t.ru_utime;
    return tv.tv_sec + (double)tv.tv_usec*1e-6;
}
void notfin(char *msg){
    fprintf(stderr, "//-----------------------------------//\n");
    fprintf(stderr, "!!! warning !!!  (not finished!: %s)\n", msg);
    fprintf(stderr, "//-----------------------------------//\n");
}
