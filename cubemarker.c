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
#define LMAX 1000000 // Maximum length of a bundle
#define FMAX 10000   // Maximum number of files
#define DMAX 8     // Maximum number of dimensions

int main(int argc,char *argv[]){
    // computation time
    double start, end, totaltime;

    setseed();
    MarkerWS ws;
    ws.lmax = LMAX; 
    ws.d = DMAX;
    ws.fmax = FMAX;

    // check arguments
    if(argc != 3 + LSET) error("usage: cmd fin, fout", argv[0]);

    // load argument
    ws.fin = argv[1];
    ws.fout = argv[2];
    if(LSET){ ws.lmax = atoi(argv[3]); }

    printf("--------------------------\n");
    fprintf(stderr, "loading data...\n");
    printf("--------------------------\n");

    // load tensor data
    LoadDB(ws.fin, FMAX, &ws.lmax, &ws.d, &ws.T);
    ws.lmax = ws.T.x->m;
    ws.d = ws.T.x->d;
    fprintf(stderr, "(d, w, n) = (%d, %d, %d)\n", ws.d, ws.T.ntrip, ws.lmax);
    printf("--------------------------\n");

    // z-normalize
    ZnormSequence(&ws.T, ws.d);

    // allocate memory
    AllocMarker(&ws);
    fprintf(stderr, "start CubeMarker...\n");
    
    // computation time 
    start = getrusage_sec();
    
    // start method
    Marker(&ws);
    
    // computation time
    end = getrusage_sec();
    totaltime = end-start;
    
    char dir[100];
    strcat(strcpy(dir, ws.fout), "final/");
    SaveMarker(&ws, dir);
    CopyFile(ws.fin, strcat(dir, "input")); 

    printf("==================================\n");
    printf("Time: %.4f sec.\n",totaltime);
    printf("Regimes: %d \n", (ws.Opt.idx));
    printf("Cost: %.0f \n", ws.costT);
    printf("==================================\n");
}