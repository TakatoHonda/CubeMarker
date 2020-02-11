CFLAGS= -g -Wall -std=c99
INCS=
CC=gcc
CODES = cubemarker
all :  $(CODES)
#------------------------------------------------------#
cubemarker: cubemarker.o marker.o cps.o tool.o nrutil.o hmm.h cubemarker.h \
	viterbi.o hmmutils.o \
	baum.o forbackward.o kmeans.o segbox.o
	$(CC) -o cubemarker \
	cubemarker.o marker.o cps.o tool.o nrutil.o \
	viterbi.o hmmutils.o \
	baum.o forbackward.o kmeans.o segbox.o \
	-lm 
#------------------------------------------------------#

#make plot: plots_all.m
#	matlab -r 'plots_all' 
#make demo: demo.sh
#	sh demo.sh

clean :
	rm -f core $(CODES) *.o *~ 
