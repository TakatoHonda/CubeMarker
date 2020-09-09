###### Version: March 27, 2018
###### Author: Takato Honda (takato@sanken.osaka-u.ac.jp)
---
# Introduction
This is a code of CubeMarker, which finds multi-aspect patterns and groups of patterns, written in C. CubeMarker is accepted as a regular paper for ICDM'19, entitled ["Multi-Aspect Mining of Complex Sensor Sequences"](https://takatohonda.github.io/paper/paper-icdm19.pdf).
![overview](http://takatohonda.github.io/assets/img/overview.png)

# Quick demo
(a) CubeMarker for a demo data
```
sh demo.sh 
```
(b) Visualization (Matlab code)
```
matlab -r 'plots_all' 
```
# Clean up all files
```
make clean
```
# CubeMarker outputs
\----------  
i|r|m|Cost  
\----------  
1 1 33 1009418  
2 2 412 1199212  
3 3 583 1161483  
4 3 583 1161483  
...

* i: # of iterations
* r: # of regimes (groups of patterns)
* m: # of segments (patterns)
* Cost: Total cost of the result

# Files and directories
* _dat/: demo data
* _list/: demo data list
* _demo/: demo result figs
* _out/: output directory
    * _out/[#]/: result for each iteration
        - _out/[#]/input: input files list
        - _out/[#]/model.[#]: HMM parameters
        - _out/[#]/segment.[#]: segments info
        - _out/[#]/segment.labels: regime info
