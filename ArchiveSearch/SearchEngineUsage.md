#Usage of Library Retrieval Search Tool
## How to start a task?
### command line
```$xslt
wulong@KEZ323 ~/data/libsearch $ ./LibraryRetrieval --config ~/bitbucket/codejam/tools/Release/libsearch.conf  --rawdata 13sep2013_yeast_longRuns_1.mzXML --libraryfile yeast_2012_04_17_qtof.sptxt --overwrite
```

### The confige file
```$xslt
~/bitbucket/codejam/tools/Release/libsearch.conf
```

## Issues
### Why it is very slow to load queries?
```$xslt
[2019-01-02 20:38:58.716] [A] [info] Loading query data file...
[13sep2013_yeast_longRuns_1.mzXML] Progress:..........10%..........20%..........30%..........40%..........50%..........60%..........70%..........80%..........90%..........100%
[Info] Empty spectra: 0

Summary: #Spec: 196299  MS1:    40039   MS2:    156260

[2019-01-02 20:39:28.221] [A] [info] Converting query data as mz file
[2019-01-02 20:39:28.221] [A] [info] Creating mz file ...
Processing ---  File 0: 13sep2013_yeast_longRuns_1.mzXML
[13sep2013_yeast_longRuns_1.mzXML] Progress:..........10%..........20%..........30%..........40%..........50%..........60%..........70%..........80%..........90%..........100%
[Info] Empty spectra: 0

Summary: #Spec: 196299  MS1:    40039   MS2:    156260

[Generating mz format] Progress:..........10%..........20%..........30%..........40%..........50%..........60%..........70%..........80%..........90%..........100%
Start to writing the file to disk
Loading compact file
filebytes: 15626000 peaknum 7813000 specnum 156260
[2019-01-02 20:40:08.166] [A] [info] mz file is created!

```

It takes 1 mins. We should make it faster!

### The EXIT without --psm option specified is Wrong

### Recall is Wrong
```$xslt
Total number of PSM above threshold: 61029
// Start recall (Type-K = 1001 means the true peptide not found in top K)
Type-K  Counts@K        Recall@K
1       10731   0.1758
2       128     0.0021
3       7       0.0001
4       3       0.0000
6       3       0.0000
7       4       0.0001
10      1       0.0000
1001    49804   0.8161
1002    348     0.0057
// End of Recall
```


### The CPU and GPU version return different neighbours
Now we know that the GPU and CPU version are different. This may not because of bug. But due to the different [implementation in floating point value calcualtion](https://github.com/facebookresearch/faiss/issues/178).


```
Total number of PSM above threshold: 42403
// Start recall (Type-K = 1001 means the true peptide not found in top K)
Type-K  Counts@K        Recall@K        ACCRecall@k
1       40532   0.9559  0.9559
2       730     0.0172  0.9731
3       156     0.0037  0.9768
4       81      0.0019  0.9787
5       40      0.0009  0.9796
6       29      0.0007  0.9803
7       16      0.0004  0.9807
8       21      0.0005  0.9812
9       10      0.0002  0.9814
10      14      0.0003  0.9817
11      5       0.0001  0.9819
12      13      0.0003  0.9822
13      8       0.0002  0.9824
14      8       0.0002  0.9825
15      4       0.0001  0.9826
16      7       0.0002  0.9828
17      6       0.0001  0.9829
18      5       0.0001  0.9831
19      5       0.0001  0.9832
20      2       0.0000  0.9832
21      2       0.0000  0.9833
22      2       0.0000  0.9833
23      8       0.0002  0.9835
24      2       0.0000  0.9836
25      1       0.0000  0.9836
26      4       0.0001  0.9837
28      1       0.0000  0.9837
29      2       0.0000  0.9838
30      1       0.0000  0.9838
31      2       0.0000  0.9838
32      3       0.0001  0.9839
33      1       0.0000  0.9839
35      1       0.0000  0.9839
39      3       0.0001  0.9840
41      1       0.0000  0.9840
42      1       0.0000  0.9841
43      1       0.0000  0.9841
45      1       0.0000  0.9841
49      1       0.0000  0.9841
1001    665     0.0157  0.9998  Correct libraries spectra NOT fount in top 50
1002    8       0.0002  1.0000  Matched to decoy protein

```
The result above come from the following command:
```bash
./LibraryRetrieval --config ~/bitbucket/codejam/tools/Release/libsearch.conf  --rawdata 13sep2013_yeast_longRuns_1.mzXML --libraryfile NIST_yeast_IT_2012_TD.sptxt  --psm interact-spectrast.pep.xml   --overwrite  --annrecall  --indexpath ivf256pq16_4indices --outputpath outputgpu --indexstrings "IVF256,PQ16+16;IVF256,PQ16+16;IVF256,PQ16+16" -g
```

So the top 5 recall is already quite good. Maybe we could try **Interaction**, instead of **Union**.


## Notes of reading binary mz file format

### How to read mz file or any binary file?
```bash
od -x 251.mzXML.mz.backup -N 100 -t u2
```

```
0000000  4aa2  5231  4cd2  38b7  4545  285e  1cb2  5008
       19106 21041 19666 14519 17733 10334  7346 20488
0000020  3a85  497c  2c6c  38ea  439a  4016  517c  22fd
       14981 18812 11372 14570 17306 16406 20860  8957
0000040  4726  4acc  3170  4d14  4e78  4816  3893  41f1
       18214 19148 12656 19732 20088 18454 14483 16881
0000060  4fbf  35a5  4854  0000  0000  0000  0000  0000
       20415 13733 18516     0     0     0     0     0
0000100  0000  0000  0000  0000  0000  0000  0000  0000
           0     0     0     0     0     0     0     0
*
0000140  0000  0000
           0     0

```


## Note for CUDA Compiler
CUDA does not support gcc higher than 5. Using gcc 6 or 7 is not possible so far. 