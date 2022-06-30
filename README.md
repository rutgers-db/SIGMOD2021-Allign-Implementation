# SIGMOD2021-Allign-Implementation
Implementation for 2021 SIGMOD Paper Allign: Aligning All-Pair Near-Duplicate Passages in Long Texts

# Content
allign.cc: main function for allign

segtree.h: implementation for segment tree 

utils.hpp: some implemented interfaces for allign

stopwords.txt: some stop words should be removed from the source and suspicious document

source-document01256.txt: one example source document

suspicious-document00001.txt: one example suspicious document

makefile: file for code compilation

# How to Run
First, generate the executable file:

```make allign```

One executable example here: 

```./allign -docFileName source-document01256.txt -queryFileName suspicious-document00001.txt -theta 0.5 -tau 50 -k 100```

For the example here, it takes around ten seconds to generate the ouptut

Parameters: 
```
-docFileName: path to source document file

-queryFileName: path to suspicious document file

-theta: given threshold for jaccard similarity

-tau: minimum range for a compact window

-k: num of hash functions