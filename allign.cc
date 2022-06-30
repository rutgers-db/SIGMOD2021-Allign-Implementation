#include <fstream>
#include <cmath>
#include <iostream>
#include <sys/time.h>
#include <vector>
#include <unordered_map>
#include <cassert>
#include <unordered_set>
#include <string>
#include <ctype.h>
#include <algorithm>
#include <set>
#include "omp.h"
#include "segtree.h"
#include <dirent.h>
#include <chrono>
#include "utils.hpp"
// #include "bbt.h"

// #define DEBUG 0
// #define RMQ 1

using namespace std;

typedef unsigned int uint;
typedef unsigned char byte;
int MAXLENGTH = 10000;//0x7fffffff; 


// placeholder range boundary and interval id
#define INVALID_RANGE -1  
#define INVALID_ID -1

enum mode {
  SMALL,
  LARGE
};

// TODO: add more code
int binaryFind(int pos, const vector<int> &ppos, mode m) {
  int l = 0;
  int r = ppos.size() - 1;
  int res = 0;
  while (l <= r) {
    int mid = (l + r) >> 1;
    if (pos >= ppos[mid]) {
      res = mid;
      l = mid + 1;
    } else {
      r = mid - 1;
    }
  }
  return m == SMALL? res: res + 1;
}


// k hash functions;
vector<pair<int, int>> ab;

#ifdef RMQ
vector<RMQSegTree *> querySegTrees;
vector<RMQSegTree *> docSegTrees;
vector<RMQSegTree *> *ptr;
#endif

// The hash value function
inline uint hval(int word, int kth_hash, int multiplicity,  const vector<vector<int>> &mulWordIds)
{
  assert(multiplicity >= 1);
  if (multiplicity == 1)
    return ab[kth_hash].first * word + ab[kth_hash].second;
  else
    return ab[kth_hash].first * mulWordIds[word][multiplicity - 2] + ab[kth_hash].second;
}


// Make string to tokens based on delimiters
// Output parameter: res
void strToTokens(string &str, const string &delimiter, vector<string> &res, vector<pair<int, int>> &offsets) {
    #ifdef NORMAL
    // Replace illegal chars
    for (int i = 0; i < str.length(); ++i) {
      str[i] = str[i] <= 0 || str[i] == '\n'? ' ': str[i];
    }
    #endif
    char *inputStr = strdup(str.c_str());
    int wordNum = 0;
    char *key = strtok(inputStr, delimiter.c_str());

    // Iterate each word
    while(key){
      // Calculate start position
      int startPos = key - inputStr;
      int endPos = startPos + strlen(key);
      wordNum ++;
      for(int p = 0; p < strlen(key); p ++)
        key[p] = tolower(key[p]);
      res.push_back(key);
      offsets.push_back({startPos, endPos});
      key = strtok(0, delimiter.c_str());
    }

    delete []inputStr;
}

void getFiles(string path, vector<string> &files)
{
    DIR *dr;
    struct dirent *en;
    string file_path;
    dr = opendir(path.c_str()); //open all directory
    if (dr)
    {
        while ((en = readdir(dr)) != NULL)
        {
            //ignore hidden files and folders
            if (en->d_name[0] != '.')
            {
                file_path = path + en->d_name;
                files.push_back(file_path);
            }
        }
        closedir(dr); //close all directory
    }
    cout << endl;
}

// Read stopwords from given file name
void readStopWords(const string &fileName, unordered_set<string> &stopWords)
{
    // one per line
    string stopWord;
    ifstream file(fileName, ios::in);
    while (getline(file, stopWord))
        stopWords.insert(stopWord);
}

// Read from file and do preprocessing
//
// Input parameter:
//  1) filename: the file name
//  2) stopwords: stopwords that don't need to be considered
//
// Output parameter:
//  1) doc: contains the word id (from word string => word id)
//  2) doc_offsets: store word offset
//  3) word2id: the mapping from word to id
//  4) id2word: the mapping from id to word
void word2int(const string &filename, vector<int> &doc, vector<pair<int, int>> &doc_offsets, unordered_map<string, int> &word2id, vector<string> &id2word, vector<int> &id2maxFreq, vector<vector<int>> &id2mulId, const unordered_set<string> &stopwords)
{
    static int mulWordId = -1;
#ifdef NORMAL
    const string delim = "\t\n\r\x0b\x0c !\"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~";
    const string period_str = ",.!?;";
//cout << "execute NORMAL here line 102" <<endl;
#endif

#ifndef NORMAL
    const string delim = " ";
    const string period_str = "\n";
//cout << "execute NORMAL here line 108" <<endl;
#endif

    // Read from file ...
    ifstream datafile(filename, ios::binary);
    datafile.seekg(0, std::ios_base::end);
    int length = datafile.tellg();
    string docstr(length + 1, '\0');
    datafile.seekg(0);
    datafile.read(&docstr[0], length);

    // Make the docstr to sentences ..
    vector<string> tokens;
    vector<pair<int, int>> tmp;
    strToTokens(docstr, delim, tokens, doc_offsets);

    unordered_map<int, int> id2freq;
    int wordcnt = 0;
    for (int i = 0; i < tokens.size(); i++)
    {

        // Skip stop words
        if (stopwords.find(tokens[i]) != stopwords.end())
            continue;
        
        if (wordcnt++ >= MAXLENGTH)
          break;

        // If a new word, add to word2id, offsets, doc
        if (word2id.find(tokens[i]) == word2id.end())
        {
            word2id[tokens[i]] = id2word.size();
            id2word.emplace_back(tokens[i]);
        }

        int wid = word2id[tokens[i]];
        doc.emplace_back(wid);
        id2freq[wid] += 1;
    }

    // Now calculate multi-id for each word that appears multiple times .. 
    id2maxFreq.resize(id2word.size(), 1);
    id2mulId.resize(id2word.size(), vector<int>());

    for (auto &entry : id2freq)
    {
        int id = entry.first;
        int freq = entry.second;
        if (id2maxFreq[id] >= freq)
            continue;

        for (int i = id2maxFreq[id] + 1; i <= freq; i++)
        {
            id2mulId[id].emplace_back(mulWordId--);
        }
        id2maxFreq[id] = freq;
    }
}


// Find the min hash from the window [l, r]
//
// Input parameters:
//  1) doc: word id vectors
//  2) hashfunc: the hash functions selected
//  3) l: left side
//  4) r: right side
//  5) mulWordIds: multiwordIds
// 
// Output parameters:
//  1) mul: the number of appear times
//  2) minHash: the minHash value
//  3) pos: the position of each appears times
void scanMinHash(const vector<int> &doc, const int hashfunc, int l, int r, const vector<vector<int>> &mulWordIds, int &mul, uint &minHash, vector<int> &pos)
{
  assert(l <= r);

  #ifdef RMQ
  if (ptr) {
    (*ptr)[hashfunc]->scanMinHash(l, r, mul, minHash, pos);
    return;
  }
  #endif

  unordered_map<int, int> freq;

  // Record the word ID with min hash
  int minId; 

  // Iterate each word in the window
  for (int x = l; x <= r; x++) {
    // Increase the frequency
    freq[doc[x]] += 1;
    // Calculate the hash value
    uint hv = hval(doc[x], hashfunc, freq[doc[x]], mulWordIds);
    // Record the minHash value
    if (x == l || hv <= minHash)
    {
      minHash = hv;
      minId = doc[x];
      mul = freq[doc[x]];
    }
  }

  // Record each appear position ..
  for (int x = l; x <= r; x++) {
    if (doc[x] == minId) {
      pos.push_back(x);
    }
  }
}

// Do multi-divide-conquer to generate compact windows
// Algorithm 10
void multiDivideConquer(const vector<int> &doc, int l, int le, int r, vector<CW> &cws, const int hashfunc, const int len_threshold, const vector<vector<int>> &mulWordIds, const vector<int> &ppos)
{
  if (r - l + 1 < len_threshold)
    return;

  int x; // multiplicity
  uint minhash;
  vector<int> minpos;
  minpos.push_back(l - 1);

  scanMinHash(doc, hashfunc, l, r, mulWordIds, x, minhash, minpos);

  int q = (int)minpos.size() - 1;

  for (int i = 0; i <= q - x; i++)
  {
    if (le >= minpos[i + 1])
    {
      // TODO:  hfunc, freq, minhash
      if (r - minpos[i] >= len_threshold && minpos[i + 1] >= minpos[i] + 1 && r >= minpos[i + x]) {
        cws.emplace_back(minpos[i] + 1, minpos[i + 1], minpos[i + x], r, hashfunc, x, minhash);
      }
      if (x == 1)
      {

        multiDivideConquer(doc, minpos[i] + 1, minpos[i + 1] - 1, minpos[i + x] - 1, cws, hashfunc, len_threshold, mulWordIds, ppos);
      }
      else
      {

        multiDivideConquer(doc, minpos[i] + 1, minpos[i + 1], minpos[i + x] - 1, cws, hashfunc, len_threshold, mulWordIds, ppos);
      }
    }
    else
    {
      if (r - minpos[i] >= len_threshold && le >= minpos[i] + 1 && r >= minpos[i + x]) {
        cws.emplace_back(minpos[i] + 1, le, minpos[i + x], r, hashfunc, x, minhash);
      }
      multiDivideConquer(doc, minpos[i] + 1, le, minpos[i + x] - 1, cws, hashfunc, len_threshold, mulWordIds, ppos);
      return;
    }
  }
  multiDivideConquer(doc, minpos[q - x + 1] + 1, le, r, cws, hashfunc, len_threshold, mulWordIds, ppos);
  return;
}


// Do multi-divide-conquer to generate compact windows
// Algorithm 9
void multiDivCompactWindows(const vector<int> &doc, const int hashfunc, vector<CW> &cws, const int minlen, const vector<vector<int>> &mulWordIds, const vector<int> &ppos)
{
  multiDivideConquer(doc, 0, doc.size() - 1, doc.size() - 1, cws, hashfunc, minlen, mulWordIds, ppos);
}

// Do multi-divide-conquer to generate compact windows
// Algorithm 10
void multiDivideConquer(const vector<int> &doc, int l, int le, int r, vector<CW> &cws, const int hashfunc, const int len_threshold, const vector<vector<int>> &mulWordIds)
{
  if (r - l + 1 < len_threshold)
    return;

  int x; // multiplicity
  uint minhash;
  vector<int> minpos;
  minpos.push_back(l - 1);

  scanMinHash(doc, hashfunc, l, r, mulWordIds, x, minhash, minpos);

  int q = (int)minpos.size() - 1;

  for (int i = 0; i <= q - x; i++)
  {
    if (le >= minpos[i + 1])
    {
      // TODO:  hfunc, freq, minhash
      if (r - minpos[i] >= len_threshold && minpos[i + 1] >= minpos[i] + 1 && r >= minpos[i + x]) {
        cws.emplace_back(minpos[i] + 1, minpos[i + 1], minpos[i + x], r, hashfunc, x, minhash);
      }
      if (x == 1)
      {
        multiDivideConquer(doc, minpos[i] + 1, minpos[i + 1] - 1, minpos[i + x] - 1, cws, hashfunc, len_threshold, mulWordIds);
      }
      else
      {
        multiDivideConquer(doc, minpos[i] + 1, minpos[i + 1], minpos[i + x] - 1, cws, hashfunc, len_threshold, mulWordIds);
      }
    }
    else
    {
      if (r - minpos[i] >= len_threshold && le >= minpos[i] + 1 && r >= minpos[i + x]) {
        cws.emplace_back(minpos[i] + 1, le, minpos[i + x], r, hashfunc, x, minhash);
      }
      multiDivideConquer(doc, minpos[i] + 1, le, minpos[i + x] - 1, cws, hashfunc, len_threshold, mulWordIds);
      return;
    }
  }
  multiDivideConquer(doc, minpos[q - x + 1] + 1, le, r, cws, hashfunc, len_threshold, mulWordIds);
  return;
}


// Do multi-divide-conquer to generate compact windows
// Algorithm 9
void multiDivCompactWindows(const vector<int> &doc, const int hashfunc, vector<CW> &cws, const int minlen, const vector<vector<int>> &mulWordIds)
{
  multiDivideConquer(doc, 0, doc.size() - 1, doc.size() - 1, cws, hashfunc, minlen, mulWordIds);
}

// Do verification on compact windows
void verifyCW(const vector<int> &doc, const vector<CW> &cws, int k, const vector<vector<int>> &mulWordIds) {

  // Verify 1: each window should be represented in CW
  for (int l = 0; l < doc.size(); ++l) {
    printf("%d/%d\n", l, doc.size());
    for (int r = l; r < doc.size(); ++r) {
      // Check [l, r] appeaks only in one windows
      int cnt = 0;
      for (auto &cw: cws) {
        cnt += l >= cw.ll && l <= cw.lr && r >= cw.rl && r <= cw.rr;
      }
      assert(cnt == k);
    }
  }

  // Verify 2: each compact window is valid
  for (int i = 0; i < cws.size(); ++i) {
    printf("%d/%d\n", i, cws.size());
    auto &cw = cws[i];
    // Check condition 1 in definiton 6.1
    assert (cw.ll >= 0);
    assert (cw.ll <= cw.lr);
    assert (cw.lr <= cw.rl);
    assert (cw.rl <= cw.rr);
    assert (cw.rr < doc.size());

    // Check the second condiiton
    for (int l = cw.ll; l <= cw.lr; ++l) {
      for (int r = cw.rl; r <= cw.rr; ++r) {
        unordered_map<int, int> freq;
        uint minHash;
        for (int x = l; x <= r; ++x) {
          freq[doc[x]] += 1;
          uint hv = hval(doc[x], cw.hfunc, freq[doc[x]], mulWordIds);
          if (x == l || hv <= minHash)
          {
            minHash = hv;
          }
        }
        assert(minHash == cw.minhash);
      }
    }

  }

}

// printing info
void wordmatchcnt(const vector<int> &doc1, const vector<int> &doc2)
{
  unordered_map<int, int> doccount;
  for (auto &entry : doc1)
    doccount[entry] += 1;
  int paircnt = 0;
  for (auto &entry : doc2)
    paircnt += doccount[entry];
  cout << "matched word pair: " << paircnt << endl;
}

// Generate random hash functions based on given seed
// ab: the output argument
void generateHashFunc(unsigned int seed, vector<pair<int, int>> &ab) {
    // TODO: knuth shuffling?
    // TODO: random_seed
    srand(seed);
    int a = 0;
    while (a == 0)
      a = rand();
    int b = rand();
    ab.emplace_back(a, b);
}

// Get subset of Wx's subset that has collision with subset of Wy
vector<CW> getCWColission(const vector<CW> &Wx, const vector<CW> &Wy) {
  vector<CW> res;
  unordered_set<uint64_t> hashAndFunc;

  for (const auto &y : Wy) {
    hashAndFunc.insert(y.hfunc << 32 | y.minhash);
  }

  for (const auto &x : Wx) {
    if (hashAndFunc.find(x.hfunc << 32 | x.minhash) != hashAndFunc.end()) {
      res.emplace_back(x);
    }
  }

  return res;
}

// Get subset of Wx's subset that has collision with subset of Wy
vector<CW> getCWColission(const vector<CW> &Wx, const vector<int> &WxId, const vector<CW> &Wy, const vector<int> &WyId) {
  vector<CW> res;
  unordered_set<uint64_t> hashAndFunc;

  for (auto yId : WyId) {
    const auto &y = Wy[yId];
    hashAndFunc.insert(y.hfunc << 32 | y.minhash);
  }

  for (auto xId : WxId) {
    const auto &x = Wx[xId];
    if (hashAndFunc.find(x.hfunc << 32 | x.minhash) != hashAndFunc.end()) {
      res.emplace_back(x);
    }
  }
  return res;
}

// Get subset of Wy that has collision with subset of Wx
vector<CW> getCWColission(const vector<CW> &Wx, const vector<int> &WxId, const vector<CW> &Wy) {
  vector<CW> res;
  unordered_set<uint64_t> hashAndFunc;

  for (const auto &xId : WxId) {
    hashAndFunc.insert(Wx[xId].hfunc << 32 | Wx[xId].minhash);
  }

  for (const auto &y : Wy) {
    if (hashAndFunc.find(y.hfunc << 32 | y.minhash) != hashAndFunc.end()) {
      res.emplace_back(y);
    }
  }

  return res;
}

// Update intersection range based on given compact window
void updateIntersection(const CW &cw, IntervalType type, Range &intersection, bool &first) {
  if (first) {
    // First time ..
    if (type == LEFT) {
      intersection.l = cw.ll;
      intersection.r = cw.lr;
    } else {
      intersection.l = cw.rl;
      intersection.r = cw.rr;
    }
    first = false;
  } else {
    if (type == LEFT) {
      intersection.l = max(intersection.l, cw.ll);
      intersection.r = min(intersection.r, cw.lr);
    } else {
      intersection.l = max(intersection.l, cw.rl);
      intersection.r = min(intersection.r, cw.rr);
    }
  }
}

// Get the intersection of compact windows 
// type == LEFT => left interval intersection
// type == RIGHT => right interval intersection
Range getIntersection(const vector<CW> &cws, IntervalType type) {
  // At least one elements ..
  assert(cws.size() > 0);
  Range intersection;
  bool first = true;
  for (const auto &cw : cws) {
    updateIntersection(cw, type, intersection, first);
  }
  return intersection;
}

// Get the intersection of compact windows 
// type == LEFT => left interval intersection
// type == RIGHT => right interval intersection
Range getIntersection(const vector<CW> &cws, const vector<int> &cwId, IntervalType type) {
  // At least one elements ..
  assert(cwId.size() > 0);
  Range intersection;
  bool first = true;
  for (auto id : cwId) {
    const auto &cw = cws[id];
    updateIntersection(cw, type, intersection, first);
  }
  return intersection;
}

// Find subsets according to certain condition
// a : the output argument
void findSubsets(ST &st, int minSize, const vector<CW> &w, vector<vector<int>> &a) {
  // Reuse st
  st.updateTimeStamp();
  // Candidates

  // First filter by left interval
  int cwId = 0;
  for (const auto & cw: w) {
    st.pushDown(cw.ll, cw.lr, cwId++, LEFT);
  }

  // Have correctly verified ..
  vector<vector<int>> aPrime;
  // Refine by first time (based on left interval)
  st.refine(minSize, w, aPrime, LEFT);
  for (auto &cws : aPrime) {
    // Reuse st
    st.updateTimeStamp();
    for (auto cwId : cws) {
      st.pushDown(w[cwId].rl, w[cwId].rr, cwId, RIGHT);
    }
    // Refine based on right interval
    st.refine(minSize, w, a, RIGHT);
  }
}



// Check whether under same sentence
bool checkPos(int l1, int r1, int l2, int r2) {
  return true;
}



vector<CW> sentenceModify(const vector<CW> &cws, const vector<int> &periodPos, int len) {
  vector<CW> res;
  for (const auto &cw : cws) {
    #ifdef NORMAL
    if (cw.lr < periodPos[0] || cw.ll > periodPos[periodPos.size() - 1]) {
      continue;
    }
    if (cw.rr < periodPos[0] || cw.rl > periodPos[periodPos.size() - 1]) {
      continue;
    }
    #endif
    int left;
    int right;
    auto pos1 = binaryFind(cw.ll, periodPos, LARGE);
    auto pos2 = binaryFind(cw.rr, periodPos, SMALL);
    if (pos1 < pos2) {
      if (periodPos[pos2] - periodPos[pos1] < len) {
        continue;
      }
      int new_ll = periodPos[pos1];
      int new_rr = periodPos[pos2] - 1;
      if (new_ll <= cw.lr && new_rr >= cw.rl) {
        CW tmp = cw;
        tmp.ll = new_ll;
        tmp.rr = new_rr;
        tmp.minSent = pos1;
        tmp.maxSent = pos2;
        res.emplace_back(tmp);
      }
    }

  }
  return res;
}

void generateRet(string docFileName, string queryFileName, const vector<ALLIGN::windowPair> &pairs){
  for (int i = 0; i < pairs.size(); i++){
    cout << "<feature name=\"detected-plagiarism\" source_reference=\"" << docFileName << " \"source_length=\"" << pairs[i].docLen;
	  cout << "\" source_offset=\"" << pairs[i].docOffset << "\" query_reference=\"" << queryFileName << "\" query_length=\"";
	  cout << pairs[i].queryLen << "\" query_offset=\"" << pairs[i].queryOffset << "\"/>" << endl;
  }
  
}


int main(int argc, char ** argv) 
{
  auto start_time = chrono::high_resolution_clock::now();
  srand(1234);
  // parameters 
  string docFileName; // data text
  string queryFileName; // query text
  double theta; // jac threshold
  int tau; // length threshold 
  int k;  // number of permutations
  string swFileName = "stopwords.txt"; // stopwords file with default name
  bool longest = true;                    // Match longest
  bool sentence = false;
  unordered_set<string> stopWords;
  if (argc < 11) {
    cerr << "Example: ./allign -docFileName data.txt -queryFileName query.txt -theta 0.5 -tau 50 -k 100" << endl;
    return -1;
  }

  // Shared data structure for two documents
  // 1) word to wordId map
  unordered_map<string, int> word2id; 
  // 2) wordId to word map
  vector<string> id2word;
  // 3) wordId to frequency map
  vector<int> id2freq;
  // 4) mulWordId for each repeated wordId
  vector<vector<int>> mulWordIds;  // the mulWordIds ids in negative of each word


  // Data structure for query 
  // 1) query word id
  vector<int> query;
  // 2) query peroid position
  vector<int> queryPPos;
  // 3) query compact windows
  vector<CW> queryCW;
  // 4) query token offsets
  vector<pair<int, int>> queryOffsets;

  // Data structre for doc
  vector<int> doc;
  vector<int> docPPos;
  vector<CW> docCW;
  // 4) data token offsets
  vector<pair<int, int>> docOffsets;
  int n = 1000;
  for (int i = 0; i < argc; i++) {
    string arg = argv[i];
    if (arg == "-docFileName")         docFileName         = string(argv[i+1]);
    if (arg == "-queryFileName")       queryFileName        = string(argv[i+1]);
    if (arg == "-theta")               theta          = atof(argv[i+1]);
    if (arg == "-tau")                 tau            = atoi(argv[i+1]);
    if (arg == "-k")                   k              = atoi(argv[i+1]);
    if (arg == "-swFileName")          swFileName  = string(argv[i+1]);    if (arg == "-longest")             longest        = atoi(argv[i+1]);
    if (arg == "-sentence")            sentence       = atoi(argv[i+1]);
    if (arg == "-n")                   n              = atoi(argv[i+1]);
  }

  timeval t_start, t_end;
  timeval t1, t2, t3, t4;
  float time_query_window = 0;
  float time_doc_window = 0;
  // Read stop words
  readStopWords(swFileName, stopWords);
  gettimeofday(&t1, NULL);
  // Preprocess documents
  gettimeofday(&t_start, NULL);
  word2int(queryFileName, query, queryOffsets, word2id, id2word, id2freq, mulWordIds, stopWords);
  
  gettimeofday(&t_end, NULL);
  time_query_window += t_end.tv_sec - t_start.tv_sec + (t_end.tv_usec - t_start.tv_usec) / 1e6;
  word2int(docFileName, doc, docOffsets, word2id, id2word, id2freq, mulWordIds, stopWords);
  
  cout << "the size of the doc file is: " << doc.size() << endl;
  cout << "the size of the query file is: " << query.size() << endl;
  

  int size1 = 0;
  int size2 = 0;
  int size3 = 0;
  // Recrod time ...
  
  #ifdef RMQ
  querySegTrees.resize(k, nullptr);
  docSegTrees.resize(k, nullptr);
  #endif
  // std::cout << query.size()
  // Generate compact windows for original test and query text
  
  for (int round = 0; round < k; round++) {
    // Generate a random universal hash function h_i
    generateHashFunc(round, ab);

    #ifdef RMQ 
    // Build new segtree
    querySegTrees[round] = new RMQSegTree(query, round, hval);
    querySegTrees[round]->build(mulWordIds, tau);

    // Build new seg tree
    docSegTrees[round] = new RMQSegTree(doc, round, hval);
    docSegTrees[round]->build(mulWordIds, tau);
    #endif

    #ifdef DEBUG
    int l = 0; 
    for (; l < query.size(); l++) {
      for (int r = l + tau; r < query.size(); r++) {
        int mul1;
        uint minHash1;
        vector<int> pos1;
        querySegTrees[round]->scanMinHash(l, r, mul1, minHash1, pos1);

        int mul2;
        uint minHash2;
        vector<int> pos2;
        scanMinHash(query, round, l, r, mulWordIds, mul2, minHash2, pos2);
        assert(mul1 == mul2);
        assert(minHash1 == minHash2);
        assert(pos1.size() == pos2.size());
        for (int i = 0; i < pos1.size(); ++i) {
          assert(pos1[i] == pos2[i]);
        }
      }
      std::cout << l << std::endl;
    }
    #endif 

    #ifdef RMQ
    ptr = &querySegTrees;
    #endif
    // Generate compact windows for query text (doc1)
    gettimeofday(&t_start, NULL);
    multiDivCompactWindows(query, round, queryCW, tau, mulWordIds);
    gettimeofday(&t_end, NULL);
    time_query_window += t_end.tv_sec - t_start.tv_sec + (t_end.tv_usec - t_start.tv_usec) / 1e6;
    #ifdef RMQ
    ptr = &docSegTrees;
    #endif
    // Generate compact windows for doc text (doc2)
    gettimeofday(&t_start, NULL);
    multiDivCompactWindows(doc, round, docCW, tau, mulWordIds);
    gettimeofday(&t_end, NULL);
    time_doc_window += t_end.tv_sec - t_start.tv_sec + (t_end.tv_usec - t_start.tv_usec) / 1e6;
  }

  // TODO: sentence level ...
  if (sentence) {
    queryCW = sentenceModify(queryCW, queryPPos, tau);
    docCW = sentenceModify(docCW, docPPos, tau);
  }


  size1 = docCW.size();
  size2 = queryCW.size();
  
  gettimeofday(&t2, NULL);
  gettimeofday(&t3, NULL);
  vector<ALLIGN::windowPair> pairs;
  cout << "info for parameters: " << endl;
  cout << "current k is: " << k << endl; 
  cout << "current threshold is: " << theta << endl;
  cout << "current n is: " << n << endl;
  cout << "current tau is: " << tau << endl;

  cout << endl;
  // Get subset of queryCW that collide with docCW
  queryCW = getCWColission(queryCW, docCW);
  size3 = queryCW.size();

  // Build a segment tree at this length (Keep the longest elements)
  ST stree(max(query.size(), doc.size()), longest);
  vector<vector<int>> Ax;
  // First find subsets for the compact windows of text query
  findSubsets(stree, ceil(theta * k), queryCW, Ax);



  std::unordered_set<int64_t> sentencePairs;
  int count = 0;
  for (const auto &WxPrime : Ax) {
    count += 1;
    // std::cout << WxPrime.size() << std::endl;
    // Get subset of docCW that has collision with WxPrime

    auto WyPrime = getCWColission(queryCW, WxPrime, docCW);

    vector<vector<int>> Ay;

    findSubsets(stree, ceil(theta * k), WyPrime, Ay);

    for (const auto &Wy2Prime : Ay) {
      auto Wx2Prime = getCWColission(queryCW, WxPrime, WyPrime, Wy2Prime);
      if (Wx2Prime.size() == 0) {
        continue;
      }
      
      // get lx range
      auto lxRange = getIntersection(Wx2Prime, LEFT);
      // get rx range
      auto rxRange = getIntersection(Wx2Prime, RIGHT);
      // get ly range
      auto lyRange = getIntersection(WyPrime, Wy2Prime, LEFT);
      // get ry range
      auto ryRange = getIntersection(WyPrime, Wy2Prime, RIGHT);

      
      // pairs.emplace_back(ALLIGN::windowPair(lxRange.l, rxRange.r - lxRange.l + 1, lyRange.l, ryRange.r - lyRange.l + 1));
      // auto docOffset = binaryFind(lyRange.l, docPPos, SMALL);
      // auto docLen = binaryFind(ryRange.r, docPPos, LARGE) - docOffset;

      // auto queryOffset = binaryFind(lxRange.l, queryPPos, SMALL);
      // auto queryLen = binaryFind(rxRange.r, queryPPos, LARGE) - queryOffset;
  
      int docOffset = docOffsets[lyRange.l].first;
      int docLen = docOffsets[ryRange.r].second - docOffset + 1;
      int queryOffset = queryOffsets[lxRange.l].first;
      int queryLen = queryOffsets[rxRange.r].second - queryOffset + 1; 
      // Insert into vector if no nested window pairs
      bool flag = true;
      for (int j = 0; j < pairs.size(); ++j) {
        if (docOffset >= pairs[j].docOffset && docOffset + docLen <= pairs[j].docOffset + pairs[j].docLen) {
          if (queryOffset >= pairs[j].queryOffset && queryOffset + queryLen <= pairs[j].queryOffset + pairs[j].queryLen) {
            flag = false;
            break;
          }
        }
      }
      if (!flag) {
        continue;
      }
      for (int j = 0; j < pairs.size(); ++j) {
        if (docOffset <= pairs[j].docOffset && docOffset + docLen >= pairs[j].docOffset + pairs[j].docLen) {
          if (queryOffset <= pairs[j].queryOffset && queryOffset + queryLen >= pairs[j].queryOffset + pairs[j].queryLen) {
            pairs[j] = ALLIGN::windowPair(docOffset, docLen, queryOffset, queryLen);
            flag = false;
            break;
          }
        }
      }
      if (flag) {
        pairs.emplace_back(ALLIGN::windowPair(docOffset, docLen, queryOffset, queryLen));
      }

    }
  }
  // TODO: need to more modification
  cout << "start reduction" << endl;
  std::vector<ALLIGN::windowPair> reduced;
  // Reduce
  for (int i = 0; i < pairs.size(); ++i) {
    int docOffset = pairs[i].docOffset;
    int docLen = pairs[i].docLen;
    int queryOffset = pairs[i].queryOffset;
    int queryLen = pairs[i].queryLen;
    bool flag = true;
    for (int j = 0; j < i; ++j) {
      if (docOffset >= pairs[j].docOffset && docOffset + docLen <= pairs[j].docOffset + pairs[j].docLen) {
        if (queryOffset >= pairs[j].queryOffset && queryOffset + queryLen <= pairs[j].queryOffset + pairs[j].queryLen) {
          flag = false;
          break;
        }
      }
    }
    if (!flag) {
      continue;
    }
    // for (int j = i + 1; j < pairs.size(); ++j) {
    //     if (docOffset <= pairs[j].docOffset && docOffset + docLen >= pairs[j].docOffset + pairs[j].docLen) {
    //       if (queryOffset <= pairs[j].queryOffset && queryOffset + queryLen >= pairs[j].queryOffset + pairs[j].queryLen) {
    //         flag = false;
    //         break;
    //       }
    //     }

    // }
    if (flag) {
      reduced.push_back(pairs[i]);
    }
  }
  gettimeofday(&t4, NULL);
  //print out result
  generateRet(docFileName, queryFileName, reduced);
  
  //ALLIGN::generateXML(docFileName.substr(docFileName.length() - 24, 24), queryFileName.substr(queryFileName.length() - 28, 28), reduced);
  // ALLIGN::generateXML(docFileName.substr(docFileName.length() - 24, 24), queryFileName.substr(queryFileName.length() - 28, 28), pairs);

  //cout << "genrate two same cws: " << t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1e6 <<endl;
  
  //cout << "Generate query compact window pairs time: " << time_query_window << "s" << endl;
  //cout << " Find window pairs " << t4.tv_sec - t3.tv_sec + (t4.tv_usec - t3.tv_usec) / 1e6 << endl;
  //cout << "the size of doc's compact window is: " << size1 <<endl;
  //cout << "the size of query's compact window is: " << size2 << endl;
  //cout << "total query time(not included build cw) is: " << t4.tv_sec - t3.tv_sec + (t4.tv_usec - t3.tv_usec) / 1e6 << "s" << endl;
  // cout << query.size() << " " << doc.size() << endl;
  auto finish_time = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::microseconds>(finish_time - start_time);
  cout << "total time is: " << duration.count() / 1000000.0 << " seconds" << endl;
  return 0;
}
