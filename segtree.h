#include <fstream>
#include <cmath>
#include <iostream>
#include <sys/time.h>
#include <vector>
#include <unordered_map>
#include <cassert>
#include <cstring>
#include <unordered_set>
#include <string>
#include <ctype.h>
#include <algorithm>
#include <functional>


using namespace std;
const int INVALID_RANGE = -1;
const int INVALID_ID = -1;
// intersection
struct Range {
  int l;
  int r;
  Range(): l(INVALID_RANGE), r(INVALID_RANGE) {}
  Range(int l, int r): l(l), r(r) {}
  Range(const pair<int, int> &p): l(p.first), r(p.second) {}
};

// Specify left interval / right interval 
enum IntervalType {
  LEFT,
  RIGHT
};

// compact window
class CW 
{
  public:
    int ll;
    int lr;
    int rl;
    int rr;
    int hfunc;
    int freq;
    uint minhash;
    // Default empty constructor
    int sentence;
    int minSent;
    int maxSent;

    CW() {}
    CW(int l1, int l2, int r1, int r2, int h, int cnt, uint hv)
      : ll(l1), lr(l2), rl(r1), rr(r2), hfunc(h), freq(cnt), minhash(hv) { }
    void display() const {
      printf("(%d, %d, %d, %d, %d, %d, %ud)\n", ll, lr, rl, rr, hfunc, freq, minhash);
    }
};

// ST Node used in findsubsets algorithm
class STNode {
  public:
  int l;  // Left interval
  int r;  // Right interval
  vector<int> set; // Store compact window ids 
  int localTimeStamp;
  STNode() : l(INVALID_RANGE), r(INVALID_RANGE) {}
};

// ST used in findsubsets algorithm
class ST {

public:
  vector<STNode> ranges;
  int leafNum;
  int root;
  int globalTimeStamp;
  bool longest;
  int nodesNum;
  // Default constructor function
  ST() {}

  // Constructor function with target size
  ST(int leafNum, bool longest): root(0), globalTimeStamp(0), leafNum(leafNum), longest(longest), ranges(4 * leafNum), nodesNum(0) {
    // Build this tree
    buildST(leafNum);
  }

  // Build the Segment Tree at Specified size
  void buildST(int leafNum) {
    buildSTHelper(0, 0, leafNum - 1);
  }

  void printRanges() {
    for (const auto &node : ranges) {
      printf("(%d, %d)\n", node.l, node.r);
    }
  }

  // Update timestamp to reuse the tree ..
  inline void updateTimeStamp() {
    globalTimeStamp++;
  }

  // Push down operation
  void pushDown(int l, int r, int cwId, IntervalType type) {
    // Notice: root is 0
    pushDownHelper(root, l, r, cwId, type);
  }

  // refine the resul
  void refine(int minSize, const vector<CW> &cws, vector<vector<int>> &A, IntervalType type) {

    // Record number of different hash functions
    unordered_set<int> numHash;
    // Get the union of set in ST[node] and all its ancestors
    vector<int> W;

    for (int node = nodesNum - 1; node >= 0; --node) {
      if (ranges[node].localTimeStamp != globalTimeStamp || ranges[node].set.size() == 0) {
        continue;
      }
      if (longest) {
        int curPos = node;
        bool flag = true;
        while (1) {
          // Move to leftmost / rightmost descents 
          if (type == LEFT) {
            curPos = (curPos << 1) + 1;
          } else {
            curPos = (curPos + 1) << 1;
          }
          // Have iterated all leftmost descents
          if (curPos >= nodesNum) {
            break;
          }
          // Check whether have add curPos to the C
          if (ranges[curPos].localTimeStamp == globalTimeStamp && ranges[curPos].set.size() > 0) {
            flag = false;
            break;
          }
        }
        // One of its descents have been added to C
        if (!flag) {
          continue;
        }
      }

      W.clear();
      auto curNode = node;

      while (1) {
        // Only merge correct time stamp node
        if (ranges[curNode].localTimeStamp == globalTimeStamp) {
          for (auto p : ranges[curNode].set) {
            W.emplace_back(p);
          }
        }
        if (curNode == root) {
          break;
        }
        // Move up
        curNode = (curNode - 1) >> 1;
      }

      if (W.size() < minSize) {
        continue;
      }

      // Sort it
      sort(W.begin(), W.end());

      // Record number of different hash functions
      numHash.clear();
      for (auto cwId: W) {
        numHash.insert(cws[cwId].hfunc);
      }

      if (numHash.size() >= minSize) {
        // De duplicate A 
        if (find(A.begin(), A.end(), W) == A.end()) {
          A.emplace_back(W);
        }
      }
    }
  }

private:

  inline bool intersect(int l1, int r1, int l2, int r2) {
    return l2 <= r1 && l1 <= r2;
  }

  void buildSTHelper(const int node, const int segSt, const int segEn) {
    // This is the leaf node
    if (segSt == segEn) {
      ranges[node].l = ranges[node].r = segSt;
      ranges[node].localTimeStamp = 0;
      // Update nodesNum
      nodesNum = max(node + 1, nodesNum);
    } else {
      // Not a leaf node
      int mid = (segSt + segEn) >> 1;
      int leftChild = (node << 1) + 1;
      int rightChild = (node << 1) + 2;

      // First build left
      buildSTHelper(leftChild, segSt, mid);
      // Next build right
      buildSTHelper(rightChild, mid + 1, segEn);

      // Assertion !
      assert(ranges[leftChild].r + 1 == ranges[rightChild].l);

      // Set left and right value
      ranges[node].l = ranges[leftChild].l;
      ranges[node].r = ranges[rightChild].r;
      ranges[node].localTimeStamp = 0;
    }
  }

  void pushDownHelper(int node, int l, int r, int cwId, IntervalType type) {
    if (ranges[node].l >= l && ranges[node].r <= r) {
      if (ranges[node].localTimeStamp != globalTimeStamp) {
        // Clear set 
        ranges[node].set.clear();
        ranges[node].localTimeStamp = globalTimeStamp;
      }
      ranges[node].set.emplace_back(cwId);
    } else {
      int leftChild = (node << 1) + 1;
      int rightChild = leftChild + 1;
      if (intersect(l, r, ranges[leftChild].l, ranges[leftChild].r)) {
        pushDownHelper(leftChild, l, r, cwId, type);
      }
      if (intersect(l, r, ranges[rightChild].l, ranges[rightChild].r)) {
        pushDownHelper(rightChild, l, r, cwId, type);
      }
    }
  }
};


// Interval used in RMQSegTree
class Interval: public CW {
  public:
  // Inherit all fileds from CW
  // Add word id field (that has the minimum hash value)
  int wid;
  // Add real appearance times in the document
  int realRank;
  Interval(int l1, int l2, int r1, int r2, int cnt, uint hv, int wid, int realRank):
    wid(wid), realRank(realRank) {
      ll = l1;
      lr = l2;
      rl = r1;
      rr = r2;
      freq = cnt;
      minhash = hv;
    }
};

class RMQSegTree2D {
  public:
  vector<int> intervalIds;  // Store the intervalIds
  int leafNum;
  int nodeNum;
  STNode *ranges;
  int *optIds;        // Store the optimal interval id (with the minimum hash value)
                      // And its minimum hash value

  ~RMQSegTree2D() {
    // if (ranges) {
    //   delete(ranges);
    // }
    // if (optIds) {
    //   delete(optIds);
    // }
  }

  // Get segments based on points set
  void getSegments(vector<int> &pointVec, vector<pair<int, int>> &segments) {
    // First sort the points
    sort(pointVec.begin(), pointVec.end());

    if (pointVec.size() > 0) {
      int prev = pointVec[0];
      segments.reserve(2 * pointVec.size() - 1);
      segments.emplace_back(prev, prev);
      int num = 1;
      for (int i = 1; i < pointVec.size(); ++i) {
        int cur = pointVec[i];
        if (cur == prev) {
          continue;
        }
        if (prev + 1 <= cur - 1) {
          segments.emplace_back(prev + 1, cur - 1);
          num++;
        }
        segments.emplace_back(cur, cur);
        num++;
        prev = cur;
      }
    }
  }


  // Build Segment tree for the end ranges
  void buildST(const vector<pair<int, int>> &segments) {
    leafNum = segments.size();
    nodeNum = 0;
    ranges = new STNode[4 * leafNum];
    buildSTHelper(0, 0, leafNum - 1, segments);
  }



  // Build segment tree based on ending positions
  void build(const vector<Interval> &intervals) {
    // Store the ending pointer of interval ends  
    vector<int> intervalEnds(intervalIds.size() * 2);
    int cnt = 0;
    for (auto intervalId: intervalIds) {
      const auto &interval = intervals[intervalId];
      intervalEnds[cnt++] =interval.rl;
      intervalEnds[cnt++] =interval.rr;
    }

    // Now, get the segments for the end positions
    vector<pair<int, int>> endSegments;
    getSegments(intervalEnds, endSegments);

    // Next, build segment tree for the ending ranges
    buildST(endSegments);

    optIds = new int[4 * leafNum];
    memset(optIds, INVALID_ID, 4 * leafNum * sizeof(int));

    // Next, insert each interval id to trees
    for (auto id: intervalIds) {
      insertInterval(0, id, intervals);
    }
  }

  // Build segment tree based on ending positions
  void build(const vector<Interval> &intervals, int ll, int lr, int tau) {
    // Store the ending pointer of interval ends  
    vector<int> intervalEnds;
    intervalEnds.reserve(intervalIds.size() * 2);
    int minL = max(lr, ll + tau - 1);
    int cnt = 0;

    intervalEnds.emplace_back(minL);
    for (auto intervalId: intervalIds) {
      const auto &interval = intervals[intervalId];
      if (interval.rl > minL) {
        intervalEnds.emplace_back(interval.rl);
      }
      if (interval.rr > minL) {
        intervalEnds.emplace_back(interval.rr);
      }
    }

    // Now, get the segments for the end positions
    vector<pair<int, int>> endSegments;
    getSegments(intervalEnds, endSegments);

    // Next, build segment tree for the ending ranges
    buildST(endSegments);

    optIds = new int[4 * leafNum];
    memset(optIds, INVALID_ID, 4 * leafNum * sizeof(int));

    // Next, insert each interval id to trees
    for (auto id: intervalIds) {
      insertInterval(0, id, intervals);
    }
  }


  // Insert interval
  void insertInterval(int node, int intervalId, const vector<Interval> &intervals) {
    const auto &interval = intervals[intervalId];
    // Record the minimum hash value point
    if (interval.rl <= ranges[node].l && ranges[node].r <= interval.rr) {
      if (optIds[node] == INVALID_ID) {
        optIds[node] = intervalId;
      } else {
        if (intervals[optIds[node]].minhash > interval.minhash) {
          optIds[node] = intervalId;
        }
      }
    } else {
      int leftChild = (node << 1) + 1;
      int rightChild = leftChild + 1;
      // Insert interval to the left child
      if (intersect(interval.rl, interval.rr, ranges[leftChild].l, ranges[leftChild].r)) {
        insertInterval(leftChild, intervalId, intervals);
      }
      // Insert interval to the right child
      if (intersect(interval.rl, interval.rr, ranges[rightChild].l, ranges[rightChild].r)) {
        insertInterval(rightChild, intervalId, intervals);
      }
    }
  }

  // Query the minimum hash interval
  int query(int r, const vector<Interval> &intervals) {
    int optId = INVALID_ID;
    if (intervalIds.size() > 0) {
      if (ranges[0].l <= r && r <= ranges[0].r) {
        queryHelper(0, r, intervals, optId);
      }
    }
    return optId;
  }

  private:

  inline bool intersect(int l1, int r1, int l2, int r2) {
    return l2 <= r1 && l1 <= r2;
  }

  // Helper function for building segment trees for ending ranges
  void buildSTHelper(int node, int l, int r, const vector<pair<int, int>> &segments) {
    if (l == r) {
      ranges[node].l = segments[l].first;
      ranges[node].r = segments[l].second;
      // Update nodeNum
      nodeNum = max(node + 1, nodeNum);
    } else {
      // Not a leaf node
      int mid = (l + r) >> 1;
      int leftChild = (node << 1) + 1;
      int rightChild = leftChild + 1;

      // First build left
      buildSTHelper(leftChild, l, mid, segments);
      // Next build right
      buildSTHelper(rightChild, mid + 1, r, segments);

      assert(ranges[leftChild].r + 1 == ranges[rightChild].l);

      // Set left and right value
      ranges[node].l = ranges[leftChild].l;
      ranges[node].r = ranges[rightChild].r;
    }
  }

  // Query the minimum hash interval
  void queryHelper(int node, int r, const vector<Interval> &intervals, int &optId) {
    int curId = optIds[node];
    if (curId != INVALID_ID) {
      if (optId == INVALID_ID || intervals[curId].minhash < intervals[optId].minhash) {
        optId = curId;
      }
    }

    int leftChild = (node << 1) + 1;
    int rightChild = leftChild + 1;
    if (leftChild < nodeNum && ranges[leftChild].l <= r && r <= ranges[leftChild].r) {
      queryHelper(leftChild, r, intervals, optId);
    }

    if (rightChild < nodeNum && ranges[rightChild].l <= r && r <= ranges[rightChild].r) {
      queryHelper(rightChild, r, intervals, optId);
    }
  }
};

// RMQSegTree used in find the minimum hash value in range [l, r]
class RMQSegTree {
  public:
  int docLen;                                                        // The length of the document
  int hashFunc;                                                      // The id of the hash function
  unordered_map<int, vector<int>> id2pos;                            // The mapping from word id to position list
  typedef function<uint(int, int, int, 
                        const vector<vector<int>> &)> hashFuncType;  // Define the function type
  hashFuncType hash;                                                 // Hash funcitons sets
  vector<Interval> intervals;                                        // Stores [ll, lr] x [rl, rr] x minHash

  vector<STNode> ranges;                                             // Segment tree for the starting ranges
  vector<RMQSegTree2D*> trees;
  int leafNum;
  int nodeNum;

  // Build RMQ SegTree based on document and hash function
  RMQSegTree(const vector<int> &doc, const int hashFunc, hashFuncType hash): docLen(doc.size()), hashFunc(hashFunc), hash(hash) {
    // Build the invert index from wordId to each apperance position
    for (int i = 0; i < doc.size(); ++i) {
      id2pos[doc[i]].emplace_back(i);
    }
  }

  ~RMQSegTree() {
    // for (int i = 0; i < ranges.size(); ++i) {
    //   if (trees[i]) {
    //     delete(trees[i]);
    //   }
    // }
  }

  // Get segments based on points set
  void getSegments(const unordered_set<int> &points, vector<pair<int, int>> &segments) {
    // First sort the points
    vector<int> pointVec;
    for (auto p: points) {
      pointVec.push_back(p);
    }
    sort(pointVec.begin(), pointVec.end());

    if (pointVec.size() > 0) {
      int prev = pointVec[0];
      segments.emplace_back(prev, prev);
      for (int i = 1; i < pointVec.size(); ++i) {
        int cur = pointVec[i];
        if (prev + 1 <= cur - 1) {
          segments.emplace_back(prev + 1, cur - 1);
        }
        segments.emplace_back(cur, cur);
        prev = cur;
      }
    }

  }

  // Build Segment tree for the start ranges
  void buildST(const vector<pair<int, int>> &segments) {
    leafNum = segments.size();
    nodeNum = 0;
    ranges.resize(4 * leafNum);
    buildSTHelper(0, 0, leafNum - 1, segments);
  }

  // Get intervals
  void getIntervals(unordered_set<int> &intervalStarts, const vector<vector<int>> &mulWordIds) {
    // Intervals start from [ll, lr] end to [rl, rr]
    for (auto &entry: id2pos) {
      // Current word id
      auto wid = entry.first;
      // Corresponding position list
      auto &posList = entry.second;
      // Record the hash value of each appearance
      vector<uint> kthHashList;
      for (int i = 0; i < posList.size(); ++i) {
        kthHashList.emplace_back(hash(wid, hashFunc, i + 1, mulWordIds));
      }

      // each pair (i, j) shows that 
      // the xth occurance has the minimum hash from [1st, (j - 1)th]
      vector<pair<int, int>> minMal;

      int i = 0;
      while (i < kthHashList.size()) {
        int j = i + 1;
        // Find the breaking point
        while (j < kthHashList.size()) {
          if (kthHashList[i] > kthHashList[j]) {
            break;
          }
          j++;
        }
        minMal.emplace_back(i + 1, j + 1);
        i = j;
      }

      for (int k = 0; k < posList.size(); ++k) {
        for (auto &entry: minMal) {
          int ith = entry.first;
          int jth = entry.second;
          // In [1st, (j - 1)th], the ith appearance hash the minimum hash value
          // We need to make sure that the k should be the ith appearance
          // So the first appearance should at position (k - i + 1)
          // And the windows should before the (k + j - i) position
          int startRank = k - ith + 1;
          if (startRank < 0) {
            continue;
          }

          // Possible range [ll, lr] for the starting positions
          int ll = 0;
          if (startRank > 0)  {
            ll = posList[startRank - 1] + 1;
          }
          int lr = posList[startRank];

          // Possible range [rl, rr] for the ending positions
          int endRank = k + jth - ith;
          int rl = posList[k];
          int rr = docLen - 1;
          if (endRank < posList.size()) {
            rr = posList[endRank] - 1;
          }
          intervals.emplace_back(Interval(ll, lr, rl, rr, ith, kthHashList[ith - 1], wid, k));
          intervalStarts.insert(ll);
          intervalStarts.insert(lr);
        }
      }
    }
  }
  
  // Insert interval to trees
  void insertInterval(int node, int intervalId) {
    const auto &interval = intervals[intervalId];
    if (interval.ll <= ranges[node].l && ranges[node].r <= interval.lr) {
      trees[node]->intervalIds.emplace_back(intervalId);
    } else {
      int leftChild = (node << 1) + 1;
      int rightChild = leftChild + 1;
      // Insert interval to the left child
      if (intersect(interval.ll, interval.lr, ranges[leftChild].l, ranges[leftChild].r)) {
        insertInterval(leftChild, intervalId);
      }
      // Insert interval to the right child
      if (intersect(interval.ll, interval.lr, ranges[rightChild].l, ranges[rightChild].r)) {
        insertInterval(rightChild, intervalId);
      }
    }
  }

  // Insert interval to trees
  void insertInterval(int node, int intervalId, int tau) {
    const auto &interval = intervals[intervalId];
    if (interval.ll <= ranges[node].l && ranges[node].r <= interval.lr) {
      if (interval.rr > ranges[node].l + tau - 1) {
        trees[node]->intervalIds.emplace_back(intervalId);
      }
    } else {
      int leftChild = (node << 1) + 1;
      int rightChild = leftChild + 1;
      // Insert interval to the left child
      if (intersect(interval.ll, interval.lr, ranges[leftChild].l, ranges[leftChild].r)) {
        insertInterval(leftChild, intervalId);
      }
      // Insert interval to the right child
      if (intersect(interval.ll, interval.lr, ranges[rightChild].l, ranges[rightChild].r)) {
        insertInterval(rightChild, intervalId);
      }
    }
  }

  // Build RMQSegTree
  void build(const vector<vector<int>> &mulWordIds) {
    // Store the starting pointer of interval starts
    unordered_set<int> intervalStarts;
    getIntervals(intervalStarts, mulWordIds);

    // Now, get the segments for start positions
    vector<pair<int, int>> startSegments;
    getSegments(intervalStarts, startSegments);

    // Next, build segment tree for the starting ranges
    buildST(startSegments);
    trees.resize(ranges.size());

    for (int i = 0; i < ranges.size(); ++i) {
      trees[i] = new RMQSegTree2D();
    }

    // Next, insert each interval id to trees that are covered by starting range
    for (int i = 0; i < intervals.size(); ++i) {
      insertInterval(0, i);
    }

    // Build segment tree for each node
    for (int i = 0; i < nodeNum; ++i) {
      if (!trees[i]->intervalIds.empty()) {
        trees[i]->build(intervals);
      }
    }
  }

  // Build RMQSegTree
  void build(const vector<vector<int>> &mulWordIds, int tau) {
    // Store the starting pointer of interval starts
    unordered_set<int> intervalStarts;
    getIntervals(intervalStarts, mulWordIds);

    // Now, get the segments for start positions
    vector<pair<int, int>> startSegments;
    getSegments(intervalStarts, startSegments);

    // Next, build segment tree for the starting ranges
    buildST(startSegments);
    trees.resize(ranges.size());

    for (int i = 0; i < ranges.size(); ++i) {
      trees[i] = new RMQSegTree2D();
    }

    // Next, insert each interval id to trees that are covered by starting range
    for (int i = 0; i < intervals.size(); ++i) {
      insertInterval(0, i, tau);
    }

    // Build segment tree for each node
    for (int i = 0; i < nodeNum; ++i) {
      if (!trees[i]->intervalIds.empty()) {
        trees[i]->build(intervals, ranges[i].l, ranges[i].r, tau);
      }
    }
  }


  void scanMinHash(int l, int r, int &mul, uint &minHash, vector<int> &pos) {
      int optId = INVALID_ID;
      queryHelper(0, l, r, optId);
      minHash = intervals[optId].minhash;
      mul = intervals[optId].freq;
      int realRank = intervals[optId].realRank;

      auto &posList = id2pos[intervals[optId].wid];
      // Make the position is the same
      for (int i = realRank - mul + 1; ; ++i) {
        if (i >= posList.size() || posList[i] > r) {
          break;
        }
        pos.emplace_back(posList[i]);
      }
  }

  private:

  inline bool intersect(int l1, int r1, int l2, int r2) {
    return l2 <= r1 && l1 <= r2;
  }

  // Helper function for building segment trees for starting ranges
  void buildSTHelper(int node, int l, int r, const vector<pair<int, int>> &segments) {
    if (l == r) {
      ranges[node].l = segments[l].first;
      ranges[node].r = segments[l].second;
      // Update nodeNum
      nodeNum = max(node + 1, nodeNum);
    } else {
      // Not a leaf node
      int mid = (l + r) >> 1;
      int leftChild = (node << 1) + 1;
      int rightChild = leftChild + 1;

      // First build left
      buildSTHelper(leftChild, l, mid, segments);
      // Next build right
      buildSTHelper(rightChild, mid + 1, r, segments);

      assert(ranges[leftChild].r + 1 == ranges[rightChild].l);

      // Set left and right value
      ranges[node].l = ranges[leftChild].l;
      ranges[node].r = ranges[rightChild].r;
    }
  }

  void queryHelper(int node, int l, int r, int &optId) {

    int curId = trees[node]->query(r, intervals);
    if (curId != INVALID_ID) {
      if (optId == INVALID_ID || intervals[curId].minhash < intervals[optId].minhash) {
        optId = curId;
      }
    }

    int leftChild = (node << 1) + 1;
    int rightChild = leftChild + 1;

    if (leftChild < nodeNum && ranges[leftChild].l <= l && l <= ranges[leftChild].r) {
      queryHelper(leftChild, l, r, optId);
    }

    if (rightChild < nodeNum && ranges[rightChild].l <= l && l <= ranges[rightChild].r) {
      queryHelper(rightChild, l, r, optId);
    }
  }

};