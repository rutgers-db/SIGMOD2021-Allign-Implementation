#pragma once

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
#include <cstring>

// #define NORMAL 0
using namespace std;
namespace ALLIGN {

struct windowPair {
  int docOffset;
  int docLen;
  int queryOffset;
  int queryLen;
  windowPair(int docOffset, int docLen, int queryOffset, int queryLen) {
    this->docOffset = docOffset;
    this->docLen = docLen;
    this->queryOffset = queryOffset;
    this->queryLen = queryLen;
  }
};

void pos2ppos(const vector<int> &docppos, vector<int> &docpos2ppos)
{
  for (auto j = 1; j < docppos.size(); j++)
  {
    int beg = docppos[j - 1];
    int end = docppos[j];
    docpos2ppos.insert(docpos2ppos.end(), end - beg, j);
  }
}

// Read stopwords from given file name
void readStopWords(const string &fileName, unordered_set<string> &stopWords) {
  // one per line
  string stopWord;
  ifstream file(fileName, ios::in);
  while (getline(file, stopWord))
    stopWords.insert(stopWord);
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

// Read from file and do preprocessing
// 
// Input parameter:
//  1) filename: the file name
//  2) stopwords: stopwords that don't need to be considered
//
// Output parameter: 
//  1) doc: contains the word id (from word string => word id)
//  2) ppos: contains peroid positions
//  3) word2id: the mapping from word to id
//  4) id2word: the mapping from id to word
//  5) id2maxFreq: the mapping from id to maximum frequency
//  6) id2mulId: if this word appears more than 1 time, it will have a unique mulId for each appearance
void word2int(const string &filename, vector<int> &doc, vector<pair<int, int>> &offsets, vector<int> &ppos, unordered_map<string, int> &word2id, vector<string> &id2word, vector<int> &id2maxFreq, vector<vector<int>> &id2mulId, const unordered_set<string> &stopwords)
{
  // Static variables accross the function
  static int mulWordId = -1;

  #ifdef NORMAL
  const string delim = "\t\n\r\x0b\x0c !\"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~";
  const string period_str = ",.!?;";
  #endif

  #ifndef NORMAL
  const string delim = " ";
  const string period_str = "\n";
  #endif

  // Read from file ...
	ifstream datafile(filename, ios::binary);
  datafile.seekg(0, std::ios_base::end);
  int length = datafile.tellg();
  string docstr(length + 1, '\0');
  datafile.seekg(0);
  datafile.read(&docstr[0], length);


  // Make the docstr to sentences ..
  vector<string> sentences;
  vector<pair<int, int>> sentenceOffsets;
  ALLIGN::strToTokens(docstr, period_str, sentences, sentenceOffsets);

  ppos.emplace_back(0);
  unordered_map<int, int> id2freq;
  for (int s = 0; s < sentences.size(); s++)
  {
    
    // get words in this sentence
    vector<string> tokens;
    vector<pair<int, int>> wordOffsets;
    ALLIGN::strToTokens(sentences[s], delim, tokens, wordOffsets);

    int offset = sentenceOffsets[s].first;
    // Iterate each word in the sentence
    for (int i = 0; i < tokens.size(); i++)
    {

      // Skip stop words
      if (stopwords.find(tokens[i]) != stopwords.end())
        continue;

      // If a new word, add to word2id
      if (word2id.find(tokens[i]) == word2id.end())
      {
        word2id[tokens[i]] = id2word.size();
        id2word.emplace_back(tokens[i]);
      }

      // Insert this word into the doc and increase word frequency
      int wid = word2id[tokens[i]];
      id2freq[wid] += 1;
      offsets.push_back({wordOffsets[i].first + offset, wordOffsets[i].second + offset});
      doc.emplace_back(wid);
    }

    // Record the period position
    if (ppos.back() != doc.size()) {
      ppos.emplace_back(doc.size());
    }
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

void generateXML(string outputAddress,string suspiciousFileName, vector<int> to, vector<int> tl, vector<string> sn, vector<int> so, vector<int> sl) {
    ofstream fout;
    string fn = outputAddress + suspiciousFileName + ".xml";
    fout.open(fn, ios::out);
    fout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
    fout << "<document reference=\"" << suspiciousFileName << ".txt" << "\">" << endl;
    if (to.size() != tl.size() || so.size() != sl.size() || to.size() != so.size() || sn.size() != tl.size()) {
        cout <<"error when generating XML file!"<< endl;
        fout.close();
        return;
    }
    for (size_t i = 0; i < to.size(); i++)
    {
        fout << "<feature name=\"plagiarism\" type=\"artificial\" obfuscation=\"low\" this_language=\"en\" this_offset=\"";
        fout << to.at(i) << "\" this_length=\"" << tl.at(i) << "\" source_reference=\"" << sn.at(i) << "\"source_language=\"en\" source_offset=\"";
        fout << so.at(i) << "\" source_length=\"" << sl.at(i) << "\" />" << endl;
    }
    fout << "</document>" << endl;
    fout.close();
}

void generateXML(string dName, string qName, const vector<windowPair> &pairs)
{
	cout << "<?xml version=\"1.0\" encoding=\"utf-8\"?><document reference=\"" << qName << "\">" << endl;
  for (int i = 0; i < pairs.size(); ++i) {
		cout << "<feature name=\"detected-plagiarism\" source_length=\"" << pairs[i].docLen;
		cout << "\" source_offset=\"" << pairs[i].docOffset << "\" source_reference=\"" << dName << "\" this_length=\"";
		cout << pairs[i].queryLen << "\" this_offset=\"" << pairs[i].queryOffset << "\"/>" << endl;
	}
	cout << "</document>" << endl;

  // #ifndef NORMAL
  // for (int i = 0; i < pairs.size(); ++i) {
  //   cout << pairs[i].docOffset << " " << pairs[i].docOffset + pairs[i].docLen - 1 << " " << pairs[i].queryOffset << " " << pairs[i].queryOffset + pairs[i].queryLen - 1 << std::endl;
  // }
  // #endif
}

void showPairs(const vector<int> &l1P, const vector<int> &r1P, const vector<int> &l2P, const vector<int> &r2P) {
  for (int i = 0; i < l1P.size(); ++i) {
    printf("[%d, %d] [%d, %d]\n", l1P[i], r1P[i], l2P[i], r2P[i]);
  }
}

}