//from http://stackoverflow.com/questions/4523220/sorting-a-vector-of-double-precision-reals-and-obtain-their-order

#ifndef PERMSORT
#define PERMSORT

#include <vector>
#include <algorithm>
#include <utility>
//#include <iostream>
#if defined CPP11 
 template<class Vals>
 void sortingPermutation(const Vals& values, std::vector<int>& v){
   int size = values.size(); 
   v.clear(); v.reserve(size);
   for(int i=0; i < size; ++i)
     v.push_back(i);
  
   std::sort(v.begin(), v.end(), [&values](int a, int b) -> bool { 
     return values[a] < values[b];
   });
 }
#else
/*template<class T>
struct CmpPairs{
  CmpPairs(const std::vector<T> &v): v_(v) {}
  std::vector<T> v_;
  bool operator()(int a, int b){ return v_[a] < v_[b]; }
};
template<class T>
CmpPairs<T> CreateCmpPairs(const std::vector<T> & v) { return CmpPairs<T>(v); };

 template<class Vals>
 void sortingPermutation(const Vals& values, std::vector<int>& v){
   int size = values.size(); 
   v.clear(); v.reserve(size);
   for(int i=0; i < size; ++i)
     v.push_back(i);
   std::sort(v.begin(), v.end(),CreateCmpPairs<Vals>(values));//<vals> dosent compile on new gcc
 }
*/
typedef std::pair<string, int> Pair;//here is string explicitly set !! 
typedef std::pair<double, int> PairD;

struct CmpPair
{
    bool operator()(const Pair& a, const Pair& b)
    { return a.first < b.first; }
};
struct CmpPairD
{
    bool operator()(const PairD& a, const PairD& b)
    { return a.first < b.first; }
};
void sortingPermutation(
    const std::vector<string>& values,
    std::vector<int>& permutation)
{       permutation.clear();
	permutation.reserve(values.size());
    std::vector<Pair> pairs;
    for (int i = 0; i < (int)values.size(); i++)
        pairs.push_back(Pair(values[i], i));

    std::sort(pairs.begin(), pairs.end(), CmpPair());
int i=0;
    typedef std::vector<Pair>::const_iterator I;
    for (I p = pairs.begin(); p != pairs.end(); ++p,++i)
        permutation[i]=(p->second);

}
void sortingPermutation(
    const std::vector<double>& values,
    std::vector<int>& permutation)
{       permutation.clear();
	permutation.reserve(values.size());
    std::vector<PairD> pairs;
    for (int i = 0; i < (int)values.size(); i++)
        pairs.push_back(PairD(values[i], i));

    std::sort(pairs.begin(), pairs.end(), CmpPairD());
int i=0;
    typedef std::vector<PairD>::const_iterator I;
    for (I p = pairs.begin(); p != pairs.end(); ++p,++i)
        permutation[i]=(p->second);

}


#endif 
#endif
