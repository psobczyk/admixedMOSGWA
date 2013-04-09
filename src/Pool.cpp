#include "Pool.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>

using namespace std;



ostream& operator << ( ostream &out, const PoolItem &p )
{
  vector <snp_index_t>::const_iterator it = p.snps.begin();
  out << "[";
  if (p.snps.size() > 0)
  {
    out << *it;
    ++it;
  }  
  for (; it != p.snps.end(); it++)
    out << ", " << *it;
  return out << "], msc = " << setprecision(10) << p.msc;
}

//-----------------------------------------------------------------------------------------------
PoolItem::PoolItem(vector <snp_index_t> &snp_, double msc)
:msc(msc)
{
   snps.insert(snps.begin(), snp_.begin(), snp_.end());
   sort(snps.begin(), snps.end());
}

//-----------------------------------------------------------------------------------------------
PoolItem::PoolItem(const PoolItem &p)
:msc(p.msc)
{
   snps.insert(snps.begin(), p.snps.begin(), p.snps.end());
}

//-----------------------------------------------------------------------------------------------
PoolItem & PoolItem::operator= (const PoolItem &p)
{
  if (this != &p)
  {
    snps.clear();
    snps.insert(snps.begin(), p.snps.begin(), p.snps.end());
    msc = p.msc;
  }
  return *this;
}

//-----------------------------------------------------------------------------------------------
bool PoolItem::operator < (const PoolItem &p) const
{
/*  if (snps.size() <= 0)
  {
    cerr << "operator <, snp.size <= 0!, Press any key" << endl;
    exit(-1);
  }
  */
  if (snps.size() == p.snps.size())
  {
    for (unsigned int i = 0; i < snps.size(); i++)
    {
      if (snps[i] != p.snps[i])
        return snps[i] < p.snps[i];
    }
    if (msc == p.msc)
      return false;
    else 
		{
     //cout << "!!!!!!!!!!" << endl; char c; cout << *this << endl << p << endl << "wciśnij "; cin >> c;
     return msc < p.msc;
		}	
  }
  return snps.size() < p.snps.size();
}
//-----------------------------------------------------------------------------------------------
bool PoolItem::findSns(snp_index_t snp)
{
  std::vector<snp_index_t>::iterator it;
  it = std::find(snps.begin(), snps.end(), snp);
  return (it != snps.end());
}