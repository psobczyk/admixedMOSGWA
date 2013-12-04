/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2012–2013, Artur Gola, Erich Dolejsi, Bernhard Bodenstorfer.	*
 *										*
 *	This program is free software; you can redistribute it and/or modify	*
 *	it under the terms of the GNU General Public License as published by	*
 *	the Free Software Foundation; either version 3 of the License, or	*
 *	(at your option) any later version.					*
 *										*
 *	This program is distributed in the hope that it will be useful,		*
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of		*
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.			*
 *	See the GNU General Public License for more details.			*
 ********************************************************************************/

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
