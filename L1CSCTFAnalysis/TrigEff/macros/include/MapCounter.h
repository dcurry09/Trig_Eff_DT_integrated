#ifndef __MAPCOUNTER_H__
#define __MAPCOUNTER_H__ (1)

#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <cmath>

class MapCounter{

 public:

  typedef unsigned long int          count_t;
  typedef std::string                name_t;
  typedef std::map <name_t, count_t> map_t;
  typedef std::vector<name_t>        vector_t;
  typedef vector_t::iterator         vector_iterator;
  typedef map_t::iterator            map_iterator;
  typedef vector_t::const_iterator   vector_const_iterator;
  typedef map_t::const_iterator      map_const_iterator;

 MapCounter():_map(){
    _map.clear();
  }
  
  // Different ways of incrementing the count at tag "name" ...

  count_t Add(const name_t& name, count_t i=1)
  { if ( _map.find(name) == _map.end() ) _ord.push_back(name);
    return _map[name]+=i; }

  count_t Add(const char *name, count_t i=1)
  { name_t tmpstr(name); return Add(tmpstr,i); }

  MapCounter& operator << (const char *rhs)
    { Add (rhs); return *this; }

  // Access individual counts by tag "name" ...

  count_t operator[](const name_t& name)
  { map_const_iterator iter = _map.find(name); 
    if (iter == _map.end()) return 0; else return iter->second;
  }

  count_t operator[](const char *name)
  { name_t tmpstr(name); return (*this)[tmpstr]; }

  // Dump all counts, see also operator << 

  std::ostream& Print(std::ostream &os=std::cout) const;
  

 private:

  // storage
  map_t     _map;  // counter storage
  vector_t  _ord;  // remember counter order
  
  // helper functions
  void   _countRange( count_t &smallest, count_t &largest ) const;
  size_t _longestStrlen() const;
  name_t _padSpace ( int len ) const;
  
};


// the other way to dump all counts, std::cout << MapCounter ...

std::ostream& operator<<(std::ostream &os, const MapCounter& mc)
{ return mc.Print(os); }


std::ostream& MapCounter::Print(std::ostream &os) const {

  vector_const_iterator iter;

  count_t loVal, hiVal;
  size_t countLong, strLong;

  _countRange( loVal, hiVal );

  countLong = floor( log10( hiVal) );
  strLong   = _longestStrlen();

  os << std::endl;
  os << "List of counted values: " << std::endl;
  os << "----------------------- " << std::endl;

  for (iter = _ord.begin(); iter != _ord.end(); iter++){

    name_t tagTmp = (*iter);
    count_t countTmp = _map.find(*iter)->second;

    os << " -- counter [" 
       << tagTmp << "]"
       << _padSpace( strLong - tagTmp.length() )
       << " == " 
       << _padSpace( countLong - floor(log10(countTmp)) )
       << countTmp 
       << std::endl;

  }

  os << "----------------------- " << std::endl;
  os << std::endl;

  return os;
  
}

void MapCounter::_countRange( MapCounter::count_t &smallest, 
			      MapCounter::count_t &largest ) const {

  map_const_iterator iter = _map.begin();

  smallest = iter -> second;
  largest  = smallest;

  for (; iter!= _map.end(); iter++){

    if ( iter -> second < smallest ) smallest = iter -> second;
    if ( iter -> second > largest  ) largest  = iter -> second;

  }

}

size_t MapCounter::_longestStrlen() const {

  vector_const_iterator iter = _ord.begin();

  size_t retVal = iter -> length();

  for (; iter != _ord.end(); iter++)
    if ( iter -> length() > retVal)
      retVal = iter -> length();

  return retVal;

}

MapCounter::name_t MapCounter::_padSpace( int len ) const {

  name_t retVal;

  for (int i = 0; i < len; i++)
    retVal += " ";

  return retVal;

}

#endif
