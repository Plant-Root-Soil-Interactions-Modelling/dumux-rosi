#include <cassert>
#include <cmath>
#include <algorithm>
#include <limits>
#include "hash.hpp"


HashClassList::~HashClassList()
{
  clear();
}


void HashClassList::insert(HashClassData *d)
{
  // insert item at the top of the list.
  HashClassListItem *i=new HashClassListItem(d);
  i->next=top;
  top=i;
}


void HashClassList::clear()
{
  HashClassListItem *i=top, *j;
  while(i) {
    j=i;
    i=i->next;
    delete j->data;
    delete j;
  }
  top=0;
}


HashClass::HashClass(): rows(1223), cols(10)
{
   initHash();
}


HashClass::HashClass(int isize): cols(10)
{
  if(isize<13) rows=13;
  else rows=nextPrime(isize);
  initHash();
}


HashClass::HashClass(int isize, int collisions): cols(collisions)
{
  if(isize<13) rows=13;
  else rows=nextPrime(isize);
  initHash();
}


HashClass::~HashClass()
{
  clear();
  // delete the hash arrays.
  delete [] key;
  delete [] colptr;
  delete [] data;
}


void HashClass::clear()
{
  // delete the data.
  for(int i=0; i<tablesize;i++) {
    if(key[i]!=empty) {
      delete data[i];
      key[i]=empty;
    }
  }
}


void HashClass::initHash()
{
  key=0; colptr=0; data=0;
  current=0;
  tablesize=rows*cols;
  key=new int[tablesize];
  assert(key);
  colptr=new int[rows];
  assert(colptr);
  data=new HashClassData*[tablesize];
  assert(data);
  empty=std::numeric_limits<int>::min();
  init(); // init hash
}


int HashClass::nextPrime(int n) {
  int i, nw;
  const int not_prime=1;
  int nicht_teilbar;
  if(n%2 == 0) n++;
  do {
    nicht_teilbar=1;
    nw=(int) sqrt((double)n);
    for(i=3;i<=nw;i+=2) {
      if( (n%i) == 0 ) {
	nicht_teilbar=0;
	break;
      }
    }
    if(nicht_teilbar) break;
    n+=2;
  } while(not_prime);
  return n;
}


void HashClass::insert(int k, HashClassData *d)
{
  int row=k%rows;
  int index=row*cols+colptr[row];
  // save old data pointer 
  if(key[index]!=empty) deleted_items.insert(data[index]);
  // save pointer to new data
  data[index]=d;
  d->key=key[index]=k;
  if(++colptr[row]==cols) colptr[row]=0;
}


HashClassData* HashClass::find(int k)
{
  // get for one value
  int row=k%rows;
  int c=colptr[row];
  int index=row*cols;
  // search backward for key, starting at position colptr-1
  for(int i=0;i<cols;i++) {
    c=prev_col(c);
    if(key[index+c]==empty) return 0;
    if(key[index+c]==k) {
      data[index+c]->used=1;
      return data[index+c];
    }
  }
  return 0;
}


void HashClass::init()
{
  std::fill(colptr, colptr+cols, 0);
  std::fill(key, key+tablesize, empty);
}


void HashClass::reset_used()
{
  int i,j,c;
  int index=0;
  for(i=0;i<rows;i++) {
    c=colptr[i];
    for(j=0;j<cols;j++) {
      c=prev_col(c);
      if(key[index+c]==empty) break;
      data[index+c]->used=0;
    }
    index+=cols;
  }
}


void HashClass::erase_unused()
{
  int i,j,c,n;
  int index=0;
  for(i=0;i<rows;i++) {
    c=colptr[i];
    n=prev_col(c);
    for(j=0;j<cols;j++) {
      c=prev_col(c);
      if(key[index+c]==empty) break;
      if(data[index+c]->used) {
	// data still used, possibly shift it.
	if(n!=c) {
	  data[index+n]=data[index+c];
	  key[index+n]=key[index+c];
	  key[index+c]=empty;
	}
	n=prev_col(n);
      }
      else {
	// data no more used.
	delete data[index+c];
	key[index+c]=empty;
      }
    }
    index+=cols;
  }
  // delete unused data.
  deleted_items.clear();
}


HashClassData* HashClass::begin()
{
  for(current=0; current<tablesize; current++)
    if(key[current]!=empty) return data[current];
  return 0;
}


HashClassData* HashClass::next()
{
  while(++current<tablesize)
    if(key[current]!=empty) return data[current];
  current=0;
  return 0;
}
