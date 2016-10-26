#ifndef HASH_CLASS_H
#define HASH_CLASS_H


class HashClassData {
public:
  int used;
  int key;
  HashClassData(): used(1), key(-1) {} 
  virtual ~HashClassData() {}
};


class HashClassListItem {
public:
  HashClassData *data;
  HashClassListItem *next;
  HashClassListItem(HashClassData *d): data(d), next(0) {} 
  virtual ~HashClassListItem() {}
};


class HashClassList {
public:
  HashClassListItem *top;
  HashClassList(): top(0) {}
  virtual ~HashClassList();
  void insert(HashClassData *d);
  void clear();
};


class HashClass {

protected:
  int rows;             // no. of rows in the hash
  int cols;             // no. of collisions
  int tablesize;        // total size of the hash
  int* key;             // [rows*cols] vector of all keys
  HashClassData** data; // [rows*cols] vector of data objects
  int *colptr;          // [rows] pointer to column in row
  int empty;            // value for free position.
  int current;          // current element in hash.
  HashClassList deleted_items; // list of deleted items.

  void initHash();
  int nextPrime(int n);
  inline int prev_col(int c) const { return (--c<0) ? cols-1 : c; }

public:
  HashClass();
  HashClass(int size);
  HashClass(int size, int collisions);
  virtual ~HashClass();

  inline int size() const { return tablesize; }
  inline int printrow() const { return current; }
  inline int printkey() const { return key[current]; }

  // routines for one value, do not use for dim>1
  void insert(int key, HashClassData *data);
  HashClassData* find(int key);
  void init();
  void clear();
  void reset_used();
  void erase_unused();
  HashClassData* begin();
  HashClassData* next();
};

#endif


