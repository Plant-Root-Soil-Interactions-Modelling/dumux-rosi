#ifndef DISTDATA_HPP
#define DISTDATA_HPP

#include <string>
#include <iostream>

class MPIWindow;
class Parallel;


/** Class ConstDataIndex.
  * This class is used for indexed access to the data. The index maintained
  * by this class can be used by more than one DistData class. The index is
  * deleted, if it is no longer used by any DistData class.
  */
class ConstDataIndex
{
public:
  int used;        ///< reference counter.
  int* data;       ///< pointer to the index data.
  MPIWindow *win;  ///< Object containing the MPI window used for the index.
  ConstDataIndex(): used(0), data(0), win(0) {}
  virtual ~ConstDataIndex() {}
  inline void set_used() { ++used; }
  inline bool deallocate() { return --used==0; } ///< Index is deleted from outside.
};


/** Class ConstantData.
  * All the data is constant.
  */
class ConstData 
{
protected:
  std::string Name;   ///< Name for the data.
  int NumPE;          ///< Number of processors.
  int myPE;           ///< Own processor number.
  // index related variables
  int NumLines;       ///< Number of lines (distributed over PEs).
  int NumCols;        ///< Numbers of data items per line.
  // index related variables, if distributed over PEs.
  int FirstLine;      ///< Number of the first line on PE.
  int LastLine;       ///< Number of the last line on PE.
  int localLines;     ///< Number of lines on local PE.
  int *FirstLineOnPE; ///< Number of the first data line for each PE.
  int LinesPerPE;     ///< estimated number of lines on each PE.
  // data related variables
  int globalDataLines;///< number of global data lines.
  int dataLines;      ///< number of local data lines.
  int dataSize;       ///< number of local data values.
  double* Data;       ///< pointer to the data.
  // buffer variables for buffered distribution.
  int distbuffersize;     ///< length of the temporary data buffer.
  int distbufferallocsize;///< allocated length of the temporary data buffer.
  int idistbuffer;        ///< current index of the temporary data buffer.
  double* distbuffer;     ///< pointer to the temporary data.

public:
  Parallel *parallel; ///< Pointer to the parallel class.
  /// constructor.
  ConstData(Parallel *parallel, const std::string &s, int nl, int nc, int alloc=1);
  virtual ~ConstData();                            ///< destructor
  virtual bool constant()    const { return true; }
  virtual bool indexed()     const { return false; }
  virtual bool distributed() const { return false; }
  inline std::string name() const { return Name; }
  inline int lines() const { return NumLines; }
  inline int cols() const { return NumCols; }
  inline int local_lines() const { return localLines; }
  inline int local_datalines() const { return dataLines; }
  inline int global_datalines() const { return globalDataLines; }
  inline int numpes() const { return NumPE; }
  inline int mype() const { return myPE; }
  inline double* address() const { return Data; }
  void init(double d);
  void calculate_distribution();
  double* allocate_buffer(int size=1024);
  double* start_buffered_distribute(int bufsize=1024);
  double* buffered_distribute(int n, int from_pe=0);
  void    end_buffered_distribute(int from_pe=0);
  double* start_buffered_sum(int bufsize=1024);
  double* buffered_sum(int n);
  void    end_buffered_sum();
  void copy(ConstData *from, int col_from=0, int col_to=0);
  void multiply(ConstData *from, int col_to=0, int col_from=0);
  void divide(ConstData *from, int col_to=0, int col_from=0);
  void divide(ConstData *from1, ConstData *from2, int col_to=0, int col_from1=0, int col_from2=0);
  virtual std::string type() const { return "constant"; }
  virtual double get(int line, int col=0);
  virtual double* get_line(int line);               ///< get line with internal buffer.
  virtual double* new_buffer();
  virtual void delete_buffer(double *buf) {}
  virtual double* get_line(int line, double* &buf); ///< get line with external buffer.
  virtual void set_data(int line, int col, double d);
  virtual void set_dataline(int line, double* d);
  virtual int onPE(int n) const { return myPE; }
  virtual int begin() const { return 0; }
  virtual int end() const { return dataLines; }
  virtual void set_dataline_from_one_PE(int line, double* d, int from_pe=0);
  virtual void distribute(int from, int to, double* d, int from_pe=0);
  /// empty methods, overloaded in the IndexData class.
  virtual bool set_index_data(int line, int idata) { return 0; }
  virtual int check_index() { return -1; }
  virtual ConstDataIndex* alloc_index() { return 0; }
  virtual int get_index(int n) { return n; }
  virtual void set_index(ConstDataIndex *ind) {}
  virtual ConstDataIndex* get_index() { return 0; }
  virtual void set_index_data_from_one_PE(int line, int idata, int from_pe=0) {}
  virtual void sum_over_pes(double *x, int from, int to);
  virtual void print(std::ostream& out);
};


/** Class LocalData.
  * All the data is completely on each prozessor.
  */
class LocalData: public ConstData
{
public:
  /// constructor
  LocalData(Parallel *parallel, const std::string &s, int nl, int nc, int alloc=1);
  virtual ~LocalData();                                         ///< destructor
  virtual std::string type() const { return "local"; }
  virtual bool constant() const { return false; }
  virtual double get(int line, int col=0);
  virtual double* get_line(int line);               ///< get line with internal buffer.
  virtual double* get_line(int line, double* &buf); ///< get line with external buffer.
  virtual void set_data(int line, int col, double d);
  virtual void set_dataline(int line, double* d);
  virtual void set_dataline_from_one_PE(int line, double* d, int from_pe=0);
  virtual void sum_over_pes(double *x, int from, int to);
};


/** Class IndexedData.
  * All the data is completely on each prozessor. The access to
  * the data is indexed. The index is also completely on each
  * processor.
  */
class IndexedData: public ConstData
{
protected:
  ConstDataIndex* index; ///< Pointer to the index of the data.
public:
  /// constructor
  IndexedData(Parallel *parallel, const std::string &s, int nl, int nc, int nd, int alloc=1);
  virtual ~IndexedData();  ///< destructor
  virtual std::string type() const { return "indexed"; }
  virtual bool constant() const { return false; }
  virtual bool indexed()  const { return true; }
  virtual double get(int line, int col=0);
  virtual double* get_line(int line);               ///< get line with internal buffer.
  virtual double* get_line(int line, double* &buf); ///< get line with external buffer.
  virtual void set_data(int idata, int col, double d);
  virtual void set_dataline(int idata, double* d);
  virtual void set_dataline_from_one_PE(int line, double* d, int from_pe=0);
  virtual bool set_index_data(int line, int idata);
  virtual int check_index();
  virtual ConstDataIndex* alloc_index();
  virtual int get_index(int n);
  virtual void set_index(ConstDataIndex *ind);
  virtual ConstDataIndex* get_index() { return index; }
  virtual void set_index_data_from_one_PE(int line, int idata, int from_pe=0);
  virtual void sum_over_pes(double *x, int from, int to);
  virtual void print(std::ostream& out);
  virtual void printindex(std::ostream& out);
};




/** Class IndexedDistData.
  * All the data is completely on each prozessor. The access to
  * the data is indexed. The index is distributed over the processors.
  */
class IndexedDistData: public IndexedData
{
protected:
  int bufbegin;      ///< start number of the last read and not local index.
  int bufend;        ///< last number of the last read and not local index.
  int *buffer;       ///< last read and not local data index.
  int buflen;        ///< length of the buffer.
public:
  /// constructor
  IndexedDistData(Parallel *parallel, const std::string &s, int nl, int nc, int nd, int alloc=1);
  virtual ~IndexedDistData();  ///< destructor
  virtual std::string type() const { return "distributed indexed"; }
  virtual bool constant() const { return false; }
  virtual bool distributed() const { return true; }
  virtual double get(int line, int col=0);
  virtual double* get_line(int line);               ///< get line with internal buffer.
  virtual double* get_line(int line, double* &buf); ///< get line with external buffer.
  virtual int get_index(int n);
  void get_remote_index(int line);
  virtual ConstDataIndex* alloc_index();
  virtual bool set_index_data(int line, int idata);
  virtual int begin() const { return FirstLine; }
  virtual int end() const { return LastLine; }
  virtual int onPE(int n) const { return n/LinesPerPE; }
  virtual void set_index_data_from_one_PE(int line, int idata, int from_pe=0);
  virtual void printindex(std::ostream& out);
};


/** class DistData.
  * The data is distributed over the prozessors.
  * 
  */
class DistData: public ConstData
{
protected:
  int bufferline;    ///< Number of the line in the buffer.
  double *buffer;    ///< Last line received from an other PE.
  MPIWindow *win;    ///< Object containing the MPI window.

public:
  /// constructor
  DistData(Parallel *parallel, const std::string &s, int nl, int nc, int alloc=1);
  virtual ~DistData();  ///< destructor
  virtual std::string type() const { return "distributed"; }
  virtual bool constant() const { return false; }
  virtual bool distributed() const { return true; }
  virtual double get(int line, int col=0);
  virtual double* get_line(int line);               ///< get line with internal buffer.
  virtual double* new_buffer();
  virtual void delete_buffer(double *buf);
  virtual double* get_line(int line, double* &buf); ///< get line with external buffer.
  virtual void set_data(int line, int col, double d);
  virtual void set_dataline(int line, double* d);
  virtual int onPE(int n) const { return n/LinesPerPE; }
  virtual int begin() const { return FirstLine; }
  virtual int end() const { return LastLine; }
  void fence(int assert=0);
  /**
   * Distribute data from line \a from to line \a to.
   * Processor \a from_pe send the data defined by the pointer \a d
   * to the corresponding processor(s).
   * Which line belongs to which processor is defined by the
   * function calculate_distribution().
   * \param from = first line number of the data.
   * \param   to = last line number of the data (exclusively).
   * \param    d = pointer to the data.
   * \param from_pe = number of the sending processor.
   */
  virtual void set_dataline_from_one_PE(int line, double* d, int from_pe=0);
  virtual void distribute(int from, int to, double* d, int from_pe=0);
  virtual void sum_over_pes(double *x, int from, int to);
  inline void reset_buffer() { bufferline=-1; }
  virtual void print(std::ostream& out);
};

std::ostream& operator<<(std::ostream& out, ConstData &i);
std::ostream& operator<<(std::ostream& out, ConstData *i);
#endif
