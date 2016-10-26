#include <iostream>
#include <algorithm>
#include "distdata.hpp"
#include "parallel.hpp"


////////////////////////////////////////////////////////////////////////////////
ConstData::ConstData(Parallel *p, const std::string &s, int nl, int nc, int alloc):
  NumLines(nl), NumCols(nc), FirstLineOnPE(0), Data(0),
  distbuffersize(0), distbufferallocsize(0), idistbuffer(0), distbuffer(0), parallel(p)
{
  Name=s;
  NumPE=parallel->cpus();
  myPE=parallel->mycpu();
  FirstLine=0;
  LastLine=localLines=NumLines;
  if(alloc) {
    dataSize=NumCols;
    Data=new double[dataSize];
    globalDataLines=dataLines=1;
  }
}


ConstData::~ConstData()
{
  delete [] FirstLineOnPE;
  delete [] Data;
  delete [] distbuffer;
}


void ConstData::calculate_distribution()
{
  LinesPerPE=NumLines/NumPE;
  if(NumLines%NumPE>0) {
    LinesPerPE++;
    if(LinesPerPE*myPE>=NumLines)
      localLines=0;
    else if(LinesPerPE*(myPE+1)>NumLines)
      localLines=NumLines-(NumPE-1)*LinesPerPE;
    else
      localLines=LinesPerPE;
  }
  else localLines=LinesPerPE;
  FirstLineOnPE=new int[NumPE+1];
  FirstLineOnPE[0]=0;
  parallel->gather(&localLines, FirstLineOnPE+1, 1);
  for(int i=2;i<=NumPE;i++) FirstLineOnPE[i]+=FirstLineOnPE[i-1];
  FirstLine=FirstLineOnPE[myPE];
  LastLine=FirstLineOnPE[myPE+1];
}


double ConstData::get(int line, int col)
{
  return Data[col];
}


double* ConstData::get_line(int line)
{
  return Data;
}
  

double* ConstData::new_buffer()
{
  return 0;
}


double* ConstData::get_line(int line, double* &buf)
{
  // buf must be allocated by call to new_buffer()
  // and removed by call to delete_buffer(buf)
  return buf=Data;
}


void ConstData::init(double d)
{
  for(int i=0;i<dataSize;i++) Data[i]=d;
}


void ConstData::set_data(int line, int col, double d)
{
  Data[col]=d;
}


void ConstData::set_dataline(int line, double* d)
{
  if(line==0) for(int i=0;i<NumCols;i++) Data[i]=d[i];
}


void ConstData::set_dataline_from_one_PE(int line, double* d, int from_pe)
{
  if(line==0) {
    if(myPE==from_pe) for(int i=0;i<NumCols;i++) Data[i]=d[i];
    parallel->broadcast(Data, NumCols, from_pe);
  }
}


void ConstData::distribute(int from, int to, double* d, int from_pe)
{
  int size=(to-from)*NumCols;
  int start=from*NumCols;
  // set local data on PE from_pe
  if(mype()==from_pe) for(int i=0;i<size;i++) Data[start+i]=d[i];
  // broadcast Data to all other PEs
  parallel->broadcast(Data+start, size, from_pe);
}


double* ConstData::allocate_buffer(int size)
{
  if(size>distbufferallocsize) {
    if(distbuffer) delete [] distbuffer;
    distbuffer=new double[size*NumCols];
    distbufferallocsize=size;
  }
  distbuffersize=size;
  return distbuffer;
}


double* ConstData::start_buffered_distribute(int bufsize)
{
  idistbuffer=0;
  return allocate_buffer(bufsize);
}


double* ConstData::buffered_distribute(int n, int from_pe)
{
  if(++idistbuffer==distbuffersize) {
    distribute(n+1-distbuffersize, n+1, distbuffer, from_pe);
    idistbuffer=0;
  }
  return distbuffer+idistbuffer*NumCols;
}


void ConstData::end_buffered_distribute(int from_pe)
{
  if(idistbuffer>0) distribute(NumLines-idistbuffer, NumLines, distbuffer, from_pe);
  idistbuffer=0;
}


double* ConstData::start_buffered_sum(int bufsize)
{
  idistbuffer=0;
  return allocate_buffer(bufsize);
}


double* ConstData::buffered_sum(int n)
{
  if(++idistbuffer==distbuffersize) {
    sum_over_pes(distbuffer, n+1-distbuffersize, n+1);
    idistbuffer=0;
  }
  return distbuffer+idistbuffer;
}


void ConstData::end_buffered_sum()
{
  if(idistbuffer>0) sum_over_pes(distbuffer, NumLines-idistbuffer, NumLines);
  idistbuffer=0;
}


void ConstData::copy(ConstData *from, int col_from, int col_to)
{
  double *from_data=from->address();
  int from_cols=from->cols();
  for(int n=0;n<dataLines;n++) {
    Data[col_to]=from_data[col_from];
    col_to+=NumCols;
    col_from+=from_cols;
  }
}


void ConstData::multiply(ConstData *from, int col_to, int col_from)
{
  for(int n=FirstLine; n<LastLine; n++,col_to+=NumCols)
    Data[col_to]*=from->get(n, col_from);
}


void ConstData::divide(ConstData *from, int col_to, int col_from)
{
  for(int n=FirstLine; n<LastLine; n++,col_to+=NumCols)
    Data[col_to]/=from->get(n, col_from);
}


void ConstData::divide(ConstData *from1, ConstData *from2,
		       int col_to, int col_from1, int col_from2)
{
  for(int n=FirstLine; n<LastLine; n++,col_to+=NumCols)
    Data[col_to]/=from1->get(n, col_from1)*from2->get(n, col_from2);
}


void ConstData::sum_over_pes(double *x, int from, int to)
{
  parallel->sum(x, Data, NumCols);
}


void ConstData::print(std::ostream& out)
{
  parallel->sync();
  if(myPE==0) {
    int i,j;
    double* pdata;
    out<<name()<<": size="<<NumLines<<"  type="<<type()<<std::endl;
    for(i=0;i<dataLines;i++) {
      pdata=get_line(i);
      out<<"i="<<i<<" data: ";
      for(j=0;j<NumCols;j++) out<<' '<<pdata[j];
      out<<std::endl<<std::flush;
    }
  }
  parallel->sync();
}


////////////////////////////////////////////////////////////////////////////////
LocalData::LocalData(Parallel *p, const std::string &s, int nl, int nc, int alloc):
  ConstData(p,s,nl,nc,0)
{
  if(alloc) {
    globalDataLines=dataLines=NumLines;
    dataSize=dataLines*NumCols;
    Data=new double[dataSize];
  }
}


LocalData::~LocalData()
{
}


double LocalData::get(int line, int col)
{
  return Data[line*NumCols+col];
}


double* LocalData::get_line(int line)
{
  return Data+line*NumCols;
}


double* LocalData::get_line(int line, double* &buf)
{
  // buf must be allocated by call to new_buffer()
  // and removed by call to delete_buffer(buf)
  return buf=Data+line*NumCols;
}


void LocalData::set_data(int line, int col, double d)
{
  Data[line*NumCols+col]=d;
}


void LocalData::set_dataline(int line, double* d)
{
  int j=line*NumCols;
  for(int i=0;i<NumCols;i++) Data[i+j]=d[i];
}


void LocalData::set_dataline_from_one_PE(int line, double* d, int from_pe)
{
  int start=line*NumCols;
  if(myPE==from_pe) for(int i=0;i<NumCols;i++) Data[i+start]=d[i];
  parallel->broadcast(Data+start, NumCols, from_pe);
}


void LocalData::sum_over_pes(double *x, int from, int to)
{
  parallel->sum(x, Data+from*NumCols, (to-from)*NumCols);
}


////////////////////////////////////////////////////////////////////////////////
// class for indexed access to the data
IndexedData::IndexedData(Parallel *p, const std::string &s,
			 int nl, int nc, int nd, int alloc):
  ConstData(p,s,nl,nc,0), index(0)
{
  if(alloc) {
    globalDataLines=dataLines=nd;
    dataSize=dataLines*NumCols;
    Data=new double[dataSize];
  }
}


IndexedData::~IndexedData()
{
  if(index && index->deallocate()) {
    if(index->win) delete index->win;
    else delete [] index->data;
    delete index;
  }
}


double IndexedData::get(int line, int col)
{
  return Data[index->data[line]*NumCols+col];
}


double* IndexedData::get_line(int line)
{
  return Data+index->data[line]*NumCols;
}


double* IndexedData::get_line(int line, double* &buf)
{
  // buf must be allocated by call to new_buffer()
  // and removed by call to delete_buffer(buf)
  return buf=Data+index->data[line]*NumCols;
}


void IndexedData::set_data(int idata, int col, double d)
{
  Data[idata*NumCols+col]=d;
}


void IndexedData::set_dataline(int idata, double* d)
{
  int j=idata*NumCols;
  for(int i=0;i<NumCols;i++) Data[i+j]=d[i];
}


void IndexedData::set_dataline_from_one_PE(int line, double* d, int from_pe)
{
  int start=line*NumCols;
  if(myPE==from_pe) for(int i=0;i<NumCols;i++) Data[i+start]=d[i];
  parallel->broadcast(Data+start, NumCols, from_pe);
}


void IndexedData::set_index(ConstDataIndex *ind)
{
  index=ind;
  index->set_used();
}


bool IndexedData::set_index_data(int line, int idata)
{
  index->data[line]=idata;
  return idata>=0 && idata<dataLines;
}


ConstDataIndex* IndexedData::alloc_index()
{
  index=new ConstDataIndex;
  index->data=new int[NumLines];
  index->set_used();
  return index;
}


int IndexedData::get_index(int n)
{
  return index->data[n];
}


int IndexedData::check_index()
{
  for(int n=0;n<NumLines;n++)
    if(index->data[n]<0 || index->data[n]>=dataLines) return n;
  return -1;
}


void IndexedData::set_index_data_from_one_PE(int line, int idata, int from_pe)
{
  if(mype()==from_pe) index->data[line]=idata;
  // broadcast Data to all other PEs
  parallel->broadcast(index->data+line, 1, from_pe);
}


void IndexedData::sum_over_pes(double *x, int from, int to)
{
  parallel->abort("operation sum_over_pes is not allowed for IndexData");
}


void IndexedData::print(std::ostream& out)
{
  parallel->sync();
  if(myPE==0) {
    int i,j;
    double* pdata=Data;
    out<<name()<<": size="<<NumLines<<"  datasize="<<dataLines<<"  type="<<type()<<std::endl;
    for(i=0;i<dataLines;i++) {
      out<<"i="<<i<<" data:";
      for(j=0;j<NumCols;j++) out<<' '<<pdata[j];
      out<<std::endl<<std::flush;
      pdata+=NumCols;
    }
    out<<"index:"<<std::endl;
  }
  printindex(out);
  parallel->sync();
}


void IndexedData::printindex(std::ostream& out)
{
  if(myPE==0) {
    int i,j;
    int indold=-1;
    j=0;
    for(i=0;i<NumLines;i++) {
      if(index->data[i]==indold) j++;
      else {
	if(j>0) {
	  if(j==1) out<<"i="<<i<<" data: "<<indold<<std::endl;
	  else out<<"i="<<i-j<<'-'<<i<<" data: "<<j<<'*'<<indold<<std::endl;
	}
	indold=index->data[i];
	j=1;
      }
    }
    if(j==1) out<<"i="<<i<<" data: "<<indold<<std::endl;
    else out<<"i="<<i-j<<'-'<<i<<" data: "<<j<<'*'<<indold<<std::endl;
  }
}


////////////////////////////////////////////////////////////////////////////////
// class for distributed indexed data
IndexedDistData::IndexedDistData(Parallel *p, const std::string &s, int nl, int nc, int nd, int alloc):
  IndexedData(p,s,nl,nc,nd,1), bufbegin(0), bufend(0), buffer(0)
{
  calculate_distribution();
  buflen=16;
  buffer=new int[buflen];
}


IndexedDistData::~IndexedDistData()
{
  delete [] buffer;
}


double IndexedDistData::get(int line, int col)
{
  if(FirstLine<=line && line<LastLine)
    return Data[index->data[line-FirstLine]*NumCols+col];
  if(bufbegin<=line && line<bufend)
    return Data[buffer[line-bufbegin]*NumCols+col];
  // line is not on PE, get line from PE p.
  get_remote_index(line);
  return Data[buffer[0]*NumCols+col];
}


double* IndexedDistData::get_line(int line)
{
  if(FirstLine<=line && line<LastLine)
    return Data + index->data[line-FirstLine]*NumCols;
  if(bufbegin<=line && line<bufend)
    return Data + buffer[line-bufbegin]*NumCols;
  // line is not on PE, get line from PE p.
  get_remote_index(line);
  return Data + buffer[0]*NumCols;
}


double* IndexedDistData::get_line(int line, double* &buf)
{
  // buf must be allocated by call to new_buffer()
  // and removed by call to delete_buffer(buf)
  return buf=get_line(line);
}


int IndexedDistData::get_index(int line)
{
  if(FirstLine<=line && line<LastLine)
    return index->data[line-FirstLine];
  if(bufbegin<=line && line<bufend)
    return buffer[line-bufbegin];
  // line is not on PE, get line from PE p.
  get_remote_index(line);
  return *buffer;
}


void IndexedDistData::get_remote_index(int line)
{
  // line is not on PE, buffered get line from PE p.
  int p=onPE(line);
  bufbegin=line;
  bufend=std::min(bufbegin+buflen, FirstLineOnPE[p+1]);
  index->win->get(buffer, line-FirstLineOnPE[p], bufend-bufbegin, p);
}


ConstDataIndex* IndexedDistData::alloc_index()
{
  // localLines indexes are stored on each PE
  index=new ConstDataIndex;
  index->win=new MPIWindow(localLines, sizeof(int));
  index->data=static_cast<int*>(index->win->data());
  index->set_used();
  return index;
}


bool IndexedDistData::set_index_data(int line, int idata)
{
  if(FirstLine<=line && line<LastLine)
    index->data[line-FirstLine]=idata;
  return idata>=0 && idata<dataLines;
}


void IndexedDistData::set_index_data_from_one_PE(int line, int idata, int from_pe)
{
  int to_pe=onPE(line);
  if(from_pe==to_pe) { // line is on local PE
    if(myPE==from_pe) index->data[line-FirstLine]=idata;
  }
  else { // line is not on PE, send line to PE.
    if(myPE==from_pe) // send data from PE from_pe to PE to_pe
      parallel->send(to_pe, &idata, 1);
    else if(myPE==to_pe) // receive data on PE to_pe from PE from_pe
      parallel->receive(from_pe, index->data+(line-FirstLine), 1);
  }
}


void IndexedDistData::printindex(std::ostream& out)
{
  int i,j, ind,indold, p;
  for(p=0;p<NumPE;p++) {
    if(myPE==p) {
      j=0;
      indold=-1;
      for(i=FirstLine;i<LastLine;i++) {
	ind=index->data[i-FirstLine];
	if(ind==indold) j++;
	else {
	  if(j>0) {
	    if(j==1) out<<"i="<<i<<" data: "<<indold<<std::endl;
	    else out<<"i="<<i-j<<'-'<<i<<" data: "<<j<<'*'<<indold<<std::endl;
	  }
	  indold=ind;
	  j=1;
	}
      }
      if(j==1) out<<"i="<<i<<" data: "<<indold<<std::endl;
      else out<<"i="<<i-j<<'-'<<i<<" data: "<<j<<'*'<<indold<<std::endl;
    }
    parallel->sync();
  }
}


////////////////////////////////////////////////////////////////////////////////
DistData::DistData(Parallel *p, const std::string &s, int nl, int nc, int alloc):
  ConstData(p,s,nl,nc,0), bufferline(-1), buffer(0)
{
  calculate_distribution();
  buffer=new double[NumCols];
  dataLines=localLines;
  globalDataLines=NumLines;
  dataSize=dataLines*NumCols;
  win=new MPIWindow(dataSize, sizeof(double));
  Data=static_cast<double*>(win->data());
}


DistData::~DistData()
{
  delete [] buffer;
  delete win;
  Data=0; // because of the destructor in the base class
}


double DistData::get(int line, int col)
{
  if(FirstLine<=line && line<LastLine)
    return Data[(line-FirstLine)*NumCols+col];
  if(line!=bufferline) {
    // line is not buffer, read line from other PE
    int p=onPE(line);
    win->get(buffer, (line-FirstLineOnPE[p])*NumCols, NumCols, p);
    bufferline=line;
  }
  return buffer[col];
}


double* DistData::get_line(int line)
{
  // Caution: this subroutine returns a pointer to the class variable buffer.
  if(FirstLine<=line && line<LastLine)
    return Data+(line-FirstLine)*NumCols;
  if(line!=bufferline) {
    // line is not buffered, read line from PE p.
    int p=onPE(line);
    win->get(buffer, (line-FirstLineOnPE[p])*NumCols, NumCols, p);
    bufferline=line;
  }
  return buffer;
}


double* DistData::new_buffer()
{
  return new double[NumCols];
}
  

void DistData::delete_buffer(double *buf)
{
  delete [] buf;
}


double* DistData::get_line(int line, double* &buf)
{
  // buf must be allocated by call to new_buffer()
  // and removed by call to delete_buffer(buf)
  if(FirstLine<=line && line<LastLine) {
    int i=(line-FirstLine)*NumCols;
    for(int n=0;n<NumCols;n++,i++) buf[n]=Data[i];
  }
  else {
    // line is not buffered, read line from PE p.
    int p=onPE(line);
    win->get(buf, (line-FirstLineOnPE[p])*NumCols, NumCols, p);
  }
  return buf;
}


void DistData::set_data(int line, int col, double d)
{
  if(FirstLine<=line && line<LastLine)
    Data[(line-FirstLine)*NumCols+col]=d;
}


void DistData::set_dataline(int line, double* d)
{
  if(FirstLine<=line && line<LastLine) {
    int j=(line-FirstLine)*NumCols;
    for(int i=0;i<NumCols;i++) Data[i+j]=d[i];
  }
}


void DistData::set_dataline_from_one_PE(int line, double* d, int from_pe)
{
  int to_pe=onPE(line);
  if(from_pe==to_pe) {
    if(myPE==from_pe) {
      // data is on local PE
      int j=(line-FirstLine)*NumCols;
      for(int i=0;i<NumCols;i++) Data[i+j]=d[i];
    }
  }
  else { // send data to remote PE
    if(myPE==from_pe) // send data from PE from_pe to PE to_pe
      parallel->send(to_pe, d, NumCols);
    else if(myPE==to_pe) // receive data on PE to_pe from PE from_pe
      parallel->receive(from_pe, Data+(line-FirstLine)*NumCols, NumCols);
  }
}


void DistData::fence(int assert)
{
  win->fence(assert);
}


void DistData::distribute(int from, int to, double* d, int from_pe)
{
  int i;
  int n1=from;
  int n2=to;
  int to_pe;
  int datapos, size;
  int dpos=0;
  while(n1<n2) {
    to_pe=onPE(n1);
    i=FirstLineOnPE[to_pe+1];
    if(i<n2) n2=i;
    datapos=(n1-FirstLineOnPE[to_pe])*NumCols;
    size=(n2-n1)*NumCols;
    if(from_pe==to_pe) {
      if(mype()==from_pe) // set local data on PE from_pe
	for(i=0;i<size;i++) Data[datapos+i]=d[dpos+i];
    }
    else {
      if(mype()==from_pe) // send data from PE from_pe to PE to_pe
	parallel->send(to_pe, d+dpos, size);
      else if(mype()==to_pe) // receive data on PE to_pe from PE from_pe
	parallel->receive(from_pe, Data+datapos, size);
    }
    // set new from and to values
    n1=n2; n2=to;
    dpos+=size;
  }
}


void DistData::sum_over_pes(double *x, int from, int to)
{
  int n1=from;
  int n2=to;
  int to_pe;
  int i, no;
  int pos=0;
  while(n1<n2) {
    to_pe=onPE(n1);
    i=FirstLineOnPE[to_pe+1];
    if(i<n2) n2=i;
    no=(n2-n1)*NumCols;
    parallel->sum(x+pos, Data+(n1-FirstLine)*NumCols, no, to_pe);
    n1=n2; n2=to;
    pos+=no;
  }
}


void DistData::print(std::ostream& out)
{
  int i,j;
  double* pdata;
  const char nl='\n';
  parallel->sync();
  if(myPE==0) {
    out<<name()<<": size="<<NumLines<<"  type="<<type()<<std::endl;
    for(i=0;i<NumPE;i++) out<<"  PE="<<i<<" start="<<FirstLineOnPE[i];
    out<<std::endl;
  }
  for(int p=0;p<NumPE;p++) {
    parallel->sync();
    if(p==myPE) {
      for(i=begin();i<end();i++) {
	pdata=get_line(i);
	out<<"i="<<i<<" data: ";
	for(j=0;j<NumCols;j++) out<<' '<<pdata[j];
	out<<nl;
      }
      out.flush();
    }
  }
  parallel->sync();
}


std::ostream& operator<<(std::ostream& out, ConstData &x)
{
  x.print(out);
  return out;
}


std::ostream& operator<<(std::ostream& out, ConstData *x)
{
  x->print(out);
  return out;
}
