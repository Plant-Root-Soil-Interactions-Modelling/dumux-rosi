#include <ctime>
#include <cmath>
#include <limits>
#include "random.hpp"


RandomClass::RandomClass(): m(2147483647e0), a(16807e0), fast_iset(1)
{ 
  time_t t;
  set_seed((double) time(&t)); 
  sgenrand((unsigned long) time(&t));
}

RandomClass::RandomClass(double s): m(2147483647e0), a(16807e0), fast_iset(1)
{ 
  set_seed(s);
  unsigned long i=(unsigned long) s;
  if(i==0) i=1234567890;
  sgenrand(i);
}

void RandomClass::set_seed(double s) 
{
  if(s<0.0) s=-s;
  // s mod m berechnen
  seed=fmod(s,m);
}

void RandomClass::sgenrand(unsigned long seed_l) /* seed should not be 0 */
{
  infinity=std::numeric_limits<double>::infinity();
  M=397; N=624;
  MATRIX_A=0x9908b0df;
  UPPER_MASK=0x80000000;
  LOWER_MASK=0x7fffffff;
  MAX32_MASK=0xffffffff;
  TEMPERING_MASK_B=0x9d2c5680;
  TEMPERING_MASK_C=0xefc60000;
  k=1;     
  max32=MAX32_MASK;
  // mag01[x] = x * MATRIX_A  for x=0,1 
  mag01[0]=0x0; mag01[1]=MATRIX_A;
  /* setting initial seeds to ptgfsr[N] using     */
  /* the generator Line 25 of Table 1 in          */
  /* [KNUTH 1981, The Art of Computer Programming */
  /*    Vol. 2 (2nd Ed.), pp102]                  */
  ptgfsr[0]= seed_l & MAX32_MASK;
  for (int i=1; i<N; i++) ptgfsr[i] = (69069 * ptgfsr[i-1]) & MAX32_MASK;
}

double RandomClass::fastdraw() 
{ 
  // Gleichverteilte ZV im Breich 0 < zv <1
  // Faktor 2-3 mal schneller all draw()
  if(k == N){ /* generate N words at one time */
    for (kk=0;kk<N-M;kk++) {
      y = (ptgfsr[kk]&UPPER_MASK)|(ptgfsr[kk+1]&LOWER_MASK);
      ptgfsr[kk] = ptgfsr[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    for (;kk<N-1;kk++) {
      y = (ptgfsr[kk]&UPPER_MASK)|(ptgfsr[kk+1]&LOWER_MASK);
      ptgfsr[kk] = ptgfsr[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    y = (ptgfsr[N-1]&UPPER_MASK)|(ptgfsr[0]&LOWER_MASK);
    ptgfsr[N-1] = ptgfsr[M-1] ^ (y >> 1) ^ mag01[y & 0x1];
    k = 0;
  }
  y = ptgfsr[k++];
  y ^= (y >> 11);
  y ^= (y << 7) & TEMPERING_MASK_B;
  y ^= (y << 15) & TEMPERING_MASK_C;
  y &= MAX32_MASK; /* you may delete this line if word size = 32 */
  y ^= (y >> 18);
  return y/max32;
}

double RandomClass::fastgauss() { 
  // Gauss-Verteilung, Mittelwert=0, Standardabweichung=1
  if(fast_iset) {
    while(1) {
      fast_v1=2*fastdraw()-1;
      fast_v2=2*fastdraw()-1;
      fast_rsq=fast_v1*fast_v1+fast_v2*fast_v2;
      if(fast_rsq>=1.0 || fast_rsq==0.0) continue;
      if( (fast_fac=sqrt(-2.0*log(fast_rsq)/fast_rsq)) != infinity )
	break;
    }
    fast_iset=0;
    return fast_v2*fast_fac;
  }
  else {
    fast_iset=1;
    return fast_v1*fast_fac;
  }
}
