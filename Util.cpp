typedef size_t ov_size;

#include <algorithm>
#include <iterator>

#include"os_predef.h"
#include"os_std.h"
#include"os_time.h"

#include"Util.h"
#include"MemoryDebug.h"
#include"Err.h"

void *UtilArrayCalloc(unsigned int *dim,ov_size ndim,ov_size atom_size)
{
  ov_size size;
  ov_size sum,product;
  ov_size chunk;
  ov_size a,b,c;
  void *result;
  char **p;
  char *q;
  
  sum = 0;
  for(a=0;a<(ndim-1);a++) {
    product = dim[0];
    for(b=1;b<=a;b++)
      product = product * dim[b];
    sum = sum + product * sizeof(void*);
  }
  size = atom_size;
  for(a=0;a<ndim;a++)
	 size = size * dim[a];
  size = size + sum;
  result = pymol::calloc<char>(size);

  if(result) {
    chunk = 1;
    p = (char**) result;
    for(c=0;c<(ndim-1);c++) {
      if(c<(ndim-2)) {
        chunk = dim[c+1] * sizeof(void*);
      } else {
        chunk = dim[c+1] * atom_size;
      }
      
      product = dim[0];
      for(b=1;b<=c;b++)
        product = product * dim[b];
      q = ((char*)p) + product * sizeof(void*); 
      for(a=0;a<product;a++) {
        *p = q;
        p++;
        q+=chunk;
      }
    }
  }
  return(result);
}
