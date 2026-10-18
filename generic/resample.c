/* Collection of procefures and functions used
   for downsampling and filtering of 3D matrices

   Jesper James Jensen 1995
   */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "polyr.h"

/* tri-linear resampling of data[size] into a new array of newsize
   points ; the new grid step (in old grid points) is size/newsize */
float *
downsample(float *data, int *size, int *newsize)
{
  float *newdata;
  float step[3];
  int x, y, z, x0, y0, z0, x1, y1, z1;
  int xsize = size[0];
  int xysize = size[0]*size[1];
  float fx, fy, fz, r;
  float c00, c01, c10, c11;

  newdata = (float *)tclsurface_alloc((size_t)newsize[0]*newsize[1]*newsize[2]
				 * sizeof(float));

  step[0] = ((float)size[0])/newsize[0];
  step[1] = ((float)size[1])/newsize[1];
  step[2] = ((float)size[2])/newsize[2];

#define D(i,j,k) data[(i) + (j)*xsize + (k)*xysize]
#define IPOL(a,b,t) ((a)*(1.f-(t)) + (b)*(t))
  for (z=0; z<newsize[2]; z++) {
    r = z*step[2]; z0 = (int)r; fz = r - z0;
    if (z0 > size[2]-1) z0 = size[2]-1;
    z1 = (z0 < size[2]-1) ? z0+1 : z0;
    for (y=0; y<newsize[1]; y++) {
      r = y*step[1]; y0 = (int)r; fy = r - y0;
      if (y0 > size[1]-1) y0 = size[1]-1;
      y1 = (y0 < size[1]-1) ? y0+1 : y0;
      for (x=0; x<newsize[0]; x++) {
	r = x*step[0]; x0 = (int)r; fx = r - x0;
	if (x0 > size[0]-1) x0 = size[0]-1;
	x1 = (x0 < size[0]-1) ? x0+1 : x0;

	c00 = IPOL(D(x0,y0,z0), D(x1,y0,z0), fx);
	c10 = IPOL(D(x0,y1,z0), D(x1,y1,z0), fx);
	c01 = IPOL(D(x0,y0,z1), D(x1,y0,z1), fx);
	c11 = IPOL(D(x0,y1,z1), D(x1,y1,z1), fx);
	newdata[x + y*newsize[0] + z*newsize[0]*newsize[1]] =
	  IPOL(IPOL(c00, c10, fy), IPOL(c01, c11, fy), fz);
      }
    }
  }
#undef D
#undef IPOL

  return newdata;
}


/* mean filters the data array with a boxfilter of size 'boxsize' */
void
boxfilter(float **data, int *size, int *boxsize) {
  float *newdata;
  int i,j,k,x,y,z,count;
  int datap,dataphelp,newdatap;
  int filterx,filtery,filterz;
  int imin,imax,jmin,jmax,kmin,kmax;
  int xsize=size[0], xysize=size[0]*size[1];
  float sum;

#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MIN(a,b) (((a)>(b)) ? (b) : (a))
  /* boxsize must be odd */
  if (((boxsize[0]-1)%2)||((boxsize[1]-1)%2)||((boxsize[2]-1)%2))
    tclsurface_fatal("filter dimensions must be odd");
  filterx=(boxsize[0]-1)/2;
  filtery=(boxsize[1]-1)/2;
  filterz=(boxsize[2]-1)/2;

  newdata = (float *)tclsurface_alloc((size_t)size[0]*size[1]*size[2]
				 * sizeof(float));
  /* the filter loop */
  for(j=0;j<size[1]; j++) {
    if (VERBOSE) {
      fprintf(stderr,"%s: Filtering: %3.0f%%\r",MY_NAME,((float)j+1)/size[1]*100);
      fflush(stderr);
    }
    jmin=MAX(j-filtery,0);
    jmax=MIN(j+filtery+1,size[1]);
    for(k=0;k<size[2]; k++) {
      newdatap=j*xsize+k*xysize;
      kmin=MAX(k-filterz,0);
      kmax=MIN(k+filterz+1,size[2]);
      for(i=0;i<size[0];i++) {
	imin=MAX(i-filterx,0);
	imax=MIN(i+filterx+1,size[0]);
	count=0;
	sum=0.0;
	for(z=kmin; z<kmax; z++) {
	  dataphelp=z*xysize+imin;
	  for(y=jmin; y<jmax; y++) {
	    datap=y*xsize+dataphelp;
	    for(x=imin; x<imax; x++) {
	      sum+=(*data)[datap];
	      datap++;
	      count++;
	    }
	  }
	}
	newdata[newdatap]=sum/count;
	newdatap++;
      }
    }
  }
#undef MAX
#undef MIN
  if (VERBOSE)
    fprintf(stderr,"%s: Filtering done   \n",MY_NAME);
  ckfree((char *)*data);
  (*data)=newdata;
}
