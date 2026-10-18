/*
 * This program implements the marching cubes surface tiler described by
 * Lorensen & Cline in the Siggraph 87 Conference Proceedings. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "polyr.h"

typedef struct {
  int nverts;
  int verts[8];
  int nedges;
  int edges[12];
  int npolys;
  int polys[30];
} CELL_ENTRY;

#include "cell_table.h"

static int XDIMYDIM;

static void calc_index(int *index, float *data, int y1, int z1, int xdim,
		       float thresh);
static void get_cell_verts(int *index, float *data, int y1, int z1, int xdim,
			   float xtrans, float ytrans, float ztrans,
			   float threshold, float *crossings);
static int get_cell_polys(int *index, int xdim, float *crossings);
static void calc_cell_normals(int *index, float *data, int y1, int z1,
			      int xdim, int ydim, int zdim, float *norms);
static void get_cell_verts_and_norms(int *index, float *data, int y1, int z1,
				     int xdim, float xtrans, float ytrans,
				     float ztrans, float thresh,
				     float *crossings, float *cnorm,
				     float *normals);
static int get_cell_polys_and_norms(int *index, int xdim, float *crossings,
				    float *normals);

int 
iso_surface(float *data, int xdim, int ydim, int zdim, float threshold)
     /* This subroutine will generate the polygonal isosurface. */
{
  
  int y, z;
  int xdim1, ydim1, zdim1;
  float xtrans, ytrans, ztrans;	/* used for centering the volume */
  int *index;
  int npolys;
  float *crossings;
  float *cell_norms;
  float *normals;
  
  if (xdim < 3 || ydim < 3 || zdim < 3)
    tclsurface_fatal("grid [%d,%d,%d] too small for isosurfacing", xdim, ydim, zdim);

  zdim1 = zdim - 1;
  ydim1 = ydim - 1;
  xdim1 = xdim - 1;

  xtrans = -((float) xdim / 2.0);
  ytrans = -((float) ydim / 2.0);
  ztrans = -((float) zdim / 2.0);

  XDIMYDIM = xdim * ydim;

  npolys = 0; /* keep count of total polygons */

  index = (int *)tclsurface_alloc(xdim1*sizeof(int));
  crossings = (float *)tclsurface_alloc(xdim1*13*3*sizeof(float));
  if (NORMAL_TYPE) {
    /* gradient normals - scary stuff kids */
    normals = (float *)tclsurface_alloc(xdim1*13*3*sizeof(float));
    cell_norms = (float *)tclsurface_alloc(xdim1*8*3*sizeof(float));
    for (z = 0; z < zdim1; z++)	/* process each cell in the volume */
      for (y = 0; y < ydim1; y++) {
	calc_index(index, data, y, z, xdim, threshold);
	calc_cell_normals(index, data, y, z, xdim, ydim, zdim, cell_norms);
	get_cell_verts_and_norms(index, 
				 data, 
				 y, 
				 z, 
				 xdim, 
				 xtrans,
				 ytrans, 
				 ztrans, 
				 threshold, 
				 crossings,
				 cell_norms, 
				 normals);
	npolys += get_cell_polys_and_norms(index, 
					   xdim, 
					   crossings, 
					   normals);
      }
    ckfree((char *)normals);
    ckfree((char *)cell_norms);
  } else {
    for (z = 0; z < zdim1; z++)	/* process each cell in the volume */
      for (y = 0; y < ydim1; y++) {
	calc_index(index, data, y, z, xdim, threshold);
	get_cell_verts(index, 
		       data, 
		       y, 
		       z, 
		       xdim, 
		       xtrans, 
		       ytrans, 
		       ztrans,
		       threshold, 
		       crossings);
	npolys += get_cell_polys(index, 
				 xdim, 
				 crossings);
      }
  }
  ckfree((char *)crossings);
  ckfree((char *)index);

  if (VERBOSE)
    fprintf(stderr, "%s: %d polygons\n", MY_NAME, npolys);
  return 0;
}


static void
calc_index(int *index, float *data, int y1, int z1, int xdim, float thresh)
     /* This subroutine calculates the index and creates some global */
     /* temporary variables (for speed). */
{
  float *tmp;
  float threshold = thresh;
  int x1;
  unsigned i = 0;

  /* first compute index of first cube */

  tmp = data + (z1 * XDIMYDIM) + (y1 * xdim) + /* x1= */ 0;

  i += (threshold <= tmp[0]);
  i += (threshold <= tmp[1]) * 2;

  tmp += xdim;
  i += (threshold <= tmp[1]) * 4;
  i += (threshold <= tmp[0]) * 8;

  tmp = tmp - xdim + XDIMYDIM;
  i += (threshold <= tmp[0]) * 16;
  i += (threshold <= tmp[1]) * 32;

  tmp += xdim;
  i += (threshold <= tmp[1]) * 64;
  i += (threshold <= tmp[0]) * 128;

  index[0] = i;

  /* now compute rest */

  tmp -= xdim + XDIMYDIM;
  for (x1 = 1; x1 < xdim-1; x1++) {

    ++tmp;

    i = ((i&0x44)<<1) | ((i&0x22)>>1);	/* resuse 4 of the bits */

    i += (threshold <= tmp[1]) * 2;
    i += (threshold <= tmp[xdim+1]) * 4;
    i += (threshold <= tmp[XDIMYDIM+1]) * 32;
    i += (threshold <= tmp[XDIMYDIM+xdim+1]) * 64;

    index[x1] = i;
  }
}


static void 
get_cell_verts(int *index, float *data, int y1, int z1, int xdim,
	       float xtrans, float ytrans, float ztrans, float threshold,
	       float *crossings)
{

  int x1, y2, z2;

#define CROSSINGS(x,a,b) crossings[x*13*3+a*3+b]
#define linterp(a1,a2,a,b1,b2) ((float)(((a-a1) * (float)(b2-b1) / (a2-a1)) + (float)b1))

  y2 = y1 + 1;
  z2 = z1 + 1;
  for (x1 = 0; x1 < xdim-1; x1++) {
    float cx, cy, cz;
    int nedges;
    int crnt_edge;
    int x2 = x1 + 1;
    int i;
    float *v1, *v4, *v5, *v8;


    if (!index[x1]) continue;

#ifdef FOO
    v1 = data + z1*XDIMYDIM + y1*xdim + x1;
    /*v2 = data + z1*XDIMYDIM + y1*xdim + x1+1;*/
    /*v3 = data + z1*XDIMYDIM + y1*xdim+xdim + x1+1;*/
    v4 = data + z1*XDIMYDIM + y1*xdim+xdim + x1;
    v5 = data + z1*XDIMYDIM+XDIMYDIM + y1*xdim + x1;
    /*v6 = data + z1*XDIMYDIM+XDIMYDIM + y1*xdim + x1+1;*/
    /*v7 = data + z1*XDIMYDIM+XDIMYDIM + y1*xdim+xdim + x1+1;*/
    v8 = data + z1*XDIMYDIM+XDIMYDIM + y1*xdim+xdim + x1;
#else
    v1 = data + z1*XDIMYDIM + y1*xdim + x1;
    v4 = v1 + xdim;
    v5 = v1 + XDIMYDIM;
    v8 = v4 + XDIMYDIM;
#endif

    nedges = cell_table[index[x1]].nedges;
    for (i = 0; i < nedges; i++) {
      crnt_edge = cell_table[index[x1]].edges[i];
      cx = xtrans; cy = ytrans; cz = ztrans;
      switch (crnt_edge) {
      case 1:
	cx += linterp(v1[0], v1[1], threshold, x1, x2);
	cy += (float) y1;
	cz += (float) z1;
	break;

      case 2:
	cy += linterp(v1[1], v4[1], threshold, y1, y2);
	cx += (float) x2;
	cz += (float) z1;
	break;

      case 3:
	cx += linterp(v4[0], v4[1], threshold, x1, x2);
	cy += (float) y2;
	cz += (float) z1;
	break;

      case 4:
	cy += linterp(v1[0], v4[0], threshold, y1, y2);
	cx += (float) x1;
	cz += (float) z1;
	break;

      case 5:
	cx += linterp(v5[0], v5[1], threshold, x1, x2);
	cy += (float) y1;
	cz += (float) z2;
	break;

      case 6:
	cy += linterp(v5[1], v8[1], threshold, y1, y2);
	cx += (float) x2;
	cz += (float) z2;
	break;

      case 7:
	cx += linterp(v8[0], v8[1], threshold, x1, x2);
	cy += (float) y2;
	cz += (float) z2;
	break;

      case 8:
	cy += linterp(v5[0], v8[0], threshold, y1, y2);
	cx += (float) x1;
	cz += (float) z2;
	break;

      case 9:
	cz += linterp(v1[0], v5[0], threshold, z1, z2);
	cy += (float) y1;
	cx += (float) x1;
	break;

      case 10:
	cz += linterp(v1[1], v5[1], threshold, z1, z2);
	cy += (float) y1;
	cx += (float) x2;
	break;

      case 11:
	cz += linterp(v4[0], v8[0], threshold, z1, z2);
	cy += (float) y2;
	cx += (float) x1;
	break;

      case 12:
	cz += linterp(v4[1], v8[1], threshold, z1, z2);
	cy += (float) y2;
	cx += (float) x2;
	break;

      } /* end switch */
      CROSSINGS(x1,crnt_edge,0) = cx;
      CROSSINGS(x1,crnt_edge,1) = cy;
      CROSSINGS(x1,crnt_edge,2) = cz;
    } /* end for */
  }
}

static int
get_cell_polys(int *index, int xdim, float *crossings)
     /* This subroutine will calculate the polygons */
{

  int num_o_polys, polys = 0;
  int poly;
  float *p1, *p2, *p3;
  float n1[3];
  int x1;

  for (x1 = 0; x1 < xdim-1; x1++) {
    if (!index[x1]) continue;
    num_o_polys = cell_table[index[x1]].npolys;

    for (poly = 0; poly < num_o_polys; poly++) {

      p1 = &CROSSINGS(x1,cell_table[index[x1]].polys[poly*3],0);
      p2 = &CROSSINGS(x1,cell_table[index[x1]].polys[poly*3 + 1],0);
      p3 = &CROSSINGS(x1,cell_table[index[x1]].polys[poly*3 + 2],0);

#define CULL
#ifdef CULL
      if ((p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2]) ||
	  (p1[0] == p3[0] && p1[1] == p3[1] && p1[2] == p3[2]) ||
	  (p2[0] == p3[0] && p2[1] == p3[1] && p2[2] == p3[2]))  {
	/*
	  printf("addpoly - degenerate triangle index = %x poly = %d x1 = %d\n", index[x1], poly, x1);
	  printf("%g %g %g\n%g %g %g\n%g %g %g\n", p1[0], p1[1], p1[2],
	  p2[0], p2[1], p2[2], p3[0], p3[1], p3[2]);
	*/
	polys--;
	continue;
      }
#endif
      calc_normal(p1, p2, p3, n1);
      add_polygon(p1, p2, p3, n1, n1, n1);
    }

    polys += num_o_polys;
  }
  return polys;
}

static void 
calc_cell_normals(int *index, float *data, int y1, int z1,
		  int xdim, int ydim, int zdim, float *norms)
     /* This subroutine calculates the vertex normals using a central difference */
{

  int i, x1;
  float *nn;
  int y11, z11;
#define D(x,y,z) data[(z)*XDIMYDIM + (y)*xdim + (x)]
#ifdef mips
# define SQRT	fsqrt
#else
# define SQRT	sqrt
#endif
  /*
    #define EPSILON		1.e-6
  */
#define EPSILON		0.

  /* clamp */
  y11 = y1+1; z11 = z1+1;
  if (y1 == 0) y1++;
  else if (y1 == ydim - 2) y11--;
  if (z1 == 0) z1++;
  else if (z1 == zdim - 2) z11--;
  nn = norms;

  if (index[0]) {
    nn[0] = D(2, y1, z1) - D(0, y1, z1);
    nn[1] = D(1, y1 + 1, z1) - D(1, y1 - 1, z1);
    nn[2] = D(1, y1, z1 + 1) - D(1, y1, z1 - 1);

    /* nn[3] = nn[0]; nn[4] = nn[1]; nn[5] = nn[2]; */

    nn[6] = D(2, y11, z1) - D(0, y11, z1);
    nn[7] = D(1, y11 + 1, z1) - D(1, y11 - 1, z1);
    nn[8] = D(1, y11, z1 + 1) - D(1, y11, z1 - 1);

    /* nn[9] = nn[6]; nn[10] = nn[7]; nn[11] = nn[8]; */

    nn[12] = D(2, y1, z11) - D(0, y1, z11);
    nn[13] = D(1, y1 + 1, z11) - D(1, y1 - 1, z11);
    nn[14] = D(1, y1, z11 + 1) - D(1, y1, z11 - 1);

    /* nn[15] = nn[12]; nn[16] = nn[13]; nn[17] = nn[14]; */

    nn[18] = D(2, y11, z11) - D(0, y11, z11);
    nn[19] = D(1, y11 + 1, z11) - D(1, y11 - 1, z11);
    nn[20] = D(1, y11, z11 + 1) - D(1, y11, z11 - 1);

    /* nn[21] = nn[18]; nn[22] = nn[19]; nn[23] = nn[20]; */

    for(i = 0; i < 8; i+=2) {
      float d = SQRT(nn[3*i+0] * nn[3*i+0] +
		     nn[3*i+1] * nn[3*i+1] +
		     nn[3*i+2] * nn[3*i+2]);
      if (d > EPSILON) {
	nn[3*i+0] /= d; nn[3*i+1] /= d; nn[3*i+2] /= d;
      }
      nn[3*i+3] = nn[3*i+0];
      nn[3*i+4] = nn[3*i+1];
      nn[3*i+5] = nn[3*i+2];
    }
  }

#ifdef CRAY
#pragma ivdep
  for (x1 = 1; x1 < xdim-2; x1++) {
    if (!index[x1]) continue;
#define NN(i)	norms[x1*24+i]
#define oNN(i)	norms[x1*24-24+i]
    if (!index[x1-1]) {
      /* no coherency */
      /* v1 */
      NN(0) = D(x1 + 1, y1, z1) - D(x1 - 1, y1, z1);
      NN(1) = D(x1, y1 + 1, z1) - D(x1, y1 - 1, z1);
      NN(2) = D(x1, y1, z1 + 1) - D(x1, y1, z1 - 1);

      /* v4 */
      NN(9) = D(x1 + 1, y11, z1) - D(x1 - 1, y11, z1);
      NN(10) = D(x1, y11 + 1, z1) - D(x1, y11 - 1, z1);
      NN(11) = D(x1, y11, z1 + 1) - D(x1, y11, z1 - 1);

      /* v5 */
      NN(12) = D(x1 + 1, y1, z11) - D(x1 - 1, y1, z11);
      NN(13) = D(x1, y1 + 1, z11) - D(x1, y1 - 1, z11);
      NN(14) = D(x1, y1, z11 + 1) - D(x1, y1, z11 - 1);

      /* v8 */
      NN(21) = D(x1 + 1, y11, z11) - D(x1 - 1, y11, z11);
      NN(22) = D(x1, y11 + 1, z11) - D(x1, y11 - 1, z11);
      NN(23) = D(x1, y11, z11 + 1) - D(x1, y11, z11 - 1);
#define NORMALIZE(i) { \
		float d = SQRT(NN(3*(i)+0) * NN(3*(i)+0) + \
			       NN(3*(i)+1) * NN(3*(i)+1) + \
			       NN(3*(i)+2) * NN(3*(i)+2)); \
		if (d > EPSILON) {	\
		    NN(3*(i)+0) /= d; NN(3*(i)+1) /= d; NN(3*(i)+2) /= d; \
		}	\
	    }
      NORMALIZE(0);	/* normalize v1, v4, v5 and v8 */
      NORMALIZE(3);
      NORMALIZE(4);
      NORMALIZE(7);
    }
    /* v2 */
    NN(3) = D(x1 + 1 + 1, y1, z1) - D(x1 + 1 - 1, y1, z1);
    NN(4) = D(x1 + 1, y1 + 1, z1) - D(x1 + 1, y1 - 1, z1);
    NN(5) = D(x1 + 1, y1, z1 + 1) - D(x1 + 1, y1, z1 - 1);

    /* v3 */
    NN(6) = D(x1 + 1 + 1, y11, z1) - D(x1 + 1 - 1, y11, z1);
    NN(7) = D(x1 + 1, y11 + 1, z1) - D(x1 + 1, y11 - 1, z1);
    NN(8) = D(x1 + 1, y11, z1 + 1) - D(x1 + 1, y11, z1 - 1);

    /* v6 */
    NN(15) = D(x1 + 1 + 1, y1, z11) - D(x1 + 1 - 1, y1, z11);
    NN(16) = D(x1 + 1, y1 + 1, z11) - D(x1 + 1, y1 - 1, z11);
    NN(17) = D(x1 + 1, y1, z11 + 1) - D(x1 + 1, y1, z11 - 1);

    /* v7 */
    NN(18) = D(x1 + 1 + 1, y11, z11) - D(x1 + 1 - 1, y11, z11);
    NN(19) = D(x1 + 1, y11 + 1, z11) - D(x1 + 1, y11 - 1, z11);
    NN(20) = D(x1 + 1, y11, z11 + 1) - D(x1 + 1, y11, z11 - 1);
    NORMALIZE(1); /* normalize v2, v3, v6 and v7 */
    NORMALIZE(2);
    NORMALIZE(5);
    NORMALIZE(6);
  }
#pragma ivdep
  for (x1 = 1; x1 < xdim-2; x1++) {
    if (!index[x1]) continue;
    if (!index[x1-1]) continue;
    /* v1 <- v2 */
    NN(0) = oNN(3); NN(1) = oNN(4); NN(2) = oNN(5);
    /* v4 <- v3 */
    NN(9) = oNN(6); NN(10) = oNN(7); NN(11) = oNN(8);
    /* v5 <- v6 */
    NN(12) = oNN(15); NN(13) = oNN(16); NN(14) = oNN(17);
    /* v8 <- v7 */
    NN(21) = oNN(18); NN(22) = oNN(19); NN(23) = oNN(20);
  }
#else
  for (x1 = 1; x1 < xdim-2; x1++) {
    if (!index[x1]) continue;
#ifdef ndef
    nn = norms+x1*24;
#endif
#define NN(i)	norms[x1*24+i]
#define oNN(i)	norms[x1*24-24+i]
    if (!index[x1-1]) {
      /* no coherency */
      /* v1 */
      NN(0) = D(x1 + 1, y1, z1) - D(x1 - 1, y1, z1);
      NN(1) = D(x1, y1 + 1, z1) - D(x1, y1 - 1, z1);
      NN(2) = D(x1, y1, z1 + 1) - D(x1, y1, z1 - 1);

      /* v4 */
      NN(9) = D(x1 + 1, y11, z1) - D(x1 - 1, y11, z1);
      NN(10) = D(x1, y11 + 1, z1) - D(x1, y11 - 1, z1);
      NN(11) = D(x1, y11, z1 + 1) - D(x1, y11, z1 - 1);

      /* v5 */
      NN(12) = D(x1 + 1, y1, z11) - D(x1 - 1, y1, z11);
      NN(13) = D(x1, y1 + 1, z11) - D(x1, y1 - 1, z11);
      NN(14) = D(x1, y1, z11 + 1) - D(x1, y1, z11 - 1);

      /* v8 */
      NN(21) = D(x1 + 1, y11, z11) - D(x1 - 1, y11, z11);
      NN(22) = D(x1, y11 + 1, z11) - D(x1, y11 - 1, z11);
      NN(23) = D(x1, y11, z11 + 1) - D(x1, y11, z11 - 1);
#define NORMALIZE(i) { \
		float d = SQRT(NN(3*(i)+0) * NN(3*(i)+0) + \
			       NN(3*(i)+1) * NN(3*(i)+1) + \
			       NN(3*(i)+2) * NN(3*(i)+2)); \
		if (d > EPSILON) {	\
		    NN(3*(i)+0) /= d; NN(3*(i)+1) /= d; NN(3*(i)+2) /= d; \
		}	\
	    }
      NORMALIZE(0);	/* normalize v1, v4, v5 and v8 */
      NORMALIZE(3);
      NORMALIZE(4);
      NORMALIZE(7);
    } else {
#ifdef ndef
      float *onn = nn - 24;
#endif
      /* v1 <- v2 */
      NN(0) = oNN(3); NN(1) = oNN(4); NN(2) = oNN(5);
      /* v4 <- v3 */
      NN(9) = oNN(6); NN(10) = oNN(7); NN(11) = oNN(8);
      /* v5 <- v6 */
      NN(12) = oNN(15); NN(13) = oNN(16); NN(14) = oNN(17);
      /* v8 <- v7 */
      NN(21) = oNN(18); NN(22) = oNN(19); NN(23) = oNN(20);
    }
    /* v2 */
    NN(3) = D(x1 + 1 + 1, y1, z1) - D(x1 + 1 - 1, y1, z1);
    NN(4) = D(x1 + 1, y1 + 1, z1) - D(x1 + 1, y1 - 1, z1);
    NN(5) = D(x1 + 1, y1, z1 + 1) - D(x1 + 1, y1, z1 - 1);

    /* v3 */
    NN(6) = D(x1 + 1 + 1, y11, z1) - D(x1 + 1 - 1, y11, z1);
    NN(7) = D(x1 + 1, y11 + 1, z1) - D(x1 + 1, y11 - 1, z1);
    NN(8) = D(x1 + 1, y11, z1 + 1) - D(x1 + 1, y11, z1 - 1);

    /* v6 */
    NN(15) = D(x1 + 1 + 1, y1, z11) - D(x1 + 1 - 1, y1, z11);
    NN(16) = D(x1 + 1, y1 + 1, z11) - D(x1 + 1, y1 - 1, z11);
    NN(17) = D(x1 + 1, y1, z11 + 1) - D(x1 + 1, y1, z11 - 1);

    /* v7 */
    NN(18) = D(x1 + 1 + 1, y11, z11) - D(x1 + 1 - 1, y11, z11);
    NN(19) = D(x1 + 1, y11 + 1, z11) - D(x1 + 1, y11 - 1, z11);
    NN(20) = D(x1 + 1, y11, z11 + 1) - D(x1 + 1, y11, z11 - 1);
    NORMALIZE(1); /* normalize v2, v3, v6 and v7 */
    NORMALIZE(2);
    NORMALIZE(5);
    NORMALIZE(6);
  }
#endif
  /* last one */
  if (index[x1]) {
    nn = norms + x1*24;
    nn[0] = D(x1 + 1, y1, z1) - D(x1 - 1, y1, z1);
    nn[1] = D(x1, y1 + 1, z1) - D(x1, y1 - 1, z1);
    nn[2] = D(x1, y1, z1 + 1) - D(x1, y1, z1 - 1);

    /* nn[3] = nn[0]; nn[4] = nn[1]; nn[5] = nn[2]; */

    nn[6] = D(x1 + 1, y11, z1) - D(x1 - 1, y11, z1);
    nn[7] = D(x1, y11 + 1, z1) - D(x1, y11 - 1, z1);
    nn[8] = D(x1, y11, z1 + 1) - D(x1, y11, z1 - 1);

    /* nn[9] = nn[6]; nn[10] = nn[7]; nn[11] = nn[8]; */

    nn[12] = D(x1 + 1, y1, z11) - D(x1 - 1, y1, z11);
    nn[13] = D(x1, y1 + 1, z11) - D(x1, y1 - 1, z11);
    nn[14] = D(x1, y1, z11 + 1) - D(x1, y1, z11 - 1);

    /* nn[15] = n[12]; nn[16] = n[13]; nn[17] = n[14]; */

    nn[18] = D(x1 + 1, y11, z11) - D(x1 - 1, y11, z11);
    nn[19] = D(x1, y11 + 1, z11) - D(x1, y11 - 1, z11);
    nn[20] = D(x1, y11, z11 + 1) - D(x1, y11, z11 - 1);

    /* nn[21] = n[18]; nn[22] = n[19]; nn[23] = n[20]; */

    for(i = 0; i < 8; i+=2) {
      float d = SQRT(nn[3*i+0] * nn[3*i+0] +
		     nn[3*i+1] * nn[3*i+1] +
		     nn[3*i+2] * nn[3*i+2]);
      if (d > EPSILON) {
	nn[3*i+0] /= d; nn[3*i+1] /= d; nn[3*i+2] /= d;
      }
      nn[3*i+3] = nn[3*i+0];
      nn[3*i+4] = nn[3*i+1];
      nn[3*i+5] = nn[3*i+2];
    }
  }
  if (NORMAL_TYPE > 1) {
    for (nn = norms, x1 = 0; x1 < xdim-1; x1++, nn+=24)
      if (index[x1])
	for(i = 0; i < 24; i++) nn[i] = -nn[i];
  }
}


static void 
get_cell_verts_and_norms(int *index, float *data, int y1, int z1, int xdim,
			 float xtrans, float ytrans, float ztrans, float thresh,
			 float *crossings, float *cnorm, float *normals)
{

  int x1, y2, z2;
  float threshold = thresh;
  float *cn;

#define NORMALS(x,a,b) normals[x*13*3+a*3+b]
#define CN(a,b) cn[a*3+b]
#define linterp(a1,a2,a,b1,b2) ((float)(((a-a1) * (float)(b2-b1) / (a2-a1)) + (float)b1))
#define lerp(a1,a2,a,b1,b2) ((a-a1)/(a2-a1)*(b2-b1) + b1)

  y2 = y1 + 1;
  z2 = z1 + 1;
  for (cn = cnorm, x1 = 0; x1 < xdim-1; x1++, cn += 24) {
    float cx, cy, cz;
    float nx = 0., ny = 0., nz = 0.;
    int nedges;
    int crnt_edge;
    int x2 = x1 + 1;
    int i;
    float *v1, *v4, *v5, *v8;


    if (!index[x1]) continue;

#ifdef FOO
    v1 = data + z1*XDIMYDIM + y1*xdim + x1;
    /*v2 = data + z1*XDIMYDIM + y1*xdim + x1+1;*/
    /*v3 = data + z1*XDIMYDIM + y1*xdim+xdim + x1+1;*/
    v4 = data + z1*XDIMYDIM + y1*xdim+xdim + x1;
    v5 = data + z1*XDIMYDIM+XDIMYDIM + y1*xdim + x1;
    /*v6 = data + z1*XDIMYDIM+XDIMYDIM + y1*xdim + x1+1;*/
    /*v7 = data + z1*XDIMYDIM+XDIMYDIM + y1*xdim+xdim + x1+1;*/
    v8 = data + z1*XDIMYDIM+XDIMYDIM + y1*xdim+xdim + x1;
#else
    v1 = data + z1*XDIMYDIM + y1*xdim + x1;
    v4 = v1 + xdim;
    v5 = v1 + XDIMYDIM;
    v8 = v4 + XDIMYDIM;
#endif

    nedges = cell_table[index[x1]].nedges;
    for (i = 0; i < nedges; i++) {
      crnt_edge = cell_table[index[x1]].edges[i];
      cx = xtrans; cy = ytrans; cz = ztrans;
      switch (crnt_edge) {
      case 1:
	cx += linterp(v1[0], v1[1], threshold, x1, x2);
	cy += (float) y1;
	cz += (float) z1;
	nx = lerp(v1[0], v1[1], threshold, CN(0,0), CN(1,0));
	ny = lerp(v1[0], v1[1], threshold, CN(0,1), CN(1,1));
	nz = lerp(v1[0], v1[1], threshold, CN(0,2), CN(1,2));
	break;

      case 2:
	cy += linterp(v1[1], v4[1], threshold, y1, y2);
	cx += (float) x2;
	cz += (float) z1;
	nx = lerp(v1[1], v4[1], threshold, CN(1,0), CN(2,0));
	ny = lerp(v1[1], v4[1], threshold, CN(1,1), CN(2,1));
	nz = lerp(v1[1], v4[1], threshold, CN(1,2), CN(2,2));
	break;

      case 3:
	cx += linterp(v4[0], v4[1], threshold, x1, x2);
	cy += (float) y2;
	cz += (float) z1;
	nx = lerp(v4[1], v4[0], threshold, CN(2,0), CN(3,0));
	ny = lerp(v4[1], v4[0], threshold, CN(2,1), CN(3,1));
	nz = lerp(v4[1], v4[0], threshold, CN(2,2), CN(3,2));
	break;

      case 4:
	cy += linterp(v1[0], v4[0], threshold, y1, y2);
	cx += (float) x1;
	cz += (float) z1;
	nx = lerp(v4[0], v1[0], threshold, CN(3,0), CN(0,0));
	ny = lerp(v4[0], v1[0], threshold, CN(3,1), CN(0,1));
	nz = lerp(v4[0], v1[0], threshold, CN(3,2), CN(0,2));
	break;

      case 5:
	cx += linterp(v5[0], v5[1], threshold, x1, x2);
	cy += (float) y1;
	cz += (float) z2;
	nx = lerp(v5[0], v5[1], threshold, CN(4,0), CN(5,0));
	ny = lerp(v5[0], v5[1], threshold, CN(4,1), CN(5,1));
	nz = lerp(v5[0], v5[1], threshold, CN(4,2), CN(5,2));
	break;

      case 6:
	cy += linterp(v5[1], v8[1], threshold, y1, y2);
	cx += (float) x2;
	cz += (float) z2;
	nx = lerp(v5[1], v8[1], threshold, CN(5,0), CN(6,0));
	ny = lerp(v5[1], v8[1], threshold, CN(5,1), CN(6,1));
	nz = lerp(v5[1], v8[1], threshold, CN(5,2), CN(6,2));
	break;

      case 7:
	cx += linterp(v8[0], v8[1], threshold, x1, x2);
	cy += (float) y2;
	cz += (float) z2;
	nx = lerp(v8[1], v8[0], threshold, CN(6,0), CN(7,0));
	ny = lerp(v8[1], v8[0], threshold, CN(6,1), CN(7,1));
	nz = lerp(v8[1], v8[0], threshold, CN(6,2), CN(7,2));
	break;

      case 8:
	cy += linterp(v5[0], v8[0], threshold, y1, y2);
	cx += (float) x1;
	cz += (float) z2;
	nx = lerp(v8[0], v5[0], threshold, CN(7,0), CN(4,0));
	ny = lerp(v8[0], v5[0], threshold, CN(7,1), CN(4,1));
	nz = lerp(v8[0], v5[0], threshold, CN(7,2), CN(4,2));
	break;

      case 9:
	cz += linterp(v1[0], v5[0], threshold, z1, z2);
	cy += (float) y1;
	cx += (float) x1;
	nx = lerp(v1[0], v5[0], threshold, CN(0,0), CN(4,0));
	ny = lerp(v1[0], v5[0], threshold, CN(0,1), CN(4,1));
	nz = lerp(v1[0], v5[0], threshold, CN(0,2), CN(4,2));
	break;

      case 10:
	cz += linterp(v1[1], v5[1], threshold, z1, z2);
	cy += (float) y1;
	cx += (float) x2;
	nx = lerp(v1[1], v5[1], threshold, CN(1,0), CN(5,0));
	ny = lerp(v1[1], v5[1], threshold, CN(1,1), CN(5,1));
	nz = lerp(v1[1], v5[1], threshold, CN(1,2), CN(5,2));
	break;

      case 11:
	cz += linterp(v4[0], v8[0], threshold, z1, z2);
	cy += (float) y2;
	cx += (float) x1;
	nx = lerp(v4[0], v8[0], threshold, CN(3,0), CN(7,0));
	ny = lerp(v4[0], v8[0], threshold, CN(3,1), CN(7,1));
	nz = lerp(v4[0], v8[0], threshold, CN(3,2), CN(7,2));
	break;

      case 12:
	cz += linterp(v4[1], v8[1], threshold, z1, z2);
	cy += (float) y2;
	cx += (float) x2;
	nx = lerp(v4[1], v8[1], threshold, CN(2,0), CN(6,0));
	ny = lerp(v4[1], v8[1], threshold, CN(2,1), CN(6,1));
	nz = lerp(v4[1], v8[1], threshold, CN(2,2), CN(6,2));
	break;

      } /* end switch */
      CROSSINGS(x1,crnt_edge,0) = cx;
      CROSSINGS(x1,crnt_edge,1) = cy;
      CROSSINGS(x1,crnt_edge,2) = cz;
      NORMALS(x1,crnt_edge,0) = nx;
      NORMALS(x1,crnt_edge,1) = ny;
      NORMALS(x1,crnt_edge,2) = nz;
    } /* end for */
  }
}

static int
get_cell_polys_and_norms(int *index, int xdim, float *crossings,
			 float *normals)
     /* This subroutine will calculate the polygons */
{

  int num_o_polys, polys = 0;
  int poly;
  float *p1, *p2, *p3;
  float *n1, *n2, *n3;
  int x1;

  for (x1 = 0; x1 < xdim-1; x1++) {
    if (!index[x1]) continue;
    num_o_polys = cell_table[index[x1]].npolys;
    for (poly = 0; poly < num_o_polys; poly++) {

      p1 = &CROSSINGS(x1,cell_table[index[x1]].polys[poly*3],0);
      p2 = &CROSSINGS(x1,cell_table[index[x1]].polys[poly*3 + 1],0);
      p3 = &CROSSINGS(x1,cell_table[index[x1]].polys[poly*3 + 2],0);
#ifdef CULL
      if ((p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2]) ||
	  (p1[0] == p3[0] && p1[1] == p3[1] && p1[2] == p3[2]) ||
	  (p2[0] == p3[0] && p2[1] == p3[1] && p2[2] == p3[2]))  {
	/*
	  printf("addpoly - degenerate triangle index = %x poly = %d x1 = %d\n", index[x1], poly, x1);
	  printf("%g %g %g\n%g %g %g\n%g %g %g\n", p1[0], p1[1], p1[2],
	  p2[0], p2[1], p2[2], p3[0], p3[1], p3[2]);
	*/
	polys--;
	continue;
      }
#endif

      n1 = &NORMALS(x1,cell_table[index[x1]].polys[poly*3],0);
      n2 = &NORMALS(x1,cell_table[index[x1]].polys[poly*3 + 1],0);
      n3 = &NORMALS(x1,cell_table[index[x1]].polys[poly*3 + 2],0);

      add_polygon(p1, p2, p3, n1, n2, n3);
    }

    polys += num_o_polys;
  }
  return polys;
}
