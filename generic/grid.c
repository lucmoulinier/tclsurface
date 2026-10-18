#include <stdio.h>
#include <string.h>
#include "polyr.h"


static float vdw[NVDW];
/* sps[n] : number of values, then the (i,j,k) offsets of the grid
   points inside the sphere of element n (NULL if no such atom) */
static int *sps[NVDW];

/* prototype functions */
static int ***expand_mask(int ***grid);
static int ***smooth_mask(int ***grid);
static int ***contract_mask(int ***grid);

#ifdef TCLSURFACE_DEBUG
static void 
show_grid(int ***grid) 
{
  int i, j, k;
  char c = '?';

  fprintf(stderr, "\nDump of Grid [%d,%d,%d] :\n", XDIM, YDIM, ZDIM);
  for(i=0;i<ZDIM;i++) {
    for(j=0;j<YDIM;j++) {
      for(k=0;k<XDIM;k++) {
	if(grid[i][j][k] == 0) 
	  c='.';
	else if (grid[i][j][k] == 1)
	  c='O';
	else 
	  fprintf(stderr, "error : %d %d %d\n", i, j, k);
	fprintf(stderr, "%c",c);
      }
      fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n\n");
  }

  return;
}


#endif

static void
init_vdw(void) 
{
  vdw[109]=0.;
  vdw[108]=0.; vdw[107]=0.; vdw[106]=0.;
  vdw[105]=0.; vdw[104]=0.; vdw[103]=0.;
  vdw[102]=0.; vdw[101]=0.; vdw[100]=0.;
  vdw[99]=0.; vdw[98]=0.; vdw[97]=0.;
 vdw[96]=0.; vdw[95]=0.; vdw[94]=0.;
 vdw[93]=0.; vdw[92]=1.86; vdw[91]=0.;
 vdw[90]=0.; vdw[89]=0.; vdw[88]=0.;
 vdw[87]=0.; vdw[86]=0.; vdw[85]=0.;
 vdw[84]=0.; vdw[83]=0.; vdw[82]=2.02;
 vdw[81]=1.96; vdw[80]=1.55; vdw[79]=1.66;
 vdw[78]=1.75; vdw[77]=0.; vdw[76]=0.;
 vdw[75]=0.; vdw[74]=0.; vdw[73]=0.;
 vdw[72]=0.; vdw[71]=0.; vdw[70]=0.;
 vdw[69]=0.; vdw[68]=0.; vdw[67]=0.;
 vdw[66]=0.; vdw[65]=0.; vdw[64]=0.;
 vdw[63]=0.; vdw[62]=0.; vdw[61]=0.;
 vdw[60]=0.; vdw[59]=0.; vdw[58]=0.;
 vdw[57]=0.; vdw[56]=0.; vdw[55]=0.;
 vdw[54]=2.16; vdw[53]=1.98; vdw[52]=2.06;
 vdw[51]=0.; vdw[50]=2.17; vdw[49]=1.93;
 vdw[48]=1.58; vdw[47]=1.72; vdw[46]=1.63;
 vdw[45]=0.; vdw[44]=0.; vdw[43]=0.;
 vdw[42]=0.; vdw[41]=0.; vdw[40]=0.;
 vdw[39]=0.; vdw[38]=0.; vdw[37]=0.;
 vdw[36]=2.02; vdw[35]=1.85; vdw[34]=1.9;
 vdw[33]=1.85; vdw[32]=0.; vdw[31]=1.87;
 vdw[30]=1.39; vdw[29]=1.4; vdw[28]=1.63;
 vdw[27]=2.23; vdw[26]=2.23; vdw[25]=2.23;
 vdw[24]=2.23; vdw[23]=2.23; vdw[22]=2.23;
 vdw[21]=2.23; vdw[20]=2.23; vdw[19]=2.75;
 vdw[18]=1.88; vdw[17]=1.75; vdw[16]=1.8;
 vdw[15]=1.8; vdw[14]=2.1; vdw[13]=1.7;
 vdw[12]=1.73; vdw[11]=2.27; vdw[10]=0.;
 vdw[9]=1.47; vdw[8]=1.52; vdw[7]=1.55;
 vdw[6]=1.70; vdw[5]=1.7; vdw[4]=1.7;
 vdw[3]=1.82; vdw[2]=1.4; vdw[1]=1.2;
 vdw[0]=0.;
}


void
free_spheres(void)
{
  int n;

  for (n=0; n<NVDW; n++) {
    if (sps[n]) {
      ckfree((char *)sps[n]);
      sps[n] = NULL;
    }
  }
}


/* build the sphere of grid points of every element present in the
   atom set and return the largest radius (in grid points) */
static int 
sample_vdw(int *na, int ndata) 
{
  int offset, r, i, j, k, n, v, a;

  free_spheres();
  offset = 0;
  for (a=0; a<ndata; a++) {
    n = na[a];
    if (vdw[n] == 0. || sps[n]) continue;
    
    r = (int)(vdw[n]*(float)SAMPLE+0.5);
    if (r > offset) offset = r;
    sps[n] = (int *)tclsurface_alloc((1 + 3*(size_t)(2*r+1)*(2*r+1)*(2*r+1))
				* sizeof(int));
    v = 1;
    for (i=-r;i<=r;i++) {
      for (j=-r;j<=r;j++) {
	for (k=-r;k<=r;k++) {
	  if ((i*i + j*j + k*k) <= r*r) {
	    sps[n][v++] = i;
	    sps[n][v++] = j;
	    sps[n][v++] = k;
	  }
	}
      }
    }
    sps[n][0] = v-1;
  }
  if (VERBOSE)
    fprintf(stderr, "%s: largest atom radius %d grid points\n", MY_NAME, offset);

  return offset;
}

/* get min and max of atom set coordinates */
static void
get_xmin_xmax(double **xyz, int ndata, float xxm[2][3]) 
{
  int i,j;
  
  /* initialise min et max with the first atom */
  for (j=0;j<3;j++)
    xxm[0][j] = xxm[1][j] = xyz[0][j];

  /* [0] = min ; [1] = max */
  for (i=1;i<ndata;i++) {
    for (j=0;j<3;j++) {
      xxm[0][j] = xyz[i][j]<xxm[0][j]?xyz[i][j]:xxm[0][j];
      xxm[1][j] = xyz[i][j]>xxm[1][j]?xyz[i][j]:xxm[1][j];
    }
  }

  if (VERBOSE) {
    fprintf(stderr, "xmin %6.3f ymin %6.3f zmin %6.3f\n", xxm[0][0],xxm[0][1],xxm[0][2]);
    fprintf(stderr, "xmax %6.3f ymax %6.3f zmax %6.3f\n", xxm[1][0],xxm[1][1],xxm[1][2]);
  }
}


/* define the grid */
int 
define_grid(int offset, float xxm[2][3], int *sx, int *sy, int *sz) 
{
  int i, gr[3];

  for (i=0;i<3;i++) {
    gr[i] = (int)(0.5+((float)xxm[1][i]-xxm[0][i])*(float)SAMPLE) + 3 + 2*offset;
    
    if(! (gr[i] % 2)) gr[i]++;
  }

  *sx = gr[0];
  *sy = gr[1];
  *sz = gr[2];

  return TCL_OK;
}


/* setup parameters from data */
int 
setup_data(double **xyz, int *na, int ndata, float xxm[2][3], int *offset) 
{
  get_xmin_xmax(xyz, ndata, xxm);
  init_vdw();
  *offset = sample_vdw(na, ndata);

  return TCL_OK;
}


static int 
***do_connolly(int ***grid) 
{
  int n, i;

  /* water radius = 1.4 A */
  n = (int)(1.4*(float)SAMPLE+0.5);
  for(i=0; i <= n; i++)
    grid = expand_mask(grid);

  return grid;
}


static int 
count_neighbours(int ***grid, int i, int j, int k, int val)
{
  int nbc;

  nbc = 
    grid[i-1][j  ][k  ] + grid[i+1][j  ][k  ] + 
    grid[i  ][j-1][k  ] + grid[i  ][j+1][k  ] + 
    grid[i  ][j  ][k-1] + grid[i  ][j  ][k+1] + 
    grid[i-1][j-1][k  ] + grid[i-1][j+1][k  ] + 
    grid[i+1][j-1][k  ] + grid[i+1][j+1][k  ] + 
    grid[i-1][j  ][k-1] + grid[i-1][j  ][k+1] + 
    grid[i+1][j  ][k-1] + grid[i+1][j  ][k+1] + 
    grid[i  ][j-1][k-1] + grid[i  ][j-1][k+1] + 
    grid[i  ][j+1][k-1] + grid[i  ][j+1][k+1] + 
    grid[i-1][j-1][k-1] + grid[i-1][j-1][k+1] + 
    grid[i-1][j+1][k-1] + grid[i-1][j+1][k+1] + 
    grid[i+1][j-1][k-1] + grid[i+1][j-1][k+1] + 
    grid[i+1][j+1][k-1] + grid[i+1][j+1][k+1];
  
  /* if we count the O's */
  if(val == 0) nbc = 26 - nbc;

  return nbc;
}


static int 
***smooth_mask(int ***grid)
{
  int i, j, k, nbn, factor, ***tmp, nchg;

  if (VERBOSE) fprintf(stderr, "smoothing mask ....\n");
  tmp = alloc_grid(XDIM, YDIM, ZDIM);
  fill_grid_sides(tmp, XDIM, YDIM, ZDIM);

  /* Mask Smoothing :
     For every non-mask point, the number of 
     neighbours which ARE IN the mask are counted; 
     if this number (0-26) is greater than or equal 
     to the parameter "FACTOR", the point is included
     in the mask 
  */

  factor = 15;
  nchg = 0;
  for(i=1; i<ZDIM-1; i++) {
    for(j=1; j<YDIM-1; j++) {
      for(k=1; k<XDIM-1; k++) {
	if(! grid[i][j][k]) {
	  nbn = count_neighbours(grid, i, j, k, 1);
	  if(nbn >= factor) {
	    tmp[i][j][k] = 1;
	    nchg++;
	  }
	  else tmp[i][j][k] = 0;
	} else {
	  tmp[i][j][k] = 1;
	}
      }
    }
  }

  dealloc_grid(grid, XDIM, YDIM, ZDIM);

  return tmp;
}


static int 
***expand_mask(int ***grid) 
{
  int i, j, k, nbn;
  int ***tmp, ***tmp2, newx, newy, newz;
  
  /* create a new extended grid */
  newx = XDIM+2; newy = YDIM+2; newz = ZDIM+2;
  if (VERBOSE)
    fprintf(stderr, "expanding to [%d,%d,%d] ....\n", newx, newy, newz);

  tmp  = alloc_grid(newx, newy, newz);
  tmp2 = alloc_grid(newx, newy, newz);

  fill_grid_sides(tmp , newx, newy, newz);
  fill_grid_sides(tmp2, newx, newy, newz);

  /* copy grid to tmp */
  for(i=1; i<newz-1; i++)
    for(j=1; j<newy-1; j++)
      for(k=1; k<newx-1; k++)
	tmp[i][j][k] = tmp2[i][j][k] = grid[i-1][j-1][k-1];

  /* do the expansion :
     For every non-mask point, the number of 
     neighbours (nbn) which ARE IN the mask are 
     counted; 
     If nbn >= 1, the point is included in the 
     mask.
  */
  for(i=1; i<newz-1; i++) {
    for(j=1; j<newy-1; j++) {
      for(k=1; k<newx-1; k++) {
	if(! tmp[i][j][k]) {
	  nbn = count_neighbours(tmp, i, j, k, 1);
	  if(nbn >= 1) tmp2[i][j][k] = 1;
	}
      }
    }
  }

  /* dealloc memory */
  dealloc_grid(grid, XDIM, YDIM, ZDIM);
  dealloc_grid(tmp, newx, newy, newz);

  XDIM = newx ; YDIM = newy; ZDIM = newz;
  GRID_SHIFT++;

  return tmp2;
}


static int 
***contract_mask(int ***grid) 
{
  int i, j, k, nbn, ***tmp, newx, newy, newz;

  /* create a new reduced grid */
  newz = ZDIM-2; newy=YDIM-2; newx=XDIM-2;
  if (newx < 3 || newy < 3 || newz < 3) {
    dealloc_grid(grid, XDIM, YDIM, ZDIM);
    tclsurface_fatal("-contract: grid too small to be contracted");
  }
  if (VERBOSE)
    fprintf(stderr, "contracting to [%d,%d,%d] ....\n", newx, newy, newz);

  tmp = alloc_grid(newx, newy, newz);
  fill_grid_sides(tmp, newx, newy, newz);

  /* 
     for each point IN the mask, remove it from the 
     mask if there is ONE neighbour which is NOT in 
     the mask
  */
  for(i=1; i<ZDIM-1; i++) {
    for(j=1; j<YDIM-1; j++) {
      for(k=1; k<XDIM-1; k++) {
	if(grid[i][j][k]) {
	  nbn = count_neighbours(grid, i, j, k, 0);
	  if(nbn) tmp[i-1][j-1][k-1] = 0;
	  else tmp[i-1][j-1][k-1] = 1;
	} else {
	  tmp[i-1][j-1][k-1] = 0;
	}
      }
    }
  }
	 
  /* dealloc memory */
  dealloc_grid(grid , XDIM, YDIM, ZDIM);
  XDIM = newx ; YDIM = newy; ZDIM = newz;
  GRID_SHIFT--;

  return tmp;
}


/* fills the grid with atoms with their vdw radius */
float 
*fill_data(double **xyz, int *na, int ndata, int offset, float xxm[2][3]) 
{
  int i, j, k, n, s, a, b, c, ***grid, ai[3];
  float *tmp;

  (void)offset;
  grid = alloc_grid(XDIM, YDIM, ZDIM);
  memset(grid[0][0], 0, (size_t)XDIM*YDIM*ZDIM*sizeof(int));
  GRID_SHIFT = 0;
  
  for(i=0; i<ndata; i++) {
    for(j=0; j<3; j++) {
      ai[j] = (int)(0.5+((float)xyz[i][j]-xxm[0][j])*(float)SAMPLE) + 1 + offset;
    }
    
    n = na[i];
    if (sps[n] == NULL) continue;	/* element without vdw radius */
    for (s=1; s <= sps[n][0]; s+=3) {
      a = ai[0] + sps[n][s];
      b = ai[1] + sps[n][s+1];
      c = ai[2] + sps[n][s+2];

      grid[c][b][a] = 1;
    }
  }
  /* show_grid(grid); */

  /* expand if required */
  for(i=0; i<NEXPAND; i++) 
    grid = expand_mask(grid);
  /* show_grid(grid); */

  if(SMOOTH) 
    grid = smooth_mask(grid);
  /* show_grid(grid); */
  
  /* contract if required */
  for(i=0; i<NCONTRACT; i++) 
    grid = contract_mask(grid);
  /* show_grid(grid); */

  if(CONNOLLY)
    grid = do_connolly(grid);

  /* fill the data array */
  tmp = (float *)attemptckalloc((size_t)XDIM*YDIM*ZDIM*sizeof(float));
  if (tmp == NULL) {
    dealloc_grid(grid, XDIM, YDIM, ZDIM);
    tclsurface_fatal("not enough memory for the data array: "
		"try a smaller -sample value");
  }
  n = 0;
  for(i=0 ; i < ZDIM ; i++) {
    for(j=0 ; j < YDIM ; j++) { 
      for(k=0 ; k < XDIM ; k++) {
	tmp[n++] = (float)grid[i][j][k] + 48.;
      }
    }
  }

  /* dealloc memory */
  dealloc_grid(grid, XDIM, YDIM, ZDIM);

  return tmp;
}


/* parses the input data : a list of {x y z atomic_number ?tag?} */
int 
parse_input(Tcl_Interp *interp, 
	    int ndata, 
	    Tcl_Obj **OData, 
	    double **xyz, 
	    int *na) 
{
  int i, j;
  Tcl_Size nval;
  Tcl_Obj **Oxyz;

  for(i=0;i<ndata;i++) {
    if (Tcl_ListObjGetElements(interp,OData[i],&nval,&Oxyz) != TCL_OK)
      goto bad;
    if (nval < 4) {
      Tcl_SetObjResult(interp, Tcl_ObjPrintf(
	  "atom %d: expected {x y z atomic_number ?tag?}, got \"%s\"",
	  i, Tcl_GetString(OData[i])));
      return TCL_ERROR;
    }
    for (j=0;j<3;j++)
      if (Tcl_GetDoubleFromObj(interp,Oxyz[j],&xyz[i][j]) != TCL_OK)
	goto bad;
    if (Tcl_GetIntFromObj(interp,Oxyz[3],&na[i]) != TCL_OK)
      goto bad;
    if (na[i] < 0 || na[i] >= NVDW) {
      Tcl_SetObjResult(interp, Tcl_ObjPrintf(
	  "atom %d: atomic number %d out of range [0,%d]", i, na[i], NVDW-1));
      return TCL_ERROR;
    }
  }

  return TCL_OK;

 bad:
  Tcl_AppendObjToErrorInfo(interp, Tcl_ObjPrintf(
      "\n    (reading atom %d of tclsurface data)", i));
  return TCL_ERROR;
}
