#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include "polyr.h"


#ifndef DLLEXPORT
#define DLLEXPORT __declspec(dllexport)
#endif

/* the computation uses global state: one tclsurface at a time per process */
TCL_DECLARE_MUTEX(tclsurfaceMutex)

/* memory owned by the running command, freed by work_cleanup() also
   after a tclsurface_fatal() */
static float *w_data = NULL;
static double **w_xyz = NULL;
static int *w_na = NULL;
static int w_ndata = 0;

static void
work_cleanup(void)
{
  int i;

  reduce_cleanup();
  release_memory();
  free_spheres();
  if (w_data) {
    ckfree((char *)w_data);
    w_data = NULL;
  }
  if (w_xyz) {
    for (i=0; i<w_ndata; i++)
      if (w_xyz[i]) ckfree((char *)w_xyz[i]);
    ckfree((char *)w_xyz);
    w_xyz = NULL;
  }
  if (w_na) {
    ckfree((char *)w_na);
    w_na = NULL;
  }
  w_ndata = 0;
}


/*
 * tclsurface ?options? atomList
 *   atomList : list of {x y z atomic_number ?tag?}
 * returns {vertices normals}, 9 values per triangle in each list
 */
static int
TclSurfaceCompute(Tcl_Interp *interp, Tcl_Size objc, Tcl_Obj *const objv[])
{
  int i;
  int xdim, ydim, zdim;	/* the dimensions of the data */
  float *newdata;
  int offset;
  int size[3], ndata;
  Tcl_Size listlen;
  int newsize[3];
  float threshold;
  float xxm[2][3];
  Tcl_Obj **OData;

  /* get rutime name */
  strncpy(MY_NAME, Tcl_GetString(objv[0]), sizeof(MY_NAME)-1);
  MY_NAME[sizeof(MY_NAME)-1] = '\0';

  if (objc < 2) {
    Tcl_WrongNumArgs(interp, 1, objv, "?options? atomList");
    return TCL_ERROR;
  }

  /* initialisations */
  first_alloc  = 1;
  NUM_VERTICES = 0;
  VERT_LIMIT   = 0;
  VERTICES     = NULL;
  NORMALS      = NULL;
  delta[0] = delta[1] = delta[2] = 1.0;

  /* get command line args */
  if (parse_options(interp, objc, objv) != TCL_OK)
    return TCL_ERROR;

  /* last arg are input data */
  if (Tcl_ListObjGetElements(interp, objv[(objc - 1)], &listlen, &OData)
      != TCL_OK)
    return TCL_ERROR;
  if (listlen > INT_MAX / (Tcl_Size)sizeof(double *)) {
    Tcl_SetObjResult(interp, Tcl_NewStringObj("too many atoms", -1));
    return TCL_ERROR;
  }
  ndata = (int)listlen;
  threshold = 48.5;

  if (VERBOSE)
    fprintf(stderr,"\n%s: <Jesper James Jensen 1995>\n",MY_NAME);

  /* any error from here on comes back through tclsurface_fatal() */
  if (setjmp(tclsurface_jmpbuf)) {
    work_cleanup();
    Tcl_SetObjResult(interp, Tcl_ObjPrintf("%s: %s", MY_NAME, tclsurface_errmsg));
    return TCL_ERROR;
  }

  if (ndata == 0) {
    /* no atom, no surface */
    Tcl_Obj *res[2];
    res[0] = Tcl_NewListObj(0, NULL);
    res[1] = Tcl_NewListObj(0, NULL);
    Tcl_SetObjResult(interp, Tcl_NewListObj(2, res));
    return TCL_OK;
  }

  /* parse data and fill table */
  w_na  = (int *)tclsurface_alloc(ndata*sizeof(int));
  w_xyz = (double **)tclsurface_alloc(ndata*sizeof(double *));
  for (i=0;i<ndata;i++)
    w_xyz[i] = NULL;
  w_ndata = ndata;
  for (i=0;i<ndata;i++)
    w_xyz[i] = (double *)tclsurface_alloc(3*sizeof(double));
  if (parse_input(interp, ndata, OData, w_xyz, w_na) != TCL_OK) {
    work_cleanup();
    return TCL_ERROR;
  }
  setup_data(w_xyz, w_na, ndata, xxm, &offset);
  define_grid(offset, xxm, &xdim, &ydim, &zdim);
  XDIM = size[0] = xdim;
  YDIM = size[1] = ydim;
  ZDIM = size[2] = zdim;

  w_data = fill_data(w_xyz, w_na, ndata, offset, xxm);

  /* as grid size may have change ... */
  xdim = size[0] = XDIM;
  ydim = size[1] = YDIM;
  zdim = size[2] = ZDIM;

  tid1=clock();

  /* filtering of the input array */
  if (FILTER[1]) {
    if (VERBOSE)
      fprintf(stderr,"%s: Mean filtering the array with a filter [%dx%dx%d]\n",
	      MY_NAME,FILTER[0],FILTER[1],FILTER[2]);
    boxfilter(&w_data,size,FILTER);
  }

  tid2=clock();

  /* Ok both downsampling and isosurfacing! */
  if ((DOWNSAMPL[0]!=1.0) ||(DOWNSAMPL[1]!=1.0) ||(DOWNSAMPL[2]!=1.0)) {
    if (VERBOSE)
      fprintf(stderr,"%s: Downsampling array with x:%.2f y:%.2f z:%.2f\n",
	      MY_NAME,DOWNSAMPL[0],DOWNSAMPL[1],DOWNSAMPL[2]);
    for (i=0; i<3; i++) {
      newsize[i] = (int)(size[i]/DOWNSAMPL[i]);
      if (newsize[i] < 3)
	tclsurface_fatal("downsampling ratio too large for a %d points grid",
		    size[i]);
      /* real step of the new grid, used to scale the vertices back */
      delta[i] = (float)size[i]/newsize[i];
    }
    xdim = newsize[0];
    ydim = newsize[1];
    zdim = newsize[2];
    newdata = downsample(w_data,size,newsize);
    ckfree((char *)w_data);
    w_data = newdata;
  }

  if (VERBOSE)
    fprintf(stderr,"%s: Doing isosurfacing\n",MY_NAME);
  iso_surface(w_data, xdim, ydim, zdim, threshold);
  ckfree((char *)w_data);
  w_data = NULL;

  /* polygon reduction, sets the result */
  reduce(interp, offset, xxm);

  if (VERBOSE)
    fflush(stderr);

  /* deallocation */
  work_cleanup();

  return TCL_OK;
}


/* Tcl 9: Tcl_ObjCmdProc2 receives objc as a Tcl_Size */
#if TCL_MAJOR_VERSION > 8
# define TCLSURFACE_CREATE_CMD Tcl_CreateObjCommand2
typedef Tcl_Size TclsurfaceObjc;
#else
# define TCLSURFACE_CREATE_CMD Tcl_CreateObjCommand
typedef int TclsurfaceObjc;
#endif

static int
TclSurfaceObjCmd(void *clientData,
	    Tcl_Interp *interp,
	    TclsurfaceObjc objc,
	    Tcl_Obj *const objv[])
{
  int code;

  (void)clientData;
  Tcl_MutexLock(&tclsurfaceMutex);
  code = TclSurfaceCompute(interp, objc, objv);
  Tcl_MutexUnlock(&tclsurfaceMutex);

  return code;
}


int DLLEXPORT Tclsurface_Init (Tcl_Interp *interp)
{
  if (Tcl_InitStubs(interp, "8.5-", 0) == NULL) {
    return TCL_ERROR;
  }

  TCLSURFACE_CREATE_CMD(interp, "tclsurface", TclSurfaceObjCmd, NULL, NULL);

  if (Tcl_PkgProvide(interp, PACKAGE_NAME, PACKAGE_VERSION) != TCL_OK) {
    return TCL_ERROR;
  }

  return TCL_OK;
}
