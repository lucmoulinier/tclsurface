#ifndef POLYR_H
#define POLYR_H

#include <tcl.h>

/* Tcl 9 : list lengths and objc are Tcl_Size (64 bits) ;
   compatibility definitions for Tcl 8.6 headers (wiki: Porting
   extensions to Tcl 9) */
#ifndef TCL_SIZE_MAX
# include <limits.h>
# define TCL_SIZE_MAX INT_MAX
# ifndef Tcl_Size
typedef int Tcl_Size;
# endif
# define TCL_SIZE_MODIFIER ""
# define Tcl_GetSizeIntFromObj Tcl_GetIntFromObj
#endif

#define FALSE 0
#define TRUE 1

/* upper bound for atomic numbers (index in the vdw radius table) */
#define NVDW 110

#ifdef MAIN

int first_alloc = 1;
int NORMAL_TYPE = 1;
char MY_NAME[80];
int VERBOSE = 0;
int STORE_POLYGONS = 0;
int SAMPLE = 1;
int SMOOTH = 0;
int CONNOLLY = 0;
int NEXPAND = 0;
int NCONTRACT = 0;
float DOWNSAMPL[3]={1.0,1.0,1.0};
int FILTER[3]={0,0,0};
int REQPOLS = 0;
int MAXPOLS = 150;
int XDIM, YDIM, ZDIM;
int GRID_SHIFT = 0;
int RAW_INPUT, ANA_INPUT, MS_INPUT, HDF_INPUT;
float delta[3]={1.0,1.0,1.0};
float tid1,tid2,tid3,tid4,tid5,tid6;
int NUM_VERTICES = 0;
int VERT_LIMIT = 0;
int VERT_INCR = 20000;
float *VERTICES = NULL;
float *NORMALS = NULL;

#else

extern int first_alloc;
extern float delta[3];
extern int NORMAL_TYPE;
extern char MY_NAME[80];
extern int REQPOLS;
extern int VERBOSE;
extern int STORE_POLYGONS;
extern int SAMPLE;
extern int SMOOTH;
extern int CONNOLLY;
extern int NEXPAND;
extern int NCONTRACT;
extern float DOWNSAMPL[3];
extern int FILTER[3];
extern int MAXPOLS;
extern int XDIM, YDIM, ZDIM;
extern int GRID_SHIFT;	/* grid points added on each side by expand/contract */
extern int RAW_INPUT, ANA_INPUT, MS_INPUT, HDF_INPUT;
extern float tid1,tid2,tid3,tid4,tid5,tid6;
extern int NUM_VERTICES;
extern int VERT_LIMIT;
extern int VERT_INCR;
extern float *VERTICES;
extern float *NORMALS;

#endif

/* misc.c
 * tclsurface_fatal() replaces the old exit() calls: it stores the message
 * and longjmp()s back to the tclsurface command, which cleans up and returns
 * TCL_ERROR instead of killing the interpreter. */
#include <setjmp.h>
extern jmp_buf tclsurface_jmpbuf;
extern char tclsurface_errmsg[256];
void tclsurface_fatal(const char *fmt, ...);
void *tclsurface_alloc(size_t size);
void fill_grid_sides(int ***tmp, int nx, int ny, int nz);
int ***alloc_grid(int nx, int ny, int nz);
int dealloc_grid(int ***tmp, int nx, int ny, int nz);
void calc_normal(float p1[3], float p2[3], float p3[3], float n[3]);
int add_polygon(float *p1, float *p2, float *p3, float *n1, float *n2, float *n3);
void release_memory(void);

/* parseopts.c */
int parse_options(Tcl_Interp *interp, Tcl_Size objc, Tcl_Obj *const objv[]);

/* grid.c */
int parse_input(Tcl_Interp *interp, int ndata, Tcl_Obj **OData,
		double **xyz, int *na);
int setup_data(double **xyz, int *na, int ndata, float xxm[2][3], int *offset);
int define_grid(int offset, float xxm[2][3], int *sx, int *sy, int *sz);
float *fill_data(double **xyz, int *na, int ndata, int offset, float xxm[2][3]);
void free_spheres(void);

/* resample.c */
float *downsample(float *data, int *size, int *newsize);
void boxfilter(float **data, int *size, int *boxsize);

/* niso.c */
int iso_surface(float *data, int xdim, int ydim, int zdim, float threshold);

/* polyreduce.c */
int reduce(Tcl_Interp *interp, int offset, float xxm[2][3]);
void reduce_cleanup(void);

#endif /* POLYR_H */
