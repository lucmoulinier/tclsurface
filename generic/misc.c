#include <stdio.h>
#include <stdarg.h>
#include <limits.h>
#include <math.h>
#include "polyr.h"

jmp_buf tclsurface_jmpbuf;
char tclsurface_errmsg[256];

/* abort the current tclsurface command with an error message */
void
tclsurface_fatal(const char *fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);
  vsnprintf(tclsurface_errmsg, sizeof(tclsurface_errmsg), fmt, ap);
  va_end(ap);
  longjmp(tclsurface_jmpbuf, 1);
}


/* ckalloc that fails with a Tcl error instead of a panic */
void *
tclsurface_alloc(size_t size)
{
  void *p;

  if (size == 0) size = 1;
  if (size > (size_t)INT_MAX
      || (p = (void *)attemptckalloc(size)) == NULL)
    tclsurface_fatal("not enough memory (%lu bytes requested): "
		"try a smaller -sample value", (unsigned long)size);
  return p;
}


void
fill_grid_sides(int ***tmp, int newx, int newy, int newz)
{
  int i, j, k;

  /* add points at 0 and newXX-1 */
  for (i=0; i < newz; i++)
    for (j=0; j < newy; j++)
      tmp[i][j][0] = tmp[i][j][newx-1] = 0;

  for (i=0; i < newz; i++)
    for (k=0; k < newx; k++)
      tmp[i][0][k] = tmp[i][newy-1][k] = 0;

  for (j=0; j < newy; j++)
    for (k=0; k < newx; k++)
      tmp[0][j][k] = tmp[newz-1][j][k] = 0;

  return;
}


/* a grid is 3 blocks: plane pointers, row pointers and the values */
int
dealloc_grid(int ***tmp, int nx, int ny, int nz)
{
  (void)nx; (void)ny; (void)nz;

  ckfree((char *)tmp[0][0]);
  ckfree((char *)tmp[0]);
  ckfree((char *)tmp);

  return 1;
}


int
***alloc_grid(int nx, int ny, int nz)
{
  int i, j, ***tmp, **rows, *vals;

  if ((double)nx*ny*nz*sizeof(int) > (double)INT_MAX)
    tclsurface_fatal("grid [%d,%d,%d] too large: try a smaller -sample value",
		nx, ny, nz);
  tmp  = (int ***)tclsurface_alloc(nz*sizeof(int **));
  rows = (int **)attemptckalloc(nz*ny*sizeof(int *));
  vals = (int *)attemptckalloc((size_t)nx*ny*nz*sizeof(int));
  if (rows == NULL || vals == NULL) {
    if (rows) ckfree((char *)rows);
    if (vals) ckfree((char *)vals);
    ckfree((char *)tmp);
    tclsurface_fatal("not enough memory for a [%d,%d,%d] grid: "
		"try a smaller -sample value", nx, ny, nz);
  }
  for (i=0; i<nz; i++) {
    tmp[i] = rows + i*ny;
    for (j=0; j<ny; j++)
      tmp[i][j] = vals + ((size_t)i*ny + j)*nx;
  }

  return tmp;
}


/* This subroutine calculates a normal from three vertices */
void
calc_normal(float p1[3], float p2[3], float p3[3], float n[3])
{
  float u[3], v[3];
  float sum, mag;

  u[0] = p3[0] - p2[0];
  u[1] = p3[1] - p2[1];
  u[2] = p3[2] - p2[2];

  v[0] = p1[0] - p2[0];
  v[1] = p1[1] - p2[1];
  v[2] = p1[2] - p2[2];

  n[0] = u[1] * v[2] - u[2] * v[1];
  n[1] = u[2] * v[0] - u[0] * v[2];
  n[2] = u[0] * v[1] - u[1] * v[0];

  sum = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
  mag = (float) sqrt((double) sum);

  if (mag == 0.0)
    mag = 1.0;

  n[0] = n[0] / mag;
  n[1] = n[1] / mag;
  n[2] = n[2] / mag;
}


/* This subroutine is for allocating memory for the vertex lists */
static void
vert_alloc(size_t size)
{
  float *v, *n;

  if (first_alloc) {
    VERTICES = (float *)tclsurface_alloc(size * 3);
    NORMALS  = (float *)tclsurface_alloc(size * 3);
    first_alloc = 0;
  } else {
    if ((v = (float *)attemptckrealloc((char *)VERTICES, size * 3)) == NULL)
      tclsurface_fatal("not enough memory to store vertices");
    VERTICES = v;
    if ((n = (float *)attemptckrealloc((char *)NORMALS, size * 3)) == NULL)
      tclsurface_fatal("not enough memory to store normals");
    NORMALS = n;
  }
}


/* This subroutine stores a polygon (triangle) in a list of vertices */
/* and connectivity.  This list can then be written out in different */
/* file formats. */
int
add_polygon(float *p1, float *p2, float *p3, float *n1, float *n2, float *n3)
{
  float *ptr;

  /* see if we have enough space to store the vertices */
  if (NUM_VERTICES >= (VERT_LIMIT - 3)) {
    /* get more space */
    VERT_LIMIT += VERT_INCR;
    vert_alloc((size_t)VERT_LIMIT * sizeof(float));
  }
  /* store the vertices */
  ptr = VERTICES + NUM_VERTICES * 3;
  *ptr++ = *p1++;		/* x of first vertex */
  *ptr++ = *p1++;		/* y of first vertex */
  *ptr++ = *p1++;		/* z of first vertex */
  *ptr++ = *p2++;		/* x of second vertex */
  *ptr++ = *p2++;		/* y of second vertex */
  *ptr++ = *p2++;		/* z of second vertex */
  *ptr++ = *p3++;		/* x of third vertex */
  *ptr++ = *p3++;		/* y of third vertex */
  *ptr++ = *p3++;		/* z of third vertex */

  ptr = NORMALS + NUM_VERTICES * 3;
  *ptr++ = *n1++; *ptr++ = *n1++; *ptr++ = *n1++;
  *ptr++ = *n2++; *ptr++ = *n2++; *ptr++ = *n2++;
  *ptr++ = *n3++; *ptr++ = *n3++; *ptr++ = *n3++;

  NUM_VERTICES += 3;

  return 0;
}


void
release_memory(void)
{
  first_alloc = 1;

  if (VERTICES != NULL) {
    ckfree((char *) VERTICES);
    VERTICES = NULL;
  }
  if (NORMALS != NULL) {
    ckfree((char *) NORMALS);
    NORMALS = NULL;
  }
  NUM_VERTICES=0;
  VERT_LIMIT=0;
}
