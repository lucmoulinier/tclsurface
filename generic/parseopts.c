#include <stdio.h>
#include "polyr.h"


/* parses the command line options : objv[1] .. objv[objc-2]
   (the last argument holds the atoms).  All options are reset to
   their default value first, so that nothing leaks from one call to
   the next. */
int
parse_options(Tcl_Interp *interp, Tcl_Size objc, Tcl_Obj *const objv[])
{
  static const char *const options[] = {
    "-sample", "-smooth", "-connolly", "-expand", "-contract",
    "-request", "-d", "-D", "-filter", "-m", "-p", NULL
  };
  enum {
    OPT_SAMPLE, OPT_SMOOTH, OPT_CONNOLLY, OPT_EXPAND, OPT_CONTRACT,
    OPT_REQUEST, OPT_DOWNS3, OPT_DOWNS, OPT_FILTER, OPT_MAXPOLS, OPT_VERBOSE
  };
  static const int nargs[] = {1, 0, 0, 1, 1, 1, 3, 1, 1, 1, 0};
  Tcl_Size i;
  int j, idx, ival, founddowns = 0;
  double dval;

  /* defaults */
  SAMPLE    = 1;
  SMOOTH    = 0;
  CONNOLLY  = 0;
  NEXPAND   = 0;
  NCONTRACT = 0;
  REQPOLS   = 0;
  MAXPOLS   = 150;
  VERBOSE   = 0;
  DOWNSAMPL[0] = DOWNSAMPL[1] = DOWNSAMPL[2] = 1.0;
  FILTER[0] = FILTER[1] = FILTER[2] = 0;

  for (i = 1; i < objc - 1; i++) {
    if (Tcl_GetIndexFromObj(interp, objv[i], options, "option", 0, &idx)
	!= TCL_OK)
      return TCL_ERROR;
    if (i + nargs[idx] >= objc - 1) {
      Tcl_SetObjResult(interp, Tcl_ObjPrintf(
	  "missing value for option \"%s\"", options[idx]));
      return TCL_ERROR;
    }

    switch (idx) {
    case OPT_SAMPLE:
    case OPT_EXPAND:
    case OPT_CONTRACT:
    case OPT_REQUEST:
    case OPT_FILTER:
    case OPT_MAXPOLS:
      if (Tcl_GetIntFromObj(interp, objv[++i], &ival) != TCL_OK)
	return TCL_ERROR;
      switch (idx) {
      case OPT_SAMPLE:
	if (ival < 1 || ival > 100) goto range;
	SAMPLE = ival;
	break;
      case OPT_EXPAND:
	if (ival < 0) goto range;
	NEXPAND = ival;
	break;
      case OPT_CONTRACT:
	if (ival < 0) goto range;
	NCONTRACT = ival;
	break;
      case OPT_REQUEST:
	if (ival <= 0 || ival >= 1000000) goto range;
	REQPOLS = ival;
	break;
      case OPT_FILTER:
	if (ival < 3 || !(ival % 2)) {
	  Tcl_SetObjResult(interp, Tcl_ObjPrintf(
	      "-filter: size must be an odd integer >= 3, got %d", ival));
	  return TCL_ERROR;
	}
	FILTER[0] = FILTER[1] = FILTER[2] = ival;
	break;
      case OPT_MAXPOLS:
	if (ival < 0) goto range;
	MAXPOLS = ival;
	break;
      }
      break;

    case OPT_SMOOTH:
      SMOOTH = 1;
      break;

    case OPT_CONNOLLY:
      CONNOLLY = 1;
      break;

    case OPT_VERBOSE:
      VERBOSE = 1;
      break;

    case OPT_DOWNS3:
    case OPT_DOWNS:
      if (founddowns) {
	Tcl_SetObjResult(interp, Tcl_NewStringObj(
	    "options -d and -D can only be given once", -1));
	return TCL_ERROR;
      }
      founddowns = 1;
      for (j = 0; j < 3; j++) {
	if (idx == OPT_DOWNS3 || j == 0) {
	  if (Tcl_GetDoubleFromObj(interp, objv[++i], &dval) != TCL_OK)
	    return TCL_ERROR;
	  if (dval < 1.0) {
	    Tcl_SetObjResult(interp, Tcl_ObjPrintf(
		"downsampling ratio must be >= 1, got %s",
		Tcl_GetString(objv[i])));
	    return TCL_ERROR;
	  }
	}
	DOWNSAMPL[j] = (float)dval;
      }
      break;
    }
  }

  return TCL_OK;

 range:
  Tcl_SetObjResult(interp, Tcl_ObjPrintf(
      "value \"%s\" out of range for option \"%s\"",
      Tcl_GetString(objv[i]), Tcl_GetString(objv[i-1])));
  return TCL_ERROR;
}
