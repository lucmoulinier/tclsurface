# tclsurface

Tcl extension that computes a triangulated molecular surface from a list of atoms.
Package `tclsurface` 0.1. It builds against Tcl 8.6 and Tcl 9 from the same sources.

```tcl
package require tclsurface
lassign [tclsurface ?options? $atoms] vertices normals
```

Each atom is `{x y z atomic_number ?tag?}`. The result is two flat lists with 9 values per triangle
(vertices in Å, unit outward normals).

| Option | Effect |
|---|---|
| `-sample n` | grid points per Å, 1..100 |
| `-smooth` | smooth the mask |
| `-connolly` | expand the mask by a 1.4 Å probe |
| `-expand n` / `-contract n` | grow / shrink the mask n times |
| `-request n` | target number of triangles |
| `-filter n` | box filter (odd n ≥ 3) |
| `-d x y z` / `-D n` | downsampling ratios |
| `-m n` | remove connected pieces with ≤ n triangles (default 150) |
| `-p` | verbose output on stderr |

## Build

TEA (autoconf) build, out of the source tree:

```sh
mkdir build-tcl8 && cd build-tcl8
../configure --with-tcl=/path/to/tcl8.6/lib && make     # -> libtclsurface0.1.so
mkdir build-tcl9 && cd build-tcl9
../configure --with-tcl=/path/to/tcl9/lib && make       # -> libtcl9tclsurface0.1.so
make install
```

On Windows, build with MinGW/MSYS2 the same way.
`configure` is already generated; run `autoconf` only after editing `configure.ac`.

## Tests

```sh
tclsh8.6 tests/cases.tcl build-tcl8 basic      # one case (basic, sphere, expand, downs, badatom...)
tclsh9.0 tests/stress.tcl build-tcl9           # every option on a 3000-atom system, repeated calls
```

## License

GPL-3.0-or-later, Copyright (c) Luc Moulinier. See [LICENSE](LICENSE).
