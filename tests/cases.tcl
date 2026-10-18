# Usage : tclsh cases.tcl <libdir> <case>
# Lance un cas de test isolé de tclsurface ; affiche un résumé du résultat.
lassign $argv libdir case
lappend auto_path $libdir
package require tclsurface

# petite molécule : quelques atomes C/N/O en chaîne, élément = numéro atomique
set mol {}
foreach {x y z na} {
    0.0 0.0 0.0 6   1.5 0.0 0.0 6   2.2 1.2 0.0 7   3.6 1.2 0.3 6
    4.3 2.4 0.3 8   4.2 0.0 0.6 6   5.6 0.0 0.9 8   0.0 -1.4 0.4 6
} { lappend mol [list $x $y $z $na 0] }
set atom1 {{10.0 20.0 30.0 6 0}}

proc summary {res} {
    lassign $res v n
    set nv [expr {[llength $v]/3}]
    set xs {}; set ys {}; set zs {}
    foreach {x y z} $v {lappend xs $x; lappend ys $y; lappend zs $z}
    if {$nv == 0} {return "0 vertices"}
    set c {}
    foreach l [list $xs $ys $zs] {
        set l [lsort -real $l]
        lappend c [format %.3f [expr {([lindex $l 0]+[lindex $l end])/2.0}]] \
                  [format %.3f [expr {([lindex $l end]-[lindex $l 0])/2.0}]]
    }
    return "[expr {$nv/3}] tris ; centre/demi-etendue x,y,z = $c"
}

set t0 [clock milliseconds]
set rc [catch {
    switch -- $case {
        basic     {tclsurface $mol}
        sphere    {tclsurface -sample 4 $atom1}
        sphere8   {tclsurface -sample 8 $atom1}
        expand    {tclsurface -sample 4 -expand 2 $atom1}
        contract  {tclsurface -sample 4 -expand 2 -contract 2 $atom1}
        connolly  {tclsurface -sample 2 -connolly $atom1}
        downs     {tclsurface -sample 4 -D 2 $atom1}
        filter4   {tclsurface -sample 2 -filter 4 $atom1}
        filter3   {tclsurface -sample 2 -filter 3 $atom1}
        request   {tclsurface -sample 4 -request 300 $mol}
        twice     {tclsurface -sample 4 -request 300 $mol; tclsurface -sample 2 $mol}
        persist   {tclsurface -sample 4 $atom1; tclsurface $atom1}
        unknown   {tclsurface -foo $atom1}
        noargs    {tclsurface}
        empty     {tclsurface {}}
        badatom   {tclsurface {{1 2 3}}}
        badnum    {tclsurface {{a b c 6 0}}}
        bigZ      {tclsurface {{0 0 0 150 0}}}
        vdw0      {tclsurface {{0 0 0 10 0}}}
        OandC     {tclsurface -sample 4 {{0 0 0 8 0} {3 0 0 6 0}}}
        faratoms  {tclsurface {{1500 0 0 6 0} {1502 0 0 6 0}}}
        default   {error "cas inconnu $case"}
    }
} res]
set dt [expr {[clock milliseconds]-$t0}]
if {$rc} { puts "$case: ERREUR ($dt ms) : [string range $res 0 150]" } else { puts "$case: OK ($dt ms) : [summary $res]" }
