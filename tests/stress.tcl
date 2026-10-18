# Usage : tclsh stress.tcl <libdir>
# Contrôles globaux : orientation des normales, appels répétés, gros système.
lappend auto_path [lindex $argv 0]
package require tclsurface

# 1) sphère : normales sortantes et unitaires
lassign [tclsurface -sample 4 {{10 20 30 6 0}}] v n
set bad 0; set nv 0
foreach {x y z} $v {nx ny nz} $n {
    incr nv
    set d [expr {($x-10)*$nx+($y-20)*$ny+($z-30)*$nz}]
    if {$d <= 0 || abs(sqrt($nx*$nx+$ny*$ny+$nz*$nz)-1) > 1e-3} {incr bad}
}
puts "sphere: $nv sommets, normales mal orientées/non unitaires : $bad"

# 2) gros système pseudo-protéique (chaîne aléatoire compacte)
expr {srand(42)}
set atoms {}; set x 0; set y 0; set z 0
set elts {6 6 6 7 8 16}
for {set i 0} {$i < 3000} {incr i} {
    set x [expr {$x + (rand()-0.5)*3}]; set y [expr {$y + (rand()-0.5)*3}]; set z [expr {$z + (rand()-0.5)*3}]
    foreach c {x y z} { if {abs([set $c]) > 18} {set $c [expr {[set $c]*0.8}]} }
    lappend atoms [list $x $y $z [lindex $elts [expr {int(rand()*6)}]] $i]
}
foreach opts {{} {-sample 2} {-sample 3 -smooth} {-sample 2 -connolly} {-sample 2 -expand 1 -contract 1} {-sample 3 -request 20000} {-sample 3 -filter 3 -D 1.5}} {
    set t [time {set r [tclsurface {*}$opts $atoms]}]
    puts [format "%-32s %7d triangles  %s" "opts {$opts}" [expr {[llength [lindex $r 0]]/9}] [lindex $t 0]µs]
}

# 3) appels répétés (fuites / état résiduel)
set r1 [tclsurface -sample 2 $atoms]
for {set i 0} {$i < 5} {incr i} { tclsurface -sample 3 -request 500 -smooth $atoms ; catch {tclsurface -bad $atoms}; catch {tclsurface {{0 0 0 200 0}}} }
set r2 [tclsurface -sample 2 $atoms]
puts "résultat identique après appels répétés : [expr {$r1 eq $r2}]"
