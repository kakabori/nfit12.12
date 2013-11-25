proc get_S {spnum} {
	global step D x0 NKparams
	set step(0) 1
	set step(1) 1
	set step(2) 1
	set step(3) 1
	DSfitting 16 $spnum
	set NKparams(5) $D
	set NKparams(10) $x0
	set paramArray(D) $D
	set paramArray(bc2b) $x0
	#setparams
	#paraset update p
}

proc dspacing {} {
	
	global lambda S x0 D
	toplevel .ds
	wm title .ds "D Spacing"
	set w .ds.p
	frame $w
	label $w.r -text R(mm)
	entry $w.er -textvariable R -width 10
	label $w.l -text lambda(A)
	entry $w.el -textvariable lambda -width 10
	checkbutton $w.a -text alpha(rad) -variable xpin(2) -borderwidth 1
	entry $w.ea -textvariable alpha -width 10
	checkbutton $w.s -text S(mm) -variable xpin(0) -borderwidth 1
	entry $w.es -textvariable S -width 10
	checkbutton $w.f -text x0(mm) -variable xpin(1) -borderwidth 1
	entry $w.ef -textvariable x0 -width 10
	checkbutton $w.d -text D(A) -variable xpin(3) -borderwidth 1
	entry $w.ed -textvariable D -width 10
	label $w.p -text "pixelsize(mm)"
	entry $w.ep -textvariable pz -width 10
	label $w.c -text "chi^2"
	entry $w.ec -textvariable chi -width 10
	grid $w.s $w.es $w.r $w.er -sticky w
	grid $w.f $w.ef $w.l $w.el -sticky w
	grid $w.d $w.ed $w.p $w.ep -sticky w
	grid $w.a $w.ea $w.c $w.ec -sticky w
	#added by NC for ref. index 3.8.04
	label $w.ln -text "ref. index n"
	entry $w.en -textvariable refn -width 10
	grid $w.ln $w.en -sticky w
	#end of addition NC 3.8.04
	
	set w .ds.ls
	frame $w
	set w .ds.ls.l
	frame $w
	pack $w
	label $w.o -text orders
	label $w.x -text sample(a)
	label $w.t -text sample(b)
	grid $w.o $w.x $w.t
	for {set i 0} {$i<16} {incr i} {
		checkbutton $w.c$i -variable boolx($i) -text $i -padx 0 -pady 0 -borderwidth 1
		#	label $w.l$i -text $i
		entry $w.t$i -textvariable xv0($i) -width 10
		entry $w.y$i -textvariable xv2($i) -width 10
		grid $w.c$i $w.t$i $w.y$i
	}

	

	set w .ds.ls.b
	frame $w
	pack $w
	button $w.d0 -text fit(a) -command "get_S 0"
	button $w.d2 -text fit(b) -command "get_S 2"
	grid $w.d0 $w.d2
	
	graph .ds.g -width 350 -height 300
	.ds.g legend configure -hide yes
	.ds.g element create x -xdata plotxDS -ydata plotyDS -pixels 3
	
	wm resizable .ds 0 0
	
	pack .ds.ls -side right
	pack .ds.p -side top
	pack .ds.g
	bind .plot <Key> { catch {set xv0(%K) $plotx} }

	set .fit.top.par.bt.lc5 $D
}


#xv1 notify always

set R 0.00001
set alpha 0
set lambda 1.1775
set S 365.1
set x0 0
set D 58.367
set pz 0.07113
set refn 0.9999979

set xpin(1) 1
set xpin(3) 1
