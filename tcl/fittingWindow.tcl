# this procedure prepares the fitting panel
proc fitpanel {} {
	global setnfitt tree node xisquare aFactor iter parameters dataset \
	dupe nocz output outputflag rotation lambda S x0 D notiff
	global paramArray freeOrFixed paramNames
	global NK_lambda NK_eta NK_xi
	global noscale
	set dataset ""
	set iter 1
	set notiff 0
	set noscale 1
	set kiyomask 0
		
	updateLambdaEta
	
	set dataset "1 520 620 150 250"
	set rotation "0 512 512"
		
	toplevel .fit
	wm geometry .fit =+570+480
	wm title .fit "Fitting"
	set w .fit.menu
	frame $w
  set m $w.file.m
  menubutton $w.file -text File -menu $m
  menu $m -tearoff 0
  $m add command -label "Import" -command {ReadPar "" }
  $m add command -label "Export" -command {ExportPar $dataset $dupe $nocz ""}
	$m add command -label "Fit All" -command {FitAll $dataset}
	#Completely stop a "Fit All" run
	$m add command -label "Stop All" -command StopAll
	#Modify certain fitting options for special purposes
	$m add command -label "Options" -command fitoptwin
	#Show the elapsed time since the beginning og the most recent fitting run
	$m add command -label "Elapsed Time" -command timewin
	#Opens up a window to set certain miscelleous variables for the theory
	#$m add command -label "Misc." -command setMiscVars
	pack $w.file -side left -anchor w
	
	set w .fit.dst
	frame $w
	label $w.dt -text "dataset: sig_back qr1 qr2 qz1 qz2" -width 50 -anchor w
	entry $w.de -textvariable dataset -width 35
	bind $w.de <Return> {setdata $dataset}

	label $w.dt1 -text "rotate: tg x0 y0" -anchor w
  entry $w.de1 -textvariable rotation -width 20
  bind $w.de1 <Return> {rotate $rotation $tree $node}
	pack $w.dt $w.de $w.dt1 $w.de1 -side top
	frame .fit.top
	set w .fit.top.par
	frame $w
	set w .fit.top.par.bt
	frame $w
	
	# Set up fields for fitting parameters
	foreach i $paramNames {
		label $w.label$i -text "$i" -width 10
		checkbutton $w.button$i -variable freeOrFixed($i)
		entry $w.entry$i -textvariable paramArray($i) -borderwidth 1 -width 10
		grid $w.label$i $w.button$i $w.entry$i
	}
	foreach i $paramNames {
		#bind $w.entry$i <Return> {setparams; paraset update p; updateLambdaEta}
		bind $w.entry$i <Return> {updateLambdaEta}
	}
	
	#setparams
	
#	bind .fit <F2> {
#      setparams
#	    paraset update p
#	    NKupdateleta $NKparams(0) $NKparams(1) $NKparams(5) $NKparams(17)
#	}
#	bind .fit <F5> {
#	    paraset update p
#	    NKupdateleta $NKparams(0) $NKparams(1) $NKparams(5) $NKparams(17)
#	}
	set w .fit.top
	pack $w.par.bt -side top -expand yes -fill both
	
	set w .fit.top.fit
	frame $w
	
	set w $w.dat
	frame $w
	label $w.e -text ""
	label $w.t -text "NFIT" -background white
	pack $w.t -side top
		
	set w .fit.top.fit.opt
	frame $w
	label $w.i -text "iterations"
	entry $w.ie -textvariable iter -width 5
	grid $w.ie $w.i
	
	label $w.ldp -text "duplicates"
	entry $w.dp -textvariable dupe -width 5
	grid $w.dp $w.ldp
	
	label $w.omega_label -text "omega (deg)"
	entry $w.omega_entry -textvariable omega -width 5
	grid $w.omega_entry $w.omega_label	
	
	label $w.lcz -text "cz"
	bind $w.lcz <1> { subcz $tree $node }
	checkbutton $w.cz -variable nocz -offvalue 1 -onvalue 0
	grid $w.cz $w.lcz
	
	label $w.noscale_label -text "noscale"
	checkbutton $w.noscale_button -variable noscale 
	grid $w.noscale_button $w.noscale_label
	
	set w .fit.top.fit.out
	frame $w
	checkbutton $w.op -indicatoron 0 -text open -borderwidth 3 -variable outputflag \
	                  -command { openOutput }
	label $w.l -text output
	button $w.sv -text close -borderwidth 3 -command { closeOutput }
	pack $w.op $w.l $w.sv -side left
	
	set w .fit.top.fit.bt
	frame $w
	checkbutton $w.fit -indicatoron 0 -text nfit -borderwidth 3 -variable startflag \
	       -command { startFit $dataset $dupe $nocz $iter $tree $node }
	entry $w.xi -textvariable xisquare -relief flat -width 6
	checkbutton $w.stop -indicatoron 0 -text stop -borderwidth 3 -variable stopflag
	label $w.e1 -text " "
	label $w.e2 -text " "
	pack $w.fit $w.e1 $w.xi $w.e2 $w.stop -side left    

		
	set w .fit.top.fit.fit
	frame $w
	label $w.f -text FIT -relief flat
	button $w.tog -text tog -borderwidth 3 -command { img tog FIT }
	button $w.togall -text "tog all" -borderwidth 3 -command { fitload $ldxdy $tree 1 }
	button $w.exp -text export -borderwidth 3 -command { exportFit FIT }
	bind $w.f <1> { img show FIT }
	pack $w.f $w.tog $w.togall $w.exp -side left
	
	set w .fit.top.fit.ref
	frame $w
	label $w.l -text refine
	entry $w.e -textvariable NKfilter -width 5
	button $w.b -text show -borderwidth 3 -command { filterdata }
	pack $w.l $w.e $w.b -side left
	
	# moved to the conf window
	set w .fit.top.fit.ldf
	frame $w
	button $w.b -text "LoadFit" -command { fitload $ldxdy $tree }
	#label $w.l -text " dx dy:"
	label $w.l -text " dx:"
	entry $w.e -textvariable ldxdy -width 5
	pack $w.b $w.l $w.e -side left
	####################
	
	set w .fit.top.fit.lbg
	frame $w
	label $w.t -text "Light Background" -background white
	frame $w.n
	label $w.n.l -text "(2 boxes) order:"
	entry $w.n.e -textvariable ord -width 5
	button $w.n.b -text fit -borderwidth 3 -command { fitbckg b1 b2 $ord 4 1 $tree $node }
	pack $w.n.l $w.n.e $w.n.b -side left
	frame $w.m
	label $w.m.l -text "box name:"
	entry $w.m.e -textvariable maskbox -width 5
	button $w.m.b -text mask -borderwidth 3 -command { maskpix $maskbox $tree $node }
	pack $w.m.l $w.m.e $w.m.b -side left
	
	pack $w.t $w.n $w.m -side top
	
	set w .fit.top.fit.stat
	frame $w
	label $w.t -text "Statistics" -background white
	frame $w.st
	label $w.st.l -text "box name:"
	entry $w.st.b -textvariable sbox -width 5
	button $w.st.stat -text stat -borderwidth 3 -command { NKstatist $sbox $tree $node  }
	pack $w.st.l $w.st.b $w.st.stat -side left
	frame $w.res
	label $w.res.l -text "sig_back:"
	entry $w.res.st -textvariable STresult -state disabled -width 8 -relief flat

  label $w.res.l1 -text "aFactor:"
  entry $w.res.st1 -textvariable aFactor -width 5
	pack $w.res.l $w.res.st $w.res.l1 $w.res.st1 -side left
	pack $w.t $w.st $w.res -side top
	
	set w .fit.top.fit
	pack $w.dat $w.opt $w.out $w.bt $w.fit $w.ref $w.lbg $w.stat -side top
	
	set w .fit.top
	pack $w.par $w.fit -side left
	pack $w.fit -side right -anchor nw

	set w .fit.leta
	frame $w
	label $w.llam -text "lambda="
	entry $w.lam -textvariable NK_lambda -state disabled -width 10 -relief flat
	label $w.leta -text "eta="
	entry $w.eta -textvariable NK_eta -state disabled -width 10 -relief flat
	label $w.lxi -text "xi="
	entry $w.xi -textvariable NK_xi -state disabled -width 10 -relief flat
	pack $w.llam $w.lam $w.leta $w.eta $w.lxi $w.xi -side left
	bind $w.llam <1> {updateLambdaEta}
	
	set w .fit.conv
	frame $w
	entry $w.pix -textvariable pixel -borderwidth 1 -width 15
	label $w.lpix -text "pixels  <==>  "
	entry $w.q -textvariable qu -borderwidth 1 -width 15
	label $w.lq -text "A-1"
	pack $w.pix $w.lpix $w.q $w.lq -side left
	#bind $w.pix <Return> { set qu [expr (4*3.1415)/$NKparams(13)* \
	#sin(0.5*atan($pixel*$NKparams(14)/$NKparams(9)))] }
	#bind $w.q <Return> { set pixel [expr $NKparams(9)/$NKparams(14)* \
	#tan(2*asin($qu*$NKparams(13)/(4*3.1415)))] }
	
	set w .fit
	pack $w.menu $w.dst $w.top $w.leta $w.conv -side top -expand yes -fill both
}

proc setparams {} {
	global paramArray
	#paraset set p $paramArray(Kc) $paramArray(B) $paramArray(Lr) $paramArray(Mz) $paramArray(D) \
	#              $paramArray(mosaic) $paramArray(edisp) $paramArray(beamFWHM) \
	#              $paramArray(sDistance) $paramArray(bc2b) $paramArray(wavelength) \
	#              $paramArray(pixelSize) $paramArray(qxzero) $paramArray(refractiveIndex) \
	#              $paramArray(T)
}

proc updateLambdaEta {} {
	global NK_lambda NK_eta NK_xi
	global paramArray
	set NK_lambda [expr 1e16 * sqrt($paramArray(Kc)/$paramArray(B)) / $paramArray(D)]
	set NK_eta [expr 2.17 * ($paramArray(T) + 273) / \
	            sqrt($paramArray(Kc) * $paramArray(B)) / $paramArray(D) / $paramArray(D)]
	set NK_xi [expr sqrt($paramArray(D) * $NK_lambda)]
}
