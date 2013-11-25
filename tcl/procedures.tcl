proc FitAll {dataset} {
	global NKparams tree node nodeMax setnfitt NKfilter parfilename NK_lambda NK_eta
	global NK_xi name savedataset fittedhash stopflag stopallflag notiff
	global outputflag output outfilename recordflag
	global dupe nocz
	set stopallflag 0
	set nfitbutton ".fit.top.fit.bt.fit"
	for {set i 1} { $i <= $nodeMax } {incr i} {
		if { [$tree exists $i] != 0} {
			set name [file rootname [$tree label $i]]
			if { [catch {set contents [glob $name*.par]}] } {
				continue
			}
			foreach item $contents {
				if { $stopallflag == 1 } {
					break
				}
				set parfilename $item
				set name $item
				set node $i
        		puts "Fitting $name"
				ReadPar $name
				if { $outputflag == 1} { 
				    puts $output "\n******************************$name******************************" 
				}
				$nfitbutton invoke
				set dataset $savedataset
				puts "\nFitted $name with $dataset\n"
				ExportPar $dataset $dupe $nocz $name
				set dir [file dirname [$tree label $i]]
				if { $notiff == 0 } {
					set fittifname [file rootname $name]
					set fittifname ${fittifname}_fit.tif
        			puts "img export FIT $fittifname"
        			img export "FIT" $fittifname
        			set fittedhash($fittifname) $dataset
        			puts "Fit Exported\n"
				}
				# Create .dat file name
				set parts [split $name .]
				set length [llength $parts]
				set index [expr $length - 2]
				set datName [join [lrange $parts 0 $index] .]
				append datName .dat
				exec mv frm.dat $datName 
			}
		}
		if { $stopallflag == 1} {
			break
		}
	}	
	puts "DONE PROCESSING\n"
	tk_messageBox -message "Done fitting" -type ok
}


#Procedure to stop a "Fit All" run as soon as possible
proc StopAll {} {
	global stopallflag stopflag
	set stopallflag 1
	set stopflag 1
}

proc updatePara {} {
	global paramNames paramArray
	foreach i $paramNames {
		updateParaObject $i $paramArray($i)
		#puts "updateParaObject $i $paramArray($i)"
	}
}

# call setNFIT c++ function if at least one parameter is set free
# otherwise, flag will be set to zero, setNFIT is not called, and 
# the procedure returns TCL_ERROR 
proc setnfit {} {
	global freeOrFixed paramNames
	set index 0
	set flag 0
	foreach i $paramNames {
		if {$freeOrFixed($i) == 1} {
			set free($i) $index
			set flag 1
		} else {
			set free($i) ""
		}
		incr index
	}
	if {$flag == 1} {
		#puts "setNFIT $free(Kc) $free(B) $free(Lr) $free(Mz) $free(D) $free(mosaic) $free(edisp) \
		 #     $free(bFWHM) $free(sdistance) $free(bc2b) $free(wavelength)\
		 #     $free(pixelSize) $free(qxzero) $free(refractiveIndex) $free(T)"
		setNFIT $free(Kc) $free(B) $free(Lr) $free(Mz) $free(D) $free(mosaic) $free(edisp) \
	          $free(bFWHM) $free(s) $free(bc2b) $free(wavelength)\
	          $free(pixelSize) $free(qxzero) $free(nindex) $free(T)
		return TCL_OK
	} else {
		return TCL_ERROR
	}
}

proc startFit {dataset dupe nocz niter tree node} {
	global output startflag flag
	set datasetsave $dataset
	set flag 1
	
	updatePara
	if {[setnfit] == "TCL_ERROR"} {
		puts "\nERROR: at least, one variable must be set free\n"
		set startflag 0
		return TCL_ERROR
	}
	
	#if { [catch { puts $output "" }] == 0 } {
  #	set imgname [file join [$tree get $node filedir] [$tree label $node]]
	#  writeout $output $imgname $parameters $dataset $dupe $nocz
	#  puts $output "niter : $niter"
	#}
  #puts "paraset disp p"      
  #paraset disp p
  setdata $dataset
  #puts "fitdata d$node mc p $niter"
  fitdata d$node mc p $niter
      
  #every 10000 update
  vwait startflag
}

proc NKupdateleta {NK_Kc NK_B NK_D NK_T} {
    global NK_lambda NK_eta NK_xi
    set NK_lambda [expr 1e16*sqrt($NK_Kc/$NK_B)/$NK_D]
    set NK_eta [expr 2.17*($NK_T+273)/sqrt($NK_Kc*$NK_B)/$NK_D/$NK_D]
    set NK_xi [expr sqrt($NK_D*$NK_lambda)]
}


# Set fitting options to newly specified values
proc setFitOpts { newKiyoMask newNoScale newNoTiff newNThreads } {
	global kiyomask noscale notiff nthreads
	set kiyomask $newKiyoMask
	set noscale $newNoScale
	set notiff $newNoTiff
	set nthreads $newNThreads
}

# Validate input giving the number of threads for the chi squared calculation
proc checkNThreads {} {
	global tempNThreads;
	if { ! [ string is int $tempNThreads ]  || ($tempNThreads <= 0)} {
		set tempNThreads 1
		tk_messageBox -message "The number of threads must be a positive integer."
	}
}

# Set up a dialogue to modify the noscale parameter, determine whether to
# supress TIFF output in "Fit All" runs, and determine the number of threads
# to use in chi squared calculations.
proc fitoptwin {} {
	global tempNThreads nthreads
	set tempNThreads $nthreads
	set tempNoScale 1
	set tempNoTiff 0
	set tempKiyoMask 0
	
	toplevel .fitopts
	wm title .fitopts "Fitting Options"

	set optframe [ frame .fitopts.optframe ]

	#Set up a label and check box to turn kiyomask on or off
	label $optframe.km -text "KiyoMask"
	checkbutton $optframe.kmask -variable tempKiyoMask
	grid $optframe.kmask $optframe.km

	#Set up a label and check box to turn scaling of theoretical values on
	#or off
	label $optframe.nsc -text "No Scale"
	checkbutton $optframe.nscale -variable tempNoScale
	grid $optframe.nscale $optframe.nsc

	#Set up a label and check box to turn suppression of TIFF output in a
	#"Fit All" run on or off
	label $optframe.tiff -text "Suppress TIFF Output"
	checkbutton $optframe.tiffbutton -variable tempNoTiff
	grid $optframe.tiffbutton $optframe.tiff

	grid $optframe -row 0

	#Set up an entry field to specify the number of threads used in a chi squared calculation
	set tframe [ frame .fitopts.tframe ]
	label $tframe.tlabel -text "Number of Threads for Chi Squared Calculation"
	entry $tframe.tentry -textvariable tempNThreads -width 4

	#A return keystroke tiggers a validation function to assess input. If
	#the input is a positive integer, accept it and update the display.
	bind $tframe.tentry <Return> {	checkNThreads }
	pack $tframe.tlabel $tframe.tentry
	grid $tframe -row 1

	#Set up a Cancel button to abort changes, and an Apply button to keep changes
	set bframe [ frame .fitopts.bframe ]
	button $bframe.cancel -text Cancel  -borderwidth 1 -command { destroy .fitopts }
	button $bframe.keep -text Apply -borderwidth 1 -command { setFitOpts $tempKiyoMask $tempNoScale $tempNoTiff $tempNThreads; destroy .fitopts }
	grid $bframe.cancel -row 0 -column 0
	grid $bframe.keep -row 0 -column 1
	grid $bframe -row 2
}

# Set up a window to display elapsed time since the beginning of the most
# recent fitting run.
proc timewin {} {
	global elapsedMins elapsedSecs

	toplevel .timer
	wm title .timer "Timer"
	label .timer.description -text "Elapsed time for the current fitting run:"
	grid .timer.description -row 0
	set infopanel [ frame .timer.ipanel ]
	label $infopanel.timer1 -textvariable elapsedMins
	label $infopanel.timer2 -text " minutes and "
	label $infopanel.timer3 -textvariable elapsedSecs
	label $infopanel.timer4 -text " seconds"
	grid $infopanel.timer1 -row 0 -column 0
	grid $infopanel.timer2 -row 0 -column 1
	grid $infopanel.timer3 -row 0 -column 2
	grid $infopanel.timer4 -row 0 -column 3
	grid $infopanel -row 1

	set bframe [ frame .timer.bframe ]
	button $bframe.close -text Close -borderwidth 2 -command { destroy .timer }
	pack $bframe.close
	grid $bframe -row 2
}

