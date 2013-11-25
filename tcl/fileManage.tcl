proc ReadPar {name} {
	global dataset dupe nocz NKfilter rotation tree node aFactor iter
	global paramArray paramNames freeOrFixed omega
	global noscale notiff kiyomask
	set types {
		{{par Files}       {.par}        }
		{{All Files}        *             }
	}

	if {$name == "" } {
		set name [tk_getOpenFile -initialdir "./" -filetypes $types]
  }
	if { $name != "" } {
		set parfilename $name 
		cd [file dirname $parfilename]
		set iname [file rootname [$tree label $node]]
		source $parfilename;
    #paraset update p
    #paraset disp p
    rotate $rotation $tree $node
    updateLambdaEta
    #NKupdateleta $NKparams(0) $NKparams(1) $NKparams(5) $NKparams(17)
    }     
}


proc ExportPar {dataset dupe nocz name } {
	global tree node NKfilter NK_lambda NK_eta NK_xi rotation iter aFactor
	global paramNames paramArray freeOrFixed
	global noscale notiff kiyomask
	if {[catch {set imgname [file join [$tree get $node filedir] [$tree label \
                                                              $node]]}] != 0 } {
		set imgname ""
  }
  set types {
		{{par Files}       {.par}        }
		{{All Files}        *             }
  }
  if {$name == "" } {
  	set name [tk_getSaveFile -initialdir "./" -filetypes $types]
  }
  set parfilename $name
  set xisquare [.fit.top.fit.bt.xi get]
  set rotname [file rootname [$tree label $node]]
	cd [file dirname $parfilename]
  if { $parfilename != "" } {
  	if { [file rootname $parfilename] == $parfilename } {
	  	set parfilename "$parfilename.par"
		}   
    set fid [open $parfilename w]
		writeout $fid $imgname $dataset $dupe $nocz
		puts $fid "# ChiSquare $xisquare"
		puts $fid "# lambda $NK_lambda"
		puts $fid "# eta $NK_eta"
		puts $fid "# xi $NK_xi"
		puts $fid "# Iteration $iter"
		puts $fid "# aFactor $aFactor"
		puts $fid "#"
		foreach i $paramNames {
			puts $fid "set paramArray($i) $paramArray($i)"
		}
		foreach i $paramNames {
    	puts $fid "set freeOrFixed($i) $freeOrFixed($i)"
    }
		puts $fid "set rotation \"$rotation\""
    puts $fid "set dataset \"$dataset\""
		puts $fid "set dupe $dupe"
		puts $fid "set nocz $nocz"
		puts $fid "set noscale $noscale"  
		puts $fid "set NKfilter $NKfilter"
		puts $fid "set iter $iter"
		puts $fid "set aFactor $aFactor"
		puts $fid "set notiff $notiff"
		close $fid
  }
}


proc writeout {fid imgname dataset dupe nocz} {
	global paramArray paramNames freeOrFixed NKfilter
	set fix { fixed free }
	set yesno { yes no }
  set noyes { no yes }
  set gap { "" \t \t \t \t \t "" \t \t \t \t \t \t "" \t "" "" }
  puts $fid "# NFIT version: nfit12.12"
  puts $fid "# [clock format [clock second]]"
  puts $fid "# image filename: $imgname"
	puts $fid "# filter = $NKfilter"
  puts $fid "# Kc\t[lindex $gap 0][lindex $fix $freeOrFixed(Kc)]\t$paramArray(Kc)\t dataset\t: $dataset"
  puts $fid "# B\t[lindex $gap 1][lindex $fix $freeOrFixed(B)]\t$paramArray(B)\t duplicates\t: $dupe"
  puts $fid "# Lr\t[lindex $gap 2][lindex $fix $freeOrFixed(Lr)]\t$paramArray(Lr)\t cz\t\t: [lindex $yesno $nocz]"
  set index 0
  foreach i $paramNames {
  	if {[string compare $i "Kc"] == 0 || [string compare $i "B"] == 0 || \
  	    [string compare $i "Lr"] == 0} {
  		# do nothing if $paramNames element is equal to Kc, Lr, or B
  	} else {
  		puts $fid "# $i\t[lindex $gap $index][lindex $fix $freeOrFixed($i)]\t$paramArray($i)"
  	}
  	incr index
  }	
}


proc openOutput { } {
	global output outputflag outfilename recordflag
    set types {
		{{otp Files}       {.otp}        }
		{{All Files}        *             }
    }
	set outfilename [tk_getSaveFile -initialdir "./" -filetypes $types]
    cd [file dirname $outfilename]
    if { $outfilename != "" } {
		if { [file rootname $outfilename] == $outfilename } {
			set outfilename $outfilename.otp
		}
		set output [open $outfilename w]
		set outputflag 1
		set recordflag 1
		#puts $output "[clock format [clock second]]"
	} else { 
		set outputflag 0
		set recordflag 0 
    }
}


proc closeOutput { } {
	global output outputflag
	if { $outputflag == 1 } { close $output }
	set outputflag 0
	set recordflag 0
}


proc recordfit { line } {
     global output outfilename flag
     puts $output $line
}
