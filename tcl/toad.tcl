package require BLT
namespace import blt::*
wm geometry . =500x430+0+0
wm title . "Toad"

set zoomMod "Control-"
set colors {#ff0000 #00ff00 #0000ff #00ffff #ff00ff #000000 #ff8800 #88ff00 #8800ff \
            #aaaaaa #0000ff #0000ff}
set bxtopt " -anchor se -padx 0 -pady 0 -bindtags edbox -fg #ffff00 -bg \"\" "
set bxlopt " -outline blue "
set descending 0
set verticalfit 0

source tcl/globalVariables.tcl

global rotatehash startflag fittedhash iter aFactor FLICAM

set startflag 0
set xyplot .plot.xy
set imgWN .imgWin

set colors {#ff0000 #00ff00 #0000ff #00ffff #ff00ff #000000 #ff8800 #88ff00 #8800ff #aaaaaa}

proc gPopt {iname} {
	if {$iname == "FIT"} {
  	return "-symbol circle -linewidth 0 -pixels 1.5"
  } else {
    return " -symbol square -linewidth 0 -pixels 1.5"
  }
}
set nodeMax 0
proc cPopt {iname} {
    global tree nodeMax colors
    set turn 0
    for {set i 1} { $i <= $nodeMax } {incr i} {
       if { [catch {$tree get $i load}]==0 } {
          set name [file rootname [$tree label $i]]
       } else {
 	  set name ""
       }
       if { $iname == $name } {
          set turn [expr ($i-1)]
       }
    }
    set turn [expr ($turn)%10]
    set c [lindex $colors $turn]
    if  { $iname == "FIT" } {
       set c #ff00aa
    }
    return " -color $c"
}


graph $imgWN -width 1024 -height 1024 -plotrelief flat -plotpadx {0 0} \
             -plotpady {0 0} -plotborderwidth 0 -plotbackground #888888

$imgWN axis create fixy -hide yes -min 0 -max 1
$imgWN axis configure y -subdivisions 5 -max 2048 -min 0 -descending $descending
$imgWN axis configure x -subdivisions 5 -max 2048 -min 0
$imgWN legend configure -hide yes
$imgWN marker create image -name bg
foreach graph  $imgWN  {
	Blt_ZoomStack $graph
	Blt_Crosshairs $graph
}

set w .menu
frame $w
set m $w.file.m
menubutton $w.file -text File -menu $m
menu $m -tearoff 0
$m add command -label "Print Plot" -command "Blt_PostScriptDialog $xyplot"
$m add command -label "Print Image" -command "Blt_PostScriptDialog $imgWN"
set m $w.tls.m
menubutton $w.tls -text Tools -menu $m
menu $m -tearoff 0
$m add command -label "Manipulation" -command {manipwin}
$m add command -label "D Spacing" -command {dspacing}
$m add command -label "Fitting" -command {fitpanel}
$m add command -label "Calculator" -command {calc}
$m add command -label "Color" -command {colorwin}
set m $w.win.m
menubutton $w.win -text Win -menu $m
menu $m -tearoff 0
$m add command -label "Configuration Window" -command {confwindow}
$m add command -label "Plot Window" -command {plotpanel}
label $w.lsw -text "swath:"
entry $w.xswath -textvariable tvar_xswath -width 5
pack $w.file $w.tls $w.win -side left
pack $w.xswath $w.lsw -side right

label $w.lx -text "x:" -fg red
label $w.xpos -textvariable tvar_xpos -width 5 -anchor w
label $w.ly -text "y:" -fg red
label $w.ypos -textvariable tvar_ypos -width 5 -anchor w
label $w.li -text "I:" -fg red
label $w.inte -textvariable tvar_inte -width 7 -anchor w
pack $w.inte $w.li $w.ypos $w.ly $w.xpos $w.lx -side right

set w .butt
frame $w
checkbutton $w.b1 -indicatoron 0 -bitmap @bmp/box.bmp -borderwidth 3 -variable editmode \
                  -command {editbox}
$w.b1 configure -selectcolor [$w.b1 cget -background]
pack $w.b1 -side top

proc miniShell { w } {
        global message
        entry $w -textvariable message -width 100
        bind $w <Return> {
                eval $message
                spincur 1
                set cmdH($curCmd) $message
                set message {}
                set selCmd $curCmd
                incr selCmd
        }
        bind $w <Key-Up> {
                spinsel -1
                set message $cmdH($selCmd)
        }
        bind $w <Key-Down> {
                spinsel 1
                set message $cmdH($selCmd)
        }
}
frame .mini
miniShell .mini.shell
pack .mini.shell -side left

pack .mini -side bottom -expand yes -fill both
pack .menu -side top -expand yes -fill both
pack .butt -side right -expand yes -fill both
pack $imgWN -side top -expand yes -fill both

$imgWN marker bind edbox <3> {
    set t [%W marker get current]
    set t [string range $t 1 end] 
    boxv del $t
}

bind pcs <1> {plotcs [%W axis invtransform x %x] [%W axis invtransform y %y]}
bind pcs <B1-Motion> {plotcs [%W axis invtransform x %x] [%W axis invtransform y %y]}

##############################
# %s  for events
# 0       <1>
# 256     <B1-Motion>
##############################
bind edbox <1> {
    boxv halo [$imgWN axis invtransform x 0] [$imgWN axis invtransform x 2] [$imgWN axis invtransform y 0] [$imgWN axis invtransform y 2]
    boxv 1st [%W axis invtransform x %x] [%W axis invtransform y %y]
}
bind edbox <B1-Motion> {
    boxv 2nd [%W axis invtransform x %x] [%W axis invtransform y %y]
}
bind edbox <ButtonRelease-1> {
    boxv 3rd
}

set temp "[bind zoom-$imgWN [lindex [bind zoom-$imgWN] 1]] \nbreak"
bind zoom-$imgWN [lindex [bind zoom-$imgWN] 1] $temp
set temp "[bind zoom-$imgWN [lindex [bind zoom-$imgWN] 0]] \nbreak"
bind zoom-$imgWN [lindex [bind zoom-$imgWN] 0] $temp
unset temp

bindtags $imgWN "pcs crosshairs-$imgWN zoom-$imgWN $imgWN"
proc editbox {} {
    global imgWN editmode
    if {$editmode} {bindtags $imgWN "zoom-$imgWN $imgWN edbox"} else {bindtags $imgWN "pcs crosshairs-$imgWN zoom-$imgWN $imgWN"}
}

proc axisshift {w d a} {
	set limits [$w axis limits $a]
	set min [lindex $limits 0]
	set max [lindex $limits 1]
	set dif [expr $max-$min]
	set max [expr $max+$dif/$d]
	set min [expr $min+$dif/$d]
	$w axis configure $a -min $min -max $max
}

proc showall {i} {
    global imgWN
    $imgWN  axis configure x -min [lindex [tv::img coords $i] 0] -max [lindex [tv::img coords $i] 2]
    $imgWN  axis configure y -min [lindex [tv::img coords $i] 1] -max [lindex [tv::img coords $i] 3]
}

proc showimg { i } {
    global imgWN
    catch "$imgWN marker delete bg"
    tv::img show $i
    $imgWN marker create image -name bg -image realimg0 -coords [tv::img coords $i] -under 1
}


namespace eval ::tv namespace export *
namespace import tv::*

proc plotxy {w x y} {
    global plotx ploty
    set ploty [$w axis invtransform y $y]
    set plotx [$w axis invtransform x $x]
}

proc swathwidth {} {
	global imgWN tvar_xswath
	set y1 [$imgWN axis transform y 0]
	set y2 [$imgWN axis transform y 2048]
	set y1 [expr int(abs(($y1-$y2)/2048.0*$tvar_xswath))]
	$imgWN crosshairs configure -linewidth $y1
}

bind $imgWN <Enter> swathwidth

proc load1k { n f } {
    img load $n $f
    img trans $n
    img flip $n x
    img show $n
}

proc xxxls {} {
    set f [glob *]
    set len [llength $f]
    for {set i 0} { $i < $len } {incr i} {
	puts [lindex $f $i]
    }
}
    
proc loadfit {name x y} {
    img load fit $name
    img shift fit $x $y
    img show fit
    img tog fit
}

set message ""
set curCmd 0
set selCmd 0

proc spincur { s } {
        global curCmd
        incr curCmd $s
        set curCmd [expr $curCmd % 50]
}
proc spinsel { s } {
        global selCmd
        incr selCmd $s
        set selCmd [expr $selCmd % 50]
}

proc ctrl {w opt} {
    miniShell $w.minis
    eval pack $w.minis $opt
}

proc phide {name} {
	global xyplot
	$xyplot element configure " $name" -hide yes
}

proc pshow {name} {
	global xyplot
	$xyplot element configure " $name" -hide no
}

proc tvopen {f} {
  global fileExtension
  set iname [file tail [file rootname $f]]
  puts "path: $f"
  puts "name: $iname"
  puts "ext: $fileExtension"
  if {$fileExtension == ".dat"} {
    readDataFromASCII $iname $f
  } else { 
    img load $iname $f
  }
}

proc SortColumn { column } {
	global tw
	set old [$tw.t sort cget -column] 
  set decreasing 0
  if { "$old" == "$column" } {
		set decreasing [$tw.t sort cget -decreasing]
		set decreasing [expr !$decreasing]
  }
  $tw.t sort configure -decreasing $decreasing -column $column -mode integer
  $tw.t configure -flat yes
  $tw.t sort auto yes

  blt::busy hold $tw.t
  update
  blt::busy release $tw.t
}

proc flip {opt tree node} {
	set iname [file rootname [$tree label $node]]
  img flip $iname $opt
  img show $iname
  if { [$tree get $node plot] } {
  	img tog $iname 1
  }
}

proc Open {tree {automate 0} {name ""} } {
  global node nodeMax
  global fileExtension
  set types { {{tif Files} {.tif}} {{ASCII Files} {.dat}} {{All Files} *} }
  if { $automate == "0" } {
    set tiffilename [tk_getOpenFile -initialdir "./" -filetypes $types]
    cd [file dirname $tiffilename]
    if { $tiffilename != "" } {
	    set name [file tail $tiffilename]
	    set fileExtension [file extension $tiffilename]
	    set directory [file dirname $tiffilename]
	    set max $nodeMax
	    set nameflag 0
	    set iname ""
	    for {set i 1} {$i <= $max} {incr i} {
	      catch { set iname [$tree label $i] }
	      if { $name==$iname } { set nameflag 1 }
      }
	    if {$nameflag == 0} {
	      incr nodeMax
	      # The following line calls img load
	      set node [$tree insert root -label $name -data "filedir $directory load 1 plot 1 note --"]
	      set iname [file rootname [$tree label $node]]
	      set fname [file join [$tree get $node filedir] [$tree label $node]]
	      #puts "img load $iname $fname"
	      treeshow $tree $node
	    } else {
	      puts "file already open"
	    }
    }
  } else {
	    set tiffilename $name
     if { $tiffilename != "" } {
        set name [file tail $tiffilename]
        set directory [file dirname $tiffilename]
        set max $nodeMax
        set nameflag 0
        set iname ""
        for {set i 1} {$i <= $max} {incr i} {
          catch { set iname [$tree label $i] }
          if { $name==$iname } { set nameflag 1 }
        }
        if {$nameflag == 0} {
           incr nodeMax
           set node [$tree insert root -label $name -data "filedir $directory load 1 plot 1 note --"]
           set iname [file rootname [$tree label $node]]
           set fname [file join [$tree get $node filedir] [$tree label $node]]
           puts "img load $iname $fname"
           treeshow $tree $node
        } else {
           puts "file already open"
        }
     }

  }


}

proc Reopen {tree node } {
	set fname [file join [$tree get $node filedir] [$tree label $node]]
	tvopen $fname
	set iname [file rootname [$tree label $node]]
	treeshow $tree $node
	set iname [file rootname [$tree label $node]]
  if { [$tree get $node plot] } {
		img tog $iname 1
  }
}

proc Close {tree node} {
     $tree set $node label ""
     $tree set $node load 0
}

proc GetAbsolutePath { dir } {
    set saved [pwd]
    cd $dir
    set path [pwd] 
    cd $saved
    return $path
}

set tree {}

proc treeshow {tree node} {
    global tw
    img show [file rootname [$tree label $node]]
    blt::tv::SetSelectionAnchor $tw.t $node
}

proc treeplot {tree node key ops} {
	if { [$tree get $node plot] } {
		if { [catch {img tog [file rootname [$tree label $node]] 1}] } {
	  	$tree set $node plot 0
		}
  } else {
		catch {img tog [file rootname [$tree label $node]] 0}
  }
}

proc treeload {tree node key ops} {
  set fname [file join [$tree get $node filedir] [$tree label $node]]
  set iname [file rootname [$tree label $node]]
  if { [$tree get $node load] == 1 } {
    tvopen $fname
  } else {
    img del $iname
    $tree del $node
  }
}


proc Relb04 {tree node } {
	set fname [file join [$tree get $node filedir] [$tree label $node]]
	tvopen $fname
  set iname [file rootname [$tree label $node]]
  img rotate $iname inf 512 512
  img flip $iname x
  img rotate $iname -0.0117 512 512
	treeshow $tree $node
  if { [$tree get $node plot] } {
		img tog $iname 1
	}
}

proc lb05 {tree} {
     global node nodeMax
     set flag $nodeMax
     Open $tree 
     if { $nodeMax>$flag } {
       set iname [file rootname [$tree label $node]]
       puts "img rotate $iname inf 512 512"
       img rotate $iname 0.015 512 512
       treeshow $tree $node
       $tree set $node plot 1
     }
}



proc confwindow {} {
    global tree tw
    set tw .conf
    toplevel $tw
    wm geometry $tw =500x250+0+480
    wm title $tw "Configuration"
    bind $tw <F3> { Open $tree }    
    bind $tw <F8> { Close $tree $node } 
    set w $tw.menu
    frame $w
    button $w.l -text Load -borderwidth 1 -command { Open $tree }
    button $w.re -text Reload -borderwidth 1 -command { Reopen $tree $node }
    button $w.un -text Unload -borderwidth 1 -command { Close $tree $node }
    set m $w.lbold.m
    menubutton $w.lbold -text LoadCHESS -menu $m
    menu $m -tearoff 0
    $m add command -label "lb04" -command {lb04 $tree}
    $m add command -label "relb04" -command {Relb04 $tree $node}
    $m add command -label "lb05" -command {lb05 $tree}
    $m add command -label "relb05" -command {Relb05 $tree $node}
    frame $w.ldf
    button $w.ldf.b -text LoadFit -borderwidth 1 -command { fitload $ldxdy $tree }
    label $w.ldf.l -text "dxy:"
    entry $w.ldf.e -textvariable ldxdy -width 8
    pack $w.ldf.b $w.ldf.l $w.ldf.e -side left
    pack $w.l $w.re $w.un $w.lbold -side left
    pack $w.ldf -side right
    pack $w -side top -anchor w
    
    set tw $tw.f
    frame $tw
    pack $tw -fill both
    set tree [tree create]
    treeview $tw.t \
	-width 0 \
	-yscrollcommand "$tw.vs set" \
	-xscrollcommand "$tw.hs set" \
	-yscrollincrement 1 \
	-selectmode single \
	-tree $tree
    $tw.t column insert 0 plot
    $tw.t column insert end note -width 200
    $tw.t column configure plot note -edit yes
    
    $tw.t bind Entry <1> {
	set node [%W index current]
	if { [catch {treeshow $tree $node}] == 0 } {
	    blt::tv::SetSelectionAnchor %W current
	}
    }
    $tw.t bind Entry <B1-Motion> {}

    scrollbar $tw.vs -orient vertical -command " $tw.t yview "
    scrollbar $tw.hs -orient horizontal -command " $tw.t xview "
      
    table $tw \
        0,0 $tw.t -fill both \
	1,0 $tw.hs -fill x \
	0,1 $tw.vs -fill y
    table configure $tw c1 r1 -resize none
          
    foreach column [$tw.t column names] {
	$tw.t column configure $column -command [list SortColumn $column]
    }
           
    $tw.t style checkbox check -onvalue 1 -offvalue 0 -showvalue no
    $tw.t column configure plot -style check
    $tree trace create all plot w treeplot
    
  $tree trace create all load w treeload 
} 

set setnfitt(11) 0
set NKfilter 0.0



proc setdata {dataset} {
     global tree node savedataset sigmab
     set iname [file rootname [$tree label $node]]
     set dname d$node
     #set datasetsave $dataset
     set datalist [split $dataset { }]
     set sigmab [lindex $datalist 0]
     set qr1 [lindex $datalist 1]
     set qr2 [lindex $datalist 2]
     set i 0
     while { [lindex $datalist [expr 2*$i+3]] != "" } {
       set qz1($i) [lindex $datalist [expr 2*$i+3]]
       set qz2($i) [lindex $datalist [expr 2*$i+4]]
       incr i
     }
     set savedataset $dataset
     puts "dataset load $dname $iname $dataset"
     if {$i==1} {dataset load $dname $iname $sigmab $qr1 $qr2 $qz1(0) $qz2(0)}
     if {$i==2} {dataset load $dname $iname $sigmab $qr1 $qr2 $qz1(0) $qz2(0) $qz1(1) $qz2(1)}
     if {$i==3} {dataset load $dname $iname $sigmab $qr1 $qr2 $qz1(0) $qz2(0) $qz1(1) $qz2(1) $qz1(2) $qz2(2)}
     if {$i==4} {dataset load $dname $iname $sigmab $qr1 $qr2 $qz1(0) $qz2(0) $qz1(1) $qz2(1) $qz1(2) $qz2(2) $qz1(3) $qz2(3)}
     if {$i==5} {dataset load $dname $iname $sigmab $qr1 $qr2 $qz1(0) $qz2(0) $qz1(1) $qz2(1) $qz1(2) $qz2(2) $qz1(3) $qz2(3) $qz1(4) $qz2(4)}
     for {set j 0} {$j<$i} {incr j} {
       puts "boxv set b$dname$j $qr1 $qr2 $qz1($j) $qz2($j)"
       boxv set b$dname$j $qr1 $qr2 $qz1($j) $qz2($j)
     }
	#set dataset $datasetsave
}

proc filterdata {} {
     global tree node
     set iname [file rootname [$tree label $node]]
     set dname d$node
     #puts "dataset filter $dname"
     dataset filter $dname
     #puts "img show FILTER"
     img show FILTER
}



proc every {ms body} {eval $body; after $ms [info level 0]}

proc exportBigFit {fitname} {   
    #set name FIT.tif
    set types {
	 {{tif Files}       {.tif}        }
	 {{All Files}        *             }
     }
     set tiffilename [tk_getSaveFile -initialdir "./" -filetypes $types]
     cd [file dirname $tiffilename]
     if { $tiffilename != "" } {
        if { [file rootname $tiffilename] == $tiffilename } {
	   set tiffilename $tiffilename.tif
	}
        puts "img export $fitname $tiffilename"
        img export $fitname $tiffilename     
     }
}

proc exportFit {fitname {automate 0}} {
    #set name FIT.tif

    global fittedhash
    set types {
         {{tif Files}       {.tif}        }
         {{All Files}        *             }
     }
     set tiffilename [tk_getSaveFile -initialdir "./" -filetypes $types]
     cd [file dirname $tiffilename]
     if { $tiffilename != "" } {
        if { [file rootname $tiffilename] == $tiffilename } {
           set tiffilename $tiffilename.tif
        }
        puts "img export FIT $tiffilename"
        img export "FIT" $tiffilename
    }
}

proc exportQzr {} {   
    set types {
	 {{dat Files}       {.dat}        }
	 {{All Files}        *             }
     }
     set qzrfilename [tk_getSaveFile -initialdir "./" -filetypes $types]
     cd [file dirname $qzrfilename]
     if { $qzrfilename != "" } {
        if { [file rootname $qzrfilename] == $qzrfilename } {
	   set qzrfilename $qzrfilename.dat
	}
        puts "qzrexport $qzrfilename"
        qzrexport $qzrfilename     
     }
}

proc NKstatist {bname tree node} {
     set iname [file rootname [$tree label $node]]
     #puts "img stat $iname $bname"
     img stat $iname $bname
}

proc fitbckg {b1 b2 ord n1 n2 tree node } {
     set iname [file rootname [$tree label $node]]
     puts "img fitlbg $iname $b1 $b2 $ord $n1 $n2"
     img fitlbg $iname $b1 $b2 $ord $n1 $n2
     #puts "img show $iname"
     img show $iname
     if { [$tree get $node plot] } {
	     #puts "img tog $iname 1"
	     img tog $iname 1
     }
}

proc fitload { ldxdy tree { automate 0}} {
     global node nodeMax fittedhash
     set ldx 0
     set ldy 0
     set ldxdy [split $ldxdy { }]
     if { [lindex $ldxdy 0] !="" } { set ldx [lindex $ldxdy 0] }
     if { [lindex $ldxdy 1] !="" } { set ldy [lindex $ldxdy 1] }
     set flag $nodeMax
    if {$automate == "0"} {
     Open $tree 
     if { $nodeMax>$flag } {
       set iname [file rootname [$tree label $node]]
       puts "img shiftfit $iname $ldx $ldy" 
       img shiftfit $iname $ldx $ldy
       treeshow $tree $node
       $tree set $node plot 1
     }
 } else {
     set name ""
     foreach name [array names fittedhash] {
     Open $tree 1 $name
     if { $nodeMax>$flag } {
       set iname [file rootname [$tree label $node]]
     
     set dataset $fittedhash($name)
     set datalist [split $dataset { }]
     set ldx [lindex $datalist 1]
     set i 0
     while { [lindex $datalist [expr 2*$i+3]] != "" } {
       set qz1($i) [lindex $datalist [expr 2*$i+3]]
       incr i
     }
     set ldy $qz1(0)

	
       puts "img shiftfit $iname $ldx $ldy"
       img shiftfit $iname $ldx $ldy
       treeshow $tree $node
       $tree set $node plot 1
     }
 }

}

}





proc maskpix { bname tree node } {
     global boxx0 boxx1 boxy0 boxy1
     boxv coords $bname
     set x0 [expr ($boxx0-5)]
     set x1 [expr ($boxx1+5)]
     boxv set tmp $x0 $x1 $boxy0 $boxy1
     fitbckg $bname tmp 2 144 1 $tree $node
     boxv del tmp
}




proc subcz {tree node} {
     set iname [file rootname [$tree label $node]]
     set dname d$node
     puts "img bgadjust $iname $dname"
     img bgadjust $iname $dname
     #puts "img show $iname"
     img show $iname
     if { [$tree get $node plot] } {
	     #puts "img tog $iname 1"
	     img tog $iname 1
     }
}




proc setbox {box} {
     set box [split $box { }]
     set name [lindex $box 0]
     set x1 [lindex $box 1]
     set x2 [lindex $box 2]
     set y1 [lindex $box 3]
     set y2 [lindex $box 4]    
     puts "boxv set $name $x1 $x2 $y1 $y2"
     boxv set $name $x1 $x2 $y1 $y2
}

proc loadfitall {} {
	global fittedhash
}
 
proc rotate {rotation tree node} {
     
     global rotatehash
     
     if {![info exists rotatehash($node)]} {
     set rotation [split $rotation { }]
     set iname [file rootname [$tree label $node]]
     set tg [lindex $rotation 0]
     set x0 [lindex $rotation 1]
     set y0 [lindex $rotation 2]
     puts "img rotate $iname $tg $x0 $y0"
     img rotate $iname $tg $x0 $y0
     #puts "img show $iname"
     img show $iname
     set rotatehash($node) 1
     if { [$tree get $node plot] } {
	     #puts "img tog $iname 1"
	     img tog $iname 1
     }
  }
}

proc addconstant {const tree node} {
     set iname [file rootname [$tree label $node]]
     puts "addconst $iname $const"
     addconst $iname $const
     #puts "img show $iname"
     img show $iname
     if { [$tree get $node plot] } {
	     #puts "img tog $iname 1"
	     img tog $iname 1
     }
}

proc multiconst {const tree node} {
     set iname [file rootname [$tree label $node]]
     puts "multi $iname $const"
     multi $iname $const
     #puts "img show $iname"
     img show $iname
     if { [$tree get $node plot] } {
	     #puts "img tog $iname 1"
	     img tog $iname 1
     }
}

proc subtract {subtr tree node} {
     set subtr [split $subtr { }]
     set iname [file rootname [$tree label $node]]
     set imgb [lindex $subtr 0]
     set b [lindex $subtr 1]
     puts "sub $iname $imgb $b"
     sub $iname $imgb $b
     #puts "img show $iname"
     img show $iname
     if { [$tree get $node plot] } {
	     #puts "img tog $iname 1"
	     img tog $iname 1
     }
}

proc shifttif {dxdy tree node} {
     set dx 0
     set dy 0
     set dxdy [split $dxdy { }]
     if { [lindex $dxdy 0] !="" } { set dx [lindex $dxdy 0] }
     if { [lindex $dxdy 1] !="" } { set dy [lindex $dxdy 1] }
     set iname [file rootname [$tree label $node]]
     puts "img shifttif $iname $dx $dy"
     img shifttif $iname $dx $dy
     #puts "img show $iname"
     img show $iname
     if { [$tree get $node plot] } {
	     #puts "img tog $iname 1"
	     img tog $iname 1
     }
}

proc shiftfit {fdxdy tree node} {
     set fdx 0
     set fdy 0
     set fdxdy [split $fdxdy { }]
     if { [lindex $fdxdy 0] !="" } { set fdx [lindex $fdxdy 0] }
     if { [lindex $fdxdy 1] !="" } { set fdy [lindex $fdxdy 1] }
     set iname [file rootname [$tree label $node]]
     puts "img shiftfit $iname $fdx $fdy"
     img shiftfit $iname $fdx $fdy
     #puts "img show $iname"
     img show $iname
     if { [$tree get $node plot] } {
	     #puts "img tog $iname 1"
	     img tog $iname 1
     }
}

proc manipwin {} {     
	global tree node
	set w .manip
	toplevel $w
	wm title $w "Manipulation"
	wm geometry $w =200x300+1200+0
#	wm geometry .fit =+570+480
	
	set w $w.flip
	frame $w
	button $w.x -text "flip x" -command { flip x $tree $node }
	button $w.y -text "flip y" -command { flip y $tree $node }
        pack $w.x $w.y -side left -expand yes -fill x
	
	set w .manip.box
	frame $w
	label $w.l -text "box: name x1 x2 y1 y2"
	entry $w.e1 -textvariable mbox1 -width 20
	entry $w.e2 -textvariable mbox2 -width 20
	pack $w.l $w.e1 $w.e2 -side top
	bind $w.e1 <Return> { setbox $mbox1 }
	bind $w.e2 <Return> { setbox $mbox2 }
	
	set w .manip.rot
	frame $w
	label $w.l -text "rotate: tg x0 y0"
	entry $w.e -textvariable rotation -width 20
	pack $w.l $w.e -side top
	bind $w.e <Return> { rotate $rotation $tree $node}
	
	set w .manip.add
	frame $w
	label $w.l -text "add: constant "
	entry $w.e -textvariable addc -width 8
	pack $w.l $w.e -side left -expand yes
	bind $w.e <Return> { addconstant $addc $tree $node}
	
	set w .manip.multi
	frame $w
	label $w.l -text "multi by: constant"
	entry $w.e -textvariable mulc -width 5
	pack $w.l $w.e -side left -fill x
	bind $w.e <Return> { multiconst $mulc $tree $node}
	
	set w .manip.shifttif
	frame $w
	label $w.l -text "shift by: dx dy"
	entry $w.e -textvariable dxdy -width 8
	pack $w.l $w.e -side left -fill x
	bind $w.e <Return> { shifttif $dxdy $tree $node}

	#moved to the conf window
	set w .manip.shiftfit
	frame $w
	label $w.l -text "shift fit by: dx dy"
	entry $w.e -textvariable fdxdy -width 8
	pack $w.l $w.e -side left -fill x
	bind $w.e <Return> { shiftfit $fdxdy $tree $node}
		
	set w .manip.sub
	frame $w
	label $w.l -text "sub: image const"
	entry $w.e -textvariable subtr -width 20
	pack $w.l $w.e -side top
	bind $w.e <Return> { subtract $subtr $tree $node}
	
	set w .manip.exp
	frame $w
	button $w.exp -text export -borderwidth 3 -command { exportBigFit [file rootname [$tree label $node]] }
	label $w.l -text image
	pack $w.exp $w.l -side left
	
	set w .manip.expq
	frame $w
	button $w.exp -text export -borderwidth 3 -command { exportQzr }
	label $w.l -text qzr-plot
	pack $w.exp $w.l -side left
	
	set w .manip
	pack $w.flip $w.box $w.rot $w.add $w.multi $w.shifttif $w.sub $w.exp $w.expq -side top
}

proc colorwin {} {
	global tree node colorSAT colorMAX colorMIN colorBOT
	#catch { set imgname [file rootname [$tree label $node]] }
	set w .color
	toplevel $w
	wm title $w "Color"
	
	colorsetup update
        
	label $w.lsat -text saturation -width 10
	entry $w.esat -textvariable colorSAT -borderwidth 1 -width 10
	label $w.lmax -text maximum -width 10
	entry $w.emax -textvariable colorMAX -borderwidth 1 -width 10
	label $w.lmin -text minimum -width 10
	entry $w.emin -textvariable colorMIN -borderwidth 1 -width 10
	label $w.lbot -text bottom -width 10
	entry $w.ebot -textvariable colorBOT -borderwidth 1 -width 10
	label $w.llog -text "log scale" -width 10
	checkbutton $w.log -variable mod -padx 0 -pady 0 -borderwidth 1 -command { colorsetup mode [expr $mod?"log":"linear"]; if { [catch { set imgname [file rootname [$tree label $node]] }] ==0 } { img show $imgname } }
	bind $w.esat <Return> { colorsetup saturation $colorSAT; if { [catch { set imgname [file rootname [$tree label $node]] }] ==0 } { img show $imgname } }
	bind $w.emin <Return> { colorsetup minmax $colorMIN $colorMAX; if { [catch { set imgname [file rootname [$tree label $node]] }] ==0 } { img show $imgname } }
	bind $w.emax <Return> { colorsetup minmax $colorMIN $colorMAX; if { [catch { set imgname [file rootname [$tree label $node]] }] ==0 } { img show $imgname } }
	bind $w.ebot <Return> { colorsetup bottom $colorBOT; if { [catch { set imgname [file rootname [$tree label $node]] }] ==0 } { img show $imgname } }

	grid $w.llog $w.log
	grid $w.lsat $w.esat
	grid $w.lmax $w.emax
	grid $w.lmin $w.emin
	grid $w.lbot $w.ebot
	bind $w <F2> {
          colorsetup update
	}
}


#option add *Graph.pen.LineWidth 0
#option add *Graph.pen.Color red
#option add *Graph.pen.Symbol square
#option add *Graph.pen.Pixels 1.5


$imgWN marker configure bg -image realimg0

source tcl/plotpanel.tcl
plotpanel
#controlpanel /data/2003feb/dopctest
wm protocol .plot WM_DELETE_WINDOW {wm withdraw .plot}

tv::colorsetup minmax 0 1000
tv::colorsetup style gray1
tv::colorsetup sat 32000

set spsummaxdx 1
set spqxmindx 0.05
set spfxmindx 0.05
set spsummindx 0.05

source tcl/tools.tcl
source tcl/tvDS.tcl
modelcc new mc
modelcc retol mc
set qagabserr 0
set qagrelerr 1e-10
set qaworelerr 1.0e-6
#ctrl
puts initialized

##### Added by NK 07/03 ####




proc lbcomp {imga imgh a h s} {
#name of LB_air LB_He & counts for air He sample
    set y [expr ($s-$h)/($a-$h)]
    set x [expr ($a-$s)/($a-$h)]
    img expr $imga *= $x
    img expr $imgh *= $y
    img expr $imga += $imgh
    return "$x $y"
}

proc multi {imga a} {
    img expr $imga *= $a
}

proc addconst {imga a} {
    img expr $imga += $a
}

proc sub {imga imgb b} {
    img expr $imgb *= $b
    img expr $imga -= $imgb
    img expr $imgb /= $b
}

proc resetcolor {min max} {
    global tree node
    colorsetup minmax $min $max
    set iname [file rootname [$tree label $node]]
    #puts "img show $iname"
    img show $iname    
}
colorsetup mode linear
confwindow

source tcl/kiyo.tcl
source tcl/fittingWindow.tcl
source tcl/fileManage.tcl
source tcl/procedures.tcl


