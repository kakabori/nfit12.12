proc plotpanel {} {
    global xyplot .plot

    destroy .plot;
    toplevel .plot
    wm geometry .plot =500x430+570+0 
    wm title .plot Plot
    set w .plot.menu
    frame $w
    set m $w.opt.m
    menubutton $w.opt -text Opt -menu $m
    menu $m -tearoff 0
    $m add command -label "xplot" -command "wm title .plot xplot; set plotopt 0; plotcs"
    $m add command -label "yplot" -command "wm title .plot yplot; set plotopt 1; plotcs"
    set m $w.scale.m
    menubutton $w.scale -text Axis -menu $m
    menu $m -tearoff 0
    $m add command -label "y_log" -command "$xyplot axis configure y -logscale yes"
    $m add command -label "y_linear" -command "$xyplot axis configure y -logscale no"
    $m add command -label "showall" -command {$xyplot axis configure y -max "" -min ""; $xyplot axis configure x -max "" -min ""}
    
    label $w.lx -text "x:" -fg red
    label $w.x -textvariable plotx -width 8 -anchor w
    label $w.ly -text "y:" -fg red
    label $w.y -textvariable ploty -width 8 -anchor w
    
    label $w.lxl -text "xlow:"
    label $w.lxh -text "xhigh:"
    entry $w.xlow -textvariable xlow -width 6
    entry $w.xhigh -textvariable xhigh -width 6
    pack $w.opt $w.scale $w.lx $w.x $w.ly $w.y -side left
    pack $w.xhigh $w.lxh $w.xlow $w.lxl -side right    
    
    pack $w -side top -anchor w -expand yes -fill both
    #ctrl .plot { -side bottom -fill x}
    graph $xyplot -width 2048 -height 2048 -plotrelief flat -plotpadx {0 0} -plotpady {0 0} -plotborderwidth 0 -bg white
    $xyplot axis configure y -subdivisions 5
    $xyplot axis configure x -subdivisions 5
    Blt_ZoomStack $xyplot
    Blt_Crosshairs $xyplot
    #$xyplot legend configure -hide yes
    $xyplot legend configure -ipadx 0 -ipady 0
    $xyplot legend bind all <1> { Insignify %W }
    $xyplot legend bind all <3> { raiselement %W }

    pack $xyplot
    bind $xyplot	<B1-Motion>	{plotxy %W %x %y}
    bind $xyplot	<1>	{plotxy %W %x %y}
    bind .plot <Key-Up> "axisshift $xyplot 10 y"
    bind .plot <Key-Down> "axisshift $xyplot -10 y"
    bind .plot <Key-Left> "axisshift $xyplot 10 x"
    bind .plot <Key-Right> "axisshift $xyplot -10 x"
    
    $xyplot pen create pen
}

proc Insignify { graph } {
    set entry [$graph legend get current]
    set active [$graph legend activate]
    if { [lsearch $active $entry] < 0 } {
	$graph legend activate $entry
	$graph element configure $entry -pixels 0
    } else {
        $graph legend deactivate $entry
	$graph element configure $entry -pixels 1.5
    }
}

proc raiselement { graph } {
    set i [$graph legend get current]
    set ori [$graph element show]
    set tgt [lsearch $ori $i]
    set ori [lreplace $ori $tgt $tgt]
    set ori [linsert $ori end $i]
    $graph element show $ori
}
