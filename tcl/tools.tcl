proc calc {} {
    toplevel .calc
    wm title .calc "Calculator"
    global nf_lambda nf_eta nf_kc nf_b nf_D nf_xi nf_T
    set w .calc
    label $w.lbT -text "T:"
    entry $w.enT -textvariable nf_T -width 30
    label $w.lb0 -text "Lambda:"
    entry $w.en0 -textvariable nf_lambda -width 30
    label $w.lb1 -text "Eta :"
    entry $w.en1 -textvariable nf_eta -width 30
    label $w.lb2 -text "D:"
    entry $w.en2 -textvariable nf_D -width 30
    label $w.lb3 -text "Kc:"
    entry $w.en3 -textvariable nf_kc -width 30
    label $w.lb4 -text "B:"
    entry $w.en4 -textvariable nf_b -width 30
    label $w.lb5 -text "xi:"
    entry $w.en5 -textvariable nf_xi -width 30
    grid $w.lbT $w.enT
    for {set i 0} {$i<6} {incr i} {
	grid $w.lb$i $w.en$i
    }
    bind $w.en0 <Return> {updatekcb}
    bind $w.en1 <Return> {updatekcb}
    bind $w.en3 <Return> {lambdaeta}
    bind $w.en4 <Return> {lambdaeta}
    bind $w.en5 <Return> {set nf_lambda [expr $nf_xi*$nf_xi/$nf_D]; updatekcb}
}

proc updatekcb {} {
    global nf_lambda nf_eta nf_kc nf_b nf_D nf_xi nf_T
    #set nf_kc [expr $nf_lambda*6.56e-14/$nf_D/$nf_eta]
    #set nf_b [expr 6.56e18/pow($nf_D, 3)/$nf_lambda/$nf_eta]
    set nf_kc [expr $nf_lambda*2.17e-16*($nf_T+273)/$nf_D/$nf_eta]
    set nf_b [expr 2.17e16*($nf_T+273)/pow($nf_D, 3)/$nf_lambda/$nf_eta]
    set nf_xi [expr sqrt($nf_D*$nf_lambda)]
}
proc lambdaeta {} {
    global nf_lambda nf_eta nf_kc nf_b nf_D nf_xi nf_T
    set nf_lambda [expr sqrt($nf_kc/$nf_b)/$nf_D*1e16]
    #set nf_eta [expr 6.56e2/sqrt($nf_kc*$nf_b)/pow($nf_D, 2)]
    set nf_eta [expr 2.17*($nf_T+273)/sqrt($nf_kc*$nf_b)/pow($nf_D, 2)]
    set nf_xi [expr sqrt($nf_D*$nf_lambda)]
}
