proc lbcomp {imga imgc a b c} {
    set y [expr ($b-$a)/($c-$a)]
    set x [expr ($c-$b)/($c-$a)]
    img expr $imga *= $x
    img expr $imgc *= $y
    img expr $imga += $imgc
    return "$x $y"
}


proc statall {b} {
    foreach m [img names] {
	puts "$m : [img stat $m $b]"
    }
}
