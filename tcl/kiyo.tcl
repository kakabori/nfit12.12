proc create_kiyomask {qrlow qrhigh qzlow qzhigh} {
	set name "kiyomask.txt"
	set fid [open $name w]
	for {set i $qzlow} {$i <= $qzhigh} {incr i} {
		puts $fid "$i $qrlow $qrhigh"
	}
	close $fid
}

proc fixedangle {} {
  cd ~/Documents/dopcdope_3to1/
}
