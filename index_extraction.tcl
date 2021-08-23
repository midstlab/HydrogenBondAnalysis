## Gets the index number of the residues
## Change the line 13 file name
## "source index_extraction.tcl" from VMD Tk console


set pro [atomselect top "protein"] # protein selection
set dhf [atomselect top "resname DHF"] # ligand selection

set ind1 [$pro get index]
set ind2 [$dhf get index]

set ind [concat $ind1 $ind2]
#set ind [lremove $ind3 0]
set outfile [open i94l-index.dat w] % !! this line needs to be changed before sourcing

set b [llength $ind]

for {set i 0} {$i < $b} {incr i} {

set a [atomselect top "index [lindex $ind $i]"]
set d [atomselect top "index [lindex $ind $i]"]
set c [atomselect top "index [lindex $ind $i]"]

 set resn [ $a get resname]
 set resi [ $d get resid]
 set nm   [ $c get name]
 
 set str [format "%1d %3s%1d %1s" [lindex $ind $i] $resn $resi $nm]
 
 join $str ""
 
 puts $outfile $str

}

close $outfile