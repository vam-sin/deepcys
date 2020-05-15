#!/usr/bin/tclsh

set TOPOFILE /home/banshee/Academics/Sem6/BioFormal/Cysteine/dataset/ThioEther/menv-server/top_all27_prot_na.inp

proc psfalias { } {

alias residue HIS HSP
alias atom ILE CD1 CD
alias atom LYS 1HZ HZ1
alias atom LYS 2HZ HZ2
alias atom LYS 3HZ HZ3
alias atom ARG 1HH1 HH11
alias atom ARG 2HH1 HH12
alias atom ARG 1HH2 HH21
alias atom ARG 2HH2 HH22
alias atom ASN 1HD2 HD21
alias atom ASN 2HD2 HD22
alias atom SER HG HG1
	
}
proc del { foo } {
	exec rm $foo
}

proc splitpdb { fname } {
	set in [open $fname r]
 	set nseg 0

	set oldres -1
	set curres -1

	foreach line [split [read -nonewline $in] \n] {

		set head [string range $line 0 5]
		if { ![string compare $head "ATOM  "] } {
			set curres [string range $line 22 25]
			set resdif [expr $curres - $oldres]
			if { $oldres != -1 && $resdif != 0 && $resdif != 1 } {
		
				puts $out END
				close $out
				set oldres -1
			}
			if { $oldres == -1 } {
				
				incr nseg
				set newname "${fname}_${nseg}.pdb"
				set out [open $newname w]
				lappend fnamelist $newname
			}
			puts $out $line
			set oldres $curres
		
		}
	}
	close $out
	close $in
	return $fnamelist
}

proc build { pdb {outname psfgen}} {
	global TOPOFILE 
	topology $TOPOFILE 
	psfalias
	set nseg 0
	foreach segfile [splitpdb $pdb] {
		incr nseg
		set segid "P${nseg}"
		segment $segid {
			pdb $segfile
		}
		coordpdb $segfile $segid
		del $segfile
	}
	guesscoord
        writepdb ${outname}.pdb 
        writepsf ${outname}.psf
}
eval build $argv
#catch {eval build $argv} 
