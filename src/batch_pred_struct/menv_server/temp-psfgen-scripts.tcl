#!/usr/bin/tclsh

exec grep -v "HETATM" 253l.pdb > 253l-protein.pdb

puts "say something"

# (3) Read topology file
topology /home/vamsi/Academics/Sem6/BioFormal/deepcys_code/src/batch_pred_struct/menv_server/top_all27_prot_na.inp


# (6) Read protein coordinates from PDB file
pdbalias residue HIS HSP
pdbalias atom ILE CD1 CD    ; # formerly "alias atom ..."
pdbalias atom ILE CD1 CD
pdbalias atom LYS 1HZ HZ1
pdbalias atom LYS 2HZ HZ2
pdbalias atom LYS 3HZ HZ3
pdbalias atom ARG 1HH1 HH11
pdbalias atom ARG 2HH1 HH12
pdbalias atom ARG 1HH2 HH21
pdbalias atom ARG 2HH2 HH22
pdbalias atom ASN 1HD2 HD21
pdbalias atom ASN 2HD2 HD22
pdbalias atom SER HG HG1

segment BPTI {
 pdb 253l-protein.pdb
}

coordpdb 253l-protein.pdb BPTI
guesscoord
writepsf new-proteins.psf
writepdb new-proteins.pdb
#END of psfgen commands
#exec ENDMOL
