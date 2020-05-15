#!/usr/bin/bash


#Takes pdb file as command line argument
#echo "Enter the pdb file:"
#read pdbfilename

filename=$(echo $1 | cut -f1 -d".")
#echo ${filename}
python3 /home/banshee/Academics/Sem6/BioFormal/Cysteine/dataset/ThioEther/menv-server/renum.py $1

#Takes the pdb file and prints unique resiudes
cat pythonresults.txt | cut -f2 | uniq > uniqueresult.pdb

count=$(wc -l uniqueresult.pdb | awk '{print $1}')

#Adds a column of numbers starting from 1
cat uniqueresult.pdb | awk '{print FNR "\t" $0}' uniqueresult.pdb > result.pdb

#Formatting the columns in the file
cat result.pdb | awk '{printf("%6s ", $1); printf("%3s %3s ", $2 ,$3); printf("%3s\n",$4);}' > uniqueresidues.txt

last_res_no=$(sed -e 1b -e '$!d' uniqueresidues.txt | awk 'NR==2{print $1}')
#echo ${last_res_no}

#Removing temporary files

#Takes menv file as command line argument
#echo "Enter the Menv file:"
#read menvfilename

#Code that takes resiude number, residue and rHpy value
python3 /home/banshee/Academics/Sem6/BioFormal/Cysteine/dataset/ThioEther/menv-server/renumCrd.py $2

#Takes only the residues with value greater than 0
cat Crdresults.txt | awk '($1 > 0)' Crdresults.txt > Menvresults.txt

last_menv_res_no=$(sed -e 1b -e '$!d' Menvresults.txt | awk 'NR==2{print $1}')
#echo ${last_menv_res_no}

if [ $last_menv_res_no -lt $last_res_no ]
then 
   echo "Residue lost in Menv file ! "
fi

#Removing temporary files

awk '{printf("%3s ",$1); printf("%3s\n",$2);}' Menvresults.txt | uniq > uniqueMenvresidues.txt

# #joining two files based on residue numbers
# join --nocheck-order -1 1 -2 1 -a1 Menvresults.txt uniqueresidues.txt > combined.txt

# #Extracting lines which dont have the matching residues in column 2 and 4
# cat combined.txt | awk '$2!=$4 {print $1,$2,$3,$4,$5,$6}' > notequal.txt

# #Extracting lines which dont have matching residues other than HSP and HIS
# awk '($2 != "HSP" && $4 != "HIS")' notequal.txt > tada.txt

# #Taking the residue no which is not there in menv file
# lost_residue_no=$(sed -e 1b -e '$!d' tada.txt | awk 'NR==1{print $1}')

# echo "Lost residue number: "${lost_residue_no}


# #Removing that particular residue in uniqueresidues file
# awk '($1 != "'$lost_residue_no'")' uniqueresidues.txt > newuniqueresidues.txt

python3 /home/banshee/Academics/Sem6/BioFormal/Cysteine/dataset/ThioEther/menv-server/residuesremoval.py

#Adding a column of numbers to this new file
cat outfile.txt | awk '{print FNR "\t" $0}' outfile.txt > totalresidues.txt
join --nocheck-order -1 1 -2 1 -a1 Menvresults.txt totalresidues.txt > combined.txt

#Joining the new file with Menvresults file
join --nocheck-order -1 1 -2 1 -a1 Menvresults.txt totalresidues.txt > finalresult.txt

cat $2 | awk '{print FNR "\t" $0}' $2 > menvfilewithnum.txt
cat finalresult.txt | awk '{print FNR "\t" $0}' finalresult.txt > finalresultwithnum.txt

last_no=$(sed -e 1b -e '$!d' finalresultwithnum.txt | awk 'NR==2{print $1}')
#echo $last_no

join --nocheck-order -1 1 -2 1 -a1 menvfilewithnum.txt finalresultwithnum.txt > output.txt

cat output.txt | awk '{printf("%5s %5s %5s %5s %5s %5s %10s %10s %10s %10s %5s %5s %10s %5s %5s %5s ",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16); printf("%5s\n", $17);}' > revisedoutput.txt


awk '($1 <= '$last_no')' revisedoutput.txt > requiredoutput.txt

awk '!($1="")!($11="")!($12="")!($13="")!($14="")!($15="")' requiredoutput.txt > finalfinal.txt
cat finalfinal.txt | awk '{printf("%5s %5s %5s %5s %5s %10s %10s %10s %10s %5s ",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10);printf("%5s\n",$11);}' > finalappended.menv

awk '{$2=""; sub(" ", " "); print}' finalappended.menv > finalcut.menv
#awk '{gsub(/\,/,"");print $1,$9,$10,$2,$3,$4,$5,$6,$7,$8}' finalcut.menv > commaremoved.menv
awk '{print $1,$9,$10,$2,$3,$4,$5,$6,$7,$8}' finalcut.menv > commaremoved.menv
cat commaremoved.menv | awk '{printf("%5s %5s %5s %5s %5s %5s %10s %10s %10s ",$1,$2,$3,$4,$5,$6,$7,$8,$9);printf("%10s\n",$10);}' > "$filename"-new.menv


# python3 /home/banshee/Academics/Sem6/BioFormal/Cysteine/dataset/ThioEther/menv-server/plotting.py "$filename"-new.menv
