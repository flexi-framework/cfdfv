#!/bin/bash
# update the references for the implementation test

cfdfv="../../$(find -name cfdfv -type f -executable)"

run_and_update() {
	ini=$(basename $1)
	echo "> Running $ini..."

	for var in $@; do ln -s ../../$var .; done
	$cfdfv $ini > "${ini%.*}.log"
	find -type l -delete

	last=$(ls ${ini%.*}_0* | sort | tail -n 1)
	mv $last reference/$last

	rm *.cgns *.log *.dem *.csv
}

# Aufgabe 0
echo "Checking Aufgabe 0"
cd check/Aufg_0 &&
	run_and_update Calc/Profil/Keilprofil/keil_coarse.ini Calc/Profil/Keilprofil/keil_coarse.msh
	run_and_update Calc/Profil/Keilprofil/keil_intermediate.ini Calc/Profil/Keilprofil/keil_intermediate.msh
	#run_and_update Calc/Profil/Keilprofil/keil_fine.ini Calc/Profil/Keilprofil/keil_fine.msh
	#run_and_update Calc/Profil/Keilprofil/keil_ultrafine.ini Calc/Profil/Keilprofil/keil_ultrafine.msh
