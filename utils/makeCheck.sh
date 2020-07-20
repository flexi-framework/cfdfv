#!/bin/bash
# test the implementation

cfdfv="../../$(find -name cfdfv -type f -executable | head -n 1)"
cgnsdiff="../../$(find -name cgnsdiff -type f -executable | head -n 1)"

run_and_check() {
	ini=$(basename $1)
	echo "> Running $ini..."

	for var in $@; do ln -s ../../$var .; done
	$cfdfv $ini > "${ini%.*}.log"
	find -type l -delete

	last=$(ls ${ini%.*}_0* | sort | tail -n 1)
	echo "> Checking $last..."
	$cgnsdiff -qd -t 0.1 $last reference/$last > "${ini%.*}.diff"

	rm *.cgns *.dem *.csv

	[ ! -s "${ini%.*}.diff" ] || (echo "ERRORS in $ini, abort!" && exit 1)
}

# Aufgabe 0
echo "Checking Aufgabe 0"
cd check/Aufg_0 &&
	run_and_check Calc/Profil/Keilprofil/keil_coarse.ini Calc/Profil/Keilprofil/keil_coarse.msh
	run_and_check Calc/Profil/Keilprofil/keil_intermediate.ini Calc/Profil/Keilprofil/keil_intermediate.msh
	#run_and_check Calc/Profil/Keilprofil/keil_fine.ini Calc/Profil/Keilprofil/keil_fine.msh
	#run_and_check Calc/Profil/Keilprofil/keil_ultrafine.ini Calc/Profil/Keilprofil/keil_ultrafine.msh
