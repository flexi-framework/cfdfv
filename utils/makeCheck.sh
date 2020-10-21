#!/bin/bash
# test the implementation

cfdfv="../../$(find -name cfdfv -type f -executable | head -n 1)"
cgnsdiff="../../$(find -name cgnsdiff -type f -executable | head -n 1)"
csvdiff="../../$(find -name csvdiff.py -type f -executable | head -n 1)"

run_and_check_cgns() {
	ini=$(basename $1)
	echo "> Running $ini..."

	# create symbolic links of the mesh files
	for var in $@; do ln -s ../../$var .; done

	# run cfdfv
	$cfdfv $ini > "${ini%.*}.log"

	# delete all the symbolic links
	find -type l -delete

	# get the last output of the cfdfv run
	last=$(ls ${ini%.*}*_0* | grep -v ex | sort | tail -n 1)

	echo "> Checking $last..."

	# run the last output through cgnsdiff
	# options are:
	# -q: pint only of there is a difference
	# -d: compare node data alsp
	# -t: tolerance for comparing floats/doubles
	$cgnsdiff -qd -t 0.1 $last reference/$last > "${ini%.*}.diff"

	# delete all remaining auxilliary files
	rm -f *.cgns *.dem *.csv

	# make sure that there was no error, if there was, then break
	[ ! -s "${ini%.*}.diff" ] || (echo "ERRORS in $ini, abort!" && exit 1)
}

run_and_check_csv() {
	ini=$(basename $1)
	echo "> Running $ini..."

	# create symbolic links of the mesh files
	for var in $@; do ln -s ../../$var .; done

	# run cfdfv
	$cfdfv $ini > "${ini%.*}.log"

	# delete all the symbolic links
	find -type l -delete

	# get the last output of the cfdfv run
	last=$(ls ${ini%.*}*_0* | grep -v ex | sort | tail -n 1)

	echo "> Checking $last..."

	# run the last output through csvdiff
	# the first argument is the float/double tolerance
	# second are third are the files that are to be checked
	$csvdiff 0.1 $last reference/$last > "${ini%.*}.diff"

	# delete all remaining auxilliary files
	rm -f *.cgns *.dem *.csv

	# make sure that there was no error, if there was, then break
	[ ! -s "${ini%.*}.diff" ] || (echo "ERRORS in $ini, abort!" && exit 1)
}

# Aufgabe 0
echo "Checking Aufgabe 0"
cd check/Aufg_0 &&
	run_and_check_cgns Calc/Profil/Keilprofil/keil_coarse.ini Calc/Profil/Keilprofil/keil_coarse.msh
	run_and_check_cgns Calc/Profil/Keilprofil/keil_intermediate.ini Calc/Profil/Keilprofil/keil_intermediate.msh
	#run_and_check_cgns Calc/Profil/Keilprofil/keil_fine.ini Calc/Profil/Keilprofil/keil_fine.msh
	#run_and_check_cgns Calc/Profil/Keilprofil/keil_ultrafine.ini Calc/Profil/Keilprofil/keil_ultrafine.msh

cd ../..

echo "Checking Aufgabe 1"
cd check/Aufg_1 &&
	run_and_check_csv Calc/RiemannProblems/sod.ini
	run_and_check_csv Calc/RiemannProblems/sod_implicit.ini
	run_and_check_csv Calc/RiemannProblems/toro1.ini
	run_and_check_csv Calc/RiemannProblems/toro2.ini
	run_and_check_csv Calc/RiemannProblems/toro3.ini
	run_and_check_csv Calc/RiemannProblems/toro4.ini
	run_and_check_csv Calc/RiemannProblems/toro5.ini
