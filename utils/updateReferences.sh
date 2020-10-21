#!/bin/bash
# update the references for the implementation test

cfdfv="../../$(find -name cfdfv -type f -executable)"

run_and_update() {
	ini=$(basename $1)
	echo "> Running $ini..."

	# create symbolic links of the mesh files
	for var in $@; do ln -s ../../$var .; done

	# run cfdfv
	$cfdfv $ini > "${ini%.*}.log"

	# delete all symbolic links
	find -type l -delete

	# update the reference
	last=$(ls ${ini%.*}*_0* | grep -v ex | sort | tail -n 1)
	echo "> Updating reference/$last..."
	mv $last reference/$last

	# remove remaining auxilliary data
	rm -f *.cgns *.log *.dem *.csv
}

# Aufgabe 0
echo "Checking BlasiusBoundaryLayer"
cd check/Aufg_0
	run_and_update Calc/Profil/Keilprofil/keil_coarse.ini Calc/Profil/Keilprofil/keil_coarse.msh
	run_and_update Calc/Profil/Keilprofil/keil_intermediate.ini Calc/Profil/Keilprofil/keil_intermediate.msh
	#run_and_update Calc/Profil/Keilprofil/keil_fine.ini Calc/Profil/Keilprofil/keil_fine.msh
	#run_and_update Calc/Profil/Keilprofil/keil_ultrafine.ini Calc/Profil/Keilprofil/keil_ultrafine.msh

cd ../..

echo "Checking Aufgabe 1"
cd check/Aufg_1
	run_and_update Calc/RiemannProblems/sod.ini
	run_and_update Calc/RiemannProblems/sod_implicit.ini
	run_and_update Calc/RiemannProblems/toro1.ini
	run_and_update Calc/RiemannProblems/toro2.ini
	run_and_update Calc/RiemannProblems/toro3.ini
	run_and_update Calc/RiemannProblems/toro4.ini
	run_and_update Calc/RiemannProblems/toro5.ini

cd ../..

echo "Checking Aufgabe 2"
cd check/Aufg_2
	#run_and_update Calc/GaussianPulse/gaussianPulse.ini
	#run_and_update Calc/GaussianPulse/gaussianPulse_circgrid.ini Calc/GaussianPulse/circgrid.msh
	#run_and_update Calc/DoubleMach/dmr.ini
	#run_and_update Calc/ForwardFacingStep/ffs.ini Calc/ForwardFacingStep/ffs.mesh
	#run_and_update Calc/ForwardFacingStep/ffs_subsonisch_pressureout.ini Calc/ForwardFacingStep/ffs_subsonisch_pressureout.mesh
	#run_and_update Calc/ForwardFacingStep/ffs_subsonisch_supersonicout.ini Calc/ForwardFacingStep/ffs_subsonisch_supersonicout.mesh

cd ../..

echo "Checking Aufgabe 3"
cd check/Aufg_3
	#run_and_update Calc/RichtmyerMeshkovInstability/rmi.ini
	#run_and_update Calc/Profil/NACA2312/NACA2312.ini Calc/Profil/NACA2312/NACA2312.mesh
	#run_and_update Calc/SineWave/1D/SineWave1D.ini
	#run_and_update Calc/SineWave/2D/SineWaveO1.ini
	#run_and_update Calc/SineWave/2D/SineWaveO2.ini

cd ../..

echo "Checking Aufgabe 4"
cd check/Aufg_4
	#run_and_update Calc/Profil/NACA0012/NACA0012.ini Calc/Profil/NACA0012/NACA0012.msh
	#run_and_update Calc/Cylinder/Cylinder.ini Calc/Cylinder/Cylinder.msh
	#run_and_update Calc/Cylinder/Cylinder_implicit.ini Calc/Cylinder/Cylinder.msh

cd ../..
