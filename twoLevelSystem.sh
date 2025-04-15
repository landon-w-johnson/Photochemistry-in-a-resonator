#!/bin/bash

########################################################################
#    This script controls the flow of the TDPT VASP add-on package.    #
#    It calls the FORTRAN scripts that do the actual calculations.     #
#    It also does some data consolidation using native VASP outputs    #
########################################################################





########################################################################
#                Read in relevant settings from PHOTCAR                #
########################################################################

CHGintStr=$(grep "CHGint" PHOTCAR)

if [[ ! "$CHGintStr" =~ [^[:space:]] ]]; then # checking for blank str
    echo "CHGint not found in PHOTCAR. Setting CHGint = 0" >> outputFile.txt
    CHGint=0
else
    #echo "CHGintStr = '$CHGintStr'" >> outputFile.txt
    CHGintVal=$(awk -F '#' '{ split($1,tempArray,"="); gsub(" ","",tempArray[2]); print tempArray[2] }' <<< "$CHGintStr")
    if [[ "$CHGintVal" =~ [^[:space:]] && ! "$CHGintVal" =~ [^[:digit:]] ]]; then
	CHGint="$CHGintVal"
	echo "CHGint = $CHGint" >> outputFile.txt
    else
	echo "CHGint not found in PHOTCAR. Setting CHGint = 0" >> outputFile.txt
	CHGint=0
    fi
fi
echo '' >> outputFile.txt





########################################################################
#                     Run initial VASP optimization                    #
########################################################################

#startTime=$(date +%s%3N)
#mpiexec -machinefile $PBS_NODEFILE -np $NUM_PROC vasp_gam >> outputFile.txt # CCAST
#srun -n 64 vasp_gam # NERSC
#if [ $? -ne 0 ]; then
#    echo "VASP crashed. Terminating twoLevelSystem.sh" >> outputFile.txt
#    exit
#fi
#endTime=$(date +%s%3N)
#runTime=$((endTime-startTime))
#echo "Initial VASP runtime:" >> outputFile.txt
#echo "$runTime ms" >> outputFile.txt
#echo '' >> outputFile.txt

##### This block selectively grabs and writes the forces on the nuclei from OUTCAR #####
#writeLine=false
#cat OUTCAR | while read line
#do
#    if [[ "$line" =~ "FORCES acting on ions" ]]; then
#	writeLine=true
#	echo "" >> outputFile.txt
#    fi
#    if [[ "$writeLine" = true ]]; then
#	if [[ ! "$line" =~ [^[:space:]] ]]; then
#	    echo "" >> outputFile.txt
#	    break
#	fi
#	echo "$line" >> outputFile.txt
#    fi
#done

##### Save a copy of CHG for visualizing electron density through time #####

if [[ ! $CHGint = 0 ]]; then
    mkdir CHG_Dir
    cp CHG CHG_Dir/CHG00000
fi





########################################################################
#                  Extract data from optimized WAVECAR                 #
#                         Then update WAVECAR                          #
########################################################################

fileName='dataFile.txt'

if [ -f "$fileName" ];
then
    echo "$fileName found, skipping initialization." >> outputFile.txt
    echo "Continuing job." >> outputFile.txt
else
    #cp WAVECAR WAVECAR_00000
    cp WAVECAR WAVECAR_OLD
    startTime=$(date +%s%3N)
    ./firstTimeStep.exe >> outputFile.txt
    if [ $? -ne 0 ]; then
	echo "firstTimeStep.exe crashed. Terminating twoLevelSystem.sh" >> outputFile.txt
	exit
    fi
    endTime=$(date +%s%3N)
    runTime=$((endTime-startTime))
    echo '' >> outputFile.txt
    echo "firstTimeStep.exe runtime:" >> outputFile.txt
    echo "$runTime ms" >> outputFile.txt
    gfortran -O3 updateWAVECAR.f -o updateWAVECAR.exe
fi





########################################################################
#                     Iterate forward through time                     #
########################################################################

for i in {1..100000}
do
    j=$(printf "%.5d" "$i")
    echo -e "\n\n\nt$j" >> outputFile.txt
    echo "-----------------------------------------------------------------------------------------------------------------------------" >> outputFile.txt
    cp WAVECAR WAVECAR_OLD
    #cp WAVECAR WAVECAR_$j
    startTime=$(date +%s%3N)
    mpiexec -machinefile $PBS_NODEFILE -np $NUM_PROC vasp_gam >> outputFile.txt # CCAST
    #srun -n 64 vasp_gam >> outputFile.txt # NERSC
    if [ $? -ne 0 ]; then
	echo "VASP crashed. Terminating twoLevelSystem.sh" >> outputFile.txt
	exit
    fi
    endTime=$(date +%s%3N)
    runTime=$((endTime-startTime))
    echo '' >> outputFile.txt
    echo "VASP runtime:" >> outputFile.txt
    echo "$runTime ms" >> outputFile.txt
    cp CONTCAR POSCAR
    grep 'F   F   T' CONTCAR >> positions.txt # This line is only useful for H2+ with one ion allowed to move along the bond axis. It is necessary to track internuclear distance vs time
    ##### This block selectively grabs and writes the forces on the nuclei from OUTCAR #####
    writeLine=false
    cat OUTCAR | while read line
    do
	if [[ "$line" =~ "FORCES acting on ions" ]]; then
	    writeLine=true
	    echo "" >> outputFile.txt
	fi
	if [[ "$writeLine" = true ]]; then
	    if [[ ! "$line" =~ [^[:space:]] ]]; then
		echo "" >> outputFile.txt
		break
	    fi
	    echo "$line" >> outputFile.txt
	fi
    done
    grep " 1 T= " OSZICAR >> outputFile.txt # Track temperature
    echo "" >> outputFile.txt
    if [[ $CHGint != 0 && $(($i%$CHGint)) = 0 ]]; then
	chgFile="CHG$j"
	cp CHG CHG_Dir/$chgFile
    fi
    startTime=$(date +%s%3N)
    ./updateWAVECAR.exe >> outputFile.txt
    if [ $? -ne 0 ]; then
       echo "updateWAVECAR.exe crashed. Terminating twoLevelSystem.sh" >> outputFile.txt
       exit
    fi
    endTime=$(date +%s%3N)
    runTime=$((endTime-startTime))
    echo '' >> outputFile.txt
    echo "updateWAVECAR.exe runtime:" >> outputFile.txt
    echo "$runTime ms" >> outputFile.txt
    echo "-----------------------------------------------------------------------------------------------------------------------------" >> outputFile.txt
done


