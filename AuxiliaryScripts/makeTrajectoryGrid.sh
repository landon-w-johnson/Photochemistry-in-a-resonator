#!/bin/bash

omega_0=4188069

for amplitude in 0.00001 #0.000001 #0.0001 0.00001 0.000001 0.0000001
do
    for ratio in 1000 #{995..1004} #{990..1010}
    do
	for dt in 0.01 #0.05 0.1 0.3 0.5
	do
	    omega=$((omega_0*ratio))
	    if [[ ratio -ge 1000 ]];
	    then
		subtractedRatio=$((ratio - 1000))
		modifiedRatio=$(printf "%.3d" "$subtractedRatio")
		dirName="E${amplitude}_wRat1.${modifiedRatio}_dt${dt}"
	    else
		dirName="E${amplitude}_wRat0.${ratio}_dt${dt}"
	    fi
	    
	    if [ -e "$dirName" ];
	    then
		echo "directory $dirName found, skipping creation"
	    else
		mkdir $dirName
		cp * $dirName    
	    fi
	    cd $dirName
	    #sed -i "s|#PBS -N V5|#PBS -N Euler_NoRWA_$dirName|g" job.pbs #CCAST
	    sed -i "s/#SBATCH --job-name=[^[:space:]]\+/#SBATCH --job-name=$dirName/g" submit-NERSC.sh #NERSC
	    wd=$(pwd)
	    scriptName=${wd}/twoLevelSystem.sh
	    #NERSC: nothing needs to be done
	    sed -i "s/E[[:space:]]*=[[:space:]]*\(\([0-9]\+[.]\?[0-9]*\)\|<tmp>\)/E = ${amplitude}/g" PHOTCAR
	    sed -i "s/omega[[:space:]]*=[[:space:]]*\(\([0-9]\+[.]\?[0-9]*\)\|<tmp>\)/omega = 0.0${omega}/g" PHOTCAR
	    sed -i "s/timeStep[[:space:]]*=[[:space:]]*\(\([0-9]\+[.]\?[0-9]*\)\|<tmp>\)/timeStep = ${dt}/g" PHOTCAR
	    sed -i "s/POTIM[[:space:]]*=[[:space:]]*\(\([0-9]\+[.]\?[0-9]*\)\|<tmp>\)/POTIM=${dt}/g" INCAR
	    gfortran -O3 firstTimeStep.f -o firstTimeStep.exe
	    #qsub job.pbs #CCAST
	    sbatch submit-NERSC.sh #NERSC
	    cd ..
	done
    done
done
