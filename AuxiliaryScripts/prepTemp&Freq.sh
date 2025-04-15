#!/bin/bash

dirArray=($(ls -d -- */))
for dir in ${dirArray[@]}; do
    echo "Entering directory $dir"
    cd $dir
    if [ -f outputFile.txt ]; then
	echo "outputFile.txt found"
	echo "Extracting times and resonance frequencies..."
	mapfile -t timeArray < <(grep "time[a-zA-Z]*[[:space:]]*=" outputFile.txt)
	times=()
	for line in "${timeArray[@]}"; do
	    timeStr=$(sed "s/[[:space:]a-zA-Z]*time[[:space:]a-zA-Z]*=[[:space:]]*[^0-9.]//g" <<< "$line")
	    times+=($timeStr)
	done
	mapfile -t omegaArray < <(grep 'omega_0' outputFile.txt)
	omegas=()
	for line in "${omegaArray[@]}"; do
	    omegaStr=$(sed "s/[[:space:]a-zA-Z]*omega_0[[:space:]a-zA-Z]*=[[:space:]]*[^0-9.]//g" <<< "$line")
	    omegas+=($omegaStr)
	done
	> resonanceVtime.txt
	cat resonanceVtime.txt
	if [ ${#times[@]} -ge ${#omegas[@]} ]; then
	    maxLen=${#omegas[@]}
	else
	    maxLen=${#times[@]}
	fi
	for ((i=0; i<$maxLen; i++)); do
	    echo "${times[$i]} ${omegas[$i]}" >> resonanceVtime.txt
	done
	mapfile -t tmpArray < <(grep ' 1 T= ' outputFile.txt)
	if [ ${#tmpArray[@]} -eq 0 ]; then
	    echo "No temperature data found"
	else
	    echo "Extracting temperatures..."
	    temperatures=()
	    for ((ind=0; ind<${#tmpArray[@]}; ind+=2)); do
		line=${tmpArray[$ind]}
		trimBeg=$(sed "s/[[:space:]]\+1[[:space:]]T=[[:space:]]\+//g" <<< "$line")
		trimEnd=$(sed "s/[[:space:]]\+E=[[:space:][:alnum:][:punct:]]\+//g" <<< "$trimBeg")
		temperatures+=("$trimEnd")
	    done
	    if [ ${#times[@]} -ge ${#temperatures[@]} ]; then
		maxLen=${#temperatures[@]}
	    else
		maxLen=${#times[@]}
	    fi
	    > tempVtime.txt
	    for ((i=0; i<$maxLen; i++)); do
		echo "${times[$i]} ${temperatures[$i]}" >> tempVtime.txt
	    done
	fi
	if [ -f positions.txt ]; then
	    echo "positions.txt found"
	    echo "Extracting positions..."
	    positions=()
	    while read line; do
		positionStr=$(sed "s/\([[:space:]]*[0-9.]\+[[:space:]]\+[0-9.]\+[[:space:]]\+\)\|\([[:space:]]\+[TFtf][[:space:]]\+[TFtf][[:space:]]\+[TFtf][[:space:]]*\)//g" <<< "$line") # THIS IS SET UP TO READ THE Z COMPONENT
		positions+=($positionStr)
	    done < positions.txt
	    if [ ${#times[@]} -ge ${#positions[@]} ]; then
		maxLen=${#positions[@]}
	    else
		maxLen=${#times[@]}
	    fi
	    > positionVtime.txt
	    for ((i=0; i<$maxLen; i++)); do
		echo "${times[$((i+1))]} ${positions[$i]}" >> positionVtime.txt
	    done
	else
	    echo "positons.txt NOT found"
	fi
    else
	echo "outputFile.txt NOT found"
    fi
    subdirArray=($(ls -d -- */ 2>/dev/null))
    if [[ ${#subdirArray[@]} -eq 0 ]]; then
	echo "no subdirectories found within $dir"
    else
	for subdir in ${subdirArray[@]}; do
	    echo "Entering subdirectory $subdir"
	    cd $subdir
	    if [ -f outputFile.txt ]; then
		echo "outputFile.txt found"
		echo "Extracting times and resonance frequencies..."
		mapfile -t timeArray < <(grep "time[a-zA-Z]*[[:space:]]*=" outputFile.txt)
		times=()
		for line in "${timeArray[@]}"; do
		    timeStr=$(sed "s/[[:space:]a-zA-Z]*time[[:space:]a-zA-Z]*=[[:space:]]*[^0-9.]//g" <<< "$line")
		    times+=($timeStr)
		done
		mapfile -t omegaArray < <(grep 'omega_0' outputFile.txt)
		omegas=()
		for line in "${omegaArray[@]}"; do
		    omegaStr=$(sed "s/[[:space:]a-zA-Z]*omega_0[[:space:]a-zA-Z]*=[[:space:]]*[^0-9.]//g" <<< "$line")
		    omegas+=($omegaStr)
		done
		> resonanceVtime.txt
		cat resonanceVtime.txt
		if [ ${#times[@]} -ge ${#omegas[@]} ]; then
		    maxLen=${#omegas[@]}
		else
		    maxLen=${#times[@]}
		fi
		for ((i=0; i<$maxLen; i++)); do
		    echo "${times[$i]} ${omegas[$i]}" >> resonanceVtime.txt
		done
		mapfile -t tmpArray < <(grep ' 1 T= ' outputFile.txt)
		if [ ${#tmpArray[@]} -eq 0 ]; then
		    echo "No temperature data found"
		else
		    echo "Extracting temperatures..."
		    temperatures=()
		    for ((ind=0; ind<${#tmpArray[@]}; ind+=2)); do
			line=${tmpArray[$ind]}
			trimBeg=$(sed "s/[[:space:]]\+1[[:space:]]T=[[:space:]]\+//g" <<< "$line")
			trimEnd=$(sed "s/[[:space:]]\+E=[[:space:][:alnum:][:punct:]]\+//g" <<< "$trimBeg")
			temperatures+=("$trimEnd")
		    done
		    if [ ${#times[@]} -ge ${#temperatures[@]} ]; then
			maxLen=${#temperatures[@]}
		    else
			maxLen=${#times[@]}
		    fi
		    > tempVtime.txt
		    for ((i=0; i<$maxLen; i++)); do
			echo "${times[$i]} ${temperatures[$i]}" >> tempVtime.txt
		    done
		fi
		if [ -f positions.txt ]; then
		    echo "positions.txt found"
		    echo "Extracting positions..."
		    positions=()
		    while read line; do
			positionStr=$(sed "s/\([[:space:]]*[0-9.]\+[[:space:]]\+[0-9.]\+[[:space:]]\+\)\|\([[:space:]]\+[TFtf][[:space:]]\+[TFtf][[:space:]]\+[TFtf][[:space:]]*\)//g" <<< "$line") # THIS IS SET UP TO READ THE Z COMPONENT
			positions+=($positionStr)
		    done < positions.txt
		    if [ ${#times[@]} -ge ${#positions[@]} ]; then
			maxLen=${#positions[@]}
		    else
			maxLen=${#times[@]}
		    fi
		    > positionVtime.txt
		    for ((i=0; i<$maxLen; i++)); do
			echo "${times[$((i+1))]} ${positions[$i]}" >> positionVtime.txt
		    done
		else
		    echo "positons.txt NOT found"
		fi
	    else
		echo "outputFile.txt NOT found"
	    fi
	    echo "Exiting subdirectory $subdir"
	    cd ..
	done
    fi
    echo "Exiting directory $dir"
    cd ..
done
