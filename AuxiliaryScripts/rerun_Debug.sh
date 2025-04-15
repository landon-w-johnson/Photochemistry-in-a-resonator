#!/bin/bash

cp updateWAVECAR_BACKUP.f updateWAVECAR.f
cp WAVECAR_BACKUP WAVECAR_OLD
clear
gfortran firstTimeStep.f -o firstTimeStep.exe
./firstTimeStep.exe
