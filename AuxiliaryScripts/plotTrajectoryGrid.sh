#!/bin/bash



for amplitude in 0.0001 0.00001 0.000001 0.0000001
do
    for omega in {810..820}
    do
	dirName="E${amplitude}_w0.$omega"
	cd $dirName
	ES_outFileName="plot_ES_${dirName}.png"
	GS_outFileName="plot_GS_${dirName}.png"
	cp ../plotter_ES.gnuplot .
	cp ../plotter_GS.gnuplot .
	sed -i "s|outFileName = |outFileName = '$ES_outFileName'|g" plotter_ES.gnuplot
	sed -i "s|outFileName = |outFileName = '$GS_outFileName'|g" plotter_GS.gnuplot
	gnuplot plotter_ES.gnuplot
	gnuplot plotter_GS.gnuplot
	cd ..
    done
done
