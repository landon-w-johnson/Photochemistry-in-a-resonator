outFileName = 
set output outFileName
set terminal png size 3000,1500
set title 'Excited State'
set ylabel 'Population'
set xlabel 'time (fs)'
stats 'dataFile.txt' u 7
plot 'dataFile.txt' u 1:7 w l title 'Population', 'dataFile.txt' u 1:(($14+1)*STATS_max/2) w l title 'Relative E-Field Amplitude'