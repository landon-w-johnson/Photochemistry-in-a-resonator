outFileName = 
set output outFileName
set terminal png size 3000,1500
set title 'Ground State'
set ylabel 'Population'
set xlabel 'time (fs)'
plot 'dataFile.txt' u 1:6 w l title 'Population'