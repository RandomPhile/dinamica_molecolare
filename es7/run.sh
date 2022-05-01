#!/bin/bash
g++ main.cpp -o main.out
./main.out
pkill -x gnuplot, gnuplot_qt
#gnuplot 'plot.plt'
i=0
while read line; do    
	param[$i]=$line
	i=$i+1
done < gnuplot.dat
gnuplot -e "N=${param[0]}" -e "N_t=${param[1]}" -e "L=${param[2]}" -e "pausa=${param[3]}" -e "skip=${param[4]}" -e "dt=${param[5]}" "animation.plt"