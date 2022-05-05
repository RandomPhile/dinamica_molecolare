set datafile sep '\t'
set grid

set term qt pos 0,300
set view equal xyz
padding = L+L/2
set xrange [-padding:L+padding]
set yrange [-padding:L+padding]
set zrange [-padding:L+padding]

set style fill transparent solid 0.5 border border lc rgb 'dark-green'
set object 1 polygon from 0,0,0 to L,0,0 to L,L,0 to 0,L,0 to 0,0,0 fillcolor rgb 'light-green'
set object 2 polygon from 0,0,0 to L,0,0 to L,0,L to 0,0,L to 0,0,0 fillcolor rgb 'light-green'
set object 3 polygon from 0,0,0 to 0,L,0 to 0,L,L to 0,0,L to 0,0,0 fillcolor rgb 'light-green'
set object 4 polygon from L,0,0 to L,L,0 to L,L,L to L,0,L to L,0,0 fillcolor rgb 'light-green'
set object 5 polygon from L,L,0 to 0,L,0 to 0,L,L to L,L,L to L,L,0 fillcolor rgb 'light-green'
set object 6 polygon from 0,0,L to L,0,L to L,L,L to 0,L,L to 0,0,L fillcolor rgb 'light-green'
T(x,d) = x+d

array Xs[8] = [0,L,L,L,L,L,L,L]
array Ys[8] = [0,0,0,0,0,0,0,0]
array Zs[8] = [0,0,0,0,0,0,0,0]

n = 0
do for [step=1:1] {
	do for [i=1:8] {
		 "coordinate.xyz" every ::n::n+N u (T($2,i)):3:4 pt 7 ps 0.5 lc rgb "red",\
	}
	n=n+N*skip
	pause pausa
}


pause-1