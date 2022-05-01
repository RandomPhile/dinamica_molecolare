set datafile sep '\t'
set grid

#N N_t L

set term qt pos 0,300
set view equal xyz
set xrange [-1:L]
set yrange [-1:L]
set zrange [-1:L]

set style fill transparent solid 0.5 border border lc rgb 'dark-green'
set object 1 polygon from 0,0,0 to L,0,0 to L,L,0 to 0,L,0 fillcolor rgb 'light-green'
set object 2 polygon from 0,0,0 to L,0,0 to L,0,L to 0,0,L fillcolor rgb 'light-green'
set object 3 polygon from 0,0,0 to 0,L,0 to 0,L,L to 0,0,L fillcolor rgb 'light-green'
set object 4 polygon from L,0,0 to L,L,0 to L,L,L to L,0,L fillcolor rgb 'light-green'
set object 5 polygon from L,L,0 to 0,L,0 to 0,L,L to L,L,L fillcolor rgb 'light-green'
set object 6 polygon from 0,0,L to L,0,L to L,L,L to 0,L,L fillcolor rgb 'light-green'


n = 1
do for [step=1:N_t] {
	#set output sprintf('png/animation%03.0f.png',step)
	splot "coordinate.xyz" every ::n::n+N-1 u 2:3:4 pt 7 ps 1 title sprintf('step = %03.0f    t = %f',step,(n-1)*dt/N),\
		  "coordinate.xyz" every ::n::n+N-1 u 2:3:4:1 w labels offset 2 title ""
	n=n+N*skip
	pause pausa
}


pause-1