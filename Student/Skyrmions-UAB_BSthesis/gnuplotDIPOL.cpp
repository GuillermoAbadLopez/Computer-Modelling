reset 
set terminal gif animate delay 1 size 1000,420 
set output "output.gif"

unset key
set xlabel "x"
set ylabel "M(x,y)" 
set xrange [0:200]
set yrange [0:10]

do for [i = 0:99]{
	plot "Sky".i.".txt" using 1:2
}



//($N*2) multiplica la longitud dels vectors per 2, every 2:2 plotea cada dos vectors, -1 is black and filled la punta plena

//SMOOOTHHH
set pm3d map
set pm3d interpolate 0,0 
splot "Sky".i.".txt" matrix
