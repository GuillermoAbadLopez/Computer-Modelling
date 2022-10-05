reset

set terminal pngcairo background rgb 'black'
set xlabel 'ylabel' tc rgb 'white'
set ylabel 'xlabel' tc rgb 'white'
set border lc rgb 'white'
set key tc rgb 'white

# png
set terminal pngcairo size 1000,420

# color definitions
for [i=0:19] set style line (i+1) lc rgb '#fde725' lt 1 lw 2 pt 7 ps 2 # --- blue


set parametric
unset key

#set xrange [-1:1.5]
#set yrange [-1.5:2]
#set zrange [-1:2]

set xlabel "x"
set ylabel "y"
set zlabel "z"
system('mkdir -p png')


n=0
do for [ii=1:126] {
    n=n+1
    set output sprintf('png/gif%03.0f.png',n)
    splot \
    	for [i=1:19] "animacion_cuerpo".i.".txt" every ::1::ii w l ls (i+1)
          
}






n=0
do for [ii=1:126] {
    n=n+1
    set output sprintf('png/gif%03.0f.png',n)
    splot 'animacion_cuerpo0.txt' every ::1::ii w l ls 1, \
          'animacion_cuerpo0.txt' every ::ii::ii w p ls 1, \
		  'animacion_cuerpo1.txt' every ::1::ii w l ls 2, \
          'animacion_cuerpo1.txt' every ::ii::ii w p ls 2         
}




