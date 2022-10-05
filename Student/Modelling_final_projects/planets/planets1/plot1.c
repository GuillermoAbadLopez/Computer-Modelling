reset

# png
set terminal pngcairo size 1000,420

# color definitions
set style line 1 lc rgb '#0c0887' lt 1 lw 2 pt 7 ps 2 # --- blue
set style line 2 lc rgb '#e56b5d' lt 1 lw 2 pt 7 ps 2 # --- red
set style line 3 lc rgb '#fde725' lt 1 lw 2 pt 7 ps 2 # --- yellow
set style line 4 lc rgb '#000004' lt 1 lw 2 pt 7 ps 2 # --- black

set parametric
unset key

set xrange [-1:1.5]
set yrange [-1:2]
set zrange [-1.0:2.0]

set xlabel "x"
set ylabel "y"
set zlabel "z"
system('mkdir -p png')


n=0
do for [ii=1:126] {
    n=n+1
    set output sprintf('png/gif%03.0f.png',n)
    splot 'animacion_cuerpo0.txt' every ::1::ii w l ls 1, \
          'animacion_cuerpo0.txt' every ::ii::ii w p ls 1, \
		  'animacion_cuerpo1.txt' every ::1::ii w l ls 2, \
          'animacion_cuerpo1.txt' every ::ii::ii w p ls 2, \
          'animacion_cuerpo2.txt' every ::1::ii w l ls 3, \
          'animacion_cuerpo2.txt' every ::ii::ii w p ls 3, \
		  'animacion_cuerpo3.txt' every ::1::ii w l ls 4, \
          'animacion_cuerpo3.txt' every ::ii::ii w p ls 4
}


