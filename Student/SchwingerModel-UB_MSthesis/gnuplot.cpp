reset 
set terminal gif animate delay 1 size 1000,420 
set output "output.gif"

unset tics
set xlabel "x"
set ylabel "Fields(t,x)"

set xrange [0:4000]
set yrange [-0.005:0.005]


do for [i = 0:99]{
	plot "Sky".i.".txt" using 1:2 t "L_{real}" ,"Sky".i.".txt" using 1:3 t "L_{imag}" , "Sky".i.".txt" using 1:4 t "|L|", "Sky".i.".txt" using 1:5 t "R_{real}", "Sky".i.".txt" using 1:6 t "R_{imag}", "Sky".i.".txt" using 1:7 t "|R|", "Sky".i.".txt" using 1:9 t "E"
}

//SPLOT ^

//----------------------------------------------------------------------------------------------------

//PLOT  v

reset 
set terminal gif animate delay 1 size 1000,420 
set output "output.gif"

set xlabel "x"
set ylabel "Fields(t,x)"
set xrange [0:4000]
set yrange [-0.03:0.03]

do for [i = 0:99]{
	plot "Sky".i.".txt" using 1:2 t "Psi_L_real" , "Sky".i.".txt" using 1:3 t "Psi_L_imag" , "Sky".i.".txt" using 1:4 t "Psi_R_real", "Sky".i.".txt" using 1:5 t "Psi_R_imag","Sky".i.".txt" using 1:7 t "E"
}


