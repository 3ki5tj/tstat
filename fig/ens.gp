unset multiplot
reset
set terminal push

set encoding cp1250 # make minus sign longer
#set encoding iso_8859_1
set terminal postscript enhanced font "Arial, 24" 
set output "ens.ps"

set style line 1 lc rgb "#ff0044" lt 1 
set style line 2 lc rgb "#22cc00" lt 2
set style line 3 lc rgb "#2200cc" lt 4


# draw a gamma distribution
n = 100;
p(x) = exp(-x-lgamma(n))*x**(n-1);
a = 10;
b = 190;
flat(x) = ((x>a && x < b) ? 1.0/(b - a) : 0.0);

unset xtics
unset ytics
set xlabel "{/Arial-Italic E}"
set ylabel "{/Symbol-Oblique r}&{/=8 .}({/Arial-Italic E}&{/=8 .})"
set samples 1000

set key spacing 1.2

set parametric
set trange [0:1]
set yrange [0:0.07]

plot 100, 0.065*t       ls 1 lw 4.0 t "Microcanonical", \
     t*200,p(t*200)     ls 2 lw 4.0 t "Canonical", \
     t*200,flat(t*200)  ls 3 lw 4.0 t "Multicanonical"

unset output

#set terminal windows enhanced
set terminal pop
reset



