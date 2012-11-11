unset multiplot
reset
set terminal push

set encoding cp1250 # make minus sign longer
#set encoding iso_8859_1
set terminal postscript enhanced font "Arial, 12" size 10, 3.5
set output "ens.ps"

set style line 1 lc rgb "#ff0044" lt 1 lw 4.0
set style line 2 lc rgb "#22cc00" lt 2 lw 4.0
set style line 3 lc rgb "#0022ff" lt 4 lw 4.0
set style line 4 lc rgb "#000000" lt 1 lw 2.0

# draw a gamma distribution
n = 100;
p(x) = exp(-x-lgamma(n))*x**(n-1);
a = 10;
b = 190;
flat(x) = ((x>a && x < b) ? 1.0/(b - a) : 0.0);

dx = 0.01
dy = 0.05

set label 1 "(a)" at screen dx, 1-dy
set label 2 "(b)" at screen .5+dx, 1-dy

set tmargin 2

set multiplot
set size 0.5, 1.0
set origin 0.0, 0.0

unset xtics
unset ytics
set xlabel "{/Arial-Italic E}"
set ylabel "{/Symbol-Oblique r}&{/=8 .}({/Arial-Italic E}&{/=8 .})"
set samples 1000

set key Left reverse spacing 1.5 width -2

set title "Energy distribution" offset 0, -0.5

set parametric
set trange [0:1]
set yrange [0:0.07]

plot 100, 0.065*t       ls 1 t "Microcanonical", \
     t*200,p(t*200)     ls 2 t "Canonical", \
     t*200,flat(t*200)  ls 3 t "Multicanonical"


set origin 0.5, 0.0

set title "Sampling weight, exp[-{/Arial-Italic f}{/=6 &.}(E{/=6 &.})]"

set ylabel "-{/Arial-Italic f}{/=6 &.}(E{/=6 &.})"

set xrange [0:1]
set yrange [-0.2:1.1]

# {\303} is the hat ^,  ~{}{}

set key width -6

plot 1 - t, t - .1    ls 2 t 'Canonical, -{/Symbol-Oblique b}{/=6 &.}{/Arial-Italic E}', \
     "ens.txt" u 1:2  w lp ls 3 t 'Multicanonical, -log{/=6 &.}~{/Symbol W}{/=8{.8\^}}{/=6 &.}({/Arial-Italic E}{/=6 &.})', \
     1 - t, t * t     ls 4 t '-log{/=6 &.}{/Symbol W}{/=6 &.}({/Arial-Italic E}{/=6 &.})'

unset multiplot

unset output

#set terminal windows enhanced
set terminal pop
reset



