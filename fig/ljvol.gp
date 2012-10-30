unset multiplot
reset

set encoding cp1250 # make minus sign longer
#set encoding iso_8859_1
set terminal postscript enhanced size 10, 4
set output "ljvol.ps"

set multiplot

set lmargin 5.5
set rmargin 1.5
set tmargin 1.5
set bmargin 3.
tcfont="Arial, 10"

N=108
set mytics 5
set ytics 0.5 font tcfont offset 0.3, 0 
set mxtics 4
set xtics .1 font tcfont offset 0, 0.3
set key spacing 1.5

dx = 0.01
dy = 0.03

set label 1 "(a)" at screen dx, 1-dy
set label 2 "(b)" at screen .5+dx, 1-dy


set size 0.5, 1.0

set origin 0.0, 0.0

#set title ""
set xlabel "{/Arial-Italic V}" offset 0, 0.5
set ylabel "{/Symbol-Oblique r}({/Arial-Italic V})" offset 2.5, 0
set key right bottom

plot [][0:2] \
  "../data/ljvol2b/avp.dat" u ($1):($5) w l t "Nose-Hoover", \
  "../data/ljvol3b/avp.dat" u ($1):($5) w l t "Langevin", \
  "../data/ljvol1b/avp.dat" u ($1):($5) w l t "Monte-Carlo"


set origin 0.5, 0.0

set ytics 0.1 
set ylabel "{/Arial-Italic p}({/Arial-Italic V})" offset 2.5, 0
set key left top
plot [][] \
  "../data/ljvol2b/avp.dat" u ($1):($2) w l t "Nose-Hoover", \
  "../data/ljvol3b/avp.dat" u ($1):($2) w l t "Langevin", \
  "../data/ljvol1b/avp.dat" u ($1):($2) w l t "Monte-Carlo"


unset multiplot
unset output

set terminal wxt
reset



