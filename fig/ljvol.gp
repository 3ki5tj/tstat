unset multiplot
reset

set encoding cp1250 # make minus sign longer
#set encoding iso_8859_1
set terminal postscript enhanced size 10, 4 lw 2.0
set output "ljvol.ps"

set multiplot

# number of particles
N = 108
# half of the bin size
hbin = 0.001


set lmargin 5.5
set rmargin 1.5
set tmargin 1.5
set bmargin 3.
tcfont="Arial, 10"

set mytics 5
set ytics 0.5 font tcfont offset 0.5, 0 
set mxtics 4
set xtics .1 font tcfont offset 0, 0.3
set key spacing 1.5

dx = 0.01
dy = 0.03

set linestyle 1 lc rgb "#ff0044" lt 1 
set linestyle 2 lc rgb "#00cc00" lt 2
set linestyle 3 lc rgb "#2200cc" lt 4


set label 1 "(a)" at screen dx, 1-dy
set label 2 "(b)" at screen .5+dx, 1-dy



set size 0.5, 1.0

set origin 0.0, 0.0

#set title ""
set xlabel "{/Arial-Italic n}" offset 0, 0.5
set ylabel "{/Symbol-Oblique r}({/Arial-Italic n})" offset 2.5, 0
set key right bottom

plot [.05:.65][0:2] \
  "../data/ljvol2b/avp.dat" u ($1 + hbin):($5) w l ls 1 t "Nose-Hoover", \
  "../data/ljvol3b/avp.dat" u ($1 + hbin):($5) w l ls 2 t "Langevin", \
  "../data/ljvol1b/avp.dat" u ($1 + hbin):($5) w l ls 3 t "Monte-Carlo"


set origin 0.5, 0.0

set ytics 0.1 
set ylabel "{/Arial-Italic p}({/Arial-Italic n})" offset 2.5, 0
set key left top
set ytics nomirror
set my2tics 5
set y2tics 5 font tcfont offset -0.7, 0
set y2label "{/Arial-Italic F}({/Arial-Italic n})" offset -3.0, 0

set rmargin 5.0

set label "{/Arial-Italic p} = 0.0696" at 0.43, 0.09
plot [.05:.65][] \
  0.0696 w l lt 4 not, \
  "../data/ljvol2b/avp.dat" u ($1 + hbin):($2) w l ls 1 t "Nose-Hoover", \
  "../data/ljvol3b/avp.dat" u ($1 + hbin):($2) w l ls 2 t "Langevin", \
  "../data/ljvol1b/avp.dat" u ($1 + hbin):($2) w l ls 3 t "Monte-Carlo", \
  "../data/ljvol2b/fe.dat"  u ($1):($2 + 0.06940*N/$1 - 139.2) axes x1y2 w l ls 1 not, \
  "../data/ljvol3b/fe.dat"  u ($1):($2 + 0.06967*N/$1 - 139.6) axes x1y2 w l ls 2 not, \
  "../data/ljvol1b/fe.dat"  u ($1):($2 + 0.06963*N/$1 - 139.5) axes x1y2 w l ls 3 not


unset multiplot
unset output

set terminal wxt
reset



