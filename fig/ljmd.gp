unset multiplot
reset

set encoding cp1250 # make minus sign longer
#set encoding iso_8859_1
set terminal postscript enhanced 
set output "ljmd.ps"

set multiplot

set lmargin 3.5
set rmargin 1.0
set tmargin 2.5
set bmargin 2.0
tcfont="Arial, 10"

N=108
set mytics 5
set ytics 0.1 font tcfont offset 0.3, 0 
set mxtics 4
set xtics 2 font tcfont offset 0, 0.3
set xlabel "{/Arial-Italic E}" offset 0, 1
set key spacing 1.5

dx = 0.01
dy = 0.03

set label 1 "(a)" at screen dx, 1-dy
set label 2 "(b)" at screen dx, .5-dy
set label 3 "(c)" at screen .3333+dx, 1-dy
set label 4 "(d)" at screen .3333+dx, .5-dy
set label 5 "(e)" at screen .6667+dx, 1-dy
set label 6 "(f)" at screen .6667+dx, .5-dy

# half of the bin size
dx2 = 0.02

set size 0.3333, 0.5

set origin 0.0, 0.5

set title "Langevin dynamics" offset 0, -0.2

plot [][0:] \
  "../data/lang.001d.001/avb.dat" u ($1/N):($5*N) w l t "{/Symbol-Oblique z} = 1", \
  "../data/lang.003d.001/avb.dat" u ($1/N):($5*N) w l t "{/Symbol-Oblique z} = 3", \
  "../data/lang.01d.001/avb.dat" u ($1/N):($5*N) w l t "{/Symbol-Oblique z} = 10", \
  "../data/langbad/avb.dat" u ($1/N):($5*N) w l t "Unmodified", \


set origin 0.0, 0.0

set title "Anderson"

plot [][0:] \
  "../data/ands.4/avb.dat" u ($1/N):($5*N) w l t "{/Symbol-Oblique b}_m = 2.5", \
  "../data/ands.1/avb.dat" u ($1/N):($5*N) w l t "{/Symbol-Oblique b}_m = 10", \
  "../data/ands0/avb.dat" u ($1/N):($5*N) w l t "Unmodified"


set origin 0.3333, 0.5

set title "Nose-Hoover"

plot [][0:] \
  "../data/nh100d.001/avb.dat" u ($1/N):($5*N) w l t "{/Arial-Italic Q} = 100", \
  "../data/nh300d.001/avb.dat" u ($1/N):($5*N) w l t "{/Arial-Italic Q} = 300", \
  "../data/nh1000d.001/avb.dat" u ($1/N):($5*N) w l t "{/Arial-Italic Q} = 1000", \
  "../data/nh300bad/avb.dat" u ($1/N):($5*N) w l t "Unmodified"

set origin 0.3333, 0.

set title "Velocity-Rescaling"

plot [][0:] \
  "../data/vr.1d.001/avb.dat" u ($1/N):($5*N)  w l t "{/Symbol-Oblique D}{/Arial-Italic t}_{/=8 VR} = 0.1", \
  "../data/vr.01d.001/avb.dat" u ($1/N):($5*N)  w l t "{/Symbol-Oblique D}{/Arial-Italic t}_{/=8 VR} = 0.01", \
  "../data/vr.001d.001/avb.dat" u ($1/N):($5*N) w l t "{/Symbol-Oblique D}{/Arial-Italic t}_{/=8 VR} = 0.001", \
  "../data/vrbad/avb.dat" u ($1/N):($5*N) w l t "Unmodified" 


set origin 0.6667, 0.5

set title "Nose-Hoover, {/Arial-Italic Q} = 30"

plot [][0:0.4] \
  "../data/nh30d.002/avb.dat" u ($1/N):($5*N)  w l t "{/Symbol-Oblique D}{/Arial-Italic t}_{/=8 MD} = 0.002", \
  "../data/nh30d.001/avb.dat" u ($1/N):($5*N)  w l t "{/Symbol-Oblique D}{/Arial-Italic t}_{/=8 MD} = 0.001", \
  "../data/nh30d.0005/avb.dat" u ($1/N):($5*N) w l t "{/Symbol-Oblique D}{/Arial-Italic t}_{/=8 MD} = 0.0005", \
  "../data/nh30d.0002/avb.dat" u ($1/N):($5*N) w l t "{/Symbol-Oblique D}{/Arial-Italic t}_{/=8 MD} = 0.0002" 

set origin 0.6667, 0.

set title "1/{/Symbol-Oblique b}({/Arial-Italic E}) = 2 {/Symbol \341}{/Arial-Italic K}{/Symbol \361} / {/Arial-Italic N_f}"

set ytics nomirror
set my2tics 5
set y2tics 10 font tcfont offset -0.5, 0
set y2range [0:50]
set rmargin 2.5

plot [][0:] \
  "../data/mcx.3rat/avb.dat" u ($1/N):($5*N)    w l t "{/Symbol-Oblique r}({/Arial-Italic E})", \
  "../data/mcx.3rat/avb.dat" u ($1/N):($5*N*$6) axes x1y2 w l t "{/Symbol-Oblique r}({/Arial-Italic E}) {/Symbol \341}{/Arial-Italic K}{/Symbol \361}"




unset multiplot
unset output

set terminal wxt
reset



