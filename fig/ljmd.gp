unset multiplot
reset

set encoding cp1250 # make minus sign longer
#set encoding iso_8859_1
set terminal postscript enhanced lw 2.0
set output "ljmd.ps"

set multiplot

# the number of particles
N=108
# half of the bin size
hbin = 0.01

# width of the three columns
x1 = 0.35
x2 = 0.325
x3 = 1 - x1 - x2

set style line 1 lc rgb "#ff0044" lt 1
set style line 2 lc rgb "#228800" lt 2
set style line 3 lc rgb "#2200cc" lt 3
set style line 4 lc rgb "#cc0088" lt 4
set style line 5 lc rgb "#cc8800" lt 5
set style line 9 lc rgb "#888888" lt 1 lw 0.5


set lmargin 3.5
set rmargin 1.0
set tmargin 2.5
set bmargin 2.0
tcfont="Arial, 10"

set mytics 5
set ytics 0.1 font tcfont offset 0.3, 0
set mxtics 4
set xtics 2 font tcfont offset 0, 0.3
set xlabel "{/Arial-Italic E}" offset 0, 1
set key Left reverse spacing 1.15 font "Arial, 11"

dx = 0.01
dy = 0.03

set label 1 "(a)" at screen dx, 1-dy
set label 2 "(b)" at screen dx, .5-dy
set label 3 "(c)" at screen x1+dx, 1-dy
set label 4 "(d)" at screen x1+dx, .5-dy
set label 5 "(e)" at screen x1+x2+dx, 1-dy
set label 6 "(f)" at screen x1+x2+dx, .5-dy

set size x3, 0.5

set origin x1+x2, 0.5

set title "Langevin dynamics" offset 0, -0.2

plot [][0:.5] \
  "../data/lang.001/avb.dat" u ($1/N + hbin):($5*N) w l ls 1 t "{/Symbol-Oblique z} = 1", \
  "../data/lang.003/avb.dat" u ($1/N + hbin):($5*N) w l ls 2 t "{/Symbol-Oblique z} = 3", \
  "../data/lang.01/avb.dat" u ($1/N + hbin):($5*N) w l ls 3 t "{/Symbol-Oblique z} = 10", \
  "../data/langbad/avb.dat" u ($1/N + hbin):($5*N) w l ls 5 t "Unmodified"


set origin x1+x2, 0.0

set title "Andersen"

plot [][0:] \
  "../data/andsx/avb.dat"   u ($1/N + hbin):($5*N) w l ls 1 t "Eq. (11)", \
  "../data/ands0a/avb.dat"  u ($1/N + hbin):($5*N) w l ls 5 t "Unmodified"


set origin x1, 0.5

set title "Velocity-rescaling"

set key width -2

plot [][0:] \
  "../data/vr.1/avb.dat" u ($1/N + hbin):($5*N)  w l ls 1 t "{/Symbol-Oblique D}{/Arial-Italic t}_{/=8 VR} = 0.1", \
  "../data/vr.01/avb.dat" u ($1/N + hbin):($5*N)  w l ls 2 t "{/Symbol-Oblique D}{/Arial-Italic t}_{/=8 VR} = 0.01", \
  "../data/vr.001/avb.dat" u ($1/N + hbin):($5*N) w l ls 3 t "{/Symbol-Oblique D}{/Arial-Italic t}_{/=8 VR} = 0.001", \
  "../data/vrbad/avb.dat" u ($1/N + hbin):($5*N) w l ls 5 t "Unmodified"



set origin x1, 0.0

set title "MC velocity-rescaling"

set ytics nomirror
set my2tics 5
set y2tics 20 font tcfont offset -0.5, 0
set y2range [0:80]
set rmargin 2.5

set key width -9

plot [][0:0.5] \
  "../data/mcx1/avb.dat"     u ($1/N + hbin):($5*N)  w l ls 1 t "{/Symbol-Oblique D} = 1.0", \
  "../data/mcx.3/avb.dat"    u ($1/N + hbin):($5*N)  w l ls 2 t "{/Symbol-Oblique D} = 0.3", \
  "../data/mcx.1/avb.dat"    u ($1/N + hbin):($5*N)  w l ls 3 t "{/Symbol-Oblique D} = 0.1", \
  "../data/mcx.3rat/avb.dat" u ($1/N + hbin):($5*N)  w l ls 5 t '~{/Symbol-Oblique b}{.3\~}{/=6 &.}({/Arial-Italic E}) = {/Arial-Italic N_f} / {/Symbol \341}{/Arial-Italic 2{/=6 &.}K{/=6 &.}}{/Symbol \361}_{/Arial-Italic=8 E}', \
  "../data/mcx.3rat/avb.dat" u ($1/N + hbin):($5*N*$6) axes x1y2 w l ls 4 \
      t '~{/Symbol-Oblique b}{.3\~}{/=6 &.}({/Arial-Italic E}), {/Symbol-Oblique r}{/=6 &.}({/Arial-Italic E}) {/Symbol \341}{/Arial-Italic K}{/Symbol \361}_{/Arial-Italic=8 E}'

unset y2tics
unset my2tics
set ytics mirror
set rmargin 1.0
# "1/{/Symbol-Oblique b}({/Arial-Italic E}) = 2 {/Symbol \341}{/Arial-Italic K}{/Symbol \361}_{/Arial-Italic=8 E} / {/Arial-Italic N_f}", \



set size x1, 0.5
set lmargin 5.5
set origin 0.0, 0.5

set title "Nose-Hoover chain"

set key width -1

set ylabel "{/Symbol-Oblique r}{/=6 &.}({/Arial-Italic E}{/=6 &.})" offset 3, 0

plot [][0.:0.35] \
  "../data/nhc100/avb.dat"    u ($1/N + hbin):($5*N) w l ls 1 t "{/Arial-Italic Q}_1 = 100", \
  "../data/nhc300/avb.dat"    u ($1/N + hbin):($5*N) w l ls 2 t "{/Arial-Italic Q}_1 = 300", \
  "../data/nhc1000/avb.dat"   u ($1/N + hbin):($5*N) w l ls 3 t "{/Arial-Italic Q}_1 = 1000", \
  "../data/nhc300bad/avb.dat" u ($1/N + hbin):($5*N) w l ls 5 t "Unmodified"


#  "../data/nh300/avb.dat" u ($1/N + hbin):($5*N) w l ls 4 t "NH, {/Arial-Italic Q} = 300"

set origin 0., 0.


set title "Nose-Hoover chain, {/Arial-Italic Q} = 10"
set key width -2

plot [][0:0.35] \
  "../data/nhc10d.002/avb.dat"  u ($1/N + hbin):($5*N) w l ls 1 t "{/Symbol-Oblique D}{/Arial-Italic t}_{/=8 MD} = 0.002", \
  "../data/nhc10d.0005/avb.dat" u ($1/N + hbin):($5*N) w l ls 2 t "{/Symbol-Oblique D}{/Arial-Italic t}_{/=8 MD} = 0.0005", \
  "../data/nhc10d.0001/avb.dat" u ($1/N + hbin):($5*N) w l ls 3 t "{/Symbol-Oblique D}{/Arial-Italic t}_{/=8 MD} = 0.0001"

unset ylabel


unset multiplot
unset output

set terminal wxt
reset



