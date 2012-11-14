unset multiplot
reset

set encoding cp1250 # make minus sign longer
#set encoding iso_8859_1
set terminal postscript enhanced size 10, 4 lw 2.0 font "Arial, 14"
set output "ljvol.ps"

set multiplot

# number of particles
N = 108
# half of the bin size
hbin = 0.001


set lmargin 6.5
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

set style line 1 lc rgb "#ff0044" lt 1 
set style line 2 lc rgb "#22cc00" lt 2
set style line 3 lc rgb "#2200cc" lt 4


set label 1 "(a)" at screen dx, 1-dy
set label 2 "(b)" at screen .5+dx, 1-dy



set size 0.5, 1.0

set origin 0.0, 0.0

#set title ""
set xlabel "{/Arial-Italic n} = {/Arial-Italic N}{/=6 &.}/{/Arial-Italic V}" offset 0, 0.5
set ylabel "Distribution {/Symbol-Oblique r}{/=6 &.}({/Arial-Italic n}{/=6 &.})" offset 2.5, 0
set key left bottom Left reverse width -0 font "Arial, 12"

plot [.05:.75][0:2] \
  "../data/ljvolnhc1/avp.dat"     u ($1 + hbin):($5) w l ls 1 t "Nose-Hoover Chain, {/Arial-Italic W} = {/Arial-Italic Q}_1 = 300", \
  "../data/ljvollangf300/avp.dat" u ($1 + hbin):($5) w l ls 2 t "Langevin, {/Arial-Italic W} = 300, {/Symbol-Oblique z} = 1", \
  "../data/ljvolmcb/avp.dat"      u ($1 + hbin):($5) w l ls 3 t "Monte-Carlo, {/Symbol-Oblique D} = 0.05"

# "../data/ljvollangb/avp.dat" u ($1 + hbin):($5) w l ls 2 t "Langevin, {/Symbol-Oblique D}t = 5{/Symbol \264}10^{-5}"

set origin 0.5, 0.0

set ytics 0.1 
set ylabel 'Pressure ~{{/Arial-Italic p}{/=8 &.}}{0.3\^}{/=10 &.}({/Arial-Italic n}{/=6 &.})' offset 2.5, 0
set yrange [:0.3]
set key left top
set ytics nomirror
set my2tics 5
set y2tics 5 font tcfont offset -0.7, 0
set y2label '~{{/Arial-Italic G}{/=12 &.}}{0.5\^}{/=10 &.}({/Arial-Italic V}{/=8 &.}) = ~{{/Arial-Italic F}{/=12 &.}}{0.5\^}{/=10 &.}({/Arial-Italic V}{/=8 &.}) + {/Arial-Italic p^*}{/=6 &.}{/Arial-Italic V}' \
    offset -2., 0 font "Arial, 12"
set y2range[:15]

set rmargin 5.5

set label "{/Arial-Italic p}^* = 0.0696" at 0.42, 0.09 font "Arial, 12"
sh = 0.2
plot [.05:.7][] \
  0.0696 w l lt 4 not, \
  "../data/ljvolnhc1/avp.dat"       u ($1 + hbin):($2) w l ls 1 t "Nose-Hoover Chain", \
  "../data/ljvollangf300/avp.dat"   u ($1 + hbin):($2) w l ls 2 t "Langevin", \
  "../data/ljvolmcb/avp.dat"        u ($1 + hbin):($2) w l ls 3 t "Monte-Carlo", \
  "../data/ljvolnhc1/fe.dat"        u ($1):($2 + 0.069625*N/$1 - 139.5 + sh)  axes x1y2 w l ls 1 lw 3 not, \
  "../data/ljvollangf300/fe.dat"    u ($1):($2 + 0.069500*N/$1 - 139.25 + sh) axes x1y2 w l ls 2 lw 3 not, \
  "../data/ljvolmcb/fe.dat"         u ($1):($2 + 0.069666*N/$1 - 139.5 + sh)  axes x1y2 w l ls 3 lw 3 not

#  "../data/ljvollangb/avp.dat" u ($1 + hbin):($2) w l ls 2 t "Langevin",
#  "../data/ljvollangb/fe.dat"   u ($1):($2 + 0.069616*N/$1 - 139.45 + sh) axes x1y2 w l ls 2 not,

unset multiplot
unset output

set terminal wxt
reset



