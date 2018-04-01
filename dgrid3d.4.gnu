set term postscript eps enhanced color
set output 'dgrid3d.4.eps'
set view 55, 290, 1, 1
unset key
set xlabel "Distance to magnet, m" offset 0,-3
set ylabel "Frequancy, Hz"
set zlabel "Power, W" offset 1,6.5
set pm3d at b
DEBUG_TERM_HTIC = 119
DEBUG_TERM_VTIC = 119
splot "3ddata.dat" using 1:2:3 with lines
