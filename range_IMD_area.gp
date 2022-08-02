reset
#set terminal epscairo enhanced font "Times-New-Roman"
#set output "image.eps"
set size ratio 0.55
#set terminal epscairo 


set key top right font ",20"
set tics scale 3
#set xlabel "Time (s)" offset 0,0 font ",30"
#set ylabel "IMD ({/Symbol m}m)" offset -0.5,0 font ",20"
set xtics out nomirror offset -2,0 font ", 30" scale 1
set ytics out nomirror offset 0,-1 font ", 30" scale 1
#set xtics out nomirror ("" 0, "" 50, "" 100, "" 150, "" 200, "" 250)
#set ytics out nomirror ("" 0, "" 0.5, "" 0.2, "" 0.3, "" 0.4, "" 0.5)
set yrange [0:0.1]
#set format x ""
#set format y ""
set xrange [1:35]
set border lw 1.5 
set label 1 "A" at 1,0.5 offset -6

p "range_clusmot_N5_D001.dat" u 1:($2-$3):($2+$3) title "" w filledcurves fc "#88808080" fs transparent,\
  "range_clusmot_N5_D001.dat" u 1:2 title "" w l lt rgb "#77000000" lw 5,\
  "range_clusmot_N5_D01.dat" u 1:($2-$3):($2+$3) title "" w filledcurves fc "#88ffd580" fs transparent,\
  "range_clusmot_N5_D01.dat" u 1:2 title "" w l lt rgb "#77ff4500" lw 5,\
  "range_clusmot_N5_D1.dat" u 1:($2-$3):($2+$3) title "" w filledcurves fc "#9990ee90" fs transparent,\
  "range_clusmot_N5_D1.dat" u 1:2 title "" w l lt rgb "#77228b22" lw 5
#  1.06*exp(x/353) w l dt 2 lw 2 lt rgb "black"

## #b5651d  = light-brown
## #964B00 = brown
## #ffd580 = light-orange
## #ff4500 = orange
## #e0ffff = light-cyan
## #00ffff = cyan

# #808080  = grey
# #000000 = black
# #ffd580 = light-orange
# #ff4500 = orange
# #90ee90 = light-green
# #228b22 = green
