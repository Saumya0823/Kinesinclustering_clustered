reset
#set terminal epscairo enhanced font "Times-New-Roman"
#set output "image.eps"
set size ratio 0.55
#set terminal epscairo 


set key bottom right font ",40"
#set xlabel "Time (s)" offset 0,0 font ",30"
#set ylabel "IMD ({/Symbol m}m)" offset -0.5,0 font ",20"
set xtics out nomirror font ", 30" scale 1
set ytics out nomirror font ", 30" scale 1
set yrange [0:0.5]
set xrange [0.2:255]
set border lw 1.5
set label 1 "A" at 1,0.5 offset -6


#set key top left font ",20"
#set tics scale 3
#set xlabel "Time (s)" offset 0,0 font ",30"
#set ylabel "IMD ({/Symbol m}m)" offset -0.5,0 font ",20"
#set xtics out nomirror offset -6,0 font ", 40" scale 1							!offset -4,-1 font ", 50" scale 1
#set ytics out nomirror offset 1,0 font ", 30" scale 1
#set ytics out nomirror ("" 0, "" 0.1, "" 0.2, "" 0.3, "" 0.4, "" 0.5)
#set yrange [0:0.5]
#set format y ""
#set xrange [0.2:8]
#set border lw 1.5 
#set label 1 "A" at 1,0.5 offset -6

p "range_clusran_D001_N5.dat" u 1:($2-$3):($2+$3) title "" w filledcurves fc "#88808080" fs transparent,\
  "range_clusran_D001_N5.dat" u 1:2 title "D= 0.01/R" w l lt rgb "#77000000" lw 5,\
  "range_clusran_D01_N5.dat" u 1:($2-$3):($2+$3) title "" w filledcurves fc "#88ffd580" fs transparent,\
  "range_clusran_D01_N5.dat" u 1:2 title "D= 0.1/R" w l lt rgb "#77ff4500" lw 5,\
  "range_clusran_D1_N5.dat" u 1:($2-$3):($2+$3) title "" w filledcurves fc "#9990ee90" fs transparent,\
  "range_clusran_D1_N5.dat" u 1:2 title "D= 1/R" w l lt rgb "#77228b22" lw 5
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
