# TER-Pyrolyse

# Gnuplot : 

set xlabel "Temps (s)"
set xtics nomirror
set xrange [0:200] #200 ici est le temps final d'étude

set ylabel "Masse volumique (kg.m-3)"
set yrange [0:1000]

set x2label "Temperature (K)"
set x2tics nomirror
set x2range [300:1000] #300 et 1000 correspondent ici à la température initiale et finale respectivement

set key outside right top

plot "densite_Euler.dat" u 6:1 w l title "{/Symbol r}_{b}", \
     "densite_Euler.dat" u 6:2 w l title "{/Symbol r}_{c}", \
     "densite_Euler.dat" u 6:3 w l title "{/Symbol r}_{g}", \
     "densite_Euler.dat" u 6:4 w l title "{/Symbol r}_{l}", \
     "densite_Euler.dat" u 6:5 w l title "{/Symbol r}_{v}", \
     "densite_Euler.dat" u 6:($1+$2+$3+$4+$5) w l title "{/Symbol r}_{totale}"