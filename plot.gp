set title "Geozentric vs Geodetic latitude"

set xlabel "time after epoch [sek]"
set ylabel "angle [deg]"

set terminal png size 1920, 1080 
set output "c_latitude.png"

plot "build/res" using 1:2 with points pointtype 7 pointsize 1 title "Data Points"

set title "Geozentric vs Geodetic longitude"

set xlabel "time after epoch [sek]"
set ylabel "angle [deg]"

set output "c_longitude.png"

plot "build/res" using 1:3 with points pointtype 7 pointsize 1 title "Data Points"

set title "Geozentric vs Geodetic hight"

set xlabel "time after epoch [sek]"
set ylabel "hight [km]"

set output "c_hight.png"

plot "build/res" using 1:4 with points pointtype 7 pointsize 1 title "Data Points"
