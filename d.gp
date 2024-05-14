set title "Sonate position over time"

set ylabel "latitude"
set xlabel "longitude"

set terminal png size 1920, 1080 
set output "d.png"
set grid
datafile = "build/d"
start_longitude = -126.624482 
end_longitude = 163.282418

set arrow 1 from start_longitude, graph 0 to start_longitude, graph 1 nohead lc rgb "red" lw 2
set arrow 2 from end_longitude, graph 0 to end_longitude, graph 1 nohead lc rgb "blue" lw 2

set arrow 3 from graph 0, first 0 to graph 1, first 0 nohead lc rgb "green" lw 2

plot "build/d" using 1:2 with points pointtype 7 pointsize 1 title "Data Points"
