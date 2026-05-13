set terminal pngcairo size 1200,800 enhanced font 'Arial,14'
set output 'three_plots_numerical_data_M=100.png'

set title "Comparison: numerical, exact, and data values"
set xlabel "X values"
set ylabel "Values"
set grid

set yrange [-0.05:1.05]

set key top left box font 'Arial,12'

plot 'data_combined_numerical.txt' using 1:2 with lines lt 1 lc rgb 'red' lw 2 title 'Numerical solution', \
     'data_combined_numerical.txt' using 1:3 with lines lt 1 lc rgb 'green' lw 2 title 'Exact (DG flag)', \
     'data_combined_numerical.txt' using 1:4 with lines lt 1 lc rgb 'blue' lw 2 title 'Data values'
