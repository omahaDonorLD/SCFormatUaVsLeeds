# mandatory to read csv
set datafile separator ","


# to use decimal separator ','
set decimalsign ','
set format "%'.0f"
set decimal locale
set title "1,000m x 1,000m Area"


## Note : following graphs only needed from figure 5
## set legends outside graph
set key right bmargin


# figure 3 : coords + one pass only results
plot 'coords.csv' u 1:2 lc rgb 'red' title "targets", 'resclassiconepass.csv' lc rgb 'blue' pt 9 title "aircraft"


# figure 4 : coords + one pass and reductions
# see above, but in place of 'resclassiconepass.csv', use 'rl.csv'


# figure 5 : stages of building a connected component
## set legends outside graph
set key right bmargin
# a : coverage
plot 'coords.csv' using 1:2 lc rgb 'red' title "targets", 'finalres.csv' lc rgb 'green' pt 5 title "aircraft", 'uavs_grounds.csv' using 1:2 w lines title "in range aircraft"
# b : connected components
plot 'coords.csv' using 1:2 lc rgb 'red' title "targets", 'finalres.csv' lc rgb 'green' pt 5 title "aircraft", 'finaluavslinkG0.csv' using 1:2 w lines title "in range aircraft"
# c : unique connected component
plot 'coords.csv' using 1:2 lc rgb 'red' title "targets", 'out/G_0_coords' using 1:2 lc rgb 'green' pt 5 title "aircraft", 'out/G_0' w lines title "in range aircraft"
# d : extension
plot 'coords.csv' using 1:2 lc rgb 'red' title "targets",\
'G_0_coords' lc rgb 'green' pt 5 notitle, 'G_0' w lines notitle,\
'G_3_coords' lc rgb '#e55964' pt 5 notitle, 'G_0_to_G_3' lc rgb '#e55964' dt 2 w lines notitle, 'G_3' lc rgb '#e55964' w lines notitle,\
'G_4_coords' lc rgb '#00008b' pt 5 notitle, 'G_0_to_G_4' lc rgb '#00008b' dt 2 w lines notitle, 'G_4' lc rgb '#00008b' w lines notitle,\
'G_5_coords' lc rgb '#b8860b' pt 5 notitle, 'G_0_to_G_5' lc rgb '#b8860b' dt 2 w lines notitle, 'G_5' lc rgb '#b8860b' w lines notitle,\
NaN w points lc rgb 'black' pt 5 title "aircraft",\
NaN w lines lc rgb 'black' title "in range aircraft",\
NaN w lines lc rgb 'black' dt 2 title "removed connections"


# figure 6 : plot coords
plot 'coords.csv' u 1:2 lc rgb 'red' title "targets"


# figure 7 : skipped


# figure 8 : final results  (!! lb, redundancy value, is set to 2)
# a : 50 targets
plot '../../../data/50_grounds_coords.csv' using 1:2 lc rgb 'red' title "targets", '../../graphs/G_0' lc rgb 'cyan' w lines title "in range aircraft", 'mip_res_50_targets.csv' lc rgb 'green' pt 5 title "cover aircraft", 'ccs_50_targets.csv' lc rgb 'orange' pt 5 title "connectivity aircraft"
# b..e : same, just change data files


# figure 9 and 10


# figure 11 
# a : see file 'UavsTargetsCovering/srcs/scripts/pythonscripts/std_deviation_size_growth.py'
# b :
# c : 1500 targets
plot '../../data/1500_grounds_coords.csv' using 1:2 lc rgb 'red' title "targets"
# d : 1500 targets connectivity uavs + cover uavs
plot '../../data/1500_grounds_coords.csv' u 1:2 lc rgb 'red' title "targets", '../../out/runtime_growth/data/mip_res_1500_targets.csv' lc rgb 'blue' pt 9 title "cover aircraft", '../../out/runtime_growth/data/ccs_50_targets.csv' lc rgb 'green' pt 2 title "connectivity aircraft"
