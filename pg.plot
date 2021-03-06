set view equal xy
set view 70,303
set xlabel 'x'
set ylabel 'y'
set zlabel 't'
set zrange [0:5]
set hidden3d
set term png
set output 'pg00.png'
splot 'pg.log' matrix i 0 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '0 deg'
set output 'pg01.png'
splot 'pg.log' matrix i 1 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '7.2 deg'
set output 'pg02.png'
splot 'pg.log' matrix i 2 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '14.4 deg'
set output 'pg03.png'
splot 'pg.log' matrix i 3 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '21.6 deg'
set output 'pg04.png'
splot 'pg.log' matrix i 4 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '28.8 deg'
set output 'pg05.png'
splot 'pg.log' matrix i 5 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '36 deg'
set output 'pg06.png'
splot 'pg.log' matrix i 6 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '43.2 deg'
set output 'pg07.png'
splot 'pg.log' matrix i 7 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '50.4 deg'
set output 'pg08.png'
splot 'pg.log' matrix i 8 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '57.6 deg'
set output 'pg09.png'
splot 'pg.log' matrix i 9 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '64.8 deg'
set output 'pg10.png'
splot 'pg.log' matrix i 10 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '72 deg'
set output 'pg11.png'
splot 'pg.log' matrix i 11 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '79.2 deg'
set output 'pg12.png'
splot 'pg.log' matrix i 12 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '86.4 deg'
set output 'pg13.png'
splot 'pg.log' matrix i 13 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '93.6 deg'
set output 'pg14.png'
splot 'pg.log' matrix i 14 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '100.8 deg'
set output 'pg15.png'
splot 'pg.log' matrix i 15 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '108 deg'
set output 'pg16.png'
splot 'pg.log' matrix i 16 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '115.2 deg'
set output 'pg17.png'
splot 'pg.log' matrix i 17 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '122.4 deg'
set output 'pg18.png'
splot 'pg.log' matrix i 18 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '129.6 deg'
set output 'pg19.png'
splot 'pg.log' matrix i 19 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '136.8 deg'
set output 'pg20.png'
splot 'pg.log' matrix i 20 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '144 deg'
set output 'pg21.png'
splot 'pg.log' matrix i 21 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '151.2 deg'
set output 'pg22.png'
splot 'pg.log' matrix i 22 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '158.4 deg'
set output 'pg23.png'
splot 'pg.log' matrix i 23 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '165.6 deg'
set output 'pg24.png'
splot 'pg.log' matrix i 24 u (0+$1*0.0416667):(0+$2*0.0416667):3 w l lc 0 t '172.8 deg'
