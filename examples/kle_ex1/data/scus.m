global tx
global ok
global ks

tx=load('tx.dat');
ok=load('ok.dat');
ks=load('ks.dat');

splini=cvt_sample(2,100,100,3,true,1011);
[ r1, seed, it_num, it_diff, energy ] = cvt ( 2, 2048, 100, 3, 3, 1000, 100, 10, 1011, splini);
