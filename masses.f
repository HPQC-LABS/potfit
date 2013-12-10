c***********************************************************************
      SUBROUTINE MASSES(IAN,IMN,NAME,GELGS,GNS,MASS,ABUND)
c***********************************************************************
c** For isotope with (input) atomic number IAN and mass number IMN,
c  return (output):  (i) as the right-adjusted 2-character variable NAME
c  the alphabetic symbol for that element,  (ii) the ground state
c  electronic degeneracy GELGS, (iii) the nuclear spin degeneracy GNS,
c  (iv) the atomic mass MASS [amu], and  (v) the natural isotopic
c  abundance ABUND [in percent].   GELGS values based on atomic states
c  in Moore's "Atomic Energy Level" tables, the isotope masses are taken
c  from the 2003 mass table [Audi, Wapstra & Thibault, Nucl.Phys. A729,
c  337-676 (2003)] and other quantities from Tables 6.2 and 6.3 of 
c  "Quantities, Units and Symbols in Physical Chemistry", by Mills et 
c  al. (Blackwell, 2'nd Edition, Oxford, 1993).
c** If the input value of IMN does not equal one of the tabulated values
c  for atomic species IAN, return the abundance-averaged standard atomic
c  weight of that atom and set GNS=-1 and ABUND=-1.
c                          COPYRIGHT 2005-2010
c** By R.J. Le Roy (with assistance from G.T. Kraemer & J.Y. Seto).
c           Last modified  25 July 2010 (added proton, d, t)
c***********************************************************************
      REAL*8 zm(123,0:10),mass,ab(123,10),abund
      INTEGER i,ian,imn,gel(123),nmn(123),mn(123,10),ns2(123,10),
     1        gelgs, gns
      CHARACTER*2 NAME,AT(123)
c
      DATA  at(1),gel(1),nmn(1),(mn(1,i),i=1,6)/' H',2,6,1,2,3,4,5,6/
      DATA  (zm(1,i),i=0,6)/1.00794d0, 1.00782503207d0, 2.0141017778d0,
     1                 3.0160492777d0,1.00727646677d0,2.013553212724d0,
     2                 3.0155007134d0/
      DATA  (ns2(1,i),i=1,3)/1,2,1/
      DATA  (ab(1,i),i=1,3)/99.985d0,0.015d0,0.d0/
c
      DATA  at(2),gel(2),nmn(2),(mn(2,i),i=1,2)/'He',1,2,3,4/
      DATA  (zm(2,i),i=0,2)/4.002602d0, 3.0160293191d0, 4.00260325415d0/
      DATA  (ns2(2,i),i=1,2)/1,0/
      DATA  (ab(2,i),i=1,2)/0.000137d0,99.999863d0/
c
      DATA  at(3),gel(3),nmn(3),(mn(3,i),i=1,10)/'Li',2,2,3,4,5,6,7,8,9,
     1                                                10,11,12/
      DATA  (zm(3,i),i=0,10)/6.941d0,3.030775,4.02719,5.01254,
     1  6.015122795d0, 7.01600455d0, 8.02248736,  9.0267895, 10.035481,
     2  11.043798, 12.05378/
      DATA  (ns2(3,i),i=1,2)/2,3/
      DATA  (ab(3,i),i=1,2)/7.5d0,92.5d0/
c
      DATA  at(4),gel(4),nmn(4),(mn(4,i),i=1,1)/'Be',1,1,9/
      DATA  (zm(4,i),i=0,1)/9.012182d0, 9.0121822d0/
      DATA  (ns2(4,i),i=1,1)/3/
      DATA  (ab(4,i),i=1,1)/100.d0/
c
      DATA at(5),gel(5),nmn(5),(mn(5,i),i=1,2)/' B',2,2,10,11/
      DATA (zm(5,i),i=0,2)/10.811d0, 10.0129370d0, 11.0093054d0/
      DATA  (ns2(5,i),i=1,2)/6,3/
      DATA  (ab(5,i),i=1,2)/19.9d0,80.1d0/
c
      DATA at(6),gel(6),nmn(6),(mn(6,i),i=1,3)/' C',1,3,12,13,14/
      DATA (zm(6,i),i=0,3)/12.011d0, 12.d0, 13.0033548378d0, 
     1                      14.003241989d0/
      DATA  (ns2(6,i),i=1,3)/0,1,0/
      DATA  (ab(6,i),i=1,3)/98.90d0,1.10d0, 0.d0/
c
      DATA at(7),gel(7),nmn(7),(mn(7,i),i=1,2)/' N',4,2,14,15/
      DATA (zm(7,i),i=0,2)/14.00674d0, 14.0030740048d0, 15.0001088982d0/
      DATA (ns2(7,i),i=1,2)/2,1/
      DATA (ab(7,i),i=1,2)/99.634d0,0.366d0/
c
      DATA at(8),gel(8),nmn(8),(mn(8,i),i=1,3)/' O',5,3,16,17,18/
      DATA (zm(8,i),i=0,3)/15.9994d0, 15.99491461956d0, 16.99913170d0,
     1                      17.9991610d0/
      DATA (ns2(8,i),i=1,3)/0,5,0/
      DATA (ab(8,i),i=1,3)/99.762d0, 0.038d0, 0.200d0/
c
      DATA at(9),gel(9),nmn(9),(mn(9,i),i=1,1)/' F',4,1,19/
      DATA (zm(9,i),i=0,1)/18.9984032d0, 18.99840322d0/
      DATA (ns2(9,i),i=1,1)/1/
      DATA (ab(9,i),i=1,1)/100.d0/
c
      DATA at(10),gel(10),nmn(10),(mn(10,i),i=1,3)/'Ne',1,3,20,21,22/
      DATA (zm(10,i),i=0,3)/20.1797d0, 19.9924401754d0, 20.99384668d0,
     1                       21.991385114d0/
      DATA (ns2(10,i),i=1,3)/0,3,0/
      DATA (ab(10,i),i=1,3)/90.48d0, 0.27d0, 9.25d0/
c
      DATA at(11),gel(11),nmn(11),(mn(11,i),i=1,1)/'Na',2,1,23/
      DATA (zm(11,i),i=0,1)/22.989768d0, 22.9897692809d0/
      DATA (ns2(11,i),i=1,1)/3/
      DATA (ab(11,i),i=1,1)/100.d0/
c
      DATA at(12),gel(12),nmn(12),(mn(12,i),i=1,3)/'Mg',1,3,24,25,26/
      DATA (zm(12,i),i=0,3)/24.3050d0, 23.985041700d0, 24.98583692d0,
     1                       25.982592929d0/
      DATA (ns2(12,i),i=1,3)/0,5,0/
      DATA (ab(12,i),i=1,3)/78.99d0, 10.00d0, 11.01d0/
c
      DATA at(13),gel(13),nmn(13),(mn(13,i),i=1,1)/'Al',2,1,27/
      DATA (zm(13,i),i=0,1)/26.981539d0, 26.98153863d0/
      DATA (ns2(13,i),i=1,1)/5/
      DATA (ab(13,i),i=1,1)/100.d0/
c
      DATA at(14),gel(14),nmn(14),(mn(14,i),i=1,3)/'Si',1,3,28,29,30/
      DATA (zm(14,i),i=0,3)/28.0855d0, 27.9769265325d0, 28.976494700d0,
     1                       29.97377017d0/
      DATA (ns2(14,i),i=1,3)/0,1,0/
      DATA (ab(14,i),i=1,3)/92.23d0, 4.67d0, 3.10d0/
 
      DATA at(15),gel(15),nmn(15),(mn(15,i),i=1,1)/' P',4,1,31/
      DATA (zm(15,i),i=0,1)/30.973762d0, 30.97376163d0/
      DATA (ns2(15,i),i=1,1)/1/
      DATA (ab(15,i),i=1,1)/100.d0/
c
      DATA at(16),gel(16),nmn(16),(mn(16,i),i=1,4)/' S',5,4,32,33,34,36/
      DATA (zm(16,i),i=0,4)/32.066d0, 31.97207100d0, 32.97145876d0,
     1                       33.96786690d0, 35.96708076d0/
      DATA (ns2(16,i),i=1,4)/0,3,0,0/
      DATA (ab(16,i),i=1,4)/95.02d0, 0.75d0, 4.21d0, 0.02d0/
c
      DATA at(17),gel(17),nmn(17),(mn(17,i),i=1,2)/'Cl',4,2,35,37/
      DATA (zm(17,i),i=0,2)/35.4527d0, 34.96885268d0, 36.96590259d0/
      DATA (ns2(17,i),i=1,2)/3,3/
      DATA (ab(17,i),i=1,2)/75.77d0, 24.23d0/
c
      DATA at(18),gel(18),nmn(18),(mn(18,i),i=1,3)/'Ar',1,3,36,38,40/
      DATA (zm(18,i),i=0,3)/39.948d0, 35.967545106d0, 37.9627324d0,
     1                       39.9623831225d0/
      DATA (ns2(18,i),i=1,3)/0,0,0/
      DATA (ab(18,i),i=1,3)/0.337d0, 0.063d0, 99.600d0/
c
      DATA at(19),gel(19),nmn(19),(mn(19,i),i=1,3)/' K',2,3,39,40,41/
      DATA (zm(19,i),i=0,3)/39.0983d0, 38.96370668d0, 39.96399848d0,
     1                       40.96182576d0/
      DATA (ns2(19,i),i=1,3)/3,8,3/
      DATA (ab(19,i),i=1,3)/93.2581d0, 0.0117d0, 6.7302d0/
 
      DATA at(20),gel(20),nmn(20),(mn(20,i),i=1,6)/'Ca',1,6,40,42,43,44,
     1                                              46,48/
      DATA (zm(20,i),i=0,6)/40.078d0, 39.96259098d0, 41.95861801d0,
     1         42.9587666d0, 43.9554818d0, 45.9536926d0, 47.952534d0/
      DATA (ns2(20,i),i=1,6)/0,0,7,0,0,0/
      DATA (ab(20,i),i=1,6)/96.941d0, 0.647d0, 0.135d0, 2.086d0,
     1                      0.004d0, 0.187d0/
c
      DATA at(21),gel(21),nmn(21),(mn(21,i),i=1,1)/'Sc',4,1,45/
      DATA (zm(21,i),i=0,1)/44.955910d0, 44.9559119d0/
      DATA (ns2(21,i),i=1,1)/7/
      DATA (ab(21,i),i=1,1)/100.d0/
c
      DATA at(22),gel(22),nmn(22),(mn(22,i),i=1,5)/'Ti',5,5,46,47,48,49,
     1                                              50/
      DATA (zm(22,i),i=0,5)/47.88d0, 45.9526316d0, 46.9517631d0,
     1         47.9479463d0, 48.9478700d0, 49.9447912d0/
      DATA (ns2(22,i),i=1,5)/0,5,0,7,0/
      DATA (ab(22,i),i=1,5)/8.0d0, 7.3d0, 73.8d0, 5.5d0, 5.4d0/
c
      DATA at(23),gel(23),nmn(23),(mn(23,i),i=1,2)/' V',4,2,50,51/
      DATA (zm(23,i),i=0,2)/50.9415d0, 49.9471585d0, 50.9439595d0/
      DATA (ns2(23,i),i=1,2)/12,7/
      DATA (ab(23,i),i=1,2)/0.250d0, 99.750d0/
c
      DATA at(24),gel(24),nmn(24),(mn(24,i),i=1,4)/'Cr',7,4,50,52,53,54/
      DATA (zm(24,i),i=0,4)/51.9961d0, 49.9460442d0, 51.9405075d0,
     1                       52.9406494d0, 53.9388804d0/
      DATA (ns2(24,i),i=1,4)/0,0,3,0/
      DATA (ab(24,i),i=1,4)/4.345d0, 83.789d0, 9.501d0, 2.365d0/
c
      DATA at(25),gel(25),nmn(25),(mn(25,i),i=1,1)/'Mn',6,1,55/
      DATA (zm(25,i),i=0,1)/54.93805d0, 54.9380451d0/
      DATA (ns2(25,i),i=1,1)/5/
      DATA (ab(25,i),i=1,1)/100.d0/
c
      DATA at(26),gel(26),nmn(26),(mn(26,i),i=1,4)/'Fe',9,4,54,56,57,58/
      DATA (zm(26,i),i=0,4)/55.847d0, 53.9396105d0, 55.9349375d0,
     1                       56.9353940d0, 57.9332756d0/
      DATA (ns2(26,i),i=1,4)/0,0,1,0/
      DATA (ab(26,i),i=1,4)/5.8d0, 91.72d0, 2.2d0, 0.28d0/
c
      DATA at(27),gel(27),nmn(27),(mn(27,i),i=1,1)/'Co',10,1,59/
      DATA (zm(27,i),i=0,1)/58.93320d0, 58.9331950d0/
      DATA (ns2(27,i),i=1,1)/7/
      DATA (ab(27,i),i=1,1)/100.d0/
c
      DATA at(28),gel(28),nmn(28),(mn(28,i),i=1,5)/'Ni',9,5,58,60,61,62,
     1                                              64/
      DATA (zm(28,i),i=0,5)/58.69d0, 57.9353429d0, 59.9307864d0,
     1         60.9310560d0, 61.9283451d0, 63.9279660d0/
      DATA (ns2(28,i),i=1,5)/0,0,3,0,0/
      DATA (ab(28,i),i=1,5)/68.077d0,26.223d0,1.140d0,3.634d0,0.926d0/
c
      DATA at(29),gel(29),nmn(29),(mn(29,i),i=1,2)/'Cu',2,2,63,65/
      DATA (zm(29,i),i=0,2)/63.546d0, 62.9295975d0,64.9277895d0/
      DATA (ns2(29,i),i=1,2)/3,3/
      DATA (ab(29,i),i=1,2)/69.17d0, 30.83d0/
c
      DATA at(30),gel(30),nmn(30),(mn(30,i),i=1,5)/'Zn',1,5,64,66,67,68,
     1                                              70/
      DATA (zm(30,i),i=0,5)/65.40d0, 63.9291422d0, 65.9260334d0,
     1         66.9271273d0, 67.9248442d0, 69.9253193d0/
      DATA (ns2(30,i),i=1,5)/0,0,5,0,0/
      DATA (ab(30,i),i=1,5)/48.6d0, 27.9d0, 4.1d0, 18.8d0, 0.6d0/
c
      DATA at(31),gel(31),nmn(31),(mn(31,i),i=1,2)/'Ga',2,2,69,71/
Coxon   DATA (zm(31,i),i=0,2)/69.723d0, 68.925581d0, 70.9247073d0/
      DATA (zm(31,i),i=0,2)/69.723d0, 68.9255736d0, 70.9247013d0/
      DATA (ns2(31,i),i=1,2)/3,3/
      DATA (ab(31,i),i=1,2)/60.108d0, 39.892d0/
c
      DATA at(32),gel(32),nmn(32),(mn(32,i),i=1,5)/'Ge',1,5,70,72,73,74,
     1                                              76/
      DATA (zm(32,i),i=0,5)/72.61d0, 69.9242474d0, 71.9220758d0,
     1         72.9234589d0, 73.9211778d0, 75.9214026d0/
      DATA (ns2(32,i),i=1,5)/0,0,9,0,0/
      DATA (ab(32,i),i=1,5)/21.23d0, 27.66d0, 7.73d0, 35.94d0, 7.44d0/
c
      DATA at(33),gel(33),nmn(33),(mn(33,i),i=1,1)/'As',4,1,75/
      DATA (zm(33,i),i=0,1)/74.92159d0, 74.9215965d0/
      DATA (ns2(33,i),i=1,1)/3/
      DATA (ab(33,i),i=1,1)/100.d0/
c
      DATA at(34),gel(34),nmn(34),(mn(34,i),i=1,6)/'Se',5,6,74,76,77,78,
     1                                              80,82/
      DATA (zm(34,i),i=0,6)/78.96d0, 73.9224764d0, 75.9192136d0,
     1         76.9199140d0, 77.9173091d0, 79.9165213d0, 81.9166994d0/
      DATA (ns2(34,i),i=1,6)/0,0,1,0,0,0/
      DATA (ab(34,i),i=1,6)/0.89d0, 9.36d0, 7.63d0, 23.78d0, 49.61d0,
     1                      8.73d0/
c
      DATA at(35),gel(35),nmn(35),(mn(35,i),i=1,2)/'Br',4,2,79,81/
      DATA (zm(35,i),i=0,2)/79.904d0, 78.9183371d0, 80.9162906d0/
      DATA (ns2(35,i),i=1,2)/3,3/
      DATA (ab(35,i),i=1,2)/50.69d0, 49.31d0/
c
      DATA at(36),gel(36),nmn(36),(mn(36,i),i=1,6)/'Kr',1,6,78,80,82,83,
     1                                              84,86/
      DATA (zm(36,i),i=0,6)/83.80d0, 77.9203648d0, 79.9163790d0,
     1          81.9134836d0, 82.914136d0, 83.911507d0, 85.91061073d0/
      DATA (ns2(36,i),i=1,6)/0,0,0,9,0,0/
      DATA (ab(36,i),i=1,6)/0.35d0, 2.25d0, 11.6d0, 11.5d0, 57.0d0,
     1                      17.3d0/
c
      DATA at(37),gel(37),nmn(37),(mn(37,i),i=1,2)/'Rb',2,2,85,87/
      DATA (zm(37,i),i=0,2)/85.4678d0, 84.911789738d0, 86.909180527d0/
      DATA (ns2(37,i),i=1,2)/5,3/
      DATA (ab(37,i),i=1,2)/72.165d0, 27.835d0/
c
      DATA at(38),gel(38),nmn(38),(mn(38,i),i=1,4)/'Sr',1,4,84,86,87,88/
      DATA (zm(38,i),i=0,4)/87.62d0, 83.913425d0, 85.9092602d0,
     1                      86.9088771d0, 87.9056121d0/
      DATA (ns2(38,i),i=1,4)/0,0,9,0/
      DATA (ab(38,i),i=1,4)/0.56d0, 9.86d0, 7.00d0, 82.58d0/
c
      DATA at(39),gel(39),nmn(39),(mn(39,i),i=1,1)/' Y',4,1,89/
      DATA (zm(39,i),i=0,1)/88.90585d0, 88.9058483d0/
      DATA (ns2(39,i),i=1,1)/1/
      DATA (ab(39,i),i=1,1)/100.d0/
c
      DATA at(40),gel(40),nmn(40),(mn(40,i),i=1,5)/'Zr',5,5,90,91,92,94,
     1                                              96/
      DATA (zm(40,i),i=0,5)/91.224d0, 89.9047044d0, 90.9056458d0,
     1                      91.9050408d0, 93.9063152d0, 95.9082734d0/
      DATA (ns2(40,i),i=1,5)/0,5,0,0,0/
      DATA (ab(40,i),i=1,5)/51.45d0, 11.22d0, 17.15d0, 17.38d0, 2.80d0/
c
      DATA at(41),gel(41),nmn(41),(mn(41,i),i=1,1)/'Nb',2,1,93/
      DATA (zm(41,i),i=0,1)/92.90638d0, 92.9063781d0/
      DATA (ns2(41,i),i=1,1)/9/
      DATA (ab(41,i),i=1,1)/100.d0/
c
      DATA at(42),gel(42),nmn(42),(mn(42,i),i=1,7)/'Mo',7,7,92,94,95,96,
     1                                              97,98,100/
      DATA (zm(42,i),i=0,7)/95.94d0, 91.906811d0, 93.9050883d0,
     1        94.9058421d0, 95.9046795d0, 96.9060215d0, 97.9054082d0,
     2        99.907477d0/
      DATA (ns2(42,i),i=1,7)/0,0,5,0,5,0,0/
      DATA (ab(42,i),i=1,7)/14.84d0, 9.25d0, 15.92d0, 16.68d0, 9.55d0,
     1                      24.13d0, 9.63d0/
c
      DATA at(43),gel(43),nmn(43),(mn(43,i),i=1,1)/'Tc',6,1,98/
      DATA (zm(43,i),i=0,1)/97.907215d0, 97.907216d0/
      DATA (ns2(43,i),i=1,1)/12/
      DATA (ab(43,i),i=1,1)/100.d0/
c
      DATA at(44),gel(44),nmn(44),(mn(44,i),i=1,7)/'Ru',11,7,96,98,99,
     1                                              100,101,102,104/
      DATA (zm(44,i),i=0,7)/101.07d0, 95.907598d0, 97.905287d0,
     1     98.9059393d0, 99.9042195d0, 100.9055821d0, 101.9043493d0,
     2     103.905433d0/
      DATA (ns2(44,i),i=1,7)/0,0,5,0,5,0,0/
      DATA (ab(44,i),i=1,7)/5.52d0, 1.88d0, 12.7d0, 12.6d0, 17.0d0,
     1                      31.6d0, 18.7d0/
c
      DATA at(45),gel(45),nmn(45),(mn(45,i),i=1,1)/'Rh',10,1,103/
      DATA (zm(45,i),i=0,1)/102.90550d0, 102.905504d0/
      DATA (ns2(45,i),i=1,1)/1/
      DATA (ab(45,i),i=1,1)/100.d0/
c
      DATA at(46),gel(46),nmn(46),(mn(46,i),i=1,6)/'Pd',1,6,102,104,105,
     1                                              106,108,110/
      DATA (zm(46,i),i=0,6)/106.42d0, 101.905609d0, 103.904036d0,
     1       104.905085d0, 105.903486d0, 107.903892d0, 109.905153d0/
      DATA (ns2(46,i),i=1,6)/0,0,5,0,0,0/
      DATA (ab(46,i),i=1,6)/1.02d0, 11.14d0, 22.33d0, 27.33d0, 26.46d0,
     1                      11.72d0/
c
      DATA at(47),gel(47),nmn(47),(mn(47,i),i=1,2)/'Ag',2,2,107,109/
      DATA (zm(47,i),i=0,2)/107.8682d0, 106.905097d0, 108.904752d0/
      DATA (ns2(47,i),i=1,2)/1,1/
      DATA (ab(47,i),i=1,2)/51.839d0, 48.161d0/
c
      DATA at(48),gel(48),nmn(48),(mn(48,i),i=1,8)/'Cd',1,8,106,108,110,
     1                                             111,112,113,114,116/ 
      DATA (zm(48,i),i=0,8)/112.411d0, 105.906459d0, 107.904184d0, 
     1       109.9030021d0, 110.9041781d0, 111.9027578d0, 112.9044017d0,
     2       113.9033585d0, 115.904756d0/
      DATA (ns2(48,i),i=1,8)/0,0,0,1,0,1,0,0/
      DATA (ab(48,i),i=1,8)/1.25d0, 0.89d0, 12.49d0, 12.80d0, 24.13d0,
     1                      12.22d0, 28.73d0, 7.49d0/
c
      DATA at(49),gel(49),nmn(49),(mn(49,i),i=1,2)/'In',2,2,113,115/
      DATA (zm(49,i),i=0,2)/114.818d0, 112.904058d0, 114.903878d0/
      DATA  (ns2(49,i),i=1,2)/9,9/
      DATA (ab(49,i),i=1,2)/4.3d0, 95.7d0/
c
      DATA at(50),gel(50),nmn(50),(mn(50,i),i=1,10)/'Sn',1,10,112,114,
     1                                 115,116,117,118,119,120,122,124/
      DATA (zm(50,i),i=0,10)/118.710d0, 111.904818d0, 113.902779d0,
     1     114.903342d0, 115.901741d0, 116.902952d0, 117.901603d0,
     2     118.903308d0, 119.9021947d0, 121.9034390d0, 123.9052739d0/
      DATA (ns2(50,i),i=1,10)/0,0,1,0,1,0,1,0,0,0/
      DATA (ab(50,i),i=1,10)/0.97d0, 0.65d0, 0.34d0, 14.53d0, 7.68d0,
     1                       24.23d0, 8.59d0, 32.59d0, 4.63d0, 5.79d0/
c
      DATA at(51),gel(51),nmn(51),(mn(51,i),i=1,2)/'Sb',4,2,121,123/
      DATA (zm(51,i),i=0,2)/121.757d0, 120.9038157d0, 122.9042140d0/
      DATA (ns2(51,i),i=1,2)/5,7/
      DATA (ab(51,i),i=1,2)/57.36d0, 42.64d0/
c
      DATA at(52),gel(52),nmn(52),(mn(52,i),i=1,8)/'Te',5,8,120,122,123,
     1                                             124,125,126,128,130/
      DATA (zm(52,i),i=0,8)/127.60d0, 119.904020d0, 121.9030439d0,
     1    122.9042700d0, 123.9028179d0, 124.9044307d0, 125.9033117d0,
     2    127.9044631d0, 129.9062244d0/
      DATA (ns2(52,i),i=1,8)/0,0,1,0,1,0,0,0/
      DATA (ab(52,i),i=1,8)/0.096d0, 2.603d0, 0.908d0, 4.816d0,
     1                      7.139d0, 18.95d0, 31.69d0, 33.80d0/
c
      DATA at(53),gel(53),nmn(53),(mn(53,i),i=1,2)/' I',4,2,127,129/
      DATA (zm(53,i),i=0,2)/126.90447d0, 126.904473d0, 128.904988d0/
      DATA (ns2(53,i),i=1,2)/5,7/
      DATA (ab(53,i),i=1,2)/100.d0,0.d0/
c
      DATA at(54),gel(54),nmn(54),(mn(54,i),i=1,9)/'Xe',1,9,124,126,128,
     1                                          129,130,131,132,134,136/
      DATA (zm(54,i),i=0,9)/131.29d0, 123.9058930d0, 125.904274d0,
     1    127.9035313d0, 128.9047794d0, 129.9035080d0, 130.9050824d0,
     2    131.9041535d0, 133.9053945d0, 135.907219d0/
      DATA (ns2(54,i),i=1,9)/0,0,0,1,0,3,0,0,0/
      DATA (ab(54,i),i=1,9)/0.10d0, 0.09d0, 1.91d0, 26.4d0, 4.1d0,
     1                      21.2d0, 26.9d0, 10.4d0, 8.9d0/
c
      DATA at(55),gel(55),nmn(55),(mn(55,i),i=1,1)/'Cs',2,1,133/
      DATA (zm(55,i),i=0,1)/132.90543d0, 132.905451933d0/
      DATA (ns2(55,i),i=1,1)/7/
      DATA (ab(55,i),i=1,1)/100.d0/
c
      DATA at(56),gel(56),nmn(56),(mn(56,i),i=1,7)/'Ba',1,7,130,132,134,
     1                                             135,136,137,138/
      DATA (zm(56,i),i=0,7)/137.327d0, 129.9063208d0, 131.9050613d0,
     1    133.9045084d0, 134.9056886d0, 135.9045759d0, 136.9058274d0,
     2    137.9052472d0/
      DATA (ns2(56,i),i=1,7)/0,0,0,3,0,3,0/
      DATA (ab(56,i),i=1,7)/0.106d0, 0.101d0, 2.417d0, 6.592d0, 
     1                      7.854d0, 11.23d0, 71.70d0/
c
      DATA at(57),gel(57),nmn(57),(mn(57,i),i=1,2)/'La',4,2,138,139/
      DATA (zm(57,i),i=0,2)/138.9055d0, 137.907112d0, 138.9063533d0/
      DATA (ns2(57,i),i=1,2)/10,7/ 
      DATA (ab(57,i),i=1,2)/0.0902d0, 99.9098d0/
c
      DATA at(58),gel(58),nmn(58),(mn(58,i),i=1,4)/'Ce',9,4,136,138,140,
     1                                             142/
      DATA (zm(58,i),i=0,4)/140.115d0, 135.907172d0, 137.905991d0,
     1    139.9054387d0, 141.909244d0/
      DATA (ns2(58,i),i=1,4)/0,0,0,0/
      DATA (ab(58,i),i=1,4)/0.19d0, 0.25d0, 88.48d0, 11.08d0/
c
      DATA at(59),gel(59),nmn(59),(mn(59,i),i=1,1)/'Pr',10,1,141/
      DATA (zm(59,i),i=0,1)/140.90765d0, 140.9076528d0/
      DATA (ns2(59,i),i=1,1)/5/
      DATA (ab(59,i),i=1,1)/100.d0/
c
      DATA at(60),gel(60),nmn(60),(mn(60,i),i=1,7)/'Nd',9,7,142,143,144,
     1                                             145,146,148,150/
      DATA (zm(60,i),i=0,7)/144.24d0, 141.9077233d0, 142.9098143d0,
     1    143.9100873d0, 144.9125736d0, 145.9131169d0, 147.916893d0,
     2    149.920891d0/
      DATA (ns2(60,i),i=1,7)/0,7,0,7,0,0,0/
      DATA (ab(60,i),i=1,7)/27.13d0, 12.18d0, 23.80d0, 8.30d0, 17.19d0,
     1                       5.76d0, 5.64d0/
c
      DATA at(61),gel(61),nmn(61),(mn(61,i),i=1,1)/'Pm',6,1,145/
      DATA (zm(61,i),i=0,1)/144.912743d0, 144.912749d0/
      DATA (ns2(61,i),i=1,1)/5/
      DATA (ab(61,i),i=1,1)/100.d0/
c
      DATA at(62),gel(62),nmn(62),(mn(62,i),i=1,7)/'Sm',1,7,144,147,148,
     1                                             149,150,152,154/
      DATA (zm(62,i),i=0,7)/150.36d0, 143.911999d0, 146.9148979d0,
     1    147.9148227d0, 148.9171847d0, 149.9172755d0, 151.9197324d0,
     2    153.9222093d0/
      DATA (ns2(62,i),i=1,7)/0,7,0,7,0,0,0/
      DATA (ab(62,i),i=1,7)/3.1d0, 15.0d0, 11.3d0, 13.8d0, 7.4d0,
     1                      26.7d0, 22.7d0/
c
      DATA at(63),gel(63),nmn(63),(mn(63,i),i=1,2)/'Eu',8,2,151,153/
      DATA (zm(63,i),i=0,2)/151.965d0, 150.9198502d0, 152.9212303d0/
      DATA (ns2(63,i),i=1,2)/5,5/
      DATA (ab(63,i),i=1,2)/47.8d0, 52.2d0/
c
      DATA at(64),gel(64),nmn(64),(mn(64,i),i=1,7)/'Gd',5,7,152,154,155,
     1                                              156,157,158,160/
      DATA (zm(64,i),i=0,7)/157.25d0, 151.9197910d0, 153.92086560d0,
     1    154.9226220d0, 155.9221227d0, 156.9239601d0, 157.9241039d0,
     2    159.9270541d0/
      DATA (ns2(64,i),i=1,7)/0,0,3,0,3,0,0/
      DATA (ab(64,i),i=1,7)/0.20d0, 2.18d0, 14.80d0, 20.47d0, 15.65d0,
     1                      24.84d0, 21.86d0/
c
      DATA at(65),gel(65),nmn(65),(mn(65,i),i=1,1)/'Tb',16,1,159/
      DATA (zm(65,i),i=0,1)/158.92534d0, 158.9253468d0/
      DATA (ns2(65,i),i=1,1)/3/
      DATA (ab(65,i),i=1,1)/100.d0/
c
      DATA at(66),gel(66),nmn(66),(mn(66,i),i=1,7)/'Dy',17,7,156,158,
     1                                           160,161,162,163,164/
      DATA (zm(66,i),i=0,7)/162.50d0, 155.924283d0, 157.924409d0,
     1    159.9251975d0, 160.9269334d0, 161.9267984d0, 162.9287312d0,
     2    163.9291748d0/
      DATA (ns2(66,i),i=1,7)/0,0,0,5,0,5,0/
      DATA (ab(66,i),i=1,7)/0.06d0, 0.10d0, 2.34d0, 18.9d0, 25.5d0,
     1                      24.9d0, 28.2d0/
c
      DATA at(67),gel(67),nmn(67),(mn(67,i),i=1,1)/'Ho',16,1,165/
      DATA (zm(67,i),i=0,1)/164.93032d0, 164.9303221d0/
      DATA (ns2(67,i),i=1,1)/7/
      DATA (ab(67,i),i=1,1)/100.d0/
     
      DATA at(68),gel(68),nmn(68),(mn(68,i),i=1,6)/'Er',13,6,162,164,
     1                                            166,167,168,170/
      DATA (zm(68,i),i=0,6)/167.26d0, 161.928778d0, 163.929200d0,
     1    165.9302931d0, 166.9320482d0, 167.9323702d0, 169.9354643d0/
      DATA (ns2(68,i),i=1,6)/0,0,0,7,0,0/
      DATA (ab(68,i),i=1,6)/0.14d0, 1.61d0, 33.6d0, 22.95d0, 26.8d0,
     1                      14.9d0/
c
      DATA at(69),gel(69),nmn(69),(mn(69,i),i=1,1)/'Tm',8,1,169/  
      DATA (zm(69,i),i=0,1)/168.93421d0, 168.9342133d0/
      DATA (ns2(69,i),i=1,1)/1/
      DATA (ab(69,i),i=1,1)/100.d0/
c
      DATA at(70),gel(70),nmn(70),(mn(70,i),i=1,7)/'Yb',1,7,168,170,171,
     1                                            172,173,174,176/
      DATA (zm(70,i),i=0,7)/173.04d0, 167.933897d0, 169.9347618d0,
     1    170.936323580d0, 171.9363815d0, 172.9382108d0, 173.9388621d0,
     2    175.9425717d0/
      DATA (ns2(70,i),i=1,7)/0,0,1,0,5,0,0/
      DATA (ab(70,i),i=1,7)/0.13d0, 3.05d0, 14.3d0, 21.9d0, 16.12d0,
     1                      31.8d0, 12.7d0/
c
      DATA at(71),gel(71),nmn(71),(mn(71,i),i=1,2)/'Lu',4,2,175,176/
      DATA (zm(71,i),i=0,2)/174.967d0, 174.9407718d0, 175.9426863d0/
      DATA (ns2(71,i),i=1,2)/7,14/
      DATA (ab(71,i),i=1,2)/97.41d0, 2.59d0/
c
      DATA at(72),gel(72),nmn(72),(mn(72,i),i=1,6)/'Hf',5,6,174,176,177,
     1                                             178,179,180/
      DATA (zm(72,i),i=0,6)/178.49d0, 173.940046d0, 175.9414086d0,
     1    176.9432207d0, 177.9436988d0, 178.9458161d0, 179.9465500d0/
      DATA (ns2(72,i),i=1,6)/0,0,7,0,9,0/
      DATA (ab(72,i),i=1,6)/0.162d0, 5.206d0, 18.606d0, 27.297d0,
     1                      13.629d0, 35.100d0/
c
      DATA at(73),gel(73),nmn(73),(mn(73,i),i=1,2)/'Ta',4,2,180,181/
      DATA (zm(73,i),i=0,2)/180.9479d0, 179.9474648d0, 180.9479958d0/
      DATA (ns2(73,i),i=1,2)/16,7/
      DATA (ab(73,i),i=1,2)/0.012d0, 99.988d0/
c
      DATA at(74),gel(74),nmn(74),(mn(74,i),i=1,5)/' W',1,5,180,182,183,
     1                                             184,186/
      DATA (zm(74,i),i=0,5)/183.84d0, 179.946704d0, 181.9482042d0,
     1    182.9502230d0, 183.9509312d0, 185.9543641d0/
      DATA (ns2(74,i),i=1,5)/0,0,1,0,0/
      DATA (ab(74,i),i=1,5)/0.13d0, 26.3d0, 14.3d0, 30.67d0, 28.6d0/
c
      DATA at(75),gel(75),nmn(75),(mn(75,i),i=1,2)/'Re',6,2,185,187/
      DATA (zm(75,i),i=0,2)/186.207d0, 184.9529550d0, 186.9557531d0/
      DATA (ns2(75,i),i=1,2)/5,5/
      DATA (ab(75,i),i=1,2)/37.40d0, 62.60d0/
c
      DATA at(76),gel(76),nmn(76),(mn(76,i),i=1,7)/'Os',9,7,184,186,187,
     1                                             188,189,190,192/
      DATA (zm(76,i),i=0,7)/190.23d0, 183.9524891d0, 185.9538382d0,
     1    186.9557505d0, 187.9558382d0, 188.9581475d0, 189.9584470d0,
     2    191.9614807d0/
      DATA (ns2(76,i),i=1,7)/0,0,1,0,3,0,0/
      DATA (ab(76,i),i=1,7)/0.02d0, 1.58d0, 1.6d0, 13.3d0, 16.1d0,
     1                      26.4d0, 41.0d0/
c
      DATA at(77),gel(77),nmn(77),(mn(77,i),i=1,2)/'Ir',10,2,191,193/
      DATA (zm(77,i),i=0,2)/192.22d0, 190.9605940d0, 192.9629264d0/
      DATA (ns2(77,i),i=1,2)/3,3/
      DATA (ab(77,i),i=1,2)/37.3d0, 62.7d0/
c
c
      DATA at(78),gel(78),nmn(78),(mn(78,i),i=1,6)/'Pt',7,6,190,192,194,
     1                                            195,196,198/
      DATA (zm(78,i),i=0,6)/195.08d0, 189.959932d0, 191.9610380d0,
     1    193.9626803d0, 194.9647911d0, 195.9649515d0, 197.967893d0/
      DATA (ns2(78,i),i=1,6)/0,0,0,1,0,0/
      DATA (ab(78,i),i=1,6)/0.01d0,0.79d0,32.9d0,33.8d0,25.3d0,7.2d0/
c
      DATA at(79),gel(79),nmn(79),(mn(79,i),i=1,1)/'Au',2,1,197/
      DATA (zm(79,i),i=0,1)/196.96654d0, 196.9665687d0/
      DATA (ns2(79,i),i=1,1)/3/
      DATA (ab(79,i),i=1,1)/100.d0/
c
      DATA at(80),gel(80),nmn(80),(mn(80,i),i=1,7)/'Hg',1,7,196,198,199,
     1                                            200,201,202,204/
      DATA (zm(80,i),i=0,7)/200.59d0, 195.965833d0, 197.9667690d0,
     1    198.9682799d0, 199.9683260d0, 200.9703023d0, 201.9706430d0,
     2    203.9734939d0/
      DATA (ns2(80,i),i=1,7)/0,0,1,0,3,0,0/
      DATA (ab(80,i),i=1,7)/0.15d0, 9.97d0, 16.87d0, 23.10d0, 13.18d0,
     1                      29.86d0, 6.87d0/
c
      DATA at(81),gel(81),nmn(81),(mn(81,i),i=1,2)/'Tl',2,2,203,205/
      DATA (zm(81,i),i=0,2)/204.3833d0, 202.9723442d0, 204.9744275d0/
      DATA (ns2(81,i),i=1,2)/1,1/
      DATA (ab(81,i),i=1,2)/29.524d0, 70.476d0/
c
      DATA at(82),gel(82),nmn(82),(mn(82,i),i=1,4)/'Pb',1,4,204,206,207,
     1                                             208/
      DATA (zm(82,i),i=0,4)/207.2d0, 203.9730436d0, 205.9744653d0,
     1    206.9758969d0, 207.9766521d0/
      DATA (ns2(82,i),i=1,4)/0,0,1,0/
      DATA (ab(82,i),i=1,4)/1.4d0, 24.1d0, 22.1d0, 52.4d0/
c
      DATA at(83),gel(83),nmn(83),(mn(83,i),i=1,1)/'Bi',4,1,209/
      DATA (zm(83,i),i=0,1)/208.98037d0, 208.9803987d0/
      DATA (ns2(83,i),i=1,1)/9/
      DATA (ab(83,i),i=1,1)/100.d0/
c
      DATA at(84),gel(84),nmn(84),(mn(84,i),i=1,1)/'Po',5,1,209/
      DATA (zm(84,i),i=0,1)/208.982404d0, 208.9824304d0/
      DATA (ns2(84,i),i=1,1)/1/
      DATA (ab(84,i),i=1,1)/100.d0/
c
      DATA at(85),gel(85),nmn(85),(mn(85,i),i=1,1)/'At',-1,1,210/
      DATA (zm(85,i),i=0,1)/209.987126d0, 209.987148d0/
      DATA (ns2(85,i),i=1,1)/10/
      DATA (ab(85,i),i=1,1)/100.d0/
c
      DATA at(86),gel(86),nmn(86),(mn(86,i),i=1,1)/'Rn',1,1,222/
      DATA (zm(86,i),i=0,1)/222.017571d0, 222.0175777d0/
      DATA (ns2(86,i),i=1,1)/0/
      DATA (ab(86,i),i=1,1)/100.d0/
c
      DATA at(87),gel(87),nmn(87),(mn(87,i),i=1,1)/'Fr',-1,1,223/
      DATA (zm(87,i),i=0,1)/223.019733d0, 223.0197359d0/
      DATA (ns2(87,i),i=1,1)/3/
      DATA (ab(87,i),i=1,1)/100.d0/
c
      DATA at(88),gel(88),nmn(88),(mn(88,i),i=1,1)/'Ra',1,1,226/
      DATA (zm(88,i),i=0,1)/226.025403d0, 226.0254098d0/
      DATA (ns2(88,i),i=1,1)/0/
      DATA (ab(88,i),i=1,1)/100.d0/
c
      DATA at(89),gel(89),nmn(89),(mn(89,i),i=1,1)/'Ac',4,1,227/
      DATA (zm(89,i),i=0,1)/227.027750d0, 227.0277521d0/
      DATA (ns2(89,i),i=1,1)/3/
      DATA (ab(89,i),i=1,1)/100.d0/
c
      DATA at(90),gel(90),nmn(90),(mn(90,i),i=1,1)/'Th',-1,1,232/
      DATA (zm(90,i),i=0,1)/232.038d0, 232.0380553d0/
      DATA (ns2(90,i),i=1,1)/0/
      DATA (ab(90,i),i=1,1)/100.d0/
c
      DATA at(91),gel(91),nmn(91),(mn(91,i),i=1,1)/'Pa',-1,1,231/
      DATA (zm(91,i),i=0,1)/231.03588d0, 231.0358840d0/
      DATA (ns2(91,i),i=1,1)/3/
      DATA (ab(91,i),i=1,1)/100.d0/
c
      DATA at(92),gel(92),nmn(92),(mn(92,i),i=1,4)/' U',-1,4,233,234,
     1                                             235,238/
      DATA (zm(92,i),i=0,4)/238.0289d0, 233.0396352d0, 234.0409521d0,
     1    235.0439299d0, 238.0507882d0/
      DATA (ns2(92,i),i=1,4)/5,0,7,0/
      DATA (ab(92,i),i=1,4)/0.d0, 0.0055d0, 0.7200d0, 99.2745d0/
c
      DATA at(93),gel(93),nmn(93),(mn(93,i),i=1,1)/'Np',-1,1,237/
      DATA (zm(93,i),i=0,1)/237.0481678d0, 237.0481734d0/
      DATA (ns2(93,i),i=1,1)/5/
      DATA (ab(93,i),i=1,1)/100.d0/
c
      DATA at(94),gel(94),nmn(94),(mn(94,i),i=1,1)/'Pu',-1,1,244/
      DATA (zm(94,i),i=0,1)/244.064199d0, 244.064204d0/
      DATA (ns2(94,i),i=1,1)/0/
      DATA (ab(94,i),i=1,1)/100.d0/
c
      DATA at(95),gel(95),nmn(95),(mn(95,i),i=1,1)/'Am',-1,1,243/
      DATA (zm(95,i),i=0,1)/243.061375d0, 243.0613811d0/
      DATA (ns2(95,i),i=1,1)/5/
      DATA (ab(95,i),i=1,1)/100.d0/
c
      DATA at(96),gel(96),nmn(96),(mn(96,i),i=1,1)/'Cm',-1,1,247/
      DATA (zm(96,i),i=0,1)/247.070347d0, 247.070354d0/
      DATA (ns2(96,i),i=1,1)/9/
      DATA (ab(96,i),i=1,1)/100.d0/
c
      DATA at(97),gel(97),nmn(97),(mn(97,i),i=1,1)/'Bk',-1,1,247/
      DATA (zm(97,i),i=0,1)/247.070300d0, 247.070307d0/
      DATA (ns2(97,i),i=1,1)/3/
      DATA (ab(97,i),i=1,1)/100.d0/
c
      DATA at(98),gel(98),nmn(98),(mn(98,i),i=1,1)/'Cf',-1,1,251/
      DATA (zm(98,i),i=0,1)/251.079580d0, 251.079587d0/
      DATA (ns2(98,i),i=1,1)/1/
      DATA (ab(98,i),i=1,1)/100.d0/
c
      DATA at(99),gel(99),nmn(99),(mn(99,i),i=1,1)/'Es',-1,1,252/
      DATA (zm(99,i),i=0,1)/252.082944d0, 252.082980d0/
      DATA (ns2(99,i),i=1,1)/10/
      DATA (ab(99,i),i=1,1)/100.d0/
c
      DATA at(100),gel(100),nmn(100),(mn(100,i),i=1,1)/'Fm',-1,1,257/
      DATA (zm(100,i),i=0,1)/257.095099d0, 257.095105d0/
      DATA (ns2(100,i),i=1,1)/9/
      DATA (ab(100,i),i=1,1)/100.d0/
c
      DATA at(101),gel(101),nmn(101),(mn(101,i),i=1,1)/'Md',-1,1,258/
      DATA (zm(101,i),i=0,1)/258.09857d0, 258.098431d0/
      DATA (ns2(101,i),i=1,1)/16/
      DATA (ab(101,i),i=1,1)/100.d0/
c
      DATA at(102),gel(102),nmn(102),(mn(102,i),i=1,1)/'No',-1,1,259/
      DATA (zm(102,i),i=0,1)/259.100931d0, 259.101030d0/
      DATA (ns2(102,i),i=1,1)/9/
      DATA (ab(102,i),i=1,1)/100.d0/
c
      DATA at(103),gel(103),nmn(103),(mn(103,i),i=1,1)/'Lr',-1,1,260/
      DATA (zm(103,i),i=0,1)/260.105320d0, 260.105500d0/
      DATA (ns2(103,i),i=1,1)/-1/
      DATA (ab(103,i),i=1,1)/100.d0/
c
      DATA at(104),gel(104),nmn(104),(mn(104,i),i=1,1)/'Rf',-1,1,261/
      DATA (zm(104,i),i=0,1)/261.10869d0, 261.108770d0/
      DATA (ns2(104,i),i=1,1)/-1/
      DATA (ab(104,i),i=1,1)/100.d0/
c
      DATA at(105),gel(105),nmn(105),(mn(105,i),i=1,1)/'Db',-1,1,262/
      DATA (zm(105,i),i=0,1)/262.11376d0, 262.114080d0/
      DATA (ns2(105,i),i=1,1)/-1/
      DATA (ab(105,i),i=1,1)/100.d0/
c
      DATA at(106),gel(106),nmn(106),(mn(106,i),i=1,1)/'Sg',-1,1,263/
      DATA (zm(106,i),i=0,1)/263.11822d0, 263.118320d0/
      DATA (ns2(106,i),i=1,1)/-1/
      DATA (ab(106,i),i=1,1)/100.d0/
c
      DATA at(107),gel(107),nmn(107),(mn(107,i),i=1,1)/'Bh',-1,1,262/
      DATA (zm(107,i),i=0,1)/262.12293d0, 262.122890d0/
      DATA (ns2(107,i),i=1,1)/-1/
      DATA (ab(107,i),i=1,1)/100.d0/
c
      DATA at(108),gel(108),nmn(108),(mn(108,i),i=1,1)/'Hs',-1,1,265/
      DATA (zm(108,i),i=0,1)/265.13016d0, 265.130090d0/
      DATA (ns2(108,i),i=1,1)/-1/
      DATA (ab(108,i),i=1,1)/100.d0/
c
      DATA at(109),gel(109),nmn(109),(mn(109,i),i=1,1)/'Mt',-1,1,266/
      DATA (zm(109,i),i=0,1)/266.13764d0, 266.137300d0/
      DATA (ns2(109,i),i=1,1)/-1/
      DATA (ab(109,i),i=1,1)/100.d0/
c
      IF((IAN.LE.0).OR.(IAN.GT.109)) THEN
          MASS= 0.d0
          NAME= 'XX'
          IMN= 0
          WRITE(6,601) IAN
          RETURN
        ELSE
          NAME= AT(IAN)
        ENDIF
      IF((IAN.EQ.1).AND.(IMN.NE.1)) THEN
c** Special case: insert common name for deuterium or tritium
          IF(IMN.EQ.2) NAME=' D'
          IF(IMN.EQ.3) NAME=' T'
          IF(IMN.EQ.4) NAME=' p'
          IF(IMN.EQ.5) NAME=' d'
          IF(IMN.EQ.6) NAME=' t'
          ENDIF
      GELGS= GEL(IAN)
      MASS= -1.d0
      GNS= -1
      ABUND = -1.d0
      DO  I= 1,NMN(IAN)
          if(i.gt.10)  write(6,606) ian,imn,nmn(ian)
          IF(IMN.EQ.MN(IAN,I)) THEN
              MASS= ZM(IAN,I)
              GNS= NS2(IAN,I)+1
              ABUND = AB(IAN,I)
              ENDIF
          ENDDO
      IF(MASS.LT.0.d0) THEN
          MASS= ZM(IAN,0)
          IF(IMN.NE.0) WRITE(6,602) AT(IAN),IMN
          IMN= 0
          ENDIF
      RETURN
  601 FORMAT(' *** MASSES Data base does not include Atomic Number=',i4)
  602 FORMAT(' *** MASSES Data base does not include ',A2,'(',i3,
     1 '), so use average atomic mass.')
  606  format(/' *** ERROR *** call MASSES for atom with  AN=',I4,
     1  '  MN=',I4,'n(MN)=',I4)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

