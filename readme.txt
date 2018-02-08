The repository contains work done on indoor location based on UWB technology (Decawave DW1000 device). the scripts re-create a field test whose results are stored in the file col.txt
Trilateration: positioning based uniquely on the intersection (חיתוך) between 3 spheres
Kalman: 1) Kalman: pure linear Kalman filter. input is synthetic (not the col.txt file)
2) EKF: Extended Kalman filter. input is synthetic not the col.txt file)
3) Trilateration_EKF: main file. positioning based on EKF whose input is the recorded file