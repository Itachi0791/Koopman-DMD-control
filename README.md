The file Direct_Encoding_Replication replicates part of thw work in the paper https://ieeexplore.ieee.org/abstract/document/10129822, and also expands on it.
The file EDMD_MPC uses EDMD based system Identification, then applies linear MPC. To use qpoases quadratic program solver in MATLAB, unzip the ./Resources/qpOASES-3.1.0.zip file, then run make.m  in ./Resources/qpOASES-3.1.0/interfaces/matlab to install qpoases for MATLAB.
The file Koopman_for_Euler finds a new approximate Koopman operator by using euler integration to simplfy analytical calculations in direct encoding.
