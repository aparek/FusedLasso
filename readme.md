I. Introduction
     The MATLAB code contained in this directory implements the algorithm
     derived in the paper [1] using Majorization-Minimization technique. 
     This file describes the details and procedures for the examples
     contained in [1]. 
     
II. Files contained and their description
     This directory contains the following files-
     
     1. demo.m
         This file demonstrates the CNC FLSA [1] for denoising a 
         sparse piecewise constant signal. This file reproduces example 1 
         in the paper [1]. 
         
     3. ECG_demo.m
         This file demonstrates the denoising of a synthetic ECG signal, 
         generated using ECGSYN (see [1] for details). 
         The synthetic ECG signal is obtained by the following commands
         
         fs = 256;
         ecg = ecgsyn(fs, 20); % 20 is the number of beats to be simulated.                                                   
      
         
     5. CNC_FLSA.m
         This function minimizes the CNC FLSA objective function
         F(x) = 0.5||y-x||_2^2 + lam0*phi(x,a0) + lam1*phi(Dx,a1)
         using the majorization-minimization technique, where phi
         is a non-convex penalty function. type `help CNC_FLSA' for more details
          
     6. soft.m
         This file implements the soft thresholding rule. 
         type `help soft.m' for more details
     
     7. tvd.c, tvd.mexmaci, tvd.mexmaci64, and tvd.mex64
         C++ implementation of TV denoising (see [1] for details)
 
For questions/comments contact: Ankit Parekh (ankit.parekh@nyu.edu)

Please cite as: 
[1] Convex fused lasso denoising with non-convex regularization and its 
     use for pulse detection. 
     Ankit Parekh and Ivan W. Selesnick, IEEE SPMB, 2015. 
