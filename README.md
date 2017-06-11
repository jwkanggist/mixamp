# mixamp
README for demonstration of 2D sparse signal separation via MixAMP and TFOCS
by Jaewook Kang (jwkkang@gist.ac.kr) 
=============================================================

Description: 
=========================================================================
 This demonstration contains methods for performing 2D sparse signal separation 
 via the MixAMP iteration and the TFOCS method.
 The solvers included in this comparison are as given below:

"	MixAMP : J. Kang, H. Jung, K. Kim, "Fast Sginal Separation of 2D Sparse Mixture 
          via Approximate Message-Passing," submitted and available at ArXiv soon.
          
"	TFOCS: S. Becker, E. J. Candes and M. Gran, “Templates for Convex Cone
         Problems with Applications to Sparse Signal Recovery,” Math. Program.
         Comput., vol. 3, num. 3, Aug. 2011. MATLAB code available
         at http://cvxr.com
         
==============================================================================

This demonstration require two external numerical solvers

1) splitBregmanROF.c (MEX version) by Tom Goldstein
	-available at http://www.ece.rice.edu/~tag7/Tom_Goldstein/Split_Bregman.html 

2) The TFOCS package 
	-available at http://cvxr.com

	-Note that users must include "\MixAMP_Framework\solvers\TFOCS-1.3.1" in the MATLAB
	search path. 
	
===============================================================================

Player Information:
=====================
"	MATLAB Version: 8.3.0.532 (R2014a) + Image Processing Toolbox + 
"	Operating System: Microsoft Windows 7 Version 6.1 (Build 7601: Service Pack 1)
"	Java Version: Java 1.7.0_11-b21 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode

Packing List: 
  *Demonstration
  
 	demo_MixAMP_direct_and_FD.m  : providing a 2D sparsity separation with direct and finite-difference sparsity mixture via the MixAMP iteration.
 	demo_MixAMP_group_and_direct.m  : providing a 2D sparsity separation with direct and group sparsity mixture via the MixAMP iteration.
 	demo_TFOCS_direct_and_FD.m   : providing a 2D sparsity separation with direct and finite-difference sparsity mixture via the TFOCS method.
 	demo_TFOCS_group_and_direct.m: providing a 2D sparsity separation with direct and group sparsity mixture via the TFOCS method.
 	
 	*MixAMP iteration solvers
 	- solve_MixAMP_direct_and_FD.m
 	- solve_MixAMP_direct_and_group.m
 	
 	* The  two external numerical solvers included
  - TFOCS-1.3.1 package
  - splitBregmanROF_mex package
  
  * some sample images for demonstration
  - 128 by 128 cameraman image (tif)
  - 128 by 128 QR code image   (png)
  - 2D shot noise with 5% sparsity  (mat)
  - 2D shot noise with 10% sparsity (mat)
  
Contact and feedback Information:
=================================
"	Office phone: +82-62-715-2264
"	E-mail      :  jwkkang@gist.ac.kr, jwkang10@gmail.com
%===============================================================================
Copyright (c) 2015, Gwangju Institute of Science and Technology.
All rights reserved.

Contributors: Jaewook Kang, Hyoyoung Jung, and Kiseon Kim
This work is suppoted by  Electronic Warfare Research Center
at Gwangju Institute of Science and Technology (GIST),
originally funded by Defense Acquisition Program Administration.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are
met:

	1. Redistributions of source code must retain the above copyright
	notice, this list of conditions and the following disclaimer.

	2. Redistributions in binary form must reproduce the above copyright
	notice, this list of conditions and the following disclaimer in the
	documentation and/or other materials provided with the distribution.

	3. Neither the names of  Gwangju Institute of Science and Technology, 
	nor the names of its contributors may
	be used to endorse or promote products derived from this software
	without specific prior written permission.

===================================================================================




      

