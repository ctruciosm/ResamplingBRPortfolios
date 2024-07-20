The package RobustVariance is a Matlab implementation of the robust
difference in variances testing based  on Ledoit & Wolf(2010).
Essentially, one tests whether the difference of two variances is 
statistically significant, taking account of the heavy tail characteristics
and autocorrelation of assets like hedge funds or mutual funds.

Before running these commands, you need to store them in your Matlab working 
directory or change  the working directory appropriately.


The main command is
				RobustVariance()
As an easy start, you can run
				RobustVariance(data)

to test the null hypothesis "H0: Difference of variances is zero", which 
will take some minutes to compute since bootstrap inference is default.
The matrix data MUST BE of size [Tx2], otherwise you get an error. 
data must contain level returns, i.e. not in logarithms already.
In this reduced form input RobustVariance(data), all other parameters of 
RobustVariance() are set to the default. Note that the computation of the 
optimal block size as in optimalblRobustVariance() is computed automatically 
as part of RobustVariance(data). The remaining time in seconds is printed in 
the console. On most machines, this should only take some minutes.
If you do not want to do bootstrap inference, type RobustVariance(data,0).

The full input options and help file can be found by typing

				help RobustVariance
in Matlab.


The second main command is
				optimalblRobustVariance()
Again, the minimal input is
				optimalblRobustVariance(data)

As said before, if no optimal block size is passed to RobustVariance(), the 
block size calibration is done automatically. Once you know the optimal 
block size for a data set, pass it to RobustVariance() so that computation 
time reduces to some seconds on most machines.

				help optimalblRobustVariance

gives details on possible inputs and outputs of the command.


All other routines are help routines.

In both RobustVariance() and optimalblRobustVariance(), one can specify the 
kernel used to compute the prewhitened kernel variance estimator. The 
two options are Gallant/Parzen or Quadratic spectral. For full compatability 
with the R implementation of RobustVariance, which is also on Michael Wolf's 
homepage, choose the Parzen/Gallant kernel. Otherwise, the Quadratic spectral 
kernel is preferable. Ledoit&Wolf report results based on the Gallant/Parzen 
kernel in their paper.

The routine andmon6cvm2.m uses parts of the gmmopt package from Mike Cliff 
version 1.1.

You are free to use this package as you see fit. Please give proper credit 
though.

Bugs or comments may be reported to dan.wunderli@econ.uzh.ch.

(c) Zurich, Institure for Empirical Research in Economics, U Zurich. Nov 2010
