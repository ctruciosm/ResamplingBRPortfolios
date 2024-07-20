These routines implement the HAC inference and the studentized
circular block bootstrap of the paper "Robust Performance Hypothesis
Testing with the Sharpe Ratio".

You can load all the routines (and the two datasets of Section 5) 
by using the following command in R:

> load("Sharpe.RData")

(Of course, you have to first put the file Sharpe.RData in your
current R working directory.)

There are basically three routines for `immediate' use, all the rest
are `help' routines:

hac.inference() computes HAC standard errors and corresponding
two-sided p-values. This routine uses the (prewhitened) Parzen
kernels instead of the (prewhitened) Quadratic Spectral kernel.
The reason is that R is somewhat slow and the former kernel has
a bounded support. Note that the two kernels yielded virtually
identical results in our simulation study, so it really should
not matter much at all. Of course, with a bit of additional work
you can modify the routines to use the (prewhitened) QS kernels 
instead.

boot.time.inference() computes a two-sided p-value based on the
studentized circular block bootstrap. Note that you have to
specify the block size b.

block.size.calibrate() implements Algorithm 3.1 for the data-dependent
choice of block size. Unfortunately, this can take quite a while 
to run in R. The problem is that R does not handle for and while
loops very efficiently. To obtain the results of the paper, I coded
up everything in C++, but these routines are not user-friendly
enough to be shared with the rest of the world.

Note: in all routines, ret stands for the T x 2 matrix of excess
returns. The routines will NOT work if you feed them a 2 x T matrix.

The dataset ret.agg contains the mutual funds data, while the dataset
ret.hedge contains the hedge funds data. Since the R routines use
the (prewhitened) Parzen kernels instead of the (prewhitened) QS
kernels, the results are slightly different from those of Table 3
if the paper. The p-values for the mutual funds data are 6.6 and 7.0,
respectively, and the p-values for the hedge funds data are 16.0 
and 24.9, respectively.

Unfortunately, I did not have the time to document or comment the
various routines in great detail. However, I have tried to make
the notation consistent with that of the paper and for the
HAC routines consistent with that of Andrews (1991) and Andrews and
Monahan (1992), respectively, as much as possible.

Last but not least, there are two separate .R functions that are 
designed to speed up the bootstrap inference by exploiting parallel
computations (which may not work on all platforms). They have been 
written and kindly made available by Mike O'Connor; see Notes.txt.
In case you have questions or comments regarding these extensions,
please contact him directly instead of me.


Zurich, Juni 2021,
Michael Wolf


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code in this distribution is released under the BSD 2-clause license.

% Copyright (c) 2014, Michael Wolf 
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
% 
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
