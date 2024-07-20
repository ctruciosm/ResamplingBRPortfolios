
These routines implement the HAC inference and the studentized
circular block bootstrap of the paper "Robust Performance Hypothesis
Testing with the Variance".

You can load all the routines by using the following command in R:

> load("Var.RData")

(Of course, you have to first put the file Var.RData in your
current R working directory.)

There are basically three routines for `immediate' use, all the rest
are `help' routines:

hac.inference.log.var() computes HAC standard errors and corresponding
two-sided p-values. This routine uses the (prewhitened) Parzen
kernels instead of the (prewhitened) Quadratic Spectral kernel.
The reason is that R is somewhat slow and the former kernel has
a bounded support. Note that the two kernels yielded virtually
identical results in our simulation study, so it really should
not matter much at all. Of course, with a bit of additional work
you can modify the routines to use the (prewhitened) QS kernels 
instead.

boot.time.inference.log.var() computes a two-sided p-value based on 
the studentized circular block bootstrap. Note that you have to
specify the block size b.

block.size.calibrate.log.var() implements Algorithm 3.1 for the data-dependent
choice of block size. Unfortunately, this can take quite a while 
to run in R. The problem is that R does not handle for and while
loops very efficiently. To obtain the results of the paper, I coded
up everything in C++, but these routines are not user-friendly
enough to be shared with the rest of the world.

Note: in all routines, ret stands for the T x 2 matrix of excess
returns. The routines will NOT work if you feed them a 2 x T matrix.

Unfortunately, I did not have the time to document or comment the
various routines in great detail. However, I have tried to make
the notation consistent with that of the paper and for the
HAC routines consistent with that of Andrews (1991) and Andrews and
Monahan (1992), respectively, as much as possible.


Zurich, September 2014,
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

