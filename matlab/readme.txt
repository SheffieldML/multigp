The multigp toolbox is a toolbox for multiple output Gaussian
processes. Version 0.1 is the first official release, we have an older
release 0.001 which was used to publish the results in our NIPS paper,
this older release is available on request.

The toolbox allows the use of multioutput Gaussian processes in
combination with the "multi" kern type, first developed for use with
the gpsim code. The aim for this toolbox is that it should be more
general than the gpsim code. In particular it should allow for sparse
approximations for multiple output GPs. 

The multioutput GPs are constructed through convolution processes. For
more details see our NIPS paper and work by Dave Higdon and the NIPS
paper by Boyle and Frean.


Version 0.13
------------

Contains updates to the code for the technical report.


Version 0.11
------------

Updated version which allows for variational outputs and includes a financial data example.

Version 0.1
-----------

First version of the software with implementation of the 2008 NIPS paper.
