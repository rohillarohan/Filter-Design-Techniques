# Filter-Design-Techniques
This projects designs a filter with same specifications but with different techniques such as: 1) Windowing Method 2) Parks-McClellan Method 3) Minimum Phase filter

In the first part I have used the Kaiser windoe to design an FIR filter. To have a liner phase from this techique I have made sure that the filter response remains symmetric.

In the second part the same filter is designed using the celebrated Parks-McClellan Algorithm. Here also the special care has been taken to towards the symmetry of the response in order to get a linear phase.

The third part uses spectral factorization method by spliting a filter into maximum and minimum phase. Then the frequency response of the filter found from this technique is compared with the frequency response obtained from the Parks-McClellan method.
