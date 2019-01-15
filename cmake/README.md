Modern lapacke provides a lapacke-config.cmake (minimum version unconfirmed, but 3.7.1 does)
so FindLAPACKE.cmake wouldn't be neccesary, (see #38)
however this isn't reliable for older/other lapacke versions, so we provide a FindModule from: https://github.com/mrirecon/bart/blob/master/cmake/FindLAPACKE.cmake.

Download with `curl -O https://raw.githubusercontent.com/mrirecon/bart/master/cmake/FindLAPACKE.cmake`

This FindLapacke handles properly the case a lapacke-config.cmake exists in the system.
