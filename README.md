# quadrature

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/5e3495d8c4004bc5a7ec7e25c65a98f8)](https://app.codacy.com/app/emsr/quadrature?utm_source=github.com&utm_medium=referral&utm_content=emsr/quadrature&utm_campaign=Badge_Grade_Dashboard)

This is a C++ quadrature library reengineered from GSL and with new things added.  This was hived off from tr29124_test.

This library began as a way to test the orthogonal polynomials that were added to the Gnu implementation of the C++ standard library as part of the implementation of TR1 Special Math Functions.  This library was dropped by the original author some years back.  Then and now, the foundation of this library is the Gnu Scientific Library (GSL) quadrature package which itself is based on the venerable QUADPACK library.

This library is not a simple translation - *that* would not add real value.

The goals and ideal of this implementation are:
- Type genericity - allow different number systems to be used with the library including multi-precision ones.  Currently, float, double, long double, and __float128 have been tested.
- Use C++ containers and algorithms throughout. The adaptive algorithms in GSL use what is esentially a priority queue to decide which chunk to work on next.  C++ obviously has a lot to offer here.
- Add new quadrature algorithms.  The GSL was focused on what I think is the middle of the spectrum: Gauss-Kronrod, Fejer and Clenshaw-Curtis, etc. This library seeks to add simple recursive midpoint, trapezoid, and Simpson rules at the "low-end" and double-exponential, sinh-tanh rules at the "high-end".
- Support contour integration in the complex plane.

Currently, this is a one-dimensional library.  Cubature is a future goal.  Monte-Carlo is another subject left out for now.  The standard C++ <random> library was developed in large part to support this use-case and it is a worthy thing to have in a C++ library.  After we get these done.
