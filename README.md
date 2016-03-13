# Golang n-vector library

This performs approximate and exact geodesic calculations using the concept of
*n-vectors*, as described in

    Gade, Kenneth, "A Non-Singular Horizontal Position Representation",
    The Journal of Navigation (2010), 63, 395-417.
    doi:10.1017/S0373463309990415

See the above or
[http://www.navlab.net/nvector/](http://www.navlab.net/nvector/) for
documentation on the methods used.

Note that in contrast to the above, this library uses a right-handed orthogonal
*XYZ* earth reference frame with *Z* aligned along the rotation axis, *X*
passing through the Prime Meridian, and *Y* passing though 90⁰E longitude.
