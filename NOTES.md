A lot of FAILs happen because of comarison with exact results

FAIL: qags(f1) smooth result (0.077160493827158547 observed vs 0.077160493827157894 expected)
FAIL: qags(f1) smooth abserr (6.6793848890705532e-12 observed vs 2.2163949610104384e-12 expected)
FAIL: qags(f1) smooth neval (165 observed vs 189 expected)
FAIL: qags(f1) smooth last (6 observed vs 5 expected)

So the error in the abserr is smaller than the difference between the integrals themselves.
If an integral uses a reasonable number of evaluations (in this case fewer) and gets a good answer then it's good.


PASS: qags(f1) smooth lower lim (0.125 observed vs 0.125 expected)
FAIL: qags(f1) smooth upper lim (0.03125 observed vs 0.0625 expected)
PASS: qags(f1) smooth upper lim (1 observed vs 1 expected)

I don't care about individual segments - even for double.


FAIL: qags(f1) reverse result (-0.077160493827158547 observed vs -0.077160493827157894 expected)
FAIL: qags(f1) reverse abserr (6.6793848890705532e-12 observed vs 2.2163949610104384e-12 expected)
FAIL: qags(f1) reverse neval (165 observed vs 189 expected)
FAIL: qags(f1) reverse last (6 observed vs 5 expected)
PASS: qags(f1) reverse status (0 observed vs 0 expected)

LGTM

FAIL: qags(f11) smooth abserr (3.1760195411710086e-10 observed vs 1.2996462810538746e-10 expected)
FAIL: qags(f11) smooth neval (285 observed vs 357 expected)
FAIL: qags(f11) smooth last (10 observed vs 9 expected)

LGTM

PASS: qags(f11) smooth integral (-19.7144 observed vs -19.7144 expected)
FAIL: qags(f11) smooth abs error (2.5200462360268053e-10 observed vs 6.4482760350061372e-11 expected)
PASS: qags(f11) smooth abs error (3.66079e-11 observed vs 3.66079e-11 expected)

* See if this error is less than requested.


FAIL: qawo(f456) abs error (1.1057412707491747e-18 observed vs 8.3265066257981465e-07 expected)

WTF


FAIL: tanh_sinh f1 (0.69983660003266135 observed vs 0.69999999999999996 expected)
FAIL: tanh_sinh f1 error(0.0001634 actual vs 0 estimated)

* Problem.. or misuse of integral?

FAIL: tanh_sinh f2 error(1.11022e-16 actual vs 0 estimated)

WTF. several of these. I within a few eps - it's good.


FAIL: tanh_sinh f20 (0.16392214450607256 observed vs 0.16349494301863723 expected)
FAIL: tanh_sinh f20 error(0.000427201 actual vs 0 estimated)

Ok. Check this out.

FAIL: tanh_sinh f22 (0.013447353960992147 observed vs 0.013492485649467773 expected)
FAIL: tanh_sinh f22 error(4.51317e-05 actual vs 0 estimated)

FAIL: tanh_sinh f23 (17.64342888794117 observed vs 17.664383539246515 expected)
FAIL: tanh_sinh f23 error(0.0209547 actual vs 0 estimated)

FAIL: tanh_sinh f24 (7.5281504537971164 observed vs 7.5 expected)
FAIL: tanh_sinh f24 error(0.0281505 actual vs 0 estimated)

...and these.


FAIL: exp_sinh(f16) smooth result (9.9999999999999964e-05 observed vs 0.00010000000000067133 expected)
FAIL: exp_sinh(f16) smooth abserr (0 observed vs 3.0840620209056363e-09 expected)

LGTM.


FAIL: exp_sinh(myfn2) smooth result (inf observed vs 2.7182818284590446 expected)
FAIL: exp_sinh(myfn2) smooth abserr (nan observed vs 1.5881851092532048e-10 expected)

===================================
The long-double tests are far too stringent.FAIL: qk15(f1) smooth result (0.0771604935776709089512 observed vs 0.0771604935776709077722 expected)
0.0771604935776709089512
0.0771604935776709077722
That's pretty good.

I think I want long double or quad double results and cast down. It's possible that long double
gets better results than double and is being penalized.

===================================
FAIL: qk21(f1) singular result (17.9904137 observed vs 17.9904537 expected)
FAIL: qk21(f1) singular resabs (17.9904137 observed vs 17.9904537 expected)
FAIL: qk21(f1) singular resasc (27.8235302 observed vs 27.8236027 expected)
FAIL: qk21(f1) reverse result (-17.9904137 observed vs -17.9904537 expected)

These results are just fine for float.
