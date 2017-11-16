
[](Tp x)
-> Tp
{ return std::exp(x); }

[](Tp x)
-> Tp
{ return Tp(x >= Tp{3} / Tp{10}); }

[](Tp x)
-> Tp
{ return std::sqrt(x); }

[](Tp x)
-> Tp
{ return (Tp{23} / Tp{25}) * std::cosh(x) - std::cos(x); }

[](Tp x)
-> Tp
{
  const auto xx = x * x;
  const auto xxxx = xx * xx;
  return Tp{1} / (xxxx + xx + Tp{9} / Tp{10});
}

[](Tp x)
-> Tp
{ return std::sqrt(x * x * x); }

[](Tp x)
-> Tp
{ return Tp{1} / std::sqrt(x); }

[](Tp x)
-> Tp
{
  const auto xx = x * x;
  const auto xxxx = xx * xx;
  return Tp{1} / (Tp{1} + xxxx);
}

[pi](Tp x)
-> Tp
{ return Tp{2} / (Tp{2} + std::sin(Tp{10} * pi * x)); }

[](Tp x)
-> Tp
{ return Tp{1} / (Tp{1} + x); }

[](Tp x)
-> Tp
{ return Tp{1} / (Tp{1} + std::exp(x)); }

[](Tp x)
-> Tp
{ return x / (std::exp(x) - Tp{1}); }

[pi](Tp x)
-> Tp
{ return std::sin(Tp{100} * pi * x) / (pi * x); }

[pi](Tp x)
-> Tp
{ return std::sqrt(Tp{50}) * std::exp(-50 * pi * x * x); }

[](Tp x)
-> Tp
{ return Tp{25} * exp(-Tp{25} * x); }

[pi](Tp x)
-> Tp
{ return Tp{50} / pi * (Tp{2500} * x * x + Tp{1}); }

[pi](Tp x)
-> Tp
{
  const auto arg = Tp{50} * pi * x;
  return Tp{50} * std::pow(std::sin(arg) / arg, Tp{2});
}

[](Tp x)
-> Tp
{
  return std::cos(std::cos(x)
	+ Tp{3} * std::sin(x)
	+ Tp{2} * std::cos(Tp{2} * x)
	+ Tp{3} * std::sin(Tp{2} * x)
	+ Tp{3} * std::cos(Tp{3} * x));
}

[](Tp x)
-> Tp
{ return std::log(x); }

[](Tp x)
-> Tp
{ return 1 / (x * x + Tp{1.005L}); }

[](Tp x)
-> Tp
{
  return 1 / std::cosh(10 * (x - Tp{0.2L}) * 2)
       + 1 / std::cosh(100 * (x - Tp{0.4L}) * 4)
       + 1 / std::cosh(1000 * (x - Tp{0.6L}) * 8);
}

[pi](Tp x)
-> Tp
{
  const auto arg = Tp{2} * pi * x;
  return 2 * pi * arg * std::sin(10 * arg) * std::cos(arg);
}

[](Tp x)
-> Tp
{
  const auto thing = Tp{230} * x - Tp{30};
  return Tp{1} / (Tp{1} + thing * thing);
}

[](Tp x)
-> Tp
{ return std::floor(std::exp(x)); }

[](Tp x)
-> Tp
{
  return (x < 1) * (x + 1)
       + (1 <= x & x <= 3) * (3 - x)
       + (x > 3) * 2;
}

limits = {
    {0.0L, 1.0L},
    {0.0L, 1.0L},
    {0.0L, 1.0L},
    {-1.0L, 1.0L},
    {-1.0L, 1.0L},
    {0.0L, 1.0L},
    {0.0L, 1.0L},
    {0.0L, 1.0L},
    {0.0L, 1.0L},
    {0.0L, 1.0L},
    {0.0L, 1.0L},
    {0.0L, 1.0L},
    {0.0L, 1.0L},
    {0.0L, 10.0L},
    {0.0L, 10.0L},
    {0.0L, 10.0L},
    {0.0L, 1.0L},
    {0.0L, pi},
    {0.0L, 1.0L},
    {-1.0L, 1.0L},
    {0.0L, 1.0L},
    {0.0L, 1.0L},
    {0.0L, 1.0L},
    {0.0L, 3.0L},
    {0.0L, 5.0L} };

f_exact = {
    1.7182818284590452354L,
    0.7L,
    2.0L / 3.0L,
    0.4794282266888016674L,
    1.5822329637296729331L,
    0.4L,
    2.0L,
    0.86697298733991103757L,
    1.1547005383792515290L,
    0.69314718055994530942L,
    0.3798854930417224753L,
    0.77750463411224827640L,
    0.49898680869304550249L,
    0.5L,
    1.0L,
    0.13263071079267703209e+08L,
    0.49898680869304550249L,
    0.83867634269442961454L,
    -1.0L,
    1.5643964440690497731L,
    0.16349494301863722618L,
    -0.63466518254339257343L,
    0.013492485649467772692L,
    17.664383539246514971L,
    7.5L };
