
gauss_legendre()
{
  mu_0 = _Tp{2};

  for (int i = 1; i <= m; ++i)
    {
      diag[i - 1] = _Tp{0};
      subd[i - 1] = _Tp(i) / std::sqrt(_Tp(4 * i * i - 1));
    }
}

/*
 * alpha = 1/2
 * beta = 1/2
 */
gauss_chebyshev_t()
{
  mu_0 = pi;

  for (int i = 0; i < m; ++i)
    diag[i] = _Tp{0};

  subd[0] = std::sqrt(_Tp{0.5L});
  for (int i = 1; i < m; ++i)
    subd[i] = _Tp{0.5L};
}

/*
 * alpha = -1/2
 * beta = -1/2
 */
gauss_chebyshev_u()
{
  mu_0 = pi / _Tp{2};

  for (int i = 0; i < m; ++i)
    {
      diag[i] = _Tp{0};
      subd[i] = _Tp{0.5L};
    }
}

gauss_gegenbauer()
{
  ab = _Tp{2} * alpha;
  mu_0 = std::pow(_Tp{2}, ab + _Tp{1})
	* std::pow(std::tgamma(alpha + _Tp{1}), _Tp{2})
	/ std::tgamma(ab + _Tp{2});

  for (int i = 0; i < m; ++i)
    diag[i] = _Tp{0};

  subd[0] = std::sqrt(_Tp{1} / (_Tp{2} * alpha + _Tp{3}));
  for (int i = 2; i <= m; ++i)
    subd[i-1] = std::sqrt(_Tp(i) * (_Tp(i) + ab)
			/ (_Tp{4} * std::pow(_Tp(i) + alpha, _Tp{2}) - _Tp{1}));
}

gauss_jacobi()
{
  ab = alpha + beta;
  abp2i = ab + _Tp{2};
  mu_0 = std::pow(_Tp{2}, ab + _Tp{1}) * std::tgamma(alpha + _Tp{1}) 
       * std::tgamma(beta + _Tp{1}) / std::tgamma(abp2i);

  diag[0] = (beta - alpha) / abp2i;
  subd[0] = _Tp{2} * std::sqrt((alpha + _Tp{1}) * (beta + _Tp{1}) 
			       / ((abp2i + _Tp{1})) / abp2i;
  a2mb2 = (beta - alpha) * (beta + alpha);
  for (int i = 1; i <= m; ++i)
    {
      abp2i += _Tp{2};
      const auto abp2ip2 = abp2i + _Tp{2};
      diag[i] = a2mb2 / abp2i / abp2ip2;
      const auto ip1 = _Tp(i + 1);
      subd[i] = std::sqrt(_Tp(4 * i) * (alpha + ip1)
			      * (beta + ip1) * (ab + ip1)
			      / (abp2ip2 * abp2ip2 - _Tp{1})) / abp2ip2;
    }
}

gauss_laguerre()
{
  mu_0 = std::tgamma(alpha + _Tp{1});

  for (int i = 1; i <= m; ++i)
    {
      diag[i - 1] = _Tp(2 * i - 1) + alpha;
      subd[i - 1] = std::sqrt(_Tp(i) * (_Tp(i) + alpha));
    }
}

gauss_hermite()
{
  mu_0 = std::tgamma((alpha + 1) / _Tp{2});

  for (int i = 0; i < m; ++i)
    diag[i] = _Tp{0};

  for (int i = 1; i <= m; ++i)
    subd[i - 1] = std::sqrt((_Tp(i) + alpha * _Tp(i % 2)) / _Tp{2});
}

/*
 * alpha > -1
 */
gauss_exponential()
{
  mu_0 = _Tp{2} / (alpha + _Tp{1});

  auto ap2i = alpha;
  for (int i = 1; i <= m; ++i)
    {
      diag[i - 1] = _Tp{0};
      ap2i += _Tp{2};
      subd[i - 1] = (_Tp(i) + _Tp(i % 2) * alpha)
		  / std::sqrt((ap2i * ap2i - _Tp{1}));
    }
}

gauss_rational()
{
  ab = alpha + beta;
  mu_0 = std::tgamma(alpha + _Tp{1}) * std::tgamma(-(ab + _Tp{1})) 
       / std::tgamma(-beta);

  ap1 = alpha + _Tp{1};
  aba = ab * ap1;
  diag[0] = -ap1 / (ab + _Tp{2});
  subd[0] = -diag[0] * (beta + _Tp{1}) / (ab + _Tp{2}) / (ab + _Tp{3});
  for (int i = 2; i <= m; ++i)
    {
      abp2i = ab + _Tp(2 * i);
      diag[i - 1] = aba + _Tp{2} * (ab + _Tp(i)) * _Tp(i - 1);
      diag[i - 1] = -diag[i-1] / abp2i / (abp2i - _Tp{2});
    }

  for (int i = 2; i <= m - 1; ++i )
    {
      abp2i = ab + _Tp(2 * i);
      subd[i-1] = _Tp(i) * (alpha + _Tp(i)) / (abp2i - _Tp{1}) * (beta + _Tp(i))
		/ (abp2i * abp2i) * (ab + _Tp(i)) / (abp2i + _Tp{1});
    }
  subd[m-1] = _Tp{0};
  for (int i = 0; i < m; ++i)
    subd[i] = std::sqrt(subd[i]);
}
