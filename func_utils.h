#ifndef FUNC_UTILS_H
#define FUNC_UTILS_H 1


/**
 * Make a function object absorbing all parameters but a single
 * real scalar argument.
 *
 * This version wraps a function pointer taking and returning the same scalar real type.
 */
template<typename _Tp, typename... _Parms>
  std::function<_Tp(_Tp)>
  make_function(_Tp(f)(_Tp, _Parms...), _Parms... p)
  { return [f, p...](_Tp x)->_Tp{ return f(x, p...); }; }

/**
 * Make a function object absorbing all parameters but a single
 * real scalar argument.
 *
 * This version wraps a more generic function object.
 *
 * Note to self: Argument deduction order matters.
 */
template<typename _Tp, typename _FuncTp, typename... _Parms,
	 typename _Ret = std::invoke_result_t<_FuncTp, _Tp, _Parms...>>
  std::function<_Ret(_Tp)>
  make_function(_FuncTp&& f, _Parms... p)
  {
    return [func = std::forward(f), p...](_Tp x)->_Ret{ return func(x, p...); };
  }

/**
 * Function wrapper to count evaluations of the target function.
 *
 * Note to self: Argument deduction order matters.
 */
template<typename _Tp, typename _FuncTp,
	 typename _Ret = std::invoke_result_t<_FuncTp, _Tp>>
  struct counted_function
  {
    //counted_function(const _FuncTp& f)
    //: m_func(f), m_neval(new int{0})
    //{ }

    //counted_function(_FuncTp&& f)
    //: m_func(std::forward(f)), m_neval(new int{0})
    //{ }

    counted_function(_FuncTp f)
    : m_func(f), m_neval(new int{0})
    { }

    _Ret
    operator()(_Tp x) const
    {
      ++(*this->m_neval);
      return this->m_func(x);
    }

    int
    num_evals() const
    { return *this->m_neval; }

    void
    num_evals(int num)
    { *this->m_neval = num; }

    void
    reset()
    { this->num_evals(0); }

  private:

    _FuncTp m_func;

    mutable std::shared_ptr<int> m_neval;
  };

/**
 * Make a counted function object.
 *
 * This version wraps a more generic function object.
 *
 * Note to self: Argument deduction order matters.
 */
template<typename _Tp, typename _FuncTp, typename... _Parms,
	 typename _Ret = std::invoke_result_t<_FuncTp, _Tp, _Parms...>>
  counted_function<_Tp, std::function<_Ret(_Tp)>>
  make_counted_function(_FuncTp&& f, _Parms... p)
  {
    return counted_function<_Tp, std::function<_Ret(_Tp)>>(make_function<_Tp>(std::forward(f), p...));
  }

#endif // FUNC_UTILS_H
