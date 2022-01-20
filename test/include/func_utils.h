#ifndef FUNC_UTILS_H
#define FUNC_UTILS_H 1

/**
 * Make a function object absorbing all parameters but a single
 * real scalar argument.
 */
template<typename Tp, typename FuncTp, typename... Parms,
	 typename Ret = std::invoke_result_t<FuncTp, Tp, Parms...>>
  std::function<Ret(Tp)>
  make_function(FuncTp f, Parms... p)
  { return [f, p...](Tp x)->Ret{ return f(x, p...); }; }

/**
 * Function wrapper to count evaluations of the target function.
 *
 * Note to self: Argument deduction order matters.
 */
template<typename Tp, typename FuncTp,
	 typename Ret = std::invoke_result_t<FuncTp, Tp>>
  struct counted_function
  {
    counted_function(FuncTp f)
    : m_func(f), m_neval(new int{0})
    { }

    Ret
    operator()(Tp x) const
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

    FuncTp m_func;

    mutable std::shared_ptr<int> m_neval;
  };

/**
 * Make a counted function object.
 *
 * This version wraps a more generic function object.
 *
 * Note to self: Argument deduction order matters.
 */
template<typename Tp, typename FuncTp, typename... Parms,
	 typename Ret = std::invoke_result_t<FuncTp, Tp, Parms...>>
  counted_function<Tp, std::function<Ret(Tp)>>
  make_counted_function(FuncTp&& f, Parms... p)
  {
    return counted_function<Tp, std::function<Ret(Tp)>>(make_function<Tp>(std::forward(f), p...));
  }

#endif // FUNC_UTILS_H
