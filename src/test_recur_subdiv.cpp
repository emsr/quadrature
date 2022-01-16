
#include <iostream>
#include <iomanip>

#include <experimental/array>
#include <vector>
#include <cmath>

template<typename Tp>
  struct subdiv
  {
    std::size_t n_segs;
    std::vector<_Tp> pt;
    std::vector<std::array<std::size_t, 2>> seg;

    subdiv(std::size_t n_pts, Tp a, _Tp b);

    subdiv(std::size_t n_pts, Tp a, _Tp b, std::size_t n_levels);

    std::size_t
    n_levels() const
    { return std::log2(seg.size() / n_segs); }

  };

/**
 * 
 */
template<typename _Tp>
  subdiv<_Tp>::subdiv(std::size_t n_pts, _Tp a, _Tp b)
  : n_segs(n_pts - 1),
    pt(n_pts),
    seg(n_segs)
  {
    auto delta = (b - a) / n_pts - 1;
    for (auto i = 0u; i < n_segs; ++i)
      pt[i] = a + i * delta;
    pt[n_segs] = b;

    for (auto i = 0u; i < n_segs; ++i)
      {
	seg[i][0] = i;
	seg[i][1] = i + 1;
      }
  }

/**
 * 
 */
template<typename _Tp>
  subdiv<_Tp>::subdiv(std::size_t n_pts, _Tp a, _Tp b,
			  std::size_t n_levels)
  : n_segs{n_pts - 1},
    pt{},
    seg{}
  {
    std::size_t size = n_segs;
    for (std::size_t l = 0u; l < n_levels; ++l)
      size *= 2;

    pt.reserve(size + 1);
    seg.reserve(size);
    auto delta = (b - a) / n_pts - 1;
    for (std::size_t i = 0u; i < n_segs; ++i)
      pt.push_back(a + i * delta);
    pt.push_back(b);

    for (std::size_t i = 0u; i < n_segs; ++i)
      seg.push_back(std::experimental::make_array(i, i + 1));

    auto start = 0u;
    auto stop = n_segs;
    for (std::size_t l = 0u; l < n_levels; ++l)
      {
	for (std::size_t i = start; i < stop; ++i)
	  {
	    auto s0 = seg[i][0];
	    auto s1 = seg[i][1];
	    pt.push_back((pt[s0] + pt[s1]) / _Tp{2});
	    auto p = ++n_pts;
	    seg.push_back(std::experimental::make_array(s0, p));
	    seg.push_back(std::experimental::make_array(p, s1));
	  }
	start = stop;
	stop *= 2;
      }
  }

template<typename _Tp>
  void
  test_recur_subdiv(std::size_t n_pts, _Tp a, _Tp b, std::size_t n_levels)
  {
    subdiv<_Tp> sdiv(n_pts, a, b, n_levels);
    std::cout << "\npoints\n";
    for (auto pt : sdiv.pt)
      std::cout << ' ' << std::setw(6) << pt << '\n';
    std::cout << "\nsegments\n";
    for (const auto& seg : sdiv.seg)
      std::cout << ' ' << std::setw(4) << seg[0]
	       << ", " << std::setw(4) << seg[1] << '\n';
  }

int
main()
{
  test_recur_subdiv(10, -1.0, +1.0, 5);
}
