/* quadrature/triangle_rules.h
 *
 * Copyright (C) 2018-2019 Free Software Foundation, Inc.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef TRIANGLE_RULES_H
#define TRIANGLE_RULES_H 1

#include <vector>
#include <array>

namespace __gnu_cxx
{
  /**
   * A generic class for a triangle integration rule.
   */
  template<typename _Tp>
    class triangle_rule
    {

    public:

      triangle_rule()
      : _M_order{0}, _M_point{}, _M_weight{}
      { }

      triangle_rule(const std::size_t __order,
		    const std::vector<_Tp>& __weight,
		    const std::vector<std::array<_Tp, 3>>& __point)
      : _M_order(__order), _M_weight(__weight), _M_point(__point)
      { }

      std::size_t
      order() const
      { return this->_M_order; }

      void
      point(const std::size_t __index, _Tp& __weight,
	    std::array<_Tp,3>& __point) const
      {
	__weight = this->_M_weight[__index];
	__point = this->_M_point[__index];
      }

      static const int _S_num_tri_rules = 6;
      static const int _S_max_tri_order = 10;

    private:

      std::size_t _M_order;
      std::vector<_Tp> _M_weight;
      std::vector<std::array<_Tp, 3>> _M_point;

      static std::size_t _S_tri_order[_S_num_tri_rules];
      static _Tp _S_tri_weight[_S_num_tri_rules][_S_max_tri_order];
      static _Tp _S_tri_point[_S_num_tri_rules][_S_max_tri_order][3];

    };

  /**
   * @todo Figure out a way to mark wich rules are interior and which
   * rules have points at vertices or edges.  One would be able easily
   * and significantly to reduce the number of function evaluations
   * over a triangulated domain.
   */

  /**
   * An array indicating the order of the canned rules.
   * This is the number of the weights and the barycenters for the rule.
   */
  template<typename _Tp>
    std::size_t
    triangle_rule<_Tp>::_S_tri_order[triangle_rule<_Tp>::_S_num_tri_rules]
    { 1, 3, 3, 4, 4, 7 };

  /**
   * An array of weight arrays for the canned rules.
   */
  template<typename _Tp>
    _Tp
    triangle_rule<_Tp>::_S_tri_weight[triangle_rule<_Tp>::_S_num_tri_rules][triangle_rule<_Tp>::_S_max_tri_order]
    {
      { _Tp{1} },
      { _Tp{1} / _Tp{3}, _Tp{1} / _Tp{3}, _Tp{1} / _Tp{3} },
      { _Tp{1} / _Tp{3}, _Tp{1} / _Tp{3}, _Tp{1} / _Tp{3} },
      { _Tp{3} / _Tp{4}, _Tp{1} / _Tp{12}, _Tp{1} / _Tp{12}, _Tp{1} / _Tp{12} },
      { -_Tp{9} / _Tp{16}, _Tp{25} / _Tp{48}, _Tp{25} / _Tp{48}, _Tp{25} / _Tp{48} },
      { _Tp{9} / _Tp{20}, _Tp{2} / _Tp{15}, _Tp{2} / _Tp{15}, _Tp{2} / _Tp{15},
	_Tp{1} / _Tp{20}, _Tp{1} / _Tp{20}, _Tp{1} / _Tp{20} }
    };

  /**
   * An array of barycentric point coordinate arrays for the canned rules.
   */
  template<typename _Tp>
    _Tp
    triangle_rule<_Tp>::_S_tri_point[triangle_rule<_Tp>::_S_num_tri_rules][triangle_rule<_Tp>::_S_max_tri_order][3]
    {
      // One point in the center.
      { { _Tp{1} / _Tp{3}, _Tp{1} / _Tp{3}, _Tp{1} / _Tp{3} } },

      // Three points near vertices.
      {
	{ _Tp{2} / _Tp{3}, _Tp{1} / _Tp{6}, _Tp{1} / _Tp{6} },
	{ _Tp{1} / _Tp{6}, _Tp{2} / _Tp{3}, _Tp{1} / _Tp{6} },
	{ _Tp{1} / _Tp{6}, _Tp{1} / _Tp{6}, _Tp{2} / _Tp{3} }
      },

      // Three points at edges.
      {
	{ _Tp{0}, _Tp{1} / _Tp{2}, _Tp{1} / _Tp{2} },
	{ _Tp{1} / _Tp{2}, _Tp{0}, _Tp{1} / _Tp{2} },
	{ _Tp{1} / _Tp{2}, _Tp{1} / _Tp{2}, _Tp{0} }
      },

      // One point at center and three at vertices.
      {
	{ _Tp{1} / _Tp{3}, _Tp{1} / _Tp{3}, _Tp{1} / _Tp{3} },
	{ _Tp{1}, _Tp{0}, _Tp{0} },
	{ _Tp{0}, _Tp{1}, _Tp{0} },
	{ _Tp{0}, _Tp{0}, _Tp{1} }
      },

      // One point at center and three near vertices.
      {
	{ _Tp{1} / _Tp{3}, _Tp{1} / _Tp{3}, _Tp{1} / _Tp{3} },
	{ _Tp{11} / _Tp{15}, _Tp{2} / _Tp{15}, _Tp{2} / _Tp{15} },
	{ _Tp{2} / _Tp{15}, _Tp{11} / _Tp{15}, _Tp{2} / _Tp{15} },
	{ _Tp{2} / _Tp{15}, _Tp{2} / _Tp{15}, _Tp{11} / _Tp{15} }
      },

      // One point at center, three points on edges, three points at vertices.
      {
	{ _Tp{1} / _Tp{3}, _Tp{1} / _Tp{3}, _Tp{1} / _Tp{3} },
	{ _Tp{0}, _Tp{1} / _Tp{2}, _Tp{1} / _Tp{2} },
	{ _Tp{1} / _Tp{2}, _Tp{0}, _Tp{1} / _Tp{2} },
	{ _Tp{1} / _Tp{2}, _Tp{1} / _Tp{2}, _Tp{0} },
	{ _Tp{1}, _Tp{0}, _Tp{0} },
	{ _Tp{0}, _Tp{1}, _Tp{0} },
	{ _Tp{0}, _Tp{0}, _Tp{1} }
      }
    };

} // namespace __gnu_cxx

//#include <ext/triangle_rules.tcc>

#endif // TRIANGLE_RULES_H
