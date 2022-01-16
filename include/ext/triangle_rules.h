/* quadrature/triangle_rules.h
 *
 * Copyright (C) 2018-2020 Free Software Foundation, Inc.
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
  template<typename Tp>
    class triangle_rule
    {

    public:

      triangle_rule()
      : m_order{0}, m_point{}, m_weight{}
      { }

      triangle_rule(const std::size_t order,
		    const std::vector<Tp>& weight,
		    const std::vector<std::array<Tp, 3>>& point)
      : m_order(order), m_weight(weight), m_point(point)
      { }

      std::size_t
      order() const
      { return this->m_order; }

      void
      point(const std::size_t index, Tp& weight,
	    std::array<Tp,3>& point) const
      {
	weight = this->m_weight[index];
	point = this->m_point[index];
      }

      static const int s_num_tri_rules = 6;
      static const int s_max_tri_order = 10;

    private:

      std::size_t m_order;
      std::vector<Tp> m_weight;
      std::vector<std::array<Tp, 3>> m_point;

      static std::size_t s_tri_order[s_num_tri_rules];
      static Tp s_tri_weight[s_num_tri_rules][s_max_tri_order];
      static Tp s_tri_point[s_num_tri_rules][s_max_tri_order][3];

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
  template<typename Tp>
    std::size_t
    triangle_rule<Tp>::s_tri_order[triangle_rule<Tp>::s_num_tri_rules]
    { 1, 3, 3, 4, 4, 7 };

  /**
   * An array of weight arrays for the canned rules.
   */
  template<typename Tp>
    Tp
    triangle_rule<Tp>::s_tri_weight[triangle_rule<Tp>::s_num_tri_rules][triangle_rule<Tp>::s_max_tri_order]
    {
      { Tp{1} },
      { Tp{1} / Tp{3}, Tp{1} / Tp{3}, Tp{1} / Tp{3} },
      { Tp{1} / Tp{3}, Tp{1} / Tp{3}, Tp{1} / Tp{3} },
      { Tp{3} / Tp{4}, Tp{1} / Tp{12}, Tp{1} / Tp{12}, Tp{1} / Tp{12} },
      { -Tp{9} / Tp{16}, Tp{25} / Tp{48}, Tp{25} / Tp{48}, Tp{25} / Tp{48} },
      { Tp{9} / Tp{20}, Tp{2} / Tp{15}, Tp{2} / Tp{15}, Tp{2} / Tp{15},
	Tp{1} / Tp{20}, Tp{1} / Tp{20}, Tp{1} / Tp{20} }
    };

  /**
   * An array of barycentric point coordinate arrays for the canned rules.
   */
  template<typename Tp>
    Tp
    triangle_rule<Tp>::s_tri_point[triangle_rule<Tp>::s_num_tri_rules][triangle_rule<Tp>::s_max_tri_order][3]
    {
      // One point in the center.
      { { Tp{1} / Tp{3}, Tp{1} / Tp{3}, Tp{1} / Tp{3} } },

      // Three points near vertices.
      {
	{ Tp{2} / Tp{3}, Tp{1} / Tp{6}, Tp{1} / Tp{6} },
	{ Tp{1} / Tp{6}, Tp{2} / Tp{3}, Tp{1} / Tp{6} },
	{ Tp{1} / Tp{6}, Tp{1} / Tp{6}, Tp{2} / Tp{3} }
      },

      // Three points at edges.
      {
	{ Tp{0}, Tp{1} / Tp{2}, Tp{1} / Tp{2} },
	{ Tp{1} / Tp{2}, Tp{0}, Tp{1} / Tp{2} },
	{ Tp{1} / Tp{2}, Tp{1} / Tp{2}, Tp{0} }
      },

      // One point at center and three at vertices.
      {
	{ Tp{1} / Tp{3}, Tp{1} / Tp{3}, Tp{1} / Tp{3} },
	{ Tp{1}, Tp{0}, Tp{0} },
	{ Tp{0}, Tp{1}, Tp{0} },
	{ Tp{0}, Tp{0}, Tp{1} }
      },

      // One point at center and three near vertices.
      {
	{ Tp{1} / Tp{3}, Tp{1} / Tp{3}, Tp{1} / Tp{3} },
	{ Tp{11} / Tp{15}, Tp{2} / Tp{15}, Tp{2} / Tp{15} },
	{ Tp{2} / Tp{15}, Tp{11} / Tp{15}, Tp{2} / Tp{15} },
	{ Tp{2} / Tp{15}, Tp{2} / Tp{15}, Tp{11} / Tp{15} }
      },

      // One point at center, three points on edges, three points at vertices.
      {
	{ Tp{1} / Tp{3}, Tp{1} / Tp{3}, Tp{1} / Tp{3} },
	{ Tp{0}, Tp{1} / Tp{2}, Tp{1} / Tp{2} },
	{ Tp{1} / Tp{2}, Tp{0}, Tp{1} / Tp{2} },
	{ Tp{1} / Tp{2}, Tp{1} / Tp{2}, Tp{0} },
	{ Tp{1}, Tp{0}, Tp{0} },
	{ Tp{0}, Tp{1}, Tp{0} },
	{ Tp{0}, Tp{0}, Tp{1} }
      }
    };

} // namespace __gnu_cxx

//#include <ext/triangle_rules.tcc>

#endif // TRIANGLE_RULES_H
