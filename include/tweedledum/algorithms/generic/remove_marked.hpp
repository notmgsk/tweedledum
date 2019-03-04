/*--------------------------------------------------------------------------------------------------
| This file is distributed under the MIT License.
| See accompanying file /LICENSE for details.
| Author(s): Bruno Schmitt
*-------------------------------------------------------------------------------------------------*/
#pragma once

#include <cassert>

namespace tweedledum {

/*! \brief Generic function to remove gates.
 *
 * **Required network functions:**
 * - `add_qubit`
 * - `emplace_gate`
 * - `foreach_cqubit`
 * - `foreach_cgate`
 * - `rewire`
 * - `rewire_map`
 */
template<class NetworkDest, class NetworkSrc>
void remove_marked(NetworkDest& dest, NetworkSrc const& src)
{
	assert(dest.size() == 0);
	src.foreach_cqubit([&](std::string const& qlabel) {
		dest.add_qubit(qlabel);
	});
	src.foreach_cgate([&](auto const& node) {
		if (src.visited(node)) {
			return;
		}
		dest.emplace_gate(node.gate);
	});
	dest.rewire(src.rewire_map());
}

/*! \brief Generic function to remove gates.
 *
 * **Required network functions:**
 * - `add_qubit`
 * - `emplace_gate`
 * - `foreach_cqubit`
 * - `foreach_cgate`
 * - `rewire`
 * - `rewire_map`
 */
template<class Network>
Network remove_marked(Network const& src)
{
	Network dest;
	src.foreach_cqubit([&](std::string const& qlabel) {
		dest.add_qubit(qlabel);
	});
	src.foreach_cgate([&](auto const& node) {
		if (src.visited(node)) {
			return;
		}
		dest.emplace_gate(node.gate);

	});
	dest.rewire(src.rewire_map());
	return dest;
}

} // namespace tweedledum
