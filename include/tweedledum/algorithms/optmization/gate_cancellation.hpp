/*--------------------------------------------------------------------------------------------------
| This file is distributed under the MIT License.
| See accompanying file /LICENSE for details.
| Author(s): Bruno Schmitt
*-------------------------------------------------------------------------------------------------*/
#pragma once

#include "../generic/remove_marked.hpp"
#include <cassert>

namespace tweedledum {

/*! \brief Cancellation of consecutive adjoint gates.
 *
 * **Required network functions:**
 * - `add_gate`
 * - `foreach_cgate`
 * - `foreach_child`
 * - `get_node`
 * - `visited`
 * - `set_visited`
 */
template<typename Network>
Network gate_cancellation(Network const& network)
{
	uint32_t num_deletions = 0u;
	network.foreach_cgate([&](auto& node) {
		network.foreach_child(node, [&](auto child_index) {
			auto& child = network.get_node(child_index);
			if (network.visited(child)) {
				return;
			}
			if (node.gate.is_adjoint(child.gate)) {
				network.set_visited(node, 1);
				network.set_visited(child, 1);
				num_deletions += 2;
				return;
			}
		});
	});
	if (num_deletions == 0) { 
		return network;
	}
	return remove_marked(network);
}

} // namespace tweedledum
