/*--------------------------------------------------------------------------------------------------
| This file is distributed under the MIT License.
| See accompanying file /LICENSE for details.
| Author(s): Mathias Soeken, Bruno Schmitt
*-------------------------------------------------------------------------------------------------*/
#pragma once

#include "../../gates/gate_set.hpp"

#include <cstdint>
#include <iostream>

namespace tweedledum {

/*! \brief Simulate a quantum circuit that has only classical gates.
 *
 * **Required gate functions:**
 * - `foreach_control`
 * - `foreach_target`
 *
 * **Required network functions:**
 * - `foreach_cgate`
 * 
 * \algtype simulation
 * \algexpects A Toffoli network
 * \algreturns The simulated pattern
 */
template<typename Network>
uint64_t simulate_pattern_classical(Network& network, uint64_t pattern)
{
	assert(network.num_qubits() <= 64);
	network.foreach_cgate([&](auto const& node) {
		auto const& gate = node.gate;
		switch (node.gate.operation()) {
		default:
			std::cerr << "[w] non-classical gate, abort simulation\n";
			pattern = 0ull;
			return false;

		case gate_set::pauli_x:
			gate.foreach_target([&](auto qid) {
				pattern ^= (1ull << qid);
			});
			break;

		case gate_set::cx:
			gate.foreach_control([&](auto control_qid) {
				auto temp_pattern = pattern;
				if (control_qid.is_complemented()) {
					temp_pattern ^= (1ull << control_qid);
				}
				gate.foreach_target([&](auto target_qid) {
					if ((temp_pattern >> control_qid) & 1ull) {
						pattern ^= (1ull << target_qid);
					}
				});
			});
			break;

		case gate_set::mcx: {
			auto control_mask = 0ull;
			auto target_mask = 0ull;
			auto temp_pattern = pattern;
			gate.foreach_control([&](auto qid) {
				control_mask |= (1ull << qid);
				if (qid.is_complemented()) {
					temp_pattern ^= (1ull << qid);
				}
			});
			gate.foreach_target([&](auto qid) {
				target_mask |= (1ull << qid);
			});
			if ((temp_pattern & control_mask) == control_mask) {
				pattern ^= target_mask;
			}
		} break;
		}
		return true;
	});
	return pattern;
}

} // namespace tweedledum
