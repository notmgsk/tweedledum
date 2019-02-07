/*--------------------------------------------------------------------------------------------------
| This file is distributed under the MIT License.
| See accompanying file /LICENSE for details.
| Author(s): Bruno Schmitt, Mathias Soeken
*-------------------------------------------------------------------------------------------------*/
#include <catch.hpp>
#include <cstdint>
#include <tweedledum/algorithms/simulation/simulate_pattern_classical.hpp>
#include <tweedledum/algorithms/synthesis/tbs.hpp>
#include <tweedledum/gates/mcmt_gate.hpp>
#include <tweedledum/networks/netlist.hpp>

TEST_CASE("Verify result of TBS", "[simulate_pattern_classical]")
{
	using namespace tweedledum;
	std::vector<uint32_t> permutation = {0, 2, 3, 5, 7, 1, 4, 6};
	const auto network = tbs<netlist<mcmt_gate>>(permutation);

	for (auto i = 0u; i < permutation.size(); ++i) {
		CHECK(simulate_pattern_classical(network, i) == permutation[i]);
	}
}
