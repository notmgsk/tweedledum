/*--------------------------------------------------------------------------------------------------
| This file is distributed under the MIT License.
| See accompanying file /LICENSE for details.
| Author(s): Bruno Schmitt
*-------------------------------------------------------------------------------------------------*/
#include <catch.hpp>
#include <sstream>
#include <vector>
#include <tweedledum/gates/mcmt_gate.hpp>
#include <tweedledum/io/write_unicode.hpp>
#include <tweedledum/networks/netlist.hpp>
#include <tweedledum/networks/qubit.hpp>

TEST_CASE("Write three parallel NOTs in Unicode", "[write_unicode]")
{
	using namespace tweedledum;
	netlist<mcmt_gate> mcmt_netlist;
	mcmt_netlist.add_qubit();
	mcmt_netlist.add_qubit();
	mcmt_netlist.add_qubit();
	std::vector<qubit_id> controls;
	std::vector<qubit_id> target = {qubit_id(0), qubit_id(1), qubit_id(2)};
	mcmt_netlist.add_gate(gate::mcx, controls, target);

	std::ostringstream os;
	write_unicode(mcmt_netlist, false, os); /* works */
	write_unicode(mcmt_netlist, true, os);  /* crashes */
}
