/*------------------------------------------------------------------------------
| This file is distributed under the MIT License.
| See accompanying file /LICENSE for details.
| Author(s): Mathias Soeken
*-----------------------------------------------------------------------------*/
#include <catch.hpp>
#include <iostream>
#include <kitty/dynamic_truth_table.hpp>
#include <sstream>
#include <tweedledum/algorithms/synthesis/single_target_gates.hpp>
#include <tweedledum/io/quil.hpp>
#include <tweedledum/io/write_projectq.hpp>
#include <tweedledum/networks/dag_path.hpp>
#include <tweedledum/networks/gates/qc_gate.hpp>

TEST_CASE("1-input function x", "stg_synthesis")
{
	using namespace tweedledum;
	dag_path<qc_gate> network;
	network.allocate_qubit();
	network.allocate_qubit();

	kitty::dynamic_truth_table func(1);
	func._bits[0] = 0x2;

	stg_from_db()(network, func, std::vector<uint8_t>{{0u, 1u}});
	std::ostringstream os;
	write_quil(network, os);
	CHECK(os.str() == "CNOT 0 1\n");
}

TEST_CASE("1-input function !x", "stg_synthesis")
{
	using namespace tweedledum;
	dag_path<qc_gate> network;
	network.allocate_qubit();
	network.allocate_qubit();

	kitty::dynamic_truth_table func(1);
	func._bits[0] = 0x1;

	stg_from_db()(network, func, std::vector<uint8_t>{{0u, 1u}});
	std::ostringstream os;
	write_quil(network, os);
	CHECK(os.str() == "CNOT 0 1\nX 1\n");
}

TEST_CASE("1-input function 0", "stg_synthesis")
{
	using namespace tweedledum;
	dag_path<qc_gate> network;
	network.allocate_qubit();
	network.allocate_qubit();

	kitty::dynamic_truth_table func(1);
	func._bits[0] = 0x0;

	stg_from_db()(network, func, std::vector<uint8_t>{{0u, 1u}});
	std::ostringstream os;
	write_quil(network, os);
	CHECK(os.str() == "");
}

TEST_CASE("1-input function 1", "stg_synthesis")
{
	using namespace tweedledum;
	dag_path<qc_gate> network;
	network.allocate_qubit();
	network.allocate_qubit();

	kitty::dynamic_truth_table func(1);
	func._bits[0] = 0x3;

	stg_from_db()(network, func, std::vector<uint8_t>{{0u, 1u}});
	std::ostringstream os;
	write_quil(network, os);
	CHECK(os.str() == "X 1\n");
}

TEST_CASE("2-input AND", "stg_synthesis")
{
	using namespace tweedledum;
	dag_path<qc_gate> network;
	network.allocate_qubit();
	network.allocate_qubit();
	network.allocate_qubit();

	kitty::dynamic_truth_table func(2);
	func._bits[0] = 0x8;

	stg_from_db()(network, func, std::vector<uint8_t>{{0u, 1u, 2u}});
	CHECK(network.num_gates() == 16u);
}

TEST_CASE("2-input OR", "stg_synthesis")
{
	using namespace tweedledum;
	dag_path<qc_gate> network;
	network.allocate_qubit();
	network.allocate_qubit();
	network.allocate_qubit();

	kitty::dynamic_truth_table func(2);
	func._bits[0] = 0xe;

	stg_from_db()(network, func, std::vector<uint8_t>{{0u, 1u, 2u}});
	CHECK(network.num_gates() == 21u);
}

TEST_CASE("2-input XOR", "stg_synthesis")
{
	using namespace tweedledum;
	dag_path<qc_gate> network;
	network.allocate_qubit();
	network.allocate_qubit();
	network.allocate_qubit();

	kitty::dynamic_truth_table func(2);
	func._bits[0] = 0x6;

	stg_from_db()(network, func, std::vector<uint8_t>{{0u, 1u, 2u}});
	CHECK(network.num_gates() == 4u);
}
