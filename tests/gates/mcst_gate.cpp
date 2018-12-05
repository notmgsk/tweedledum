/*-------------------------------------------------------------------------------------------------
| This file is distributed under the MIT License.
| See accompanying file /LICENSE for details.
| Author(s): Bruno Schmitt
*------------------------------------------------------------------------------------------------*/
#include <algorithm>
#include <catch.hpp>
#include <random>
#include <tweedledum/gates/gate_set.hpp>
#include <tweedledum/gates/mcst_gate.hpp>
#include <tweedledum/gates/gate_base.hpp>
#include <tweedledum/utils/angle.hpp>
#include <vector>

TEST_CASE("MCST gate constructor", "[mcst_gate]")
{
	using namespace tweedledum;
	// Generate random qubit ids
	std::random_device rd;
	auto qid_target = rd();
	auto qid_control0 = rd();
	auto qid_control1 = rd();
	INFO("Target qubit id is " << qid_target);
	INFO("Control qubits ids are " << qid_control0 << " and " << qid_control1);

	SECTION("Single-qubit gate")
	{
		mcst_gate h_gate(gate::hadamard, qid_target);
		CHECK(h_gate.operation() == gate_set::hadamard);
		CHECK(h_gate.num_controls() == 0u);
		CHECK(h_gate.num_targets() == 1u);
		CHECK(h_gate.rotation_angle() == symbolic_angles::one_half);
	}

	SECTION("Controlled gate")
	{
		mcst_gate cx_gate(gate::cx, qid_control0, qid_target);
		CHECK(cx_gate.operation() == gate_set::cx);
		CHECK(cx_gate.num_controls() == 1u);
		CHECK(cx_gate.num_targets() == 1u);
		CHECK(cx_gate.rotation_angle() == symbolic_angles::one_half);

		// Using vectors to define a single controlled gate
		auto qids_control = std::vector<uint32_t>({qid_control0});
		auto qids_target = std::vector<uint32_t>({qid_target});
		mcst_gate cx_gate2(gate::cx, qids_control, qids_target);
		CHECK(cx_gate2.operation() == gate_set::cx);
		CHECK(cx_gate2.num_controls() == 1u);
		CHECK(cx_gate2.num_targets() == 1u);
		CHECK(cx_gate2.rotation_angle() == symbolic_angles::one_half);
	}

	SECTION("Multiple controlled gate")
	{
		auto qids_control = std::vector<uint32_t>({qid_control0, qid_control1});
		auto qids_target = std::vector<uint32_t>({qid_target});
		mcst_gate mcx_gate(gate::mcx, qids_control, qids_target);
		CHECK(mcx_gate.operation() == gate_set::mcx);
		CHECK(mcx_gate.num_controls() == 2u);
		CHECK(mcx_gate.num_targets() == 1u);
		CHECK(mcx_gate.rotation_angle() == symbolic_angles::one_half);
	}
}

TEST_CASE("MCST gate iterators", "[mcst_gate]")
{
	using namespace tweedledum;
	// Generate random qubit ids
	std::random_device rd;
	auto qid_target = rd();
	auto qid_control0 = rd();
	auto qid_control1 = rd();
	if (qid_control0 > qid_control1) {
		std::swap(qid_control0, qid_control1);
	}
	INFO("Target qubit id is " << qid_target);
	INFO("Control qubits ids are " << qid_control0 << " and " << qid_control1);

	SECTION("Single-qubit gate")
	{
		mcst_gate h_gate(gate::hadamard, qid_target);
		h_gate.foreach_target([&qid_target](auto qid) { CHECK(qid_target == qid); });
		h_gate.foreach_control([](auto qid) {
			// This function should not be called
			(void) qid;
			CHECK(0);
		});
	}

	SECTION("Controlled gate")
	{
		auto control = std::vector<uint32_t>({qid_control0});
		auto qids_target = std::vector<uint32_t>({qid_target});
		mcst_gate cx_gate(gate::cx, control, qids_target);
		cx_gate.foreach_target([&qid_target](auto qid) { CHECK(qid_target == qid); });
		cx_gate.foreach_control([&qid_control0](auto qid) { CHECK(qid_control0 == qid); });
	}

	SECTION("Multiple controlled gate")
	{
		auto qids_control = std::vector<uint32_t>({qid_control0, qid_control1});
		auto qids_target = std::vector<uint32_t>({qid_target});
		mcst_gate mcx_gate(gate::mcx, qids_control, qids_target);
		mcx_gate.foreach_target([&qid_target](auto qid) { CHECK(qid_target == qid); });
		auto i = 0u;
		mcx_gate.foreach_control([&i, &qid_control0, &qid_control1](auto qid) {
			CHECK(((i == 0 && qid_control0 == qid) || (i == 1 && qid_control1 == qid)) == true);
			++i;
		});
	}
}