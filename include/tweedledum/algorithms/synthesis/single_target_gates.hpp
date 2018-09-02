/*------------------------------------------------------------------------------
| This file is distributed under the MIT License.
| See accompanying file /LICENSE for details.
| Author(s): Mathias Soeken
*-----------------------------------------------------------------------------*/
#pragma once

#include "../../io/dotqc.hpp"
#include "gray_synth.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/esop.hpp>
#include <kitty/operations.hpp>
#include <kitty/print.hpp>
#include <kitty/spectral.hpp>
#include <map>
#include <numeric>
#include <sstream>
#include <tweedledee/dotqc/dotqc.hpp>
#include <vector>

namespace tweedledum {

struct stg_from_pprm {
	template<class Network>
	void operator()(Network& net, kitty::dynamic_truth_table const& function,
	                std::vector<uint8_t> const& qubit_map)
	{
		const auto num_controls = function.num_vars();
		assert(qubit_map.size()
		       == static_cast<std::size_t>(num_controls) + 1u);

		for (auto const& cube : esop_from_pprm(function)) {
			assert(cube._bits == cube._mask); /* PPRM property */
			std::vector<uint32_t> controls;
			auto bits = cube._bits;
			for (auto v = 0; v < num_controls; ++v) {
				if (bits & 1) {
					controls.push_back(qubit_map[v]);
				}
				bits >>= 1;
			}
			net.add_toffoli(controls, {qubit_map.back()});
		}
	}
};

struct stg_from_pkrm {
	template<class Network>
	void operator()(Network& net, kitty::dynamic_truth_table const& function,
	                std::vector<uint8_t> const& qubit_map)
	{
		const auto num_controls = function.num_vars();
		assert(qubit_map.size()
		       == static_cast<std::size_t>(num_controls) + 1u);

		for (auto const& cube : esop_from_optimum_pkrm(function)) {
			std::vector<uint32_t> controls, negations;
			auto bits = cube._bits;
			auto mask = cube._mask;
			for (auto v = 0; v < num_controls; ++v) {
				if (mask & 1) {
					controls.push_back(qubit_map[v]);
					if (!(bits & 1))
						negations.push_back(qubit_map[v]);
				}
				bits >>= 1;
				mask >>= 1;
			}
			for (auto n : negations)
				net.add_gate(gate_kinds_t::cx, n);
			net.add_toffoli(controls, {qubit_map.back()});
			for (auto n : negations)
				net.add_gate(gate_kinds_t::cx, n);
		}
	}
};

struct stg_from_spectrum {
	inline double pi()
	{
		static double _pi = std::atan(1) * 4;
		return _pi;
	}

	template<class Network>
	void operator()(Network& net, kitty::dynamic_truth_table const& function,
	                std::vector<uint8_t> const& qubit_map)
	{
		const auto num_controls = function.num_vars();
		assert(qubit_map.size() == num_controls + 1u);

		auto g = kitty::extend_to(function, num_controls + 1);
		auto xt = g.construct();
		kitty::create_nth_var(xt, num_controls);
		g &= xt;

		std::vector<uint32_t> parities;
		std::vector<float> angles;

		float nom = pi();
		nom /= (1 << g.num_vars());

		const auto spectrum = kitty::rademacher_walsh_spectrum(g);
		for (auto i = 1u; i < spectrum.size(); ++i) {
			if (spectrum[i] == 0)
				continue;
			parities.push_back(i);
			angles.push_back(nom * spectrum[i]);
		}

		net.add_gate(gate_kinds_t::hadamard, qubit_map.back());
		gray_synth(net, parities, angles, qubit_map);
		net.add_gate(gate_kinds_t::hadamard, qubit_map.back());
	}
};

#define STG_FROM_DB_VERBOSE 0

struct stg_from_db {
	template<class Network>
	void operator()(Network& net, kitty::dynamic_truth_table const& function,
	                std::vector<uint8_t> const& qubit_map) const
	{
		const auto cls = kitty::get_spectral_class(function);
		auto db_repr = function.construct();
		db_repr._bits[0] = representatives[function.num_vars()][cls];

		std::vector<kitty::detail::spectral_operation> trans;
		const auto r1 = kitty::exact_spectral_canonization(
		    db_repr, [&trans](auto const& ops) {
			    std::copy(ops.begin(), ops.end(),
			              std::back_inserter(trans));
		    });
		const auto r2 = kitty::exact_spectral_canonization(
		    function, [&trans](auto const& ops) {
			    std::copy(ops.rbegin(), ops.rend(),
			              std::back_inserter(trans));
		    });

#if STG_FROM_DB_VERBOSE
		std::cout << "function: " << kitty::to_hex(function) << "\n";
		std::cout << "db repr.: " << kitty::to_hex(db_repr) << "\n";
		std::cout << "class: " << cls << "\n";
		std::cout << "representative(s): " << kitty::to_hex(r1) << " "
		          << kitty::to_hex(r2) << "\n";
		std::cout << "transformations:\n";
		for (auto const& t : trans) {
			std::cout << static_cast<int>(t._kind) << " " << t._var1
			          << " " << t._var2 << "\n";
		}
#endif

		auto [mat_r, vec_r]
		    = to_affine_operation(trans, function.num_vars() + 1);
#if STG_FROM_DB_VERBOSE
		std::cout << "transformations for RHS:\n";
		for (auto const& t : trans) {
			std::cout << static_cast<int>(t._kind) << " " << t._var1
			          << " " << t._var2 << "\n";
		}
#endif

		std::reverse(trans.begin(), trans.end());
		trans.erase(
		    std::remove_if(
		        trans.begin(), trans.end(),
		        [](auto const& t) {
			        return (t._kind
			                == kitty::detail::spectral_operation::kind::output_negation)
			               || (t._kind
			                   == kitty::detail::spectral_operation::
			                          kind::disjoint_translation);
		        }),
		    trans.end());
		auto [mat_l, vec_l]
		    = to_affine_operation(trans, function.num_vars());

#if STG_FROM_DB_VERBOSE
		std::cout << "transformations for LHS:\n";
		for (auto const& t : trans) {
			std::cout << static_cast<int>(t._kind) << " " << t._var1
			          << " " << t._var2 << "\n";
		}
#endif

#if STG_FROM_DB_VERBOSE
		std::cout << "matrix for LHS:\n";
		for (auto const& r : mat_l)
			std::cout << r << "\n";
		std::cout << "vector for LHS: " << vec_l << "\n";
#endif

		cnot_patel(net, mat_l, mat_l.size(), qubit_map);
		for (auto i = 0u; i < mat_l.size(); ++i)
			if ((vec_l >> i) & 1)
				net.add_gate(gate_kinds_t::pauli_x, qubit_map[i]);

		// check if ancilla is needed
		auto map_copy = qubit_map;
		if ((ancilla_needed[function.num_vars()] >> cls) & 1) {
			int32_t anc{-1};
			for (auto q = 0u; q < net.num_qubits(); ++q)
				if (std::find(qubit_map.begin(),
				              qubit_map.end(), q)
				    == qubit_map.end())
					anc = q;
			if (anc != -1)
				map_copy.push_back(anc);
			else {
				map_copy.push_back(net.num_qubits());
				net.allocate_qubit();
			}
		}
		std::istringstream buffer(circuits[function.num_vars()][cls]);
		opt_dotqc_reader reader(net, map_copy);
		tweedledee::dotqc_read(buffer, reader, identify_gate_kind());

#if STG_FROM_DB_VERBOSE
		std::cout << "matrix for RHS:\n";
		for (auto const& r : mat_r)
			std::cout << r << "\n";
		std::cout << "vector for RHS: " << vec_r << "\n";
#endif
		cnot_patel(net, mat_r, mat_r.size(), qubit_map);
		for (auto i = 0u; i < mat_r.size(); ++i)
			if ((vec_r >> i) & 1)
				net.add_gate(gate_kinds_t::pauli_x, qubit_map[i]);
	}

private:
	template<class Network>
	class opt_dotqc_reader : public tweedledee::dotqc_reader<gate_kinds_t> {
	public:
		opt_dotqc_reader(Network& net,
		                 std::vector<uint8_t> const& qubit_map)
		    : net(net)
		    , qubit_map(qubit_map)
		{}

		void on_gate(gate_kinds_t kind, std::string l) override
		{
			net.add_gate(kind, qubit_map[l[0] - 'a']);
		}

		void on_two_qubit_gate(gate_kinds_t kind, std::string l1,
		                       std::string l2) override
		{
			switch (kind) {
			case gate_kinds_t::pauli_x:
				kind = gate_kinds_t::cx;
				break;
			case gate_kinds_t::pauli_z:
				kind = gate_kinds_t::cz;
				break;
			default:
				break;
			}
			net.add_controlled_gate(kind, qubit_map[l1[0] - 'a'],
			                        qubit_map[l2[0] - 'a']);
		}

	private:
		Network& net;
		std::vector<uint8_t> const& qubit_map;
	};

	inline uint32_t swap_bits(uint32_t n, uint32_t i, uint32_t j) const
	{
		auto x = ((n >> i) & 1) ^ ((n >> j) & 1);
		return n ^ ((x << i) | (x << j));
	}

	inline uint32_t log2(uint32_t n) const
	{
		return _log2.find(n)->second;
	}

	inline std::pair<std::vector<uint32_t>, uint32_t> to_affine_operation(
	    std::vector<kitty::detail::spectral_operation> const& trans,
	    uint32_t nvars) const
	{
		/* initialize matrix and vector */
		std::vector<uint32_t> mat(nvars);
		uint32_t vec{0};
		uint32_t c{0};
		std::generate(mat.begin(), mat.end(),
		              [&c]() { return 1 << c++; });

		for (auto const& t : trans) {
			switch (t._kind) {
			default:
				assert(false);
			case kitty::detail::spectral_operation::kind::permutation: {
				const auto v1 = log2(t._var1);
				const auto v2 = log2(t._var2);
				std::swap(mat[v1], mat[v2]);
				vec = swap_bits(vec, v1, v2);
			} break;
			case kitty::detail::spectral_operation::kind::input_negation:
				vec ^= t._var1;
				break;
			case kitty::detail::spectral_operation::kind::output_negation:
				vec ^= 1 << (nvars - 1);
				break;
			case kitty::detail::spectral_operation::kind::spectral_translation: {
				const auto v1 = log2(t._var1);
				const auto v2 = log2(t._var2);
				mat[v1] ^= mat[v2];
				vec ^= ((vec >> v2) & 1) << v1;
			} break;
			case kitty::detail::spectral_operation::kind::disjoint_translation: {
				const auto v1 = log2(t._var1);
				mat[nvars - 1] ^= mat[v1];
				vec ^= ((vec >> v1) & 1) << (nvars - 1);
			} break;
			}
		}

		return {mat, vec};
	}

	static inline const std::map<uint32_t, uint32_t> _log2
	    = {{1, 0},  {2, 1},  {4, 2},  {8, 3},
	       {16, 4}, {32, 5}, {64, 6}, {128, 7}};
	static inline const std::vector<std::vector<uint32_t>> representatives
	    = {{0x0},
	       {0x0},
	       {0x0, 0x8},
	       {0x00, 0x80, 0x88},
	       {0x0000, 0x8000, 0x8080, 0x0888, 0x8888, 0x7080, 0x7880, 0x7888},
	       {0x00000000, 0x80000000, 0x80008000, 0x0a000010, 0x80808080,
	        0x88800800, 0x80088820, 0x88808080, 0x2a808080, 0x70080088,
	        0xf0800080, 0xc0c8c0c8, 0xf780b880, 0x9ba00000, 0xe8080808,
	        0x8808a808, 0xc8888888, 0x88888888, 0xd5808080, 0x70807080,
	        0xe1808880, 0xea808080, 0xcc808880, 0xd8808880, 0x7f008000,
	        0xe820c088, 0x4e144404, 0xe8a08880, 0xf8808880, 0xe2222220,
	        0xa038a028, 0xe6804c80, 0x7f808080, 0x0231da51, 0xa002bc88,
	        0xf8087888, 0xeca08088, 0xf0888888, 0x8a80cac0, 0x78807880,
	        0xf8284888, 0xfca08880, 0xdac08a80, 0x38887888, 0x78887888,
	        0xa6cc60a0, 0x62c8ea40, 0x6ac8e240}};
	static inline const uint64_t ancilla_needed[]
	    = {0b0, 0b0, 0b00, 0b010, 0b01001010, 0b000011111000011001010};
	// clang-format off
	static inline const std::vector<std::vector<std::string>> circuits
	    = {{""},
	       {""},
	       {"", "H c\nT a\nT b\nT c\ntof b a\ntof c b\ntof a c\nT* b\ntof a b\nT* a\nT* b\nT c\ntof c b\ntof a c\ntof b a\nH c\n"},
	       {"",
	        "H e\nT e\ntof b e\nT* e\ntof a e\nT e\ntof b e\nT* e\nH e\nH d\ntof d e\nT* e\ntof c e\nT e\ntof d e\nT* e\ntof c e\nT e\nH e\nT e\ntof b e\nT* e\ntof a e\nT e\ntof b e\nT* e\nH e\nT* e\ntof c e\nT e\ntof d e\nT* e\ntof c e\nT e\ntof d e\nH d\n",
	        "H d\ntof a d\nT a\nT b\nT* d\ntof d a\ntof b d\ntof a d\ntof a b\nT a\nT* b\nT* d\ntof d a\ntof d b\ntof a d\ntof a b\nT a\ntof b a\ntof d a\nH d\n"},
	       {"",
	        "H f\nT f\ntof c f\nT* f\nH f\ntof a f\nT f\ntof b f\nT* f\ntof a f\nT f\ntof b f\nT* f\nH f\nT f\ntof c f\nT* f\nH f\nH e\ntof e f\nT* f\ntof d f\nT f\ntof e f\nT* f\ntof d f\nT f\nH f\nT f\ntof c f\nT* f\nH f\nT f\ntof b f\nT* f\ntof a f\nT f\ntof b f\nT* f\ntof a f\nH f\nT f\ntof c f\nT* f\nH f\nT* f\ntof d f\nT f\ntof e f\nT* f\ntof d f\nT f\ntof e f\nH e\n",
	        "H d\nT d\ntof b d\nT* d\ntof a d\nT d\ntof b d\nT* d\nH d\nH e\ntof e d\nT* d\ntof c d\nT d\ntof e d\nT* d\ntof c d\nT d\nH d\nT d\ntof b d\nT* d\ntof a d\nT d\ntof b d\nT* d\nH d\nT* d\ntof c d\nT d\ntof e d\nT* d\ntof c d\nT d\ntof e d\nH e\n",
	        "H e\ntof b f\ntof a f\nT a\nT b\nT e\nT* f\ntof e a\ntof e b\ntof b f\ntof b e\ntof a f\ntof a e\nT* b\nT* a\nT e\ntof e b\ntof e a\ntof a e\ntof b e\nH f\ntof f c\nT f\nT* c\ntof c f\nH c\ntof c b\ntof a c\ntof b a\nT b\nT* a\nT c\ntof a b\ntof c b\ntof c a\ntof b c\nT* b\nH b\ntof b f\nT b\nT* f\ntof f b\nH f\ntof e d\ntof f d\ntof f e\ntof e f\ntof d f\nT* e\nT d\nT* f\ntof d e\ntof f e\ntof f d\ntof e f\nT e\nH e\ntof e b\nT e\nT* b\ntof b e\nH b\ntof c a\ntof b a\ntof b c\nT b\nT* c\nT a\ntof c b\ntof a b\ntof a c\ntof b c\ntof b a\nT* b\ntof c b\nH c\ntof c e\nT c\nT* e\ntof e c\nH e\ntof e f\ntof f d\nT* e\nT f\nT* d\ntof f e\ntof d e\ntof d f\ntof e f\ntof e d\nT e\ntof f e\nH e\ntof b a\ntof a b\ntof b a\n",
	        "H e\ntof a e\nT a\nT b\nT* e\ntof e a\ntof b e\ntof a e\ntof a b\nT a\nT* b\nT* e\ntof e a\ntof e b\ntof a e\ntof a b\nT a\ntof b a\ntof e a\nH e\n",
	        "H d\ntof b a\ntof d a\ntof d b\nT d\nT* b\nT a\nT c\ntof b d\ntof a d\ntof a b\ntof d b\ntof d a\nT* d\ntof b d\ntof d b\nH b\nH e\ntof b c\ntof e b\ntof b e\nT* b\nT* c\nT e\ntof c b\ntof e b\ntof b e\ntof b c\nT b\nT c\nT* e\ntof e b\ntof e c\ntof c e\nH b\ntof a d\ntof b d\ntof b a\nT b\nT* a\nT d\ntof a b\ntof d b\ntof d a\ntof b a\ntof b d\nT* b\ntof a b\ntof b a\nH a\nH e\ntof b d\ntof d b\ntof b d\ntof d a\ntof a d\ntof d a\n",
	        "H e\nT a\nT b\nT e\ntof b a\ntof e b\ntof a e\nT* b\ntof a b\nT* a\nT* b\nT e\ntof e b\ntof a e\ntof b a\nH e\nX c\nX d\nH f\nT f\ntof c f\nT* f\nH f\ntof a f\nT f\ntof b f\nT* f\ntof a f\nT f\ntof b f\nT* f\nH f\nT f\ntof c f\nT* f\nH f\nH e\ntof e f\nT* f\ntof d f\nT f\ntof e f\nT* f\ntof d f\nT f\nH f\nT f\ntof c f\nT* f\nH f\nT f\ntof b f\nT* f\ntof a f\nT f\ntof b f\nT* f\ntof a f\nH f\nT f\ntof c f\nT* f\nH f\nT* f\ntof d f\nT f\ntof e f\nT* f\ntof d f\nT f\ntof e f\nH e\nX c\nX d\nH e\nT c\nT d\nT e\ntof d c\ntof e d\ntof c e\nT* d\ntof c d\nT* c\nT* d\nT e\ntof e d\ntof c e\ntof d c\nH e\n",
	        "H e\ntof a e\nT a\nT b\nT* e\nT c\nT d\ntof b a\ntof e a\ntof e b\ntof d c\ntof e d\ntof b d\ntof a d\ntof a e\nT* a\nT* e\nT b\nT* d\nT* c\ntof b a\ntof b e\ntof c e\ntof c b\ntof c d\ntof e c\ntof e b\ntof a b\nT* d\nT e\ntof c d\ntof c e\nS c\nH c\ntof e d\ntof d e\ntof e c\ntof c e\ntof e c\n"},
	       {"",
	        "H g\ntof a b\nT a\nT g\nT* b\ntof g a\ntof b a\ntof b g\ntof a b\nT* a\nT* b\nT g\ntof g a\ntof g b\ntof b g\ntof a g\nH g\nS b\nH a\ntof a e\nT a\nT* e\ntof e a\ntof a e\ntof e a\ntof a e\nH a\ntof c a\ntof d c\ntof a c\ntof a d\nT a\nT* d\nT c\ntof d a\ntof c a\ntof c d\ntof a c\nT* a\nH a\ntof a e\nT a\nT* e\ntof e a\ntof a e\ntof e a\ntof a e\nH a\nH f\ntof f a\ntof g f\ntof f g\ntof a g\ntof a f\nT* a\nT* f\nT g\ntof f a\ntof g a\ntof a f\nT a\ntof g a\ntof f a\nH g\ntof g e\nT g\nT* e\ntof e g\ntof g e\ntof e g\ntof g e\nH g\ntof c d\ntof g d\ntof g c\nT g\nT* c\nT d\ntof c g\ntof d g\ntof d c\ntof g c\ntof g d\nT* g\ntof c g\nH c\ntof c e\nT c\nT* e\ntof e c\ntof c e\ntof e c\ntof c e\nH c\ntof a c\ntof a b\ntof f a\ntof c b\ntof c a\ntof c f\nT c\nT* f\nT a\nT* b\ntof f c\ntof a c\ntof b c\ntof a f\ntof b a\ntof c b\ntof c a\nH f\ntof f b\ntof f c\nT f\nT* c\nT* b\ntof c f\ntof b f\ntof c b\ntof f b\ntof f c\nT f\ntof c f\ntof b f\nH b\nH c\ntof c e\nT c\nT* e\ntof e c\ntof c e\ntof e c\ntof c e\nH c\ntof g c\ntof d g\ntof c g\ntof c d\nT c\nT* d\nT g\ntof d c\ntof g c\ntof g d\ntof c g\nT* c\nH c\ntof c e\nT c\nT* e\ntof e c\ntof c e\ntof e c\ntof c e\nH c\ntof a c\ntof b a\ntof c a\ntof c b\nT* c\nT b\nT* a\ntof b c\ntof a c\ntof a b\ntof c a\nT c\nH c\ntof c e\nT c\nT* e\ntof e c\ntof c e\ntof e c\ntof c e\nH c\ntof g d\ntof c d\ntof c g\nT c\nT* g\nT d\ntof g c\ntof d c\ntof d g\ntof c g\ntof c d\nT* c\ntof g c\nH g\ntof g e\nT g\nT* e\ntof e g\ntof g e\ntof e g\ntof g e\nH g\ntof a b\ntof g b\ntof g a\nT* g\nT a\nT* b\ntof a g\ntof b g\ntof b a\ntof g a\ntof g b\nT g\ntof a g\nH g\ntof f g\ntof g f\ntof f g\ntof g b\ntof b g\ntof g b\n",
	        "H e\nT e\ntof c e\nT* e\nH e\ntof a e\nT e\ntof b e\nT* e\ntof a e\nT e\ntof b e\nT* e\nH e\nT e\ntof c e\nT* e\nH e\nH f\ntof f e\nT* e\ntof d e\nT e\ntof f e\nT* e\ntof d e\nT e\nH e\nT e\ntof c e\nT* e\nH e\nT e\ntof b e\nT* e\ntof a e\nT e\ntof b e\nT* e\ntof a e\nH e\nT e\ntof c e\nT* e\nH e\nT* e\ntof d e\nT e\ntof f e\nT* e\ntof d e\nT e\ntof f e\nH f\n",
	        "H b\ntof d b\ntof b d\nT d\nT* b\nT* a\ntof b d\ntof d b\ntof b d\ntof d b\nX c\nH d\ntof d a\ntof a c\nX c\nT a\nT* d\nT c\ntof d a\ntof c a\ntof a c\nT* a\ntof d a\nH d\ntof b d\ntof d b\nT b\nT* d\ntof d b\ntof b d\ntof d b\ntof b d\nH b\nH f\ntof e b\ntof e f\ntof b f\ntof b e\nT* f\nT* b\nT e\ntof b f\ntof e f\ntof f b\nT f\ntof e f\ntof b f\nH e\ntof d e\ntof e d\nT d\nT* e\ntof e d\ntof d e\ntof e d\ntof d e\nX a\nH d\ntof d c\ntof d a\nX a\nT d\nT* c\nT a\ntof c d\ntof a d\ntof c a\ntof d a\ntof d c\nT* d\ntof c d\ntof a d\nH a\ntof e a\ntof a e\nT e\nT* a\ntof a e\ntof e a\ntof a e\ntof e a\nH e\ntof f e\ntof e f\ntof e b\nP* f\nT* f\nT e\nT* b\ntof e f\ntof b f\ntof b e\ntof d e\ntof e d\ntof f d\ntof f b\ntof f e\nT f\nT* e\ntof d f\ntof d e\nX e\nX a\nH g\ntof g e\ntof g d\nX d\nT g\nT e\nT d\ntof e g\ntof d g\ntof e d\ntof g d\ntof g e\nT g\ntof d g\nX g\nX b\nH d\nH g\ntof g b\nX b\nT g\nT b\ntof b g\ntof g b\ntof b g\ntof g b\nX g\nX b\nH g\ntof c g\ntof g a\ntof g c\nX a\nT g\nT* c\nT a\ntof c g\ntof a g\ntof g a\nT* g\ntof c g\nX g\nH c\ntof c b\nX b\nT c\nT b\ntof b c\ntof c b\ntof b c\ntof c b\nX c\nH c\ntof f c\ntof d f\ntof c f\ntof c d\nT* c\nT d\nT* f\ntof d c\ntof f c\ntof f d\ntof c f\nT c\nX b\nH c\ntof c b\nX b\nT c\nT b\ntof b c\ntof c b\ntof b c\ntof c b\nX c\nX b\nH c\ntof c g\ntof c a\nX g\nT c\nT* a\nT g\ntof a c\ntof g c\ntof a g\ntof c g\ntof c a\nT* c\ntof a c\ntof g c\nX a\nH g\ntof g b\nX b\nT g\nT b\ntof b g\ntof g b\ntof b g\ntof g b\nX g\nX e\nH g\ntof f g\ntof f e\ntof d f\ntof g e\ntof g f\ntof g d\nX e\nT g\nT* d\nT f\nT e\ntof d g\ntof f d\ntof e f\ntof f e\ntof g e\ntof g f\nX f\nH d\ntof d e\ntof d f\nX f\nT d\nT* e\nT f\ntof e d\ntof f d\ntof f e\ntof d e\ntof d f\nT* d\ntof f d\ntof e d\nX b\nH e\nH f\ntof f b\nX b\nT f\nT b\ntof b f\ntof f b\ntof b f\ntof f b\nX f\nX b\nH f\ntof c f\ntof f a\ntof f c\nX a\nT f\nT* c\nT a\ntof c f\ntof a f\ntof f a\nT* f\ntof c f\nX f\nH c\ntof c b\nX b\nT c\nT b\ntof b c\ntof c b\ntof b c\ntof c b\nX c\nH c\ntof g c\ntof e g\ntof c g\ntof c e\nT* c\nT e\nT* g\ntof e c\ntof g c\ntof g e\ntof c g\nT c\nX b\nH c\ntof c b\nX b\nT c\nT b\ntof b c\ntof c b\ntof b c\ntof c b\nX c\nX b\nH c\ntof c f\ntof c a\nX f\nT c\nT* a\nT f\ntof a c\ntof f c\ntof a f\ntof c f\ntof c a\nT* c\ntof a c\ntof f c\nH f\ntof f b\nX b\nT f\nT b\ntof b f\ntof f b\ntof b f\ntof f b\nX f\nH f\ntof g e\ntof f e\ntof f g\nT* f\nT g\nT* e\ntof g f\ntof e f\ntof e g\ntof f g\ntof f e\nT f\ntof g f\nH f\ntof e g\ntof g e\ntof e b\ntof b e\ntof e b\ntof b d\ntof d b\ntof b d\ntof d a\ntof a d\ntof d a\n",
	        "H d\nT d\ntof b d\nT* d\ntof a d\nT d\ntof b d\nT* d\nH d\nH f\ntof f d\nT* d\ntof c d\nT d\ntof f d\nT* d\ntof c d\nT d\nH d\nT d\ntof b d\nT* d\ntof a d\nT d\ntof b d\nT* d\nH d\nT* d\ntof c d\nT d\ntof f d\nT* d\ntof c d\nT d\ntof f d\nH f\n",
	        "tof d e\nH d\ntof e c\ntof c e\nT e\nT d\nT* c\ntof d e\ntof c e\ntof c d\ntof e c\nT* e\nT* c\nT d\ntof d e\ntof d c\ntof c d\ntof e d\nH d\nS c\nH e\ntof b a\ntof e a\ntof e b\nT e\nT* b\nT a\ntof b e\ntof a e\ntof a b\ntof e b\ntof e a\nT* e\ntof b e\ntof e b\nH b\nH f\ntof b d\ntof b f\nT* f\nT* d\nT b\ntof d f\ntof b f\ntof f d\nT f\ntof b f\ntof d f\nH b\ntof a e\ntof b e\ntof b a\nT b\nT* a\nT e\ntof a b\ntof e b\ntof e a\ntof b a\ntof b e\nT* b\ntof a b\ntof b a\nH a\ntof d a\ntof d f\ntof d c\ntof f d\ntof a d\ntof a f\ntof a c\nT f\nT* d\nT a\nT* c\ntof d f\ntof a f\ntof c f\ntof a d\ntof c a\ntof f c\ntof f a\nH a\nH d\ntof d c\ntof d f\nT d\nT* f\nT* c\ntof f d\ntof c d\ntof f c\ntof d c\ntof d f\nT d\ntof f d\ntof c d\nH c\ntof c d\ntof b e\ntof e b\ntof b e\ntof e d\ntof d e\ntof e d\ntof d c\ntof c d\ntof d c\ntof c f\ntof f c\ntof c f\ntof f a\ntof a f\ntof f a\n",
	        "H d\ntof a b\nX c\nX e\nT d\nT a\nT* c\nT e\nT b\ntof c d\ntof e d\ntof a b\ntof d c\ntof d e\nT* d\nT e\nT* c\ntof c d\ntof c e\ntof e c\ntof d c\nS e\nP* b\nH c\nX b\nX e\nH g\ntof g a\ntof g b\nX b\nT g\nT* a\nT b\ntof a g\ntof b g\ntof a b\ntof g b\ntof g a\nT* g\ntof a g\ntof b g\nX c\nH b\nH g\ntof g e\nX e\nT g\nT e\ntof e g\ntof g e\ntof e g\ntof g e\nX g\nX e\nH g\ntof d g\ntof g c\ntof g d\nX c\nT g\nT* d\nT c\ntof d g\ntof c g\ntof g c\nT* g\ntof d g\nH d\ntof d e\nX e\nT d\nT e\ntof e d\ntof d e\ntof e d\ntof d e\nX d\nX e\nH d\nH f\ntof f d\ntof f b\ntof d f\ntof d b\nT* d\nT* b\nT f\ntof b d\ntof f d\ntof d b\nT d\ntof f d\ntof b d\nX g\nH f\ntof f e\nX e\nT f\nT e\ntof e f\ntof f e\ntof e f\ntof f e\nX f\nX e\nH f\ntof f g\ntof f c\nX g\nT f\nT* c\nT g\ntof c f\ntof g f\ntof c g\ntof f g\ntof f c\nT* f\ntof c f\ntof g f\nX c\nH g\ntof g e\nX e\nT g\nT e\ntof e g\ntof g e\ntof e g\ntof g e\nX g\nX a\nH g\ntof d g\ntof d a\ntof d b\ntof b d\ntof g a\ntof g b\ntof g d\nX a\nT g\nT* d\nT b\nT a\ntof d g\ntof b g\ntof a g\ntof b d\ntof a b\ntof g a\ntof g b\nX g\nH d\ntof d a\ntof d g\nX g\nT d\nT* a\nT g\ntof a d\ntof g d\ntof g a\ntof d a\ntof d g\nT* d\ntof g d\ntof a d\nX e\nH a\nH g\ntof g e\nX e\nT g\nT e\ntof e g\ntof g e\ntof e g\ntof g e\nX g\nX e\nH g\ntof f g\ntof g c\ntof g f\nX c\nT g\nT* f\nT c\ntof f g\ntof c g\ntof g c\nT* g\ntof f g\nH f\ntof f e\nX e\nT f\nT e\ntof e f\ntof f e\ntof e f\ntof f e\nX f\nX e\nH f\ntof b f\ntof a b\ntof f b\ntof f a\nT* f\nT a\nT* b\ntof a f\ntof b f\ntof b a\ntof f b\nT f\nX g\nH f\ntof f e\nX e\nT f\nT e\ntof e f\ntof f e\ntof e f\ntof f e\nX f\nX e\nH f\ntof f g\ntof f c\nX g\nT f\nT* c\nT g\ntof c f\ntof g f\ntof c g\ntof f g\ntof f c\nT* f\ntof c f\ntof g f\nX c\nH g\ntof g e\nX e\nT g\nT e\ntof e g\ntof g e\ntof e g\ntof g e\nX g\nH g\ntof b a\ntof g a\ntof g b\nT* g\nT b\nT* a\ntof b g\ntof a g\ntof a b\ntof g b\ntof g a\nT g\ntof b g\nH f\ntof d b\ntof f b\ntof f d\nT f\nT* d\nT b\ntof d f\ntof b f\ntof b d\ntof f d\ntof f b\nT* f\ntof f d\nH f\ntof g c\ntof f c\ntof f g\nT* c\nT g\nT* f\ntof g c\ntof f c\ntof f g\ntof c f\nT c\nH c\ntof b d\ntof c d\ntof c b\nT c\nT* b\nT d\ntof b c\ntof d c\ntof d b\ntof c b\ntof c d\nT* c\ntof c b\nH c\ntof g c\ntof g e\ntof f g\ntof c g\ntof c f\ntof c e\nT g\nT* f\nT c\nT e\ntof f g\ntof c g\ntof e g\ntof c f\ntof e c\ntof g e\ntof g c\nP* e\nH c\nX e\nH f\ntof f e\ntof g f\ntof f g\nX e\nT g\nT* f\nT e\ntof f g\ntof e g\ntof f e\ntof g e\ntof g f\nT* g\ntof f g\ntof e g\nH e\ntof b d\ntof d b\ntof b d\ntof d e\ntof e d\ntof d e\ntof e g\ntof g e\ntof e g\ntof g a\ntof a g\ntof g a\ntof f c\ntof c f\ntof f c\n",
	        "H g\ntof a b\nT a\nT g\nT* b\ntof g a\ntof b a\ntof b g\ntof a b\nT* a\nT* b\nT g\ntof g a\ntof g b\ntof b g\ntof a g\nH g\nS b\nH a\ntof a e\nT a\nT* e\ntof e a\ntof a e\ntof e a\ntof a e\nX c\nH a\ntof d a\ntof d c\ntof a c\ntof a d\nX c\nT a\nT* d\nT* c\ntof d a\ntof c a\ntof a c\nT a\ntof c a\ntof d a\nH d\ntof d e\nT d\nT* e\ntof e d\ntof d e\ntof e d\ntof d e\nH d\nH f\ntof f d\ntof g f\ntof f g\ntof d g\ntof d f\nT* d\nT* f\nT g\ntof f d\ntof g d\ntof d f\nT d\ntof g d\ntof f d\nH g\ntof g e\nT g\nT* e\ntof e g\ntof g e\ntof e g\ntof g e\nX a\nH g\ntof c a\ntof g a\ntof g c\nX a\nT g\nT* c\nT* a\ntof c g\ntof a g\ntof a c\ntof g c\ntof g a\nT g\ntof c g\nH c\ntof c e\nT c\nT* e\ntof e c\ntof c e\ntof e c\ntof c e\nH c\ntof d c\ntof d b\ntof f d\ntof c b\ntof c d\ntof c f\nT c\nT* f\nT d\nT* b\ntof f c\ntof d c\ntof b c\ntof d f\ntof b d\ntof c b\ntof c d\nH f\ntof f b\ntof f c\nT f\nT* c\nT* b\ntof c f\ntof b f\ntof c b\ntof f b\ntof f c\nT f\ntof c f\ntof b f\nH b\nH c\ntof c e\nT c\nT* e\ntof e c\ntof c e\ntof e c\ntof c e\nX g\nH c\ntof a c\ntof a g\ntof c g\ntof c a\nX g\nT c\nT* a\nT* g\ntof a c\ntof g c\ntof c g\nT c\ntof g c\ntof a c\nH a\ntof a e\nT a\nT* e\ntof e a\ntof a e\ntof e a\ntof a e\nH a\ntof d a\ntof b d\ntof a d\ntof a b\nT* a\nT b\nT* d\ntof b a\ntof d a\ntof d b\ntof a d\nT a\nH a\ntof a e\nT a\nT* e\ntof e a\ntof a e\ntof e a\ntof a e\nX c\nH a\ntof g c\ntof a c\ntof a g\nX c\nT a\nT* g\nT* c\ntof g a\ntof c a\ntof c g\ntof a g\ntof a c\nT a\ntof g a\nH g\ntof g e\nT g\nT* e\ntof e g\ntof g e\ntof e g\ntof g e\nH g\ntof d b\ntof g b\ntof g d\nT* g\nT d\nT* b\ntof d g\ntof b g\ntof b d\ntof g d\ntof g b\nT g\ntof d g\nH c\ntof f d\ntof c d\ntof c f\nT c\nT* f\nT d\ntof f c\ntof d c\ntof d f\ntof c f\ntof c d\nT* c\ntof c f\nH c\ntof c a\ntof g c\ntof a g\nT* c\nT g\nT* a\ntof g c\ntof a c\ntof a g\ntof c a\nT c\nH c\ntof d f\ntof c f\ntof c d\nT c\nT* d\nT f\ntof d c\ntof f c\ntof f d\ntof c d\ntof c f\nT* c\ntof c d\nH c\ntof c a\ntof g a\ntof a g\nT* c\nT g\nT* a\ntof g c\ntof a c\ntof a g\ntof c g\ntof c a\nT c\ntof g c\nH c\ntof d g\ntof g d\ntof d g\ntof g b\ntof b g\ntof g b\ntof b f\ntof f b\ntof b f\ntof f c\ntof c f\ntof f c\ntof c a\ntof a c\ntof c a\n",
	        "H d\nT d\ntof b d\nT* d\ntof a d\nT d\ntof b d\nT* d\nH d\nH f\ntof f d\nT* d\ntof c d\nT d\ntof f d\nT* d\ntof c d\nT d\nH d\nT d\ntof b d\nT* d\ntof a d\nT d\ntof b d\nT* d\nH d\nT* d\ntof c d\nT d\ntof f d\nT* d\ntof c d\nT d\ntof f d\nH f\nH b\nT b\ntof d b\nT* b\ntof a b\nT b\ntof d b\nT* b\nH b\nH f\ntof f b\nT* b\ntof e b\nT b\ntof f b\nT* b\ntof e b\nT b\nH b\nT b\ntof d b\nT* b\ntof a b\nT b\ntof d b\nT* b\nH b\nT* b\ntof e b\nT b\ntof f b\nT* b\ntof e b\nT b\ntof f b\nH f\n",
	        "H d\ntof a b\nT a\nT d\nT* b\ntof d a\ntof b a\ntof b d\ntof a b\nT* a\nT* b\nT d\ntof d a\ntof d b\ntof b d\ntof a d\nH d\nH c\ntof b a\ntof c a\ntof c b\nT c\nT* b\nT a\ntof b c\ntof a c\ntof a b\ntof c b\ntof c a\nT* c\ntof b c\ntof c b\nH b\nH f\ntof b d\ntof b f\nT* f\nT* d\nT b\ntof d f\ntof b f\ntof f d\nT f\ntof b f\ntof d f\nH b\ntof a c\ntof b c\ntof b a\nT b\nT* a\nT c\ntof a b\ntof c b\ntof c a\ntof b a\ntof b c\nT* b\ntof a b\ntof b a\nH a\ntof f a\ntof d a\ntof d f\ntof a d\ntof a f\nT* f\nT d\nT* a\ntof d f\ntof a f\ntof a d\ntof f d\ntof f a\nT f\ntof d f\nS c\nH b\ntof a d\ntof b d\ntof b a\nT b\nT* a\nT d\ntof a b\ntof d b\ntof d a\ntof b a\ntof b d\nT* b\ntof b a\nH b\ntof f b\ntof e f\ntof b f\ntof b e\nT* b\nT e\nT* f\ntof e b\ntof f b\ntof f e\ntof b f\nT b\nH b\ntof d a\ntof b a\ntof b d\nT b\nT* d\nT a\ntof d b\ntof a b\ntof a d\ntof b d\ntof b a\nT* b\ntof b d\nH b\ntof f b\ntof f c\ntof e f\ntof b c\ntof b f\ntof b e\nT b\nT* e\nT f\nT* c\ntof e b\ntof f b\ntof c b\ntof f e\ntof c f\ntof b c\ntof b f\nH f\nH a\ntof a c\ntof a b\nT a\nT* b\nT* c\ntof b a\ntof c a\ntof b c\ntof a c\ntof a b\nT a\ntof b a\ntof c a\nH c\ntof b a\ntof a b\ntof b a\ntof d c\ntof c d\ntof d c\n",

	        "X d\nH e\nT e\ntof c e\nT* e\nH e\ntof a e\nT e\ntof b e\nT* e\ntof a e\nT e\ntof b e\nT* e\nH e\nT e\ntof c e\nT* e\nH e\nH f\ntof f e\nT* e\ntof d e\nT e\ntof f e\nT* e\ntof d e\nT e\nH e\nT e\ntof c e\nT* e\nH e\nT e\ntof b e\nT* e\ntof a e\nT e\ntof b e\nT* e\ntof a e\nH e\nT e\ntof c e\nT* e\nH e\nT* e\ntof d e\nT e\ntof f e\nT* e\ntof d e\nT e\ntof f e\nH f\nX d\nH a\nT a\ntof d a\nT* a\ntof c a\nT a\ntof d a\nT* a\nH a\nH f\ntof f a\nT* a\ntof e a\nT a\ntof f a\nT* a\ntof e a\nT a\nH a\nT a\ntof d a\nT* a\ntof c a\nT a\ntof d a\nT* a\nH a\nT* a\ntof e a\nT a\ntof f a\nT* a\ntof e a\nT a\ntof f a\nH f\n",
	        "X d\nH e\ntof e c\nX b\nX c\nT e\nT* b\nT c\ntof b e\ntof c e\ntof b c\ntof e c\ntof e b\ntof b c\nX c\nH c\ntof c a\nT a\nT e\nT b\ntof e a\ntof c e\nT* a\nT e\nT* c\ntof e a\ntof c e\nX b\nH c\ntof c b\nX b\nT c\nT b\ntof b c\nX b\nH b\nH f\ntof f e\ntof b d\ntof f b\ntof b f\nX d\nT* b\nT f\nT d\ntof f b\ntof b d\ntof b c\nT b\nT* c\nT* d\ntof e c\ntof d e\ntof f e\ntof f d\ntof e d\ntof b d\ntof b e\ntof b c\nX c\nH f\ntof f c\nX c\nT f\nT c\ntof c f\nX f\nX c\nH c\ntof b a\ntof c a\ntof c b\nT c\nT* b\nT a\ntof b c\ntof a c\ntof a b\ntof c b\ntof c a\nT* c\ntof b c\nH b\ntof b f\nX f\nT b\nT f\ntof f b\nX e\nX f\nH f\ntof f e\ntof d f\ntof f d\nX e\nT* d\nT f\nT* e\ntof f d\ntof e d\ntof e f\ntof d f\ntof d e\nT d\ntof e d\ntof f d\nH e\nX c\nZ c\nS* c\nX c\nS c\nX a\nZ a\nS* a\nX a\nS a\ntof c b\ntof b c\ntof c b\ntof b a\ntof a b\ntof b a\ntof f e\ntof e f\ntof f e\n",
	        "H e\ntof e a\ntof e b\nT e\nT* b\nT* a\ntof b e\ntof a e\ntof a b\ntof e b\ntof e a\nT e\ntof a e\ntof b e\nH b\nH d\ntof a e\ntof d e\ntof d a\nT d\nT* a\nT e\nT b\ntof a d\ntof e d\ntof e a\ntof d a\ntof d e\nT* d\ntof a d\ntof d a\nH a\nH f\ntof a c\ntof f a\ntof a f\nT* a\nT* c\nT f\ntof c a\ntof f a\ntof c b\ntof a b\ntof a c\nT a\nT c\nT* b\ntof c a\ntof f a\ntof b c\ntof c b\nH f\ntof e d\ntof f d\ntof f e\nT f\nT* e\nT d\ntof e f\ntof d f\ntof d e\ntof f e\ntof f d\nT* f\ntof e f\ntof f e\nH e\ntof e a\ntof c e\ntof e c\ntof a b\nT a\nT* b\nT* e\ntof c a\ntof c b\ntof c e\ntof a b\nH g\ntof g f\ntof g d\ntof e c\ntof b c\ntof b e\ntof c e\nT g\nT* d\nT* f\nT c\ntof d g\ntof f g\ntof f d\ntof b c\ntof e c\ntof g d\ntof g f\nT g\ntof f g\ntof d g\nH d\nS g\nZ f\ntof f g\nP* g\ntof f g\nH g\ntof g c\nT g\nT* c\ntof c g\ntof g c\ntof c g\ntof g c\nH g\ntof a g\ntof e a\ntof g a\ntof g e\nT g\nT* e\nT a\ntof e g\ntof a g\ntof a e\ntof g a\nT* g\nH g\ntof g c\nT g\nT* c\ntof c g\ntof g c\ntof c g\ntof g c\nH g\ntof b g\ntof d b\ntof g b\ntof g d\nT* g\nT d\nT* b\ntof d g\ntof b g\ntof b d\ntof g b\nT g\nH g\ntof g c\nT g\nT* c\ntof c g\ntof g c\ntof c g\ntof g c\nH g\ntof a e\ntof g e\ntof g a\nT g\nT* a\nT e\ntof a g\ntof e g\ntof e a\ntof g a\ntof g e\nT* g\ntof a g\nH a\ntof a c\nT a\nT* c\ntof c a\ntof a c\ntof c a\ntof a c\nH a\ntof b a\ntof b f\ntof b d\ntof d b\ntof a f\ntof a d\ntof a b\nT a\nT* b\nT d\nT* f\ntof b a\ntof d a\ntof f a\ntof d b\ntof f d\ntof a f\ntof a d\nH b\ntof b f\ntof b a\nT b\nT* a\nT* f\ntof a b\ntof f b\ntof a f\ntof b f\ntof b a\nT b\ntof a b\ntof f b\nH f\nH a\ntof a c\nT a\nT* c\ntof c a\ntof a c\ntof c a\ntof a c\nH a\ntof g a\ntof e g\ntof a g\ntof a e\nT a\nT* e\nT g\ntof e a\ntof g a\ntof g e\ntof a g\nT* a\nH a\ntof a c\nT a\nT* c\ntof c a\ntof a c\ntof c a\ntof a c\nH a\ntof d a\ntof f d\ntof a d\ntof a f\nT* a\nT f\nT* d\ntof f a\ntof d a\ntof d f\ntof a d\nT a\nH a\ntof a c\nT a\nT* c\ntof c a\ntof a c\ntof c a\ntof a c\nH a\ntof g e\ntof a e\ntof a g\nT a\nT* g\nT e\ntof g a\ntof e a\ntof e g\ntof a g\ntof a e\nT* a\ntof g a\nH g\ntof g c\nT g\nT* c\ntof c g\ntof g c\ntof c g\ntof g c\nH g\ntof d f\ntof g f\ntof g d\nT* g\nT d\nT* f\ntof d g\ntof f g\ntof f d\ntof g d\ntof g f\nT g\ntof d g\nX c\nH d\ntof d e\nT d\nT* e\ntof e d\ntof d e\ntof e d\ntof d e\nX b\nH d\ntof a d\ntof a b\ntof d b\ntof d a\nX b\nT d\nT* a\nT* b\ntof a d\ntof b d\ntof d b\nT d\ntof b d\ntof a d\nH a\ntof a e\nT a\nT* e\ntof e a\ntof a e\ntof e a\ntof a e\nH a\ntof g a\ntof a c\ntof a g\nX c\nT* a\nT g\nT* c\ntof g a\ntof c a\ntof a c\nT a\ntof g a\nX a\nH g\ntof g e\nT g\nT* e\ntof e g\ntof g e\ntof e g\ntof g e\nX d\nH g\ntof b d\ntof g d\ntof g b\nX d\nT g\nT* b\nT* d\ntof b g\ntof d g\ntof d b\ntof g b\ntof g d\nT g\ntof b g\nH b\ntof b e\nT b\nT* e\ntof e b\ntof b e\ntof e b\ntof b e\nX g\nH b\ntof c b\ntof c g\ntof c a\ntof a c\ntof b c\ntof b a\ntof b g\nX g\nX a\nX c\nT b\nT* g\nT* a\nT c\ntof g b\ntof a b\ntof c b\ntof a g\ntof c a\ntof a c\ntof g c\ntof b c\ntof b g\nH a\nH g\ntof g c\ntof g b\nT g\nT* b\nT* c\ntof b g\ntof c g\ntof b c\ntof g c\ntof g b\nT g\ntof b g\ntof c g\nH c\ntof b g\ntof g b\ntof b g\ntof g f\ntof f g\ntof g f\ntof f a\ntof a f\ntof f a\ntof d e\ntof e d\ntof d e\ntof e c\ntof c e\ntof e c\n",
	        "H a\ntof d b\ntof a b\ntof a d\nT a\nT* d\nT b\nT e\ntof d a\ntof b a\ntof b d\ntof a d\ntof a b\nT* a\ntof a d\nH a\nH f\ntof f a\ntof f e\ntof a f\ntof a e\nT* a\nT* e\nT f\ntof e a\ntof f a\ntof a f\ntof a e\nT a\nT e\nT* f\ntof f a\ntof f e\ntof e f\nH a\ntof b d\ntof a d\ntof a b\nT a\nT* b\nT d\ntof b a\ntof d a\ntof d b\ntof a b\ntof a d\nT* a\ntof a b\nH a\ntof e a\ntof e b\ntof f e\ntof a b\ntof a e\ntof a f\nT a\nT* f\nT e\nT* b\ntof f a\ntof e a\ntof b a\ntof e f\ntof b e\ntof a b\ntof a e\nH g\ntof g b\ntof g a\nT g\nT* a\nT* b\ntof a g\ntof b g\ntof a b\ntof g b\ntof g a\nT g\ntof a g\ntof b g\nS g\nS c\nH b\nH a\ntof a f\nT a\nT* f\ntof f a\ntof a f\ntof f a\ntof a f\nX c\nH a\ntof d a\ntof d c\ntof a c\ntof a d\nX c\nT a\nT* d\nT* c\ntof d a\ntof c a\ntof a c\nT a\ntof c a\ntof d a\nH d\ntof d f\nT d\nT* f\ntof f d\ntof d f\ntof f d\ntof d f\nH d\ntof e d\ntof b e\ntof d e\ntof d b\nT* d\nT b\nT* e\ntof b d\ntof e d\ntof e b\ntof d e\nT d\nH d\ntof d f\nT d\nT* f\ntof f d\ntof d f\ntof f d\ntof d f\nX a\nH d\ntof c a\ntof d a\ntof d c\nX a\nT d\nT* c\nT* a\ntof c d\ntof a d\ntof a c\ntof d c\ntof d a\nT d\ntof c d\nH c\ntof c f\nT c\nT* f\ntof f c\ntof c f\ntof f c\ntof c f\nH c\ntof e c\ntof e g\ntof b e\ntof c g\ntof c e\ntof c b\nT c\nT* b\nT e\nT* g\ntof b c\ntof e c\ntof g c\ntof e b\ntof g e\ntof c g\ntof c e\nH b\ntof b g\ntof b c\nT b\nT* c\nT* g\ntof c b\ntof g b\ntof c g\ntof b g\ntof b c\nT b\ntof c b\ntof g b\nH g\nH c\ntof c f\nT c\nT* f\ntof f c\ntof c f\ntof f c\ntof c f\nX d\nH c\ntof a c\ntof a d\ntof c d\ntof c a\nX d\nT c\nT* a\nT* d\ntof a c\ntof d c\ntof c d\nT c\ntof d c\ntof a c\nH a\ntof a f\nT a\nT* f\ntof f a\ntof a f\ntof f a\ntof a f\nH a\ntof e a\ntof g e\ntof a e\ntof a g\nT* a\nT g\nT* e\ntof g a\ntof e a\ntof e g\ntof a e\nT a\nH a\ntof a f\nT a\nT* f\ntof f a\ntof a f\ntof f a\ntof a f\nX c\nH a\ntof d c\ntof a c\ntof a d\nX c\nT a\nT* d\nT* c\ntof d a\ntof c a\ntof c d\ntof a d\ntof a c\nT a\ntof d a\nH d\ntof d f\nT d\nT* f\ntof f d\ntof d f\ntof f d\ntof d f\nH d\ntof e g\ntof d g\ntof d e\nT d\nT e\nT* g\ntof e d\ntof g d\ntof g e\ntof d e\ntof d g\nT d\ntof e d\nH c\ntof c e\ntof c a\nT c\nT* a\nT* e\ntof a c\ntof e c\ntof a e\ntof c e\ntof c a\nT c\ntof a c\ntof e c\ntof a c\nP* c\ntof a c\nH e\ntof f e\ntof d f\ntof f d\ntof e d\ntof e f\nT d\nT* e\nT* f\ntof e d\ntof f d\ntof e f\ntof d f\ntof d e\nT d\ntof e d\ntof f d\nH e\nH f\ntof f a\ntof f c\nT f\nT* c\nT* a\ntof c f\ntof a f\ntof c a\ntof f a\ntof f c\nT f\ntof c f\ntof a f\nH a\ntof c f\ntof f c\ntof c f\ntof f e\ntof e f\ntof f e\ntof e d\ntof d e\ntof e d\ntof d a\ntof a d\ntof d a\n",
	        "H g\nX b\nT g\nT* a\nT* b\ntof a g\ntof b a\ntof g a\ntof g b\nT g\nT b\nT a\ntof b g\ntof a b\ntof g a\ntof g b\nH a\nP* b\nX g\nH g\ntof g e\nT g\nT* e\ntof e g\ntof g e\ntof e g\ntof g e\nH g\ntof c g\ntof d c\ntof g c\ntof g d\nT g\nT* d\nT c\ntof d g\ntof c g\ntof c d\ntof g c\nT* g\nH g\ntof g e\nT g\nT* e\ntof e g\ntof g e\ntof e g\ntof g e\nH g\nH f\ntof f g\ntof a f\ntof f a\ntof g a\ntof g f\nT* g\nT* f\nT a\ntof f g\ntof a g\ntof g f\nT g\ntof a g\ntof f g\nH a\ntof a e\nT a\nT* e\ntof e a\ntof a e\ntof e a\ntof a e\nH a\ntof c d\ntof a d\ntof a c\nT a\nT* c\nT d\ntof c a\ntof d a\ntof d c\ntof a c\ntof a d\nT* a\ntof c a\nH c\ntof c e\nT c\nT* e\ntof e c\ntof c e\ntof e c\ntof c e\nX b\nH c\ntof g c\ntof g b\ntof f g\ntof c b\ntof c g\ntof c f\nX b\nT c\nT* f\nT g\nT b\ntof f c\ntof g c\ntof b c\ntof g f\ntof b g\ntof c b\ntof c g\nX c\nH f\ntof f b\ntof f c\nX c\nT f\nT* b\nT c\ntof b f\ntof c f\ntof c b\ntof f b\ntof f c\nT* f\ntof c f\ntof b f\nH b\nH c\ntof c e\nT c\nT* e\ntof e c\ntof c e\ntof e c\ntof c e\nH c\ntof a c\ntof d a\ntof c a\ntof c d\nT c\nT* d\nT a\ntof d c\ntof a c\ntof a d\ntof c a\nT* c\nH c\ntof c e\nT c\nT* e\ntof e c\ntof c e\ntof e c\ntof c e\nH c\ntof g c\ntof b g\ntof c g\ntof c b\nT* c\nT b\nT* g\ntof b c\ntof g c\ntof g b\ntof c g\nT c\nH c\ntof c e\nT c\nT* e\ntof e c\ntof c e\ntof e c\ntof c e\nH c\ntof a d\ntof c d\ntof c a\nT c\nT* a\nT d\ntof a c\ntof d c\ntof d a\ntof c a\ntof c d\nT* c\ntof a c\nH a\ntof a e\nT a\nT* e\ntof e a\ntof a e\ntof e a\ntof a e\nH a\ntof g b\ntof a b\ntof a g\nT* a\nT g\nT* b\ntof g a\ntof b a\ntof b g\ntof a g\ntof a b\nT a\ntof g a\nX g\nH g\ntof d c\ntof g c\ntof g d\nT g\nT* d\nT c\ntof d g\ntof c g\ntof c d\ntof g d\ntof g c\nT* g\ntof g d\nH g\ntof a g\ntof a e\ntof e a\ntof g e\ntof g a\nT* g\nT a\nT* e\ntof a g\ntof e g\ntof e a\ntof g e\nT g\nH g\ntof c d\ntof g d\ntof g c\nT g\nT* c\nT d\ntof c g\ntof d g\ntof d c\ntof g c\ntof g d\nT* g\ntof g c\nH g\ntof a e\ntof e a\ntof g e\ntof g a\nT* g\nT a\nT* e\ntof a g\ntof e g\ntof e a\ntof g a\ntof g e\nT g\ntof a g\nH d\ntof f a\ntof d a\ntof d f\nT d\nT* f\nT a\ntof f d\ntof a d\ntof a f\ntof d f\ntof d a\nT* d\ntof d f\nX c\nH d\ntof d c\ntof g c\ntof d g\nX c\nT* g\nT d\nT* c\ntof d g\ntof c g\ntof g c\nT g\ntof d g\nH d\ntof a f\ntof d f\ntof d a\nT d\nT* a\nT f\ntof a d\ntof f d\ntof f a\ntof d a\ntof d f\nT* d\ntof d a\nX g\nH d\ntof d g\ntof c d\ntof d c\nX g\nT* c\nT d\nT* g\ntof d c\ntof g c\ntof g d\ntof c d\ntof c g\nT c\ntof g c\ntof d c\nH g\ntof f g\ntof g f\ntof f g\ntof g b\ntof b g\ntof g b\n",
	        "H f\ntof b a\ntof a d\nT b\nT f\nT* a\nT* d\ntof f b\ntof a f\ntof a d\ntof b a\nT* b\nT* a\nT f\ntof f b\ntof f a\ntof a f\ntof b f\ntof d a\nH g\ntof g b\ntof g a\nT g\nT* a\nT* b\ntof a g\ntof b g\ntof b a\ntof g a\ntof g b\nT g\ntof b g\ntof a g\nH a\nS g\nS b\nH g\ntof g e\nT g\nT* e\ntof e g\ntof g e\ntof e g\ntof g e\nH g\ntof c g\ntof d c\ntof g c\ntof g d\nT g\nT* d\nT c\ntof d g\ntof c g\ntof c d\ntof g c\nT* g\nH g\ntof g e\nT g\nT* e\ntof e g\ntof g e\ntof e g\ntof g e\nH g\ntof f g\ntof a f\ntof g f\ntof g a\nT* g\nT a\nT* f\ntof a g\ntof f g\ntof f a\ntof g f\nT g\nH g\ntof g e\nT g\nT* e\ntof e g\ntof g e\ntof e g\ntof g e\nH g\ntof c d\ntof g d\ntof g c\nT g\nT* c\nT d\ntof c g\ntof d g\ntof d c\ntof g c\ntof g d\nT* g\ntof c g\nH c\ntof c e\nT c\nT* e\ntof e c\ntof c e\ntof e c\ntof c e\nH c\ntof f c\ntof f b\ntof a f\ntof c b\ntof c f\ntof c a\nT c\nT* a\nT f\nT* b\ntof a c\ntof f c\ntof b c\ntof f a\ntof b f\ntof c b\ntof c f\nH a\ntof a b\ntof a c\nT a\nT* c\nT* b\ntof c a\ntof b a\ntof c b\ntof a b\ntof a c\nT a\ntof c a\ntof b a\nH b\nH c\ntof c e\nT c\nT* e\ntof e c\ntof c e\ntof e c\ntof c e\nH c\ntof g c\ntof d g\ntof c g\ntof c d\nT c\nT* d\nT g\ntof d c\ntof g c\ntof g d\ntof c g\nT* c\nH c\ntof c e\nT c\nT* e\ntof e c\ntof c e\ntof e c\ntof c e\nH c\ntof f c\ntof b f\ntof c f\ntof c b\nT* c\nT b\nT* f\ntof b c\ntof f c\ntof f b\ntof c f\nT c\nH c\ntof c e\nT c\nT* e\ntof e c\ntof c e\ntof e c\ntof c e\nH c\ntof g d\ntof c d\ntof c g\nT c\nT* g\nT d\ntof g c\ntof d c\ntof d g\ntof c g\ntof c d\nT* c\ntof g c\nH g\ntof g e\nT g\nT* e\ntof e g\ntof g e\ntof e g\ntof g e\nH g\ntof f b\ntof g b\ntof g f\nT* g\nT f\nT* b\ntof f g\ntof b g\ntof b f\ntof g f\ntof g b\nT g\ntof f g\nH d\ntof a f\ntof d f\ntof d a\nT d\nT* a\nT f\ntof a d\ntof f d\ntof f a\ntof d a\ntof d f\nT* d\ntof d a\nH d\ntof d c\ntof g d\ntof c g\nT* d\nT g\nT* c\ntof g d\ntof c d\ntof c g\ntof d c\nT d\nH d\ntof f a\ntof d a\ntof d f\nT d\nT* f\nT a\ntof f d\ntof a d\ntof a f\ntof d f\ntof d a\nT* d\ntof d f\nH d\ntof d c\ntof g c\ntof c g\nT* d\nT g\nT* c\ntof g d\ntof c d\ntof c g\ntof d g\ntof d c\nT d\ntof g d\nH d\ntof g a\ntof f d\ntof d f\ntof f d\ntof d g\ntof g d\ntof d g\ntof g b\ntof b g\ntof g b\ntof b a\ntof a b\ntof b a\n",
	        "H g\nX a\nX g\nS b\nT b\nT* a\nT g\ntof a b\ntof g b\ntof g a\ntof b g\nT a\nT* b\nT g\ntof g a\ntof g b\ntof a g\ntof a b\nT* a\ntof b a\ntof g a\nX a\nH g\nH a\ntof a e\nT a\nT* e\ntof e a\ntof a e\ntof e a\ntof a e\nH a\ntof c a\ntof d c\ntof a c\ntof a d\nT a\nT* d\nT c\ntof d a\ntof c a\ntof c d\ntof a c\nT* a\nH a\ntof a e\nT a\nT* e\ntof e a\ntof a e\ntof e a\ntof a e\nH a\nH f\ntof f a\ntof g f\ntof f g\ntof a g\ntof a f\nT* a\nT* f\nT g\ntof f a\ntof g a\ntof f b\ntof a b\ntof a f\nT a\nT f\nT* b\ntof f a\ntof g a\ntof b f\ntof f b\nH g\ntof g e\nT g\nT* e\ntof e g\ntof g e\ntof e g\ntof g e\nH g\ntof c d\ntof g d\ntof g c\nT g\nT* c\nT d\ntof c g\ntof d g\ntof d c\ntof g c\ntof g d\nT* g\ntof c g\nH c\ntof c e\nT c\nT* e\ntof e c\ntof c e\ntof e c\ntof c e\nH c\ntof a c\ntof a f\ntof b a\ntof c f\ntof c a\ntof c b\nT c\nT* b\nT a\nT* f\ntof b c\ntof a c\ntof f c\ntof a b\ntof f a\ntof c f\ntof c a\nH b\ntof b f\ntof b c\nT b\nT* c\nT* f\ntof c b\ntof f b\ntof c f\ntof b f\ntof b c\nT b\ntof c b\ntof f b\nH f\nH c\ntof c e\nT c\nT* e\ntof e c\ntof c e\ntof e c\ntof c e\nH c\ntof g c\ntof d g\ntof c g\ntof c d\nT c\nT* d\nT g\ntof d c\ntof g c\ntof g d\ntof c g\nT* c\nH c\ntof c e\nT c\nT* e\ntof e c\ntof c e\ntof e c\ntof c e\nH c\ntof a c\ntof f a\ntof c a\ntof c f\nT* c\nT f\nT* a\ntof f c\ntof a c\ntof a f\ntof c a\nT c\nH c\ntof c e\nT c\nT* e\ntof e c\ntof c e\ntof e c\ntof c e\nH c\ntof g d\ntof c d\ntof c g\nT c\nT* g\nT d\ntof g c\ntof d c\ntof d g\ntof c g\ntof c d\nT* c\ntof g c\nH g\ntof g e\nT g\nT* e\ntof e g\ntof g e\ntof e g\ntof g e\nH g\ntof a g\ntof a b\ntof g a\ntof g b\ntof g f\nT g\nT* f\nT b\ntof f g\ntof b g\ntof b f\ntof a f\ntof g f\ntof g b\nT* g\ntof b g\ntof a g\ntof a b\nP* a\nS b\ntof a b\nX a\nH b\ntof g f\ntof f g\ntof g f\ntof f b\ntof b f\ntof f b\n",
	        "H f\nT a\nT b\nT f\ntof b a\ntof f b\ntof a f\nT* b\ntof a b\nT* a\nT* b\nT f\ntof f b\ntof a f\ntof b a\nH f\n",
	        "H d\nT d\ntof b d\nT* d\ntof a d\nT d\ntof b d\nT* d\nH d\nH f\ntof f d\nT* d\ntof c d\nT d\ntof f d\nT* d\ntof c d\nT d\nH d\nT d\ntof b d\nT* d\ntof a d\nT d\ntof b d\nT* d\nH d\nT* d\ntof c d\nT d\ntof f d\nT* d\ntof c d\nT d\ntof f d\nH f\nX a\nH b\nT b\ntof d b\nT* b\ntof a b\nT b\ntof d b\nT* b\nH b\nH f\ntof f b\nT* b\ntof e b\nT b\ntof f b\nT* b\ntof e b\nT b\nH b\nT b\ntof d b\nT* b\ntof a b\nT b\ntof d b\nT* b\nH b\nT* b\ntof e b\nT b\ntof f b\nT* b\ntof e b\nT b\ntof f b\nH f\nX a\n",
	        "H d\ntof b a\ntof d a\ntof d b\nT d\nT* b\nT a\nT c\ntof b d\ntof a d\ntof a b\ntof d b\ntof d a\nT* d\ntof b d\ntof d b\nH b\nH f\ntof b c\ntof f b\ntof b f\nT* b\nT* c\nT f\ntof c b\ntof f b\ntof b f\ntof b c\nT b\nT c\nT* f\ntof f b\ntof f c\ntof c f\nH b\ntof a d\ntof b d\ntof b a\nT b\nT* a\nT d\ntof a b\ntof d b\ntof d a\ntof b a\ntof b d\nT* b\ntof a b\ntof b a\nH a\nH f\ntof b d\ntof d b\ntof b d\ntof d a\ntof a d\ntof d a\n",

	        "X e\nH c\nT c\ntof d c\nT* c\nH c\ntof a c\nT c\ntof b c\nT* c\ntof a c\nT c\ntof b c\nT* c\nH c\nT c\ntof d c\nT* c\nH c\nH f\ntof f c\nT* c\ntof e c\nT c\ntof f c\nT* c\ntof e c\nT c\nH c\nT c\ntof d c\nT* c\nH c\nT c\ntof b c\nT* c\ntof a c\nT c\ntof b c\nT* c\ntof a c\nH c\nT c\ntof d c\nT* c\nH c\nT* c\ntof e c\nT c\ntof f c\nT* c\ntof e c\nT c\ntof f c\nH f\nX e\nX d\nH e\nT e\ntof c e\nT* e\nH e\ntof a e\nT e\ntof b e\nT* e\ntof a e\nT e\ntof b e\nT* e\nH e\nT e\ntof c e\nT* e\nH e\nH f\ntof f e\nT* e\ntof d e\nT e\ntof f e\nT* e\ntof d e\nT e\nH e\nT e\ntof c e\nT* e\nH e\nT e\ntof b e\nT* e\ntof a e\nT e\ntof b e\nT* e\ntof a e\nH e\nT e\ntof c e\nT* e\nH e\nT* e\ntof d e\nT e\ntof f e\nT* e\ntof d e\nT e\ntof f e\nH f\nX d\nX a\nX b\nH c\nT c\ntof d c\nT* c\nH c\ntof a c\nT c\ntof b c\nT* c\ntof a c\nT c\ntof b c\nT* c\nH c\nT c\ntof d c\nT* c\nH c\nH f\ntof f c\nT* c\ntof e c\nT c\ntof f c\nT* c\ntof e c\nT c\nH c\nT c\ntof d c\nT* c\nH c\nT c\ntof b c\nT* c\ntof a c\nT c\ntof b c\nT* c\ntof a c\nH c\nT c\ntof d c\nT* c\nH c\nT* c\ntof e c\nT c\ntof f c\nT* c\ntof e c\nT c\ntof f c\nH f\nX a\nX b\nH a\nT a\ntof d a\nT* a\ntof c a\nT a\ntof d a\nT* a\nH a\nH f\ntof f a\nT* a\ntof e a\nT a\ntof f a\nT* a\ntof e a\nT a\nH a\nT a\ntof d a\nT* a\ntof c a\nT a\ntof d a\nT* a\nH a\nT* a\ntof e a\nT a\ntof f a\nT* a\ntof e a\nT a\ntof f a\nH f\n"
	       }};
	// clang-format on
};

} /* namespace tweedledum */
