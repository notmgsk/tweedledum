/*--------------------------------------------------------------------------------------------------
| This file is distributed under the MIT License.
| See accompanying file /LICENSE for details.
| Author(s): Bruno Schmitt
*-------------------------------------------------------------------------------------------------*/
#pragma once

#include "../gates/gate_set.hpp"
#include "../networks/qubit.hpp"

#include <lorina/pla.hpp>
#include <utility>

namespace tweedledum {

/*! \brief Lorina reader callback for PLA files.
 *
 * **Required network functions:**
 * - `create_pi`
 * - `create_po`
 * - `create_not`
 * - `create_nary_and`
 * - `create_nary_or`
 * - `create_nary_xor`
 *
   \verbatim embed:rst
   Example
   .. code-block:: c++
      aig_network aig;
      lorina::read_pla( "file.pla", pla_reader( aig ) );
      mig_network mig;
      lorina::read_pla( "file.pla", pla_reader( mig ) );
   \endverbatim
 */
template<typename Network>
class pla_reader : public lorina::pla_reader {
public:
	explicit pla_reader(Network& ntk)
	    : network_(ntk)
	{}

	void on_number_of_inputs(std::size_t number_of_inputs) const override
	{
		input_qubits_.resize(number_of_inputs);
		std::generate(input_qubits_.begin(), input_qubits_.end(), [this]() { return network_.add_qubit(); });
	}

	void on_number_of_outputs(std::size_t number_of_outputs) const override
	{
		output_qubits_.resize(number_of_outputs);
		std::generate(output_qubits_.begin(), output_qubits_.end(), [this]() { return network_.add_qubit(); });
	}

	bool on_keyword(const std::string& keyword, const std::string& value) const override
	{
		if (keyword == "type" && value == "esop") {
			is_esop_ = true;
			return true;
		}
		return false;
	}

	void on_end() const override
	{
		for (auto gate : gates_) {
			network_.add_gate(gate_set::mcx, gate.first, gate.second);
		}
	}

	void on_term(const std::string& term, const std::string& out) const override
	{
		std::vector<qubit_id> controls;
		for (auto i = 0u; i < term.size(); ++i) {
			switch (term[i]) {
			default:
				std::cerr << "[w] unknown character '" << term[i]
				          << "' in PLA input term, treat as don't care\n";
			case '-':
				break;

			case '0':
				controls.push_back(!input_qubits_[i]);
				break;

			case '1':
				controls.push_back(input_qubits_[i]);
				break;
			}
		}

		std::vector<qubit_id> targets;
		for (auto i = 0u; i < out.size(); ++i) {
			switch (out[i]) {
			default:
				std::cerr << "[w] unknown character '" << out[i]
				          << "' in PLA output term, treat is 0\n";
			case '0':
				break;

			case '1':
				targets.push_back(output_qubits_[i]);
				break;
			}
		}
		gates_.emplace_back(controls, targets);
	}

private:
	Network& network_;
	mutable std::vector<qubit_id> input_qubits_;
	mutable std::vector<qubit_id> output_qubits_;
	mutable std::vector<std::pair<std::vector<qubit_id>, std::vector<qubit_id>>> gates_;
	mutable bool is_esop_ = false;
};

} /* namespace tweedledum */