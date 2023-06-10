/**
  * FileName: config.h
  * Author: Zeng Zhenxiang
  * Version: 1.0
  * Date:
  * Description:
***/

#include "config.h"

namespace poisson_reconstruction {

	Config::Config(const std::string& input, const std::string& output, uint8_t max_d, bool bin) :
		input_filename(input),
		output_filename(output),
		binary(bin),
		max_depth(max_d),
		kernel_depth(max_depth-3),
		average_samples_per_node(1.0)
	{ }

	std::ostream& operator<<(std::ostream& os, const Config& config) {
		std::cout << "Poisson reconstruction configuration: " << std::endl;
		std::cout << '\t' << "input  file : " << config.input_filename << std::endl;
		std::cout << '\t' << "output file : " << config.output_filename << std::endl;
		std::cout << '\t' << "binary      : ";
		if (config.binary) std::cout << "true" << std::endl;
		else std::cout << "false" << std::endl;
		std::cout << '\t' << "max depth   : " << (int)config.max_depth << std::endl;
		std::cout << '\t' << "kernel depth: " << (int)config.kernel_depth << std::endl;
		std::cout << '\t' << "scale ratio : " << config.get_sacle_ratio() << std::endl;
		std::cout << '\t' << "average samples per node : " << config.average_samples_per_node << std::endl;
		return os;
	}

}