/**
  * FileName: config.h
  * Author: Zeng Zhenxiang
  * Version: 1.0
  * Date:
  * Description:
***/

#ifndef POISSON_RECONSTRUCTION_CONFIG_H
#define POISSON_RECONSTRUCTION_CONFIG_H

#include <string>
#include <iostream>

namespace poisson_reconstruction {

	// Poisson surface reconstruction configuration option
	class Config {
	public:
		const std::string input_filename;  // contains oriented points data
		const std::string output_filename;
		bool binary;			// whether the data in input_file is stored in binary format
		uint8_t max_depth;		// max depth for the octree
		uint8_t kernel_depth;	// a kernel density estimator need 
		double average_samples_per_node;
		Config(const std::string& input, const std::string& output, uint8_t max_d, bool bin);
		Config(const Config& c) = default;

		double get_sacle_ratio() const { return scale_ratio; }

	private:
		// here the default setting of the ratio between the [0, 1]^3 cube and the sample cube is set to: 1.25
		double scale_ratio{ 1.25 };
	};

	std::ostream& operator<<(std::ostream& os, const Config& config);
}

#endif // !POISSON_RECONSTRUCTION_CONFIG_H
