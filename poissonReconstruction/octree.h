/**
  * FileName: octree.h
  * Author: Zeng Zhenxiang
  * Version: 1.0
  * Date: 
  * Description: A brief implementation of octree.
***/

#ifndef POISSON_RECONSTRUCTION_OCTREE_H
#define POISSON_RECONSTRUCTION_OCTREE_H

#include <iostream>
#include <cstring>
#include <cmath>
#include <initializer_list>
#include <fstream>
#include <stdexcept>
#include <limits>
#include <vector>
#include <utility>
#include <unordered_map>
#include <algorithm>
#include "node.h"
#include "octdata.h"
#include "gaussian.h"
#include "config.h"
#include "matrix.h"
#include "marching_cube.h"
#include "mesh.h"
#include "magic_cube.h"

namespace poisson_reconstruction {

	struct Samples_property {
		double center[3];
		double max_stretch;
		int valid_samples_num{0};
		double mean_sample_density{ 0 };
	};

	//---------------------------------------------------------------------------------------
	// Some file manipulations 
	//---------------------------------------------------------------------------------------

	inline void open_file(std::ifstream& file, const std::string& file_name, bool binary) {
		if (binary) {
			file.open(file_name, std::ifstream::in | std::ifstream::binary);
			if (!file) {
				throw std::runtime_error("can't open input file in binary model: " + file_name);
			}
		}
		else {
			file.open(file_name, std::ifstream::in);
			if (!file) {
				throw std::runtime_error("can't open input file: " + file_name);
			}
		}
	}

	// read data: each row has 6 coordinate values
	inline void read_one_line_from_file(std::ifstream& file, float coor[6], bool binary) {
		if (binary) { file.read((char*)coor, sizeof(float) * 6); }
		else {
			for (uint8_t i = 0; i < 6; ++i) file >> coor[i];
		}
	}

	inline void seek_to_begin(std::ifstream& file) { file.clear(); file.seekg(0, file.beg); }

	inline void close_file(std::ifstream& file) { file.close(); }

	//---------------------------------------------------------------------------------------
	// class Octree
	//---------------------------------------------------------------------------------------

	// An octree that uses a 3-degree box filter convolution and has data accuracy of double
	template <uint16_t N = 3, typename T = double>
	class Octree {	
		static const double epsilon;
	public:
		Octree(const Config& c): config(c), octdata(c.max_depth), root(new Node){ }
		~Octree() = default;

		const Samples_property& get_sample_property() { return samples_property; }
		Octdata<N, T>& get_octdata() { return octdata; }
		const Mesh_data& get_mesh_data() { return mesh_data; }
		const Node* get_root() const { return root; }

		/* 1-th: adaptive construct octree according to sample set
		 */
		int adaptive_construct_octree_from_samples();

		/* 2-th step of poisson reconstruction
		 * The following functions are used to estimate local sampling density
		 * i.e a kernel density estimator 
		 */
		void _set_weight_alpha(float coor[3], Magic_cube& curr_magic, Symmetric_poly<N - 1, T>& gaussian_appro);
		void set_node_weight_alpha();
		// Obtain the local sampling density of position pos at kernel depth
		double get_local_sampling_density(const Position<double>& pos);

		/* 3-th step of poisson reconstruction
		 * The following functions are used to compute the vector field
		 */
		// the average sampling density over all of the samples
		double get_average_sampling_density_W();
		uint8_t get_adaptive_depth_for_sample(const Position<double>& pos, double local_density);
		// i.e for node, set 1/W^D(s.p) * alpha_o,s * s.N
		void set_node_vector_field_factor();

		/* 4-th step of poisson reconstruction
		 * Because we use adaptive depth (in the 2-th step) for the vector field contribution of each sample point
		 * so, there may be some empty leaf node, i.e vector_field_factor[3] = {0.0, 0.0, 0.0}
		 * we clip these nodes
		 */
		void clip_empty_leaf_nodes() { Node::clip_tree(root); }
		 
		/* 5-th step of poisson reconstruction
		 * Sets the projection of the divergence of the vector field onto each node of the octree
		 */
		struct Projection_of_divergence_node2node {
			Octree<N, T>* octree;
			const std::vector<std::vector<T>>* ff_dot;
			const std::vector<std::vector<T>>* df_dot;
			// calculate the projection of node onto curr_proc_node
			void operator()(Node* node, Node* curr_proc_node);
		};
		void set_projection_of_divergence_vector_field();

		/* 6-th step of poisson reconstruction
		 * Solving linear system, i.e solve large linear equations
		 */
		struct Projection_of_indicator_node2node {
			Octree<N, T>* octree;
			const std::vector<std::vector<T>>* ff_dot;
			const std::vector<std::vector<T>>* dd_dot;
			std::vector<Entry<T>> row_elements;
			void operator()(Node* node, Node* curr_proc_node);
		};
		void solve_linear_system();

		/* 7-th step of poisson reconstruction
		 * Selecting an isovalue 
		 */
		struct calculate_indicator_func_value {
			double func_value{ 0.0 };
			Position<double> center;
			double dx, dy, dz;
			void operator()(Node* node, const Position<double>& point);
		};
		double select_isovalue();
		struct calculate_center_indicator_value {
			const std::vector<std::vector<T>>* value_table;
			double func_value{ 0.0 };
			uint16_t center_index[3];
			void operator()(Node* node, const Position<double>& point);
		};
		double select_isovalue_scheme2();

		/* 8-th step of poisson reconstruction
		 * Extraction isosurface 
		 */
		struct calculate_corner_indicator_value {
			const std::vector<std::vector<T>>* value_table;
			double func_value{ 0.0 };
			uint16_t point_index[3];
			void operator()(Node* node, const Position<double>& point);
		};
		void set_cube_corners_value(double isovalue);
		// these functions are used to refine coarse nodes
		int edge_root_count(Node* curr_node, uint8_t edge, double isovalue);
		int face_root_count(Node* curr_node, uint8_t face, double isovalue);
		void divide_node(Node* curr_node, double isovalue);
		void coarse_node_refinement(Node* curr_node, double isovalue);
		// these functions are used to set the intersection
		bool search_edge_root_key(Node* curr_node, uint8_t edge, double isovalue, long long& key);
		bool calculate_edge_root(Node* curr_node, uint8_t edge, double isovalue, Position<float>& intersect_pos);
		void set_intersection(double isovalue);
		// get triangles
		void get_triangles(double isovalue);

		/* 9-th step(final step) of poisson reconstruction
		 * Output ply data files
		 */
		// void ouput_ply_file();

	private:
		Node* root{nullptr};
		// Configuration of Poisson reconstruction
		Config config;
		// Properties of the sample set
		Samples_property samples_property;
		// octree function data
		Octdata<N, T> octdata;
		// mesh data
		Mesh_data mesh_data;
		// mapping: edge key -> index of container mesh_data.intersection_points
		std::unordered_map<long long, int> edge_key_to_index;
	};

	template<uint16_t N, typename T>
	const double Octree<N, T>::epsilon = 1e-6;

	//---------------------------------------------------------------------------------------
	// Octree member functions
	//---------------------------------------------------------------------------------------

	template<uint16_t N, typename T>
	inline int Octree<N, T>::adaptive_construct_octree_from_samples() {
		int total_samples_num = 0;
		// open input file
		std::ifstream file;
		open_file(file, config.input_filename, config.binary);

		// Scale and shift the sample in order to fit into a [0, 1]^3 cube
		double scale_ratio = config.get_sacle_ratio();  // the default setting is: 1.25
		
		// Record the maximum and minimum coordinates of the sample in the x y z direction
		double xyz_min[3];
		xyz_min[0] = xyz_min[1] = xyz_min[2] = std::numeric_limits<double>::max();
		double xyz_max[3];
		xyz_max[0] = xyz_max[1] = xyz_max[2] = std::numeric_limits<double>::min();

		double xyz_stretch[3], max_stretch = 0;
		double samples_center[3];
		float coor[6];  // must be float
		int nums = 0;
		while (!file.eof()) {
			// read data: Each row has 6 coordinate values
			if (config.binary) { file.read((char*)coor, sizeof(float) * 6); }
			else {
				for (uint8_t i = 0; i < 6; ++i) file >> coor[i];
			}
			for (uint8_t i = 0; i < 3; ++i) {
				if (coor[i] < xyz_min[i]) xyz_min[i] = coor[i];
				if (coor[i] > xyz_max[i]) xyz_max[i] = coor[i];
				//std::cout << coor[i] << " ";
			}
			++nums;
			//std::cout << '\n';
		}
		std::cout << "nums: " << nums << '\n';

		for (uint8_t i = 0; i < 3; ++i) {
			xyz_stretch[i] = xyz_max[i] - xyz_min[i];
			max_stretch = std::max(max_stretch, xyz_stretch[i]);
			samples_center[i] = (xyz_min[i] + xyz_max[i]) / 2.0;
		}

		std::cout << "xyz_stretch: \n";
		for (uint8_t i = 0; i < 3; ++i) {
			std::cout << xyz_stretch[i] << " ";
		}
		std::cout << '\n';
		std::cout << "samples_center: \n";
		for (uint8_t i = 0; i < 3; ++i) {
			std::cout << samples_center[i] << " ";
		}
		std::cout << '\n';
		std::cout << "max_stretch: " << max_stretch << '\n';

		for (uint8_t i = 0; i < 3; ++i) {
			samples_property.center[i] = samples_center[i];
		}
		samples_property.max_stretch = max_stretch;

		// return to file start position
		seek_to_begin(file);

		// Sample by sample processing
		while (!file.eof()) {
			// read data: each row has 6 coordinate values
			read_one_line_from_file(file, coor, config.binary);

			// transform to [0, 1]^3 range
			for (uint8_t dir = 0; dir < 3; ++dir) {
				coor[dir] = (coor[dir] - samples_center[dir]) / (max_stretch * scale_ratio) + 0.5;
				//std::cout << coor[dir] << " ";
			}
			//std::cout << '\n';

			bool satisfy = true;
			for (uint8_t i = 0; i < 3; ++i) {
				if (coor[i] < 0 || coor[i] > 1) { satisfy = false; break; }
			}
			if (!satisfy) continue;
			

			++total_samples_num;

			// Octree is built adaptively according to the current location of sample points
			// The maximum depth of self-adaptation is the kernel depth set in config
			uint8_t kernel_depth = config.kernel_depth;

			uint8_t curr_depth = 0;
			double curr_center[3] = { 0.5, 0.5, 0.5 };
			double curr_width = 1.0;
			uint8_t node_number = 0;
			Node* curr_node = root;
			while (curr_depth < kernel_depth) {
				if (curr_node->children == nullptr) curr_node->init_eight_children();

				// Calculate the child node's number, center and width of the sample point coor[3]
				curr_width /= 2.0;
				node_number = 0;
				for (uint8_t dir = 0; dir < 3; ++dir) {
					if (coor[dir] < curr_center[dir]) {
						curr_center[dir] -= curr_width / 2.0;
					}
					else {
						curr_center[dir] += curr_width / 2.0;
						node_number |= (1 << dir);
					}
				}
				curr_node = &curr_node->children[node_number];
				++curr_depth;
			}
		}

		close_file(file);
		samples_property.valid_samples_num = total_samples_num;
		return total_samples_num;
	}

	template<uint16_t N, typename T>
	inline void Octree<N, T>::_set_weight_alpha(float coor[3], Magic_cube& curr_magic, Symmetric_poly<N - 1, T>& gaussian_appro) {
		// set weight alpha
		double dx, dy, dz;
		Node* temp;
		Position<double> center;
		for (uint8_t x = 0; x < 3; ++x)
			for (uint8_t y = 0; y < 3; ++y)
				for (uint8_t z = 0; z < 3; ++z) {
					temp = curr_magic(Position<uint8_t>(x, y, z));
					if (temp != nullptr) {
						center = temp->get_center();
						dx = coor[0] - center.x;
						dy = coor[1] - center.y;
						dz = coor[2] - center.z;
						temp->weight_alpha += gaussian_appro(dx) * gaussian_appro(dy) * gaussian_appro(dz);
					}
				}
	}

	template<uint16_t N, typename T>
	inline void Octree<N, T>::set_node_weight_alpha() {

		// Here the weight alpha is set using the n-degree convolution of the box function to simulate Gaussian function
		Symmetric_poly<N-1, T> gaussian_appro = Gaussian_approximation<N, T>(0.5);

		// open input file
		std::ifstream file;
		open_file(file, config.input_filename, config.binary);
		
		// sample by sample processing
		// set node weight alpha layer by layer: from 0 to kernel_depth
		static const uint8_t kernel_depth = config.kernel_depth;
		static const double scale_ratio = config.get_sacle_ratio();
		float coor[6];  // must be float
		while (!file.eof()) {

			read_one_line_from_file(file, coor, config.binary);
			
			// coordinate transform: transform to [0, 1]^3 range
			for (uint8_t dir = 0; dir < 3; ++dir) {
				coor[dir] = (coor[dir] - samples_property.center[dir]) / (samples_property.max_stretch * scale_ratio) + 0.5;
			}

			bool satisfy = true;
			for (uint8_t i = 0; i < 3; ++i) {
				if (coor[i] < 0 || coor[i] > 1) { satisfy = false; break; }
			}
			if (!satisfy) continue;

			uint8_t curr_depth = 0;
			double curr_center[3] = { 0.5, 0.5, 0.5 };
			double curr_width = 1.0;
			Node* curr_node = root;
			Magic_cube curr_magic;
			curr_magic.construct(curr_node);
			// current node's magic cube structure
			while (curr_depth <= kernel_depth) {

				auto gaussian = gaussian_appro.scale(curr_width);

				// set weight alpha
				_set_weight_alpha(coor, curr_magic, gaussian);
				
				// Prepare for the next layer
				++curr_depth;
				//std::cout << "curr_depth: " << (int)curr_depth << std::endl;
				if (curr_depth > kernel_depth) break;

				curr_width /= 2.0;
				uint8_t node_number = 0;
				for (uint8_t dir = 0; dir < 3; ++dir) {
					if (coor[dir] < curr_center[dir]) {
						curr_center[dir] -= curr_width / 2.0;
					}
					else {
						curr_center[dir] += curr_width / 2.0;
						node_number |= (1 << dir);
					}
				}
				curr_node = &curr_node->children[node_number];
				curr_magic.construct(curr_node);
			}
		}
		
		close_file(file);
	}

	// The position pos here is the converted coordinate
	template<uint16_t N, typename T>
	inline double poisson_reconstruction::Octree<N, T>::get_local_sampling_density(const Position<double>& pos)
	{
		double ans_weight = 0.0;

		double coor[3] = { pos.x, pos.y, pos.z };

		static uint8_t kernel_depth = config.kernel_depth;

		uint8_t curr_depth = 0;
		double curr_center[3] = { 0.5, 0.5, 0.5 };
		double curr_width = 1.0;
		Node* curr_node = root;
		uint8_t node_number = 0;
		// Find the node at kernel depth where the pos is located
		while (curr_depth < kernel_depth && curr_node->children) {
			curr_width /= 2.0;
			node_number = 0;
			for (uint8_t dir = 0; dir < 3; ++dir) {
				if (coor[dir] < curr_center[dir]) {
					curr_center[dir] -= curr_width / 2.0;
				}
				else {
					curr_center[dir] += curr_width / 2.0;
					node_number |= (1 << dir);
				}
			}
			++curr_depth;
			curr_node = &curr_node->children[node_number];
		}

		// This is a part of the node function centered at the origin
		// i.e Box*N(t/width)
		Symmetric_poly<N - 1, T> gaussian_appro = Gaussian_approximation<N, T>(0.5).scale(curr_width);
		gaussian_appro = gaussian_appro / gaussian_appro(0);

		Node* temp;
		Position<double> center;
		double dx, dy, dz;
		// Build a Rubik's Cube structure centered on curr_node
		Magic_cube curr_node_magic;
		curr_node_magic.construct(curr_node);
		for (uint8_t x = 0; x < 3; ++x)
			for (uint8_t y = 0; y < 3; ++y)
				for (uint8_t z = 0; z < 3; ++z) {
					temp = curr_node_magic(Position<uint8_t>(x, y, z));
					if (temp != nullptr) {
						center = temp->get_center();
						dx = coor[0] - center.x;
						dy = coor[1] - center.y;
						dz = coor[2] - center.z;
						ans_weight += temp->weight_alpha * gaussian_appro(dx) * gaussian_appro(dy) * gaussian_appro(dz);
					}
				}
		return ans_weight * (1.0 / std::pow(curr_width, 3.0));
	}

	template<uint16_t N, typename T>
	inline double Octree<N, T>::get_average_sampling_density_W() {
		double W = 0;
		
		// open output file
		std::ofstream output;
		std::string fname("sample density");
		output.open(fname, output.out | output.trunc);
		if (!output) {
			throw std::runtime_error("can't open output file: " + fname);
		}

		// open input file
		std::ifstream file;
		open_file(file, config.input_filename, config.binary);

		static const uint8_t kernel_depth = config.kernel_depth;
		static const double scale_ratio = config.get_sacle_ratio();
		float coor[6];  // must be float
		double w_temp = 0;
		while (!file.eof()) {
			read_one_line_from_file(file, coor, config.binary);

			// coordinate transform: transform to [0, 1]^3 range
			for (uint8_t dir = 0; dir < 3; ++dir) {
				coor[dir] = (coor[dir] - samples_property.center[dir]) / (samples_property.max_stretch * scale_ratio) + 0.5;
			}

			bool satisfy = true;
			for (uint8_t i = 0; i < 3; ++i) {
				if (coor[i] < 0 || coor[i] > 1) { satisfy = false; break; }
			}
			if (!satisfy) continue;

			w_temp = get_local_sampling_density(Position<double>(coor[0], coor[1], coor[2]));
			W += w_temp;

			output << w_temp << std::endl;
		}

		close_file(file);
		output.close();

		samples_property.mean_sample_density = W / samples_property.valid_samples_num;
		return samples_property.mean_sample_density;
	}

	// The position pos here is the converted coordinate
	template<uint16_t N, typename T>
	inline uint8_t poisson_reconstruction::Octree<N, T>::get_adaptive_depth_for_sample(const Position<double>& pos, double local_density) {
		double ans_adaptive_depth = 0.0;

		static double kernel_depth = (double)config.kernel_depth;
		static double max_depth = (double)config.max_depth;

		// local_density = get_local_sampling_density(pos);
		

		//ans_adaptive_depth = kernel_depth + std::log(local_density / samples_property.mean_sample_density) / std::log(4.0);
		//if (ans_adaptive_depth < kernel_depth) ans_adaptive_depth = std::floor(ans_adaptive_depth);
		//else ans_adaptive_depth = std::ceil(ans_adaptive_depth);

		ans_adaptive_depth = max_depth + std::log(local_density / samples_property.mean_sample_density) / std::log(4.0);

		if (ans_adaptive_depth < 0) ans_adaptive_depth = 0;
		if (ans_adaptive_depth > max_depth) ans_adaptive_depth = max_depth;

		return (uint8_t)ans_adaptive_depth;
	}

	template<uint16_t N, typename T>
	inline void poisson_reconstruction::Octree<N, T>::set_node_vector_field_factor() {
		// open sample weight file
		std::ifstream sample_density_file;
		open_file(sample_density_file, "sample density", false);

		// open input file
		std::ifstream file;
		open_file(file, config.input_filename, config.binary);

		static const double scale_ratio = config.get_sacle_ratio();

		// sample by sample processing
		float coor[6];  // must be float
		double local_density = 0;
		while (!file.eof()) {
			read_one_line_from_file(file, coor, config.binary);

			// coordinate transform: transform to [0, 1]^3 range
			for (uint8_t dir = 0; dir < 3; ++dir) {
				coor[dir] = (coor[dir] - samples_property.center[dir]) / (samples_property.max_stretch * scale_ratio) + 0.5;
			}

			bool satisfy = true;
			for (uint8_t i = 0; i < 3; ++i) {
				if (coor[i] < 0 || coor[i] > 1) { satisfy = false; break; }
			}
			if (!satisfy) continue;

			// normal vector normalization
			double normal_len = std::sqrt(coor[3] * coor[3] + coor[4] * coor[4] + coor[5] * coor[5]);
			for (uint8_t dir = 0; dir < 3; ++dir) {
				if (std::fabs(normal_len) < epsilon) coor[dir + 3] = 0;
				else coor[dir + 3] = coor[dir + 3] / normal_len;
			}

			sample_density_file >> local_density;
			// adaptive depth for this sample
			uint8_t adaptive_depth = get_adaptive_depth_for_sample(Position<double>(coor[0], coor[1], coor[2]), local_density);
			// std::cout << "adaptive depth : " << (int)adaptive_depth << std::endl;
			double reciprocal_local_density = 1.0 / local_density;

			// Determine the adaptive_depth layer node where the sample point is located
			// May be need to initialize children nodes
			uint8_t curr_depth = 0;
			double curr_center[3] = { 0.5, 0.5, 0.5 };
			double curr_width = 1.0;
			Node* curr_node = root;
			uint8_t node_number = 0;
			while (curr_depth < adaptive_depth) {

				if (curr_node->children == nullptr) curr_node->init_eight_children();
				
				++curr_depth;
				curr_width /= 2.0;
				node_number = 0;
				for (uint8_t dir = 0; dir < 3; ++dir) {
					if (coor[dir] < curr_center[dir]) {
						curr_center[dir] -= curr_width / 2.0;
					}
					else {
						curr_center[dir] += curr_width / 2.0;
						node_number |= (1 << dir);
					}
				}
				curr_node = &curr_node->children[node_number];
			}

			// Here the weight alpha is set using the n-degree convolution of the box function to simulate Gaussian function
			Symmetric_poly<N-1, T> gaussian_appro = Gaussian_approximation<N, T>(0.5);
			gaussian_appro.scale(curr_width);

			// Build a Rubik's Cube structure centered on curr_node
			Magic_cube curr_node_magic;
			curr_node_magic.construct(curr_node);
			Node* temp;
			Position<double> center;
			double dx, dy, dz;
			double vector_field_weight_alpha = 0;
			// set vector field factor
			for (uint8_t x = 0; x < 3; ++x)
				for (uint8_t y = 0; y < 3; ++y)
					for (uint8_t z = 0; z < 3; ++z) {
						temp = curr_node_magic(Position<uint8_t>(x, y, z));
						if (temp != nullptr) {
							center = temp->get_center();
							dx = coor[0] - center.x;
							dy = coor[1] - center.y;
							dz = coor[2] - center.z;
							// alpha
							vector_field_weight_alpha = gaussian_appro(dx) * gaussian_appro(dy) * gaussian_appro(dz);
							
							for (uint8_t dir = 0; dir < 3; ++dir) {
								temp->vector_field_factor[dir] += reciprocal_local_density * vector_field_weight_alpha * coor[dir+3];
							}
							temp->set_vector_field = true;
						}
					}
		}

		close_file(file);
		sample_density_file.close();
	}

	template<uint16_t N, typename T>
	inline void Octree<N, T>::Projection_of_divergence_node2node::operator()(Node* node, Node* curr_proc_node) {
		Node::Offset node_offset = node->offset;
		uint16_t node_index[3] = { 0, 0, 0 };
		node_index[0] = node_offset.x + (1 << node->depth) - 1;
		node_index[1] = node_offset.y + (1 << node->depth) - 1;
		node_index[2] = node_offset.z + (1 << node->depth) - 1;

		Node::Offset curr_proc_node_offset = curr_proc_node->offset;
		uint16_t curr_proc_node_index[3] = { 0, 0, 0 };
		curr_proc_node_index[0] = curr_proc_node_offset.x + (1 << curr_proc_node->depth) - 1;
		curr_proc_node_index[1] = curr_proc_node_offset.y + (1 << curr_proc_node->depth) - 1;
		curr_proc_node_index[2] = curr_proc_node_offset.z + (1 << curr_proc_node->depth) - 1;

		double x_dot, y_dot, z_dot;
		x_dot = (*ff_dot)[node_index[0]][curr_proc_node_index[0]];
		y_dot = (*ff_dot)[node_index[1]][curr_proc_node_index[1]];
		z_dot = (*ff_dot)[node_index[2]][curr_proc_node_index[2]];

		double dx_dot, dy_dot, dz_dot;
		dx_dot = (*df_dot)[node_index[0]][curr_proc_node_index[0]];
		dy_dot = (*df_dot)[node_index[1]][curr_proc_node_index[1]];
		dz_dot = (*df_dot)[node_index[2]][curr_proc_node_index[2]];

		double factor = (1 / std::pow(node->get_width(), 3)) * (1 / std::pow(curr_proc_node->get_width(), 3));

		curr_proc_node->proj_div_vec_field += factor * (
			node->vector_field_factor[0] * dx_dot * y_dot * z_dot +
			node->vector_field_factor[1] * x_dot * dy_dot * z_dot +
			node->vector_field_factor[2] * x_dot * y_dot * dz_dot
			);

	}

	template<uint16_t N, typename T>
	inline void Octree<N, T>::set_projection_of_divergence_vector_field() {
		octdata.set_ff_dot();
		octdata.set_df_dot();

		Projection_of_divergence_node2node pf;
		pf.octree = this;
		pf.ff_dot = &octdata.get_ff_dot();
		pf.df_dot = &octdata.get_df_dot();

		Node* curr_node = root->get_next_node();
		while (curr_node != nullptr) {
			Node::top_down_process_node_to_node(root, N * 0.5, curr_node, N * 0.5, pf);
			curr_node = root->get_next_node(curr_node);
		}

	}

	template<uint16_t N, typename T>
	inline void Octree<N, T>::Projection_of_indicator_node2node::operator()(Node* node, Node* curr_proc_node) {
		Node::Offset node_offset = node->offset;
		uint16_t node_index[3] = { 0, 0, 0 };
		node_index[0] = node_offset.x + (1 << node->depth) - 1;
		node_index[1] = node_offset.y + (1 << node->depth) - 1;
		node_index[2] = node_offset.z + (1 << node->depth) - 1;

		Node::Offset curr_proc_node_offset = curr_proc_node->offset;
		uint16_t curr_proc_node_index[3] = { 0, 0, 0 };
		curr_proc_node_index[0] = curr_proc_node_offset.x + (1 << curr_proc_node->depth) - 1;
		curr_proc_node_index[1] = curr_proc_node_offset.y + (1 << curr_proc_node->depth) - 1;
		curr_proc_node_index[2] = curr_proc_node_offset.z + (1 << curr_proc_node->depth) - 1;

		double x_dot, y_dot, z_dot;
		x_dot = (*ff_dot)[node_index[0]][curr_proc_node_index[0]];
		y_dot = (*ff_dot)[node_index[1]][curr_proc_node_index[1]];
		z_dot = (*ff_dot)[node_index[2]][curr_proc_node_index[2]];

		double ddx_dot, ddy_dot, ddz_dot;
		ddx_dot = (*dd_dot)[node_index[0]][curr_proc_node_index[0]];
		ddy_dot = (*dd_dot)[node_index[1]][curr_proc_node_index[1]];
		ddz_dot = (*dd_dot)[node_index[2]][curr_proc_node_index[2]];

		double factor = (1 / std::pow(node->get_width(), 3)) * (1 / std::pow(curr_proc_node->get_width(), 3));

		double value = factor * (
			ddx_dot * y_dot * z_dot +
			x_dot * ddy_dot * z_dot +
			x_dot * y_dot * ddz_dot
			);

		if (std::fabs(value) > epsilon) {
			Entry<T> entry;
			entry.value = value;
			entry.col = node->index;
			row_elements.push_back(entry);
		}
	}

	template<uint16_t N, typename T>
	inline void poisson_reconstruction::Octree<N, T>::solve_linear_system() {
		octdata.set_dd_dot();
		Projection_of_indicator_node2node pf;
		pf.octree = this;
		pf.ff_dot = &octdata.get_ff_dot();
		pf.dd_dot = &octdata.get_dd_dot();

		// All nodes on the octree are pushed into the container in depth-first search order
		std::vector<Node*> octree_total_nodes;
		Node* temp = root->get_next_node();
		uint32_t num = 0;
		while (temp != nullptr) {
			temp->index = num;	// Number each node
			octree_total_nodes.push_back(temp);
			temp = root->get_next_node(temp);
			++num;
		}

		uint32_t total_nodes_num = octree_total_nodes.size();
		Sparse_sym_matrix<T> matrix;
		std::vector<T> b(total_nodes_num);
		std::vector<T> solution(total_nodes_num);

		// process each node
		for (uint32_t i = 0; i < total_nodes_num; ++i) {
			b[i] = octree_total_nodes[i]->proj_div_vec_field;

			pf.row_elements.clear();
			Node::top_down_process_node_to_node(root, N * 0.5, octree_total_nodes[i], N * 0.5, pf);
			
			matrix.push_back(std::move(pf.row_elements));
		}

		//std::cout << matrix << std::endl;
		//std::cout << b << std::endl;

		// solve
		CG(matrix, b, solution, uint32_t(std::pow(total_nodes_num, 0.7)));

		// set the solution to each node
		//std::cout << "solution: \n";
		for (uint32_t i = 0; i < total_nodes_num; ++i) {
			octree_total_nodes[i]->solution_component = solution[i];
			//std::cout << solution[i] << std::endl;
		}

	}

	template<uint16_t N, typename T>
	inline void Octree<N, T>::calculate_indicator_func_value::operator()(Node* node, const Position<double>& point) {
		static Symmetric_poly<N - 1, T> gaussian_appro = Gaussian_approximation<N, T>(0.5);
		gaussian_appro = gaussian_appro / gaussian_appro(0);

		center = node->get_center();
		dx = point.x - center.x;
		dy = point.y - center.y;
		dz = point.z - center.z;

		double width = node->get_width();
		auto gaussian_scale = gaussian_appro.scale(width);

		func_value += node->solution_component * gaussian_scale(dx) * gaussian_scale(dy) * gaussian_scale(dz) * (1 / std::pow(width, 3.0));
	}

	template<uint16_t N, typename T>
	inline double Octree<N, T>::select_isovalue() {
		double ans_isovalue = 0;

		calculate_indicator_func_value pf;

		// open sample weight file
		std::ifstream sample_density_file;
		open_file(sample_density_file, "sample density", false);

		// open input file
		std::ifstream file;
		open_file(file, config.input_filename, config.binary);

		static const double scale_ratio = config.get_sacle_ratio();

		double weight_sum = 0;
		// sample by sample processing
		float coor[6];  // must be float
		double local_density = 0;
		double reciprocal_local_density = 0;
		while (!file.eof()) {
			read_one_line_from_file(file, coor, config.binary);

			// coordinate transform: transform to [0, 1]^3 range
			for (uint8_t dir = 0; dir < 3; ++dir) {
				coor[dir] = (coor[dir] - samples_property.center[dir]) / (samples_property.max_stretch * scale_ratio) + 0.5;
			}

			bool satisfy = true;
			for (uint8_t i = 0; i < 3; ++i) {
				if (coor[i] < 0 || coor[i] > 1) { satisfy = false; break; }
			}
			if (!satisfy) continue;

			Position<double> sample(coor[0], coor[1], coor[2]);
			
			sample_density_file >> local_density;
			reciprocal_local_density = 1.0 / local_density;
			weight_sum += reciprocal_local_density;

			pf.func_value = 0;
			Node::top_down_process_node_to_point(root, N * 0.5, sample, pf);
			//std::cout << "pf.func_value: " << pf.func_value << std::endl;

			ans_isovalue += pf.func_value * reciprocal_local_density;

			//std::cout << "weight_sum: " << weight_sum << " , ans_isovalue: " << ans_isovalue << std::endl;
		}

		close_file(file);
		sample_density_file.close();
		return ans_isovalue / weight_sum;
	}

	template<uint16_t N, typename T>
	inline void Octree<N, T>::calculate_center_indicator_value::operator()(Node* node, const Position<double>& point) {
		Node::Offset node_offset = node->offset;
		uint16_t node_index[3] = { 0, 0, 0 };
		node_index[0] = node_offset.x + (1 << node->depth) - 1;
		node_index[1] = node_offset.y + (1 << node->depth) - 1;
		node_index[2] = node_offset.z + (1 << node->depth) - 1;

		func_value += node->solution_component * (*value_table)[node_index[0]][center_index[0]] *
			(*value_table)[node_index[1]][center_index[1]] *
			(*value_table)[node_index[2]][center_index[2]] *
			(1.0 / std::pow(node->get_width(), 3.0));
	}

	template<uint16_t N, typename T>
	inline double Octree<N, T>::select_isovalue_scheme2() {
		calculate_center_indicator_value pf;
		octdata.set_values();
		pf.value_table = &octdata.get_values();

		double ans_isovalue = 0;
		double weight_sum = 0;
		double weight_temp = 0;
		double vec[3];
		uint16_t center_index[3];
		Node* curr_node = root->get_next_leaf_node();
		while (curr_node) {
			std::memcpy(vec, curr_node->vector_field_factor, sizeof(double) * 3);
			weight_temp = std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
			weight_sum += weight_temp;
			
			pf.func_value = 0;
			Node::get_center_index(curr_node, center_index, config.max_depth);
			std::memcpy(pf.center_index, center_index, sizeof(uint16_t) * 3);
			Node::top_down_process_node_to_point(root, N * 0.5, curr_node->get_center(), pf);
			ans_isovalue += weight_temp * pf.func_value;

			curr_node = root->get_next_leaf_node(curr_node);
		}
		return ans_isovalue / weight_sum;
	}

	template<uint16_t N, typename T>
	inline void Octree<N, T>::calculate_corner_indicator_value::operator()(Node* node, const Position<double>& point) {
		Node::Offset node_offset = node->offset;
		uint16_t node_index[3] = { 0, 0, 0 };
		node_index[0] = node_offset.x + (1 << node->depth) - 1;
		node_index[1] = node_offset.y + (1 << node->depth) - 1;
		node_index[2] = node_offset.z + (1 << node->depth) - 1;

		func_value += node->solution_component * (*value_table)[node_index[0]][point_index[0]] *
					  (*value_table)[node_index[1]][point_index[1]] *
					  (*value_table)[node_index[2]][point_index[2]] *
					  (1.0 / std::pow(node->get_width(), 3.0));
	}

	template<uint16_t N, typename T>
	inline void Octree<N, T>::set_cube_corners_value(double isovalue) {
		octdata.set_values();
		calculate_corner_indicator_value pf;
		pf.value_table = &octdata.get_values();
		// Avoid double counting
		std::unordered_map<long long, double> map_corners_value;

		// Gets all the leaf nodes
		std::vector<Node*> total_leaf_nodes;
		Node* temp = root->get_next_leaf_node();
		while (temp != nullptr) {
			total_leaf_nodes.push_back(temp);
			temp = root->get_next_leaf_node(temp);
		}

		uint32_t total_num = total_leaf_nodes.size();
		Node* curr_node;
		long long key;
		double corner_values[8];
		Position<double> corner_position;
		for (uint32_t i = 0; i < total_num; ++i) {
			curr_node = total_leaf_nodes[i];
			// 8 corners
			for (uint8_t corner = 0; corner < 8; ++corner) {
				key = Node::get_corner_index(curr_node, corner, pf.point_index, config.max_depth);
				Node::get_corner_coordinate(curr_node, corner, corner_position);
				pf.func_value = 0;
				if (map_corners_value.find(key) == map_corners_value.end()) {
					Node::top_down_process_node_to_point(root, N * 0.5, corner_position, pf);
					map_corners_value[key] = pf.func_value;
				}
				corner_values[corner] = map_corners_value[key];
			}
			// get cube index and set 8 corners's value
			curr_node->cube_index = Marching_cube::get_cube_index(corner_values, isovalue);
			//for (int ii = 0; ii < 8; ++ii) {
			//	std::cout << corner_values[ii] << " ";
			//}
			//std::cout << std::endl;
			//std::cout << "cube index: " << (int)curr_node->cube_index << std::endl;
			std::memcpy(curr_node->corner_values, corner_values, sizeof(double) * 8);
		}

		std::sort(total_leaf_nodes.begin(), total_leaf_nodes.end(), [](const Node* lhs, const Node* rhs) {
			return lhs->depth > rhs->depth;
			});
		for (uint32_t i = 0; i < total_num; ++i) {
			coarse_node_refinement(total_leaf_nodes[i], isovalue);
		}
	}

	template<uint16_t N, typename T>
	inline int Octree<N, T>::edge_root_count(Node* curr_node, uint8_t edge, double isovalue) {
		// Gets the two faces adjacent to this edge
		uint8_t face1, face2;
		Node::faces_adjacent_to_edge(edge, face1, face2);

		uint8_t search_edge = edge;
		Node* search_node = curr_node;
		
		Node* temp_node;
		// This is the only case where finer nodes may exist around this node
		if (curr_node->depth < config.max_depth) {
			temp_node = Node::node_adjacent_to_face(curr_node, face1);
			// There are more finer nodes
			if (temp_node && temp_node->children) {
				search_edge = edge ^ 1;
				search_node = temp_node;
			}
			else {
				temp_node = Node::node_adjacent_to_face(curr_node, face2);
				if (temp_node && temp_node->children) {
					search_edge = edge ^ 1;
					search_node = temp_node;
				}
				else {
					temp_node = Node::node_adjacent_to_edge(curr_node, edge);
					if (temp_node && temp_node->children) {
						search_edge = edge ^ 3;  // 3 <-> 0011
						search_node = temp_node;
					}
				}
			}
		}

		uint8_t c1, c2;
		Node::corner_adjacent_to_edge(search_edge, c1, c2);

		if (search_node->children) {
			return edge_root_count(&search_node->children[c1], search_edge, isovalue)
				+ edge_root_count(&search_node->children[c2], search_edge, isovalue);
		}
		else {
			// search_node is the finest node
			if (search_node->corner_values[c1] <= isovalue && search_node->corner_values[c2] <= isovalue ||
				search_node->corner_values[c1] >= isovalue && search_node->corner_values[c2] >= isovalue)
				return 0;
			else return 1;
		}

		// never reach
		return 0;
	}

	template<uint16_t N, typename T>
	inline int Octree<N, T>::face_root_count(Node* curr_node, uint8_t face, double isovalue) {
		static const uint8_t search[6][4] = {
			// left-top    right-bottom
			{    4, 10,      6, 8     },  // face 0
			{    5, 11,      7, 9     },  // face 1
			{    0,  9,      2, 8     },  // face 2
			{    1, 11,      1,10     },  // face 3
			{    0,  5,      1, 4     },  // face 4
			{    3,  7,      2, 6     },  // face 5
		};
		Node* search_node = curr_node; 
		uint8_t search_face = face;

		Node* adjacent_node = Node::node_adjacent_to_face(curr_node, face);
		if (adjacent_node && adjacent_node->children) {
			search_face = face ^ 1;
			search_node = adjacent_node;
		}

		if (search_node->children) {
			uint8_t c1, c2, c3, c4;
			Node::corner_adjacent_to_face(search_face, c1, c2, c3, c4);
			int ans = 0;
			ans += edge_root_count(&search_node->children[c1], search[search_face][0], isovalue) +
				edge_root_count(&search_node->children[c1], search[search_face][1], isovalue);

			ans += edge_root_count(&search_node->children[c4], search[search_face][2], isovalue) +
				edge_root_count(&search_node->children[c4], search[search_face][3], isovalue);

			ans += face_root_count(&search_node->children[c1], search_face, isovalue) +
				face_root_count(&search_node->children[c2], search_face, isovalue) +
				face_root_count(&search_node->children[c3], search_face, isovalue) +
				face_root_count(&search_node->children[c4], search_face, isovalue);
			return ans;
		}
		
		return 0;
	}

	template<uint16_t N, typename T>
	inline void Octree<N, T>::divide_node(Node* curr_node, double isovalue) {
		calculate_corner_indicator_value pf;
		pf.value_table = &octdata.get_values();

		double children_corner_values[8][8];
		curr_node->init_eight_children();

		// 8 of 27 corners set
		for (uint8_t i = 0; i < 8; ++i) {
			children_corner_values[i][i] = curr_node->corner_values[i];
		}

		uint16_t center_index[3];
		Node::get_center_index(curr_node, center_index, config.max_depth);
		Node::Offset offset = curr_node->offset;
		Position<double> center_coor = curr_node->get_center();
		double center_width = curr_node->get_width();

		pf.func_value = 0;
		std::memcpy(pf.point_index, center_index, sizeof(uint16_t) * 3);
		Node::top_down_process_node_to_point(root, N * 0.5, center_coor, pf);
		// 9 of 27 set
		for (uint8_t i = 0; i < 8; ++i) {
			children_corner_values[i][i ^ 7] = pf.func_value;  // 7 <-> 111
		}

		uint8_t corners[4];
		// 15 of 27 set
		for (uint8_t face = 0; face < 6; ++face) {
			std::memcpy(pf.point_index, center_index, sizeof(uint16_t) * 3);
			Position<double> temp_coor(center_coor);
			uint8_t dir = (face & 6) >> 1;
			switch (dir) {
			case 0:  // x-direction
				pf.point_index[dir] = (offset.x + (face & 1)) << (config.max_depth + 1 - curr_node->depth);
				if (face & 1) temp_coor.x += center_width * 0.5;
				else temp_coor.x -= center_width * 0.5;
				break;
			case 1:  // y-direction
				pf.point_index[dir] = (offset.y + (face & 1)) << (config.max_depth + 1 - curr_node->depth);
				if (face & 1) temp_coor.y += center_width * 0.5;
				else temp_coor.y -= center_width * 0.5;
				break;
			case 2:  // z-direction
				pf.point_index[dir] = (offset.z + (face & 1)) << (config.max_depth + 1 - curr_node->depth);
				if (face & 1) temp_coor.z += center_width * 0.5;
				else temp_coor.z -= center_width * 0.5;
				break;
			}
			pf.func_value = 0;
			Node::top_down_process_node_to_point(root, N * 0.5, temp_coor, pf);
			
			Node::corner_adjacent_to_face(face, corners[0], corners[1], corners[2], corners[3]);
			for (uint8_t i = 0; i < 4; ++i) {
				switch (dir) {
				case 0:
					children_corner_values[corners[i]][corners[i] ^ 6] = pf.func_value;  // Negate y and z,  6 <-> 110
					break;
				case 1: 
					children_corner_values[corners[i]][corners[i] ^ 5] = pf.func_value;  // Negate x and z,  5 <-> 101
					break;
				case 2:
					children_corner_values[corners[i]][corners[i] ^ 3] = pf.func_value;  // Negate x and y,  3 <-> 011
					break;
				}
			}
		}

		uint8_t c1, c2;
		// 27 of 27 set
		for (uint8_t edge = 0; edge < 12; ++edge) {
			std::memcpy(pf.point_index, center_index, sizeof(uint16_t) * 3);
			Position<double> temp_coor(center_coor);
			uint8_t dir = (edge & 12) >> 2;  // 12<->1100
			switch (dir) {
			case 0:  // x-direction

				if (edge & 1) temp_coor.y += center_width * 0.5;
				else temp_coor.y -= center_width * 0.5;
				if (edge & 2) temp_coor.z += center_width * 0.5;
				else temp_coor.z -= center_width * 0.5;

				pf.point_index[1] = (offset.y + (edge & 1)) << (config.max_depth + 1 - curr_node->depth);
				pf.point_index[2] = (offset.z + ((edge & 2) >> 1)) << (config.max_depth + 1 - curr_node->depth);

				break;

			case 1:  // y-direction

				if (edge & 1) temp_coor.x += center_width * 0.5;
				else temp_coor.x -= center_width * 0.5;
				if (edge & 2) temp_coor.z += center_width * 0.5;
				else temp_coor.z -= center_width * 0.5;

				pf.point_index[0] = (offset.x + (edge & 1)) << (config.max_depth + 1 - curr_node->depth);
				pf.point_index[2] = (offset.z + ((edge & 2) >> 1)) << (config.max_depth + 1 - curr_node->depth);

				break;
			case 2:  // z-direction

				if (edge & 1) temp_coor.x += center_width * 0.5;
				else temp_coor.x -= center_width * 0.5;
				if (edge & 2) temp_coor.y += center_width * 0.5;
				else temp_coor.y -= center_width * 0.5;

				pf.point_index[0] = (offset.x + (edge & 1)) << (config.max_depth + 1 - curr_node->depth);
				pf.point_index[1] = (offset.y + ((edge & 2)>>1)) << (config.max_depth + 1 - curr_node->depth);

				break;
			}
			pf.func_value = 0;
			Node::top_down_process_node_to_point(root, N * 0.5, temp_coor, pf);

			Node::corner_adjacent_to_edge(edge, c1, c2);
			children_corner_values[c1][c2] = pf.func_value;
			children_corner_values[c2][c1] = pf.func_value;
		}

		// for 8 children
		// set 8 corners's value ans set corresponding cube index
		for (uint8_t i = 0; i < 8; ++i) {
			std::memcpy(curr_node->children[i].corner_values, children_corner_values[i], sizeof(double) * 8);
			curr_node->children[i].cube_index = Marching_cube::get_cube_index(curr_node->children[i].corner_values, isovalue);
		}

	}

	template<uint16_t N, typename T>
	inline void Octree<N, T>::coarse_node_refinement(Node* curr_node, double isovalue) {
		if (curr_node->depth >= config.max_depth) return;

		// Determine whether the curr_node needs to be divided
		bool divide = false;

		// Check if the node has edges with more than one root
		// 检查边i与等值面的相交点数, 若相交数目大于1，将sub设置为1
		for (uint8_t edge = 0; edge < 12 && !divide; ++edge) {
			if (edge_root_count(curr_node, edge, isovalue) > 1) {
				divide = true;
			} 
		}

		for (uint8_t face = 0; face < 6 && !divide; ++face) {
			if (face_root_count(curr_node, face, isovalue) > 0) {
				divide = true;
			}
		}

		if (divide) {
			divide_node(curr_node, isovalue);
			Node* temp_node;
			for (uint8_t face = 0; face < 6; ++face) {
				temp_node = Node::node_adjacent_to_face(curr_node, face);
				if (temp_node && !temp_node->children) { 
					coarse_node_refinement(temp_node, isovalue); 
				}
			}
			for (uint8_t edge = 0; edge < 12; ++edge) {
				temp_node = Node::node_adjacent_to_edge(curr_node, edge);
				if (temp_node && !temp_node->children) { 
					coarse_node_refinement(temp_node, isovalue);
				}
			}
			for (uint8_t i = 0; i < 8; ++i) {
				coarse_node_refinement(&curr_node->children[i], isovalue);
			}
		}
	}

	template<uint16_t N, typename T>
	inline bool Octree<N, T>::search_edge_root_key(Node* curr_node, uint8_t edge, double isovalue, long long& key) {
		// Gets the two faces adjacent to this edge
		uint8_t face1, face2;
		Node::faces_adjacent_to_edge(edge, face1, face2);

		uint8_t search_edge = edge;
		Node* search_node = curr_node;

		Node* temp_node;
		// This is the only case where finer nodes may exist around this node
		if (curr_node->depth < config.max_depth) {
			temp_node = Node::node_adjacent_to_face(curr_node, face1);
			// There are more finer nodes
			if (temp_node && temp_node->children) {
				search_edge = edge ^ 1;
				search_node = temp_node;
			}
			else {
				temp_node = Node::node_adjacent_to_face(curr_node, face2);
				if (temp_node && temp_node->children) {
					search_edge = edge ^ 1;
					search_node = temp_node;
				}
				else {
					temp_node = Node::node_adjacent_to_edge(curr_node, edge);
					if (temp_node && temp_node->children) {
						search_edge = edge ^ 3;  // 3 <-> 0011
						search_node = temp_node;
					}
				}
			}
		}

		uint8_t c1, c2;
		Node::corner_adjacent_to_edge(search_edge, c1, c2);
		if (search_node->children) {
			if (search_edge_root_key(&search_node->children[c1], search_edge, isovalue, key)) return true;
			else if (search_edge_root_key(&search_node->children[c2], search_edge, isovalue, key)) return true;
			else return false;
		}
		else {
			// search_node is the finest node
			if (search_node->corner_values[c1] <= isovalue && search_node->corner_values[c2] <= isovalue ||
				search_node->corner_values[c1] >= isovalue && search_node->corner_values[c2] >= isovalue)
				return false;
			
			Node::Offset offset = search_node->offset;
			uint8_t dir = (search_edge & 12) >> 2;
			// the center of edge
			uint16_t index[3];
			switch (dir) {
			case 0: 
				index[1] = (offset.y + (search_edge & 1)) << (config.max_depth + 1 - search_node->depth);
				index[2] = (offset.z + ((search_edge & 2) >> 1)) << (config.max_depth + 1 - search_node->depth);
				index[0] = ((offset.x << 1) + 1) << (config.max_depth - search_node->depth);
				break;
			case 1:
				index[0] = (offset.x + (search_edge & 1)) << (config.max_depth + 1 - search_node->depth);
				index[2] = (offset.z + ((search_edge & 2) >> 1)) << (config.max_depth + 1 - search_node->depth);
				index[1] = ((offset.y << 1) + 1) << (config.max_depth - search_node->depth);
				break;
			case 2:
				index[0] = (offset.x + (search_edge & 1)) << (config.max_depth + 1 - search_node->depth);
				index[1] = (offset.y + ((search_edge & 2) >> 1)) << (config.max_depth + 1 - search_node->depth);
				index[2] = ((offset.z << 1) + 1) << (config.max_depth - search_node->depth);
				break;
			}
			key = (long long)index[0] | ((long long)index[1] << 15) | ((long long)index[2] << 30);

			return true;
		}

		// never reach
		return false;
	}

	template<uint16_t N, typename T>
	inline bool Octree<N, T>::calculate_edge_root(Node* curr_node, uint8_t edge, double isovalue, Position<float>& intersect_pos) {
		// Gets the two faces adjacent to this edge
		uint8_t face1, face2;
		Node::faces_adjacent_to_edge(edge, face1, face2);

		uint8_t search_edge = edge;
		Node* search_node = curr_node;

		Node* temp_node;
		// This is the only case where finer nodes may exist around this node
		if (curr_node->depth < config.max_depth) {
			temp_node = Node::node_adjacent_to_face(curr_node, face1);
			// There are more finer nodes
			if (temp_node && temp_node->children) {
				search_edge = edge ^ 1;
				search_node = temp_node;
			}
			else {
				temp_node = Node::node_adjacent_to_face(curr_node, face2);
				if (temp_node && temp_node->children) {
					search_edge = edge ^ 1;
					search_node = temp_node;
				}
				else {
					temp_node = Node::node_adjacent_to_edge(curr_node, edge);
					if (temp_node && temp_node->children) {
						search_edge = edge ^ 3;  // 3 <-> 0011
						search_node = temp_node;
					}
				}
			}
		}

		uint8_t c1, c2;
		Node::corner_adjacent_to_edge(search_edge, c1, c2);
		if (search_node->children) {
			if (calculate_edge_root(&search_node->children[c1], search_edge, isovalue, intersect_pos)) return true;
			else if (calculate_edge_root(&search_node->children[c2], search_edge, isovalue, intersect_pos)) return true;
			else return false;
		}
		else {
			// search_node is the finest node
			if (search_node->corner_values[c1] <= isovalue && search_node->corner_values[c2] <= isovalue ||
				search_node->corner_values[c1] >= isovalue && search_node->corner_values[c2] >= isovalue)
				return false;

			Position<double> corner1_pos, corner2_pos;
			Node::get_corner_coordinate(search_node, c1, corner1_pos);
			Node::get_corner_coordinate(search_node, c2, corner2_pos);
			double corner1_value = search_node->corner_values[c1];
			double corner2_value = search_node->corner_values[c2];
			uint8_t dir = (search_edge & 12) >> 2;
			switch (dir) {
			case 0:
				intersect_pos.y = corner1_pos.y;
				intersect_pos.z = corner1_pos.z;
				intersect_pos.x = corner1_pos.x + (isovalue - corner1_value) * (corner2_pos.x - corner1_pos.x) / (corner2_value - corner1_value);
				break;
			case 1:
				intersect_pos.x = corner1_pos.x;
				intersect_pos.z = corner1_pos.z;
				intersect_pos.y = corner1_pos.y + (isovalue - corner1_value) * (corner2_pos.y - corner1_pos.y) / (corner2_value - corner1_value);
				break;
			case 2:
				intersect_pos.x = corner1_pos.x;
				intersect_pos.y = corner1_pos.y;
				intersect_pos.z = corner1_pos.z + (isovalue - corner1_value) * (corner2_pos.z - corner1_pos.z) / (corner2_value - corner1_value);
				break;
			}

			return true;
		}

		// never reach
		return false;
	}

	template<uint16_t N, typename T>
	inline void Octree<N, T>::set_intersection(double isovalue) {
		edge_key_to_index.clear();
		std::vector<Node*> total_leaf_nodes;
		Node* temp = root->get_next_leaf_node();
		while (temp) {
			total_leaf_nodes.push_back(temp);
			temp = root->get_next_leaf_node(temp);
		}

		std::sort(total_leaf_nodes.begin(), total_leaf_nodes.end(), [](const Node* lhs, const Node* rhs) {
			return lhs->depth > rhs->depth;
			});

		uint32_t total_num = total_leaf_nodes.size();
		Node* curr;
		long long key;
		Position<float> intersect_pos;
		for (uint32_t i = 0; i < total_num; ++i) {
			curr = total_leaf_nodes[i];
			if (curr->cube_index != 0 && curr->cube_index != 255) {
				// edge by edge processing
				for (uint8_t edge = 0; edge < 12; ++edge) {
					if (search_edge_root_key(curr, edge, isovalue, key)) {
						if (edge_key_to_index.find(key) == edge_key_to_index.end()) {
							calculate_edge_root(curr, edge, isovalue, intersect_pos);
							mesh_data.intersection_points.push_back(intersect_pos);
							edge_key_to_index[key] = mesh_data.intersection_points.size() - 1;
						}
					}
				}
			}

		}
	}

	template<uint16_t N, typename T>
	inline void Octree<N, T>::get_triangles(double isovalue) {
		int num = 0;

		long long key;
		int index;
		std::vector<int> tris;
		std::vector<int> tri_data;
		uint8_t tri_num = 0;
		Node* curr_node = root->get_next_leaf_node();
		while (curr_node) {
			tri_num = Marching_cube::add_triangles(curr_node->cube_index, tris);

			if (tri_num > 0) {
				++num;
				//std::cout << "tri_num: " << (int)tri_num << " , process num: " << num << std::endl;
				int real_tri_num = 0;
				for (uint8_t i = 0; i < tri_num; ++i) {
					for (uint8_t j = 0; j < 3; ++j) {
						if (search_edge_root_key(curr_node, tris[i * 3 + j], isovalue, key)) {
							if (edge_key_to_index.find(key) != edge_key_to_index.end()) {
								index = edge_key_to_index[key];
								tri_data.push_back(index);
							}
							else {
								std::cout << "The expected key was not found in the edge mapping table, key: " << key << std::endl;
							}
						}
						else {
							std::cout << "No intersection on the expected edge: " << tris[i * 3 + j] << std::endl;
						}
					}
					if (tri_data.size() == 3) {
						mesh_data.triangles.push_back(tri_data);
						++real_tri_num;
					}
					else {
						std::cout << "The expected triangle was not found, just find " << tri_data.size() << " vertex" << std::endl;
					}
					tri_data.clear();
				}
				//std::cout << "real triangle num / expected triangle num: " << real_tri_num << " / " << (int)tri_num << std::endl;
			}

			tris.clear();
			curr_node = root->get_next_leaf_node(curr_node);
		}
	}

}

#endif // !POISSON_RECONSTRUCTION_OCTREE_H
