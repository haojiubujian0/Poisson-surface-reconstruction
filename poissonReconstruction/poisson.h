/**
  * FileName: poisson.h
  * Author: Zeng Zhenxiang
  * Version: 1.0
  * Date:
  * Description:
***/

#ifndef POISSON_RECONSTRUCTION_POISSON_H
#define POISSON_RECONSTRUCTION_POISSON_H

#include <iostream>
#include <string>
#include "octree.h"
#include "config.h"
#include "mesh.h"
#include "ply.h"
#include "node.h"

namespace poisson_reconstruction {

	// Starting function of Poisson reconstruction
	// Using a 3-degree box filter convolution and has data accuracy of double
	template<uint8_t N = 3, typename T = double>	
	void run_poisson(const Config& config) {
		const Node* temp;
		std::cout << config;

		std::cout << std::endl;
		std::cout << "Octree initialization..." << std::endl;
		Octree<N, T> octree(config);
		std::cout << "Octree initialization finish." << std::endl;

		std::cout << std::endl;
		std::cout << "***************************************************" << std::endl;
		std::cout << "1-th stage: adaptive construct octree from samples" << std::endl;
		std::cout << "***************************************************" << std::endl;
		int total_samples_num = octree.adaptive_construct_octree_from_samples();
		std::cout << "total sample size: " << total_samples_num << std::endl;
		temp = octree.get_root();
		std::cout << "total nodes / leaf nodes: " 
				  << Node::total_nodes_count(temp) << " / " << Node::total_leaf_node_count(temp) << std::endl;
		std::cout << "1-th stage finish." << std::endl;

		std::cout << std::endl;
		std::cout << "***************************************************" << std::endl;
		std::cout << "2-th stage: set node weight alpha" << std::endl;
		std::cout << "***************************************************" << std::endl;
		octree.set_node_weight_alpha();
		std::cout << "2-th stage finish." << std::endl;

		std::cout << std::endl;
		std::cout << "***************************************************" << std::endl;
		std::cout << "3-th stage: set node vector field factor" << std::endl;
		std::cout << "***************************************************" << std::endl;
		double w = octree.get_average_sampling_density_W();
		std::cout << " the average sampling density W over all of the samples: " << w << std::endl;
		octree.set_node_vector_field_factor();
		temp = octree.get_root();
		std::cout << "total nodes / leaf nodes: "
			<< Node::total_nodes_count(temp) << " / " << Node::total_leaf_node_count(temp) << std::endl;
		std::cout << "3-th stage finish." << std::endl;

		std::cout << std::endl;
		std::cout << "***************************************************" << std::endl;
		std::cout << "4-th stage: clip empty leaf nodes" << std::endl;
		std::cout << "***************************************************" << std::endl;
		octree.clip_empty_leaf_nodes();
		temp = octree.get_root();
		std::cout << "total nodes / leaf nodes: "
			<< Node::total_nodes_count(temp) << " / " << Node::total_leaf_node_count(temp) << std::endl;
		std::cout << "4-th stage finish." << std::endl;

		std::cout << std::endl;
		std::cout << "***************************************************" << std::endl;
		std::cout << "5-th stage: set the projection of the divergence of vector field" << std::endl;
		std::cout << "***************************************************" << std::endl;
		octree.set_projection_of_divergence_vector_field();
		std::cout << "5-th stage finish." << std::endl;

		std::cout << std::endl;
		std::cout << "***************************************************" << std::endl;
		std::cout << "6-th stage: solve linear system" << std::endl;
		std::cout << "***************************************************" << std::endl;
		octree.solve_linear_system();
		std::cout << "6-th stage finish." << std::endl;

		std::cout << std::endl;
		std::cout << "***************************************************" << std::endl;
		std::cout << "7-th stage: select an isovalue" << std::endl;
		std::cout << "***************************************************" << std::endl;
		double isovalue = octree.select_isovalue_scheme2();
		std::cout << "isovalue is: " << isovalue << std::endl;
		std::cout << "7-th stage finish." << std::endl;

		// run MC
		std::cout << std::endl;
		std::cout << "***************************************************" << std::endl;
		std::cout << "8-th stage: extract isosurface " << std::endl;
		std::cout << "***************************************************" << std::endl;

		std::cout << "8-1: set cube corners value and refinement" << std::endl;
		octree.set_cube_corners_value(isovalue);
		temp = octree.get_root();
		std::cout << "total nodes / leaf nodes: "
			<< Node::total_nodes_count(temp) << " / " << Node::total_leaf_node_count(temp) << std::endl;
		std::cout << "8-1 finish." << std::endl;
		
		std::cout << "8-2: set intersection" << std::endl;
		octree.set_intersection(isovalue);
		std::cout << "8-2 finish." << std::endl;

		std::cout << "8-3: get iso-triangles" << std::endl;
		octree.get_triangles(isovalue);
		std::cout << "8-3 finish." << std::endl;

		std::cout << "8-th stage finish." << std::endl;

		Mesh_data mesh = octree.get_mesh_data();
		Samples_property sample_pro = octree.get_sample_property();
		std::cout << std::endl;
		std::cout << "total vertices : " << mesh.intersection_points.size() << std::endl;
		std::cout << "total triangles: " << mesh.triangles.size() << std::endl;

		// output .ply file
		std::cout << std::endl;
		std::cout << "***************************************************" << std::endl;
		std::cout << "9-th stage: output .ply triangles " << std::endl;
		std::cout << "***************************************************" << std::endl;
		PlyWriteTriangles((char*)config.output_filename.c_str(), mesh, PLY_BINARY_NATIVE, 
			sample_pro.center, sample_pro.max_stretch, config.get_sacle_ratio(), nullptr, 0);
		std::cout << "9-th stage finish." << std::endl;
	}

}

#endif // !POISSON_RECONSTRUCTION_POISSON_H
