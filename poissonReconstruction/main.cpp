/**
  * FileName: main.cpp
  * Author: Zeng Zhenxiang
  * Version: 1.0
  * Date:
  * Description:
***/
#include <iostream>
#include <string>
#include <time.h>
#include "node.h"
#include "poisson.h"
#include "gaussian.h"
#include "octdata.h"
#include <windows.h>
#include <psapi.h>
#pragma comment(lib,"psapi.lib")
namespace poisson = poisson_reconstruction;

void show_memory_info()
{
    HANDLE handle = GetCurrentProcess();
    PROCESS_MEMORY_COUNTERS pmc;
    double _MB = 1024.0 * 1024;
    int ret = GetProcessMemoryInfo(handle, &pmc, sizeof(pmc));
    if (ret != 0) {
        std::cout << "memory usage: (MB)\n"
            << "WorkingSetSize:\t" << pmc.WorkingSetSize / _MB << ",\n"
            << "PeakWorkingSetSize:\t" << pmc.PeakWorkingSetSize / _MB << ",\n"
            << "PagefileUsage:\t" << pmc.PagefileUsage / _MB << ",\n"
            << "PeakPagefileUsage:\t" << pmc.PeakPagefileUsage / _MB << ".\n";
    }
    else {
        std::cerr << "call function show_memory_info failed.\n";
        throw std::runtime_error("call function show_memory_info error.");
    }
}

int main() {
	time_t start, end;
	time(&start);
    
	//std::string input_file("bunny.bnpts");
	//std::string output_file("bunny_d9_kd9_cg0.5.ply");
    std::string input_file("horse.npts.txt");
    std::string output_file("horse_d8_kd5_cg0.7.ply");
    //std::string input_file("gargo_v150000.txt");
    //std::string output_file("gargo_d9_kd9_cg0.5.ply");
    //std::string input_file("horse_v100000_gaussian_noise_0.01.txt");
    //std::string output_file("horse_v100000_gaussian_noise_0.01_d8_kd5_cg0.6.ply");
	uint8_t max_depth = 8;
	poisson::Config config(input_file, output_file, max_depth, false);
	poisson::run_poisson<3, double>(config);

	time(&end);
	std::cout << "total time: " << end - start << "sec\n";
    show_memory_info();
	return 0;
}