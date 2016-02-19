// fgt — fast Gauss transforms
// Copyright (C) 2016 Peter J. Gadomski <pete.gadomski@gmail.com>
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA

#include <chrono>
#include <iostream>
#include <random>

#include "fgt.hpp"

const char* USAGE =
    "USAGE: bench <mode> <dataset> <rows> <cols> <bandwidth>\n\n"
    "  where:\n"
    "    mode = direct, direct_tree, ifgt, ifgt_adaptive\n"
    "    dataset = random\n"
    "    rows = any number\n"
    "    cols = any number\n"
    "    bandwidth = any number\n";
const double DEFAULT_EPSILON = 1e-4;

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cout << "Incorrect number of arguments.\n" << USAGE;
        return 1;
    }

    double epsilon = DEFAULT_EPSILON;

    std::string mode = argv[1];
    std::string dataset = argv[2];
    size_t rows = std::stoi(argv[3]);
    size_t cols = std::stoi(argv[4]);
    double bandwidth = std::stod(argv[5]);

    if (dataset != "random") {
        std::cout << "Invalid dataset: " << dataset << "\n";
        return 1;
    }

    std::vector<double> source(rows * cols);
    std::vector<double> target(rows * cols);
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            source[i * cols + j] = distribution(generator);
            target[i * cols + j] = distribution(generator);
        }
    }

    std::vector<double> result;
    auto tic = std::chrono::high_resolution_clock::now();
    if (mode == "direct") {
        fgt::direct(source.data(), rows, target.data(), rows, cols, bandwidth);
    } else if (mode == "direct_tree") {
        fgt::direct_tree(source.data(), rows, target.data(), rows, cols,
                         bandwidth, epsilon);
    } else if (mode == "ifgt") {
        fgt::ifgt(source.data(), rows, target.data(), rows, cols, bandwidth,
                  epsilon);
    } else {
        std::cout << "Invalid mode: " << mode << "\n";
        return 1;
    }
    auto toc = std::chrono::high_resolution_clock::now();

    auto runtime = toc - tic;
    std::cout << double(std::chrono::duration_cast<std::chrono::microseconds>(
                            runtime)
                            .count()) *
                     1e-6f
              << "\n";
    return 0;
}
