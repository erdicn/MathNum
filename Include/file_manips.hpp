#ifndef FILE_MANIPS_HPP
#define FILE_MANIPS_HPP

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdint>

template <typename T>
void saveArrayToFile(T* arr, uint64_t arr_len, const char* filename) {
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error opening file " << filename 
                  << " in code file " << __FILE__ 
                  << " line " << __LINE__ << "\n";
        return;
    }

    // Force 16 decimal places for floating-point types (ignored by ints)
    file << std::fixed << std::setprecision(16); // TODO might need to optimise or adapt this if it is float  sizeof(T) <= 4 ? 8 : 16

    for (uint64_t i = 0; i < arr_len; ++i) {
        file << arr[i] << "\n"; 
    }

    file.close();
}

#endif /* FILE_MANIPS_HPP */
