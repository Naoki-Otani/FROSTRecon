#pragma once

#include <vector>

// Force nested STL types to be visible to rootcling.
using VecInt = std::vector<int>;
using VecVecInt = std::vector<std::vector<int>>;
using VecVecVecInt = std::vector<std::vector<std::vector<int>>>;

using VecDouble = std::vector<double>;
using VecVecDouble = std::vector<std::vector<double>>;
using VecVecVecDouble = std::vector<std::vector<std::vector<double>>>;
