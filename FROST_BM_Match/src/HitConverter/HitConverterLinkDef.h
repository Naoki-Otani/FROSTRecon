#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class std::vector<int>+;
#pragma link C++ class std::vector<std::vector<int>>+;
#pragma link C++ class std::vector<std::vector<std::vector<int>>>+;

#pragma link C++ class std::vector<double>+;
#pragma link C++ class std::vector<std::vector<double>>+;
#pragma link C++ class std::vector<std::vector<std::vector<double>>>+;

#pragma link C++ typedef VecInt;
#pragma link C++ typedef VecVecInt;
#pragma link C++ typedef VecVecVecInt;
#pragma link C++ typedef VecDouble;
#pragma link C++ typedef VecVecDouble;
#pragma link C++ typedef VecVecVecDouble;
#endif
