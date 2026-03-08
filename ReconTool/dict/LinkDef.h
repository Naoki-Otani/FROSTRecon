#ifndef FROST_LINKDEF_H
#define FROST_LINKDEF_H

// LinkDef.h is used only by rootcling to generate ROOT dictionaries.
// Do NOT include this file from normal C++ source files.

#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class std::vector<double>+;
#pragma link C++ class std::vector<int>+;
#pragma link C++ class std::vector<std::vector<double>>+;
#pragma link C++ class std::vector<std::vector<int>>+;
#pragma link C++ class std::vector<std::vector<std::vector<int>>>+;

#endif

#endif
