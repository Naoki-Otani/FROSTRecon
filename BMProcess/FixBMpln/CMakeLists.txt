cmake_minimum_required(VERSION 3.16)
project(FixBMpln VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ROOT REQUIRED COMPONENTS Core Hist Tree RIO)
include(${ROOT_USE_FILE})

include(GNUInstallDirs)

set(BM_DATA_CLASSES
  BMBasicRecon
  BMBeaminfo
  PMRecon
  BMBSD
)

set(BM_DATA_SOURCES "")
set(BM_DATA_DICT_SOURCES "")

foreach(CLASS ${BM_DATA_CLASSES})

  message(STATUS "=== Generating dictionary for ${CLASS} ===")

  ROOT_GENERATE_DICTIONARY(${CLASS}Dict
    ${PROJECT_SOURCE_DIR}/include/project/${CLASS}.hpp
    LINKDEF ${PROJECT_SOURCE_DIR}/src/${CLASS}LinkDef.hh
    OPTIONS
      -I${PROJECT_SOURCE_DIR}/include
      -I${PROJECT_SOURCE_DIR}/include/project
      -I${PROJECT_SOURCE_DIR}/src
  )

  list(APPEND BM_DATA_SOURCES
    ${PROJECT_SOURCE_DIR}/src/${CLASS}.cpp
  )

  list(APPEND BM_DATA_DICT_SOURCES
    ${CMAKE_CURRENT_BINARY_DIR}/${CLASS}Dict.cxx
  )

endforeach()

add_library(bm_ana_data SHARED
  ${BM_DATA_SOURCES}
  ${BM_DATA_DICT_SOURCES}
)

target_include_directories(bm_ana_data
  PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/include/project
    ${PROJECT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(bm_ana_data
  PUBLIC
    ROOT::Core
    ROOT::Hist
    ROOT::Tree
    ROOT::RIO
)

add_executable(fix_bmbasicrecon_class
  ${PROJECT_SOURCE_DIR}/app/fix_bmbasicrecon_class.cpp
)

target_include_directories(fix_bmbasicrecon_class
  PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/include/project
    ${PROJECT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(fix_bmbasicrecon_class
  PRIVATE
    bm_ana_data
    ROOT::Core
    ROOT::Hist
    ROOT::Tree
    ROOT::RIO
)

add_executable(remove_deadpln_topview
  ${PROJECT_SOURCE_DIR}/app/remove_deadpln_topview.cpp
)

target_include_directories(remove_deadpln_topview
  PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/include/project
    ${PROJECT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(remove_deadpln_topview
  PRIVATE
    bm_ana_data
    ROOT::Core
    ROOT::Hist
    ROOT::Tree
    ROOT::RIO
)

install(
  TARGETS bm_ana_data fix_bmbasicrecon_class remove_deadpln_topview
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

foreach(CLASS ${BM_DATA_CLASSES})
  install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/lib${CLASS}Dict_rdict.pcm
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    OPTIONAL
  )
endforeach()
