# Generate ROOT dictionary
root_generate_dictionary(G__NTBMDict
	NTBMSummary.hh
	LINKDEF NTBMLinkDef.hh)

# libNTBM.so shared library
add_library(libNTBM SHARED NTBMSummary.hh NTBMSummary.cc G__NTBMDict.cxx)

# rename liblibNTBM into libNTBM
set_target_properties(libNTBM PROPERTIES OUTPUT_NAME NTBM)

# include the current directory
target_include_directories(libNTBM PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

# Link with ...
target_link_libraries(libNTBM PUBLIC
			      ${ROOT_LIBRARIES}
			      ${Boost_LIBRARIES}
			      Boost::system
			      Boost::filesystem
			      Boost::log)

# list all target headers
file(GLOB NTBM_LIB_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/*.hh")
set_target_properties(libNTBM PROPERTIES PUBLIC_HEADER "${NTBM_LIB_INCLUDES}")

# install the library in the lib folder
install(TARGETS libNTBM
	LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
	PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCDIR}")

# install the ROOT dictionsry in the lib folder
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libNTBMDict_rdict.pcm DESTINATION "${CMAKE_INSTALL_LIBDIR}")
