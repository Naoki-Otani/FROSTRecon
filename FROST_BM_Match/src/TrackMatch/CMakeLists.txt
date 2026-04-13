message (STATUS "TrackMatch...")

add_executable(TrackMatch
	TrackMatch.cpp
	TrackMatch.hpp)

target_link_libraries(TrackMatch
	${ROOT_LIBRARIES}
	${Geant4_LIBRARIES}
	${Boost_LIBRARIES}
	${B2MC_LIBRARY}
	libNTBM
)

# Make TrackMatch find libNTBM.so at runtime
# while running from the build tree and after installation.
set_target_properties(TrackMatch PROPERTIES
	BUILD_RPATH "$ORIGIN/../../lib/ninja/recon"
	INSTALL_RPATH "$ORIGIN/../../lib/ninja/recon"
	INSTALL_RPATH_USE_LINK_PATH TRUE
)

# install the execute in the bin folder
install(TARGETS TrackMatch DESTINATION "${CMAKE_INSTALL_BINDIR}/TrackMatch")
