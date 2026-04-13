message(STATUS "HitConverter...")

add_executable(HitConverter
    HitConverter.cc
)

target_compile_features(HitConverter PUBLIC cxx_std_17)

target_link_libraries(HitConverter
    ${ROOT_LIBRARIES}
    ${Geant4_LIBRARIES}
    ${Boost_LIBRARIES}
    ${B2MC_LIBRARY}
    libNTBM
)

# install the executable in the HitConverter bin folder
install(TARGETS HitConverter DESTINATION "${CMAKE_INSTALL_BINDIR}/HitConverter")
