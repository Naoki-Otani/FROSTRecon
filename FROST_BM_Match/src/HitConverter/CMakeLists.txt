message(STATUS "HitConverter...")

ROOT_GENERATE_DICTIONARY(G__HitConverterDict
    ${CMAKE_CURRENT_SOURCE_DIR}/HitConverterDict.hh
    LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/HitConverterLinkDef.h
)

add_executable(HitConverter
    HitConverter.cc
    G__HitConverterDict.cxx
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
