# Unit test configuration for Millepede 

add_executable(testMilleC testMilleC.cc)
target_link_libraries(testMilleC PUBLIC MilleC)    
install (
        TARGETS testMilleC
        DESTINATION bin
    )
add_test(NAME TestMilleC COMMAND testMilleC)

if (build_ROOT_comps)
add_executable(testMilleROOT testMilleROOT.cc)
target_link_libraries(testMilleROOT PUBLIC MilleC)    
install (
        TARGETS testMilleROOT
        DESTINATION bin
    )
add_test(NAME TestMilleROOT COMMAND testMilleROOT)
endif()

if (build_ROOT_comps)
add_executable(testConvertMilleFormats testConvertMilleFormats.cc)
target_link_libraries(testConvertMilleFormats PUBLIC MilleC)    
install (
        TARGETS testConvertMilleFormats
        DESTINATION bin
    )
add_test(NAME TestConvertMilleFormats COMMAND testConvertMilleFormats)
endif()
