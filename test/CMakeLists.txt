# Unit test configuration for Millepede 

add_executable(testMilleC testMilleC.cc)
target_link_libraries(testMilleC PUBLIC Mille readc)    
install (
        TARGETS testMilleC
        DESTINATION bin
    )
add_test(NAME IOfromCbinary COMMAND testMilleC)

if (build_ROOT_comps)
add_executable(testMilleROOT testMilleROOT.cc)
target_link_libraries(testMilleROOT PUBLIC Mille readROOT)    
install (
        TARGETS testMilleROOT
        DESTINATION bin
    )
add_test(NAME IOfromROOTbinary COMMAND testMilleROOT)
endif()

if (build_ROOT_comps)
add_executable(testCtoROOT testCtoROOT.cc)
target_link_libraries(testCtoROOT PUBLIC Mille readROOT)    
install (
        TARGETS testCtoROOT
        DESTINATION bin
    )
add_test(NAME IOviaCtoROOT COMMAND testCtoROOT)
endif()
