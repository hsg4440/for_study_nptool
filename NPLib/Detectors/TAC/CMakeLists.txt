add_custom_command(OUTPUT TTACPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TTACPhysics.h TTACPhysicsDict.cxx TTACPhysics.rootmap libNPTAC.dylib DEPENDS TTACPhysics.h)
add_custom_command(OUTPUT TTACDataDict.cxx COMMAND ../../scripts/build_dict.sh TTACData.h TTACDataDict.cxx TTACData.rootmap libNPTAC.dylib DEPENDS TTACData.h)
add_library(NPTAC SHARED TTACSpectra.cxx TTACData.cxx TTACPhysics.cxx TTACDataDict.cxx TTACPhysicsDict.cxx )
target_link_libraries(NPTAC ${ROOT_LIBRARIES} NPCore) 
install(FILES TTACData.h TTACPhysics.h TTACSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

