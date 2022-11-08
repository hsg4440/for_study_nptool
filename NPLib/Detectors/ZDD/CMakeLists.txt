add_custom_command(OUTPUT TZDDPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TZDDPhysics.h TZDDPhysicsDict.cxx TZDDPhysics.rootmap libNPZDD.dylib DEPENDS TZDDPhysics.h)
add_custom_command(OUTPUT TZDDDataDict.cxx COMMAND ../../scripts/build_dict.sh TZDDData.h TZDDDataDict.cxx TZDDData.rootmap libNPZDD.dylib DEPENDS TZDDData.h)
add_library(NPZDD SHARED TZDDSpectra.cxx TZDDData.cxx TZDDPhysics.cxx TZDDDataDict.cxx TZDDPhysicsDict.cxx )
target_link_libraries(NPZDD ${ROOT_LIBRARIES} NPCore) 
install(FILES TZDDData.h TZDDPhysics.h TZDDSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

