add_custom_command(OUTPUT TPISTAPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TPISTAPhysics.h TPISTAPhysicsDict.cxx TPISTAPhysics.rootmap libNPPISTA.dylib DEPENDS TPISTAPhysics.h)

add_custom_command(OUTPUT TPISTADataDict.cxx COMMAND ../../scripts/build_dict.sh TPISTAData.h TPISTADataDict.cxx TPISTAData.rootmap libNPPISTA.dylib DEPENDS TPISTAData.h)

add_custom_command(OUTPUT TFPMWDataDict.cxx COMMAND ../../scripts/build_dict.sh TFPMWData.h TFPMWDataDict.cxx TFPMWData.rootmap libNPPISTA.dylib DEPENDS TFPMWData.h)

add_custom_command(OUTPUT TICDataDict.cxx COMMAND ../../scripts/build_dict.sh TICData.h TICDataDict.cxx TICData.rootmap libNPPISTA.dylib DEPENDS TICData.h)

add_library(NPPISTA SHARED TPISTASpectra.cxx TPISTAData.cxx TPISTAPhysics.cxx TPISTADataDict.cxx TPISTAPhysicsDict.cxx TFPMWData.cxx TFPMWDataDict.cxx TICData.cxx TICDataDict.cxx)
target_link_libraries(NPPISTA ${ROOT_LIBRARIES} NPCore) 
install(FILES TPISTAData.h TPISTAPhysics.h TPISTASpectra.h TFPMWData.h TICData.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

