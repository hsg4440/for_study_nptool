cmake_minimum_required (VERSION 2.8) 
# Setting the policy to match Cmake version
cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
## write compile_commands.json to clangd
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# include the default NPAnalysis cmake file
include("../../NPLib/ressources/CMake/NPAnalysis.cmake")
