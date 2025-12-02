# CMakeList.txt : CMake project for PulsarAnalyser

cmake_minimum_required (VERSION 3.10)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("PulsarAnalyser")

# Use C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Main executables
add_executable(PulsarAnalyser
    PulsarAnalyser.cpp
    TimeFrequency.cpp
    TimeFrequency.h
)

# Test executable
add_executable(UnitTests
    UnitTests.cpp
    TimeFrequency.cpp
    TimeFrequency.h
)
