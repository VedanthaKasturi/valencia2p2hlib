cmake_minimum_required(VERSION 3.12)
project(valencia2p2hlib CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find ROOT package
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

# Find GSL package
find_package(GSL REQUIRED)

# Add include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)

# Define source files
set(SOURCES
    main.cc
        src/MatrixDir.cc
        src/generatormt.cpp
        src/1p1h1pion.cpp
        src/hadron_tensor_2p2h.cc
        src/hadron_tensor_2p2h_2b.cpp
        src/hadron_tensor_2p2h_rho.cpp
        src/2p2h_MC.cpp
        src/3p3h_MC.cpp
        src/utils.cpp
        src/event.cpp
        src/event_3p3h.cpp
)

# Generate ROOT dictionaries
set(EVENT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/lib/event.h)
set(EVENT_LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/lib/eventLinkDef.h)
ROOT_GENERATE_DICTIONARY(eventdict ${EVENT_HEADERS} LINKDEF ${EVENT_LINKDEF})

set(EVENT_3P3H_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/lib/event_3p3h.h)
set(EVENT_3P3H_LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/lib/event_3p3hLinkDef.h)
ROOT_GENERATE_DICTIONARY(event_3p3hdict ${EVENT_3P3H_HEADERS} LINKDEF ${EVENT_3P3H_LINKDEF})

# Add the generated dictionary sources to the SOURCES list
list(APPEND SOURCES 
    ${CMAKE_CURRENT_BINARY_DIR}/eventdict.cxx
    ${CMAKE_CURRENT_BINARY_DIR}/event_3p3hdict.cxx
)

# Create main executable
add_executable(program ${SOURCES})

# Create xsec_batch executable (batch d(sigma)/dx^10 evaluator for NuGenFlow)
add_executable(xsec_batch
    xsec_batch.cc
        src/MatrixDir.cc
        src/generatormt.cpp
        src/1p1h1pion.cpp
        src/hadron_tensor_2p2h.cc
        src/hadron_tensor_2p2h_2b.cpp
        src/hadron_tensor_2p2h_rho.cpp
        src/2p2h_MC.cpp
        src/3p3h_MC.cpp
        src/utils.cpp
        src/event.cpp
        src/event_3p3h.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/eventdict.cxx
    ${CMAKE_CURRENT_BINARY_DIR}/event_3p3hdict.cxx
)

# Create myroot executable
add_executable(myroot
        src/event.cpp
        src/generatormt.cpp
    myroot.cpp
        src/event_3p3h.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/eventdict.cxx
    ${CMAKE_CURRENT_BINARY_DIR}/event_3p3hdict.cxx
)

# Set compile options
target_compile_options(program PRIVATE -Wall -g -O2)
target_compile_options(myroot PRIVATE -Wall -g -O2)
target_compile_options(xsec_batch PRIVATE -Wall -g -O2)

# Link libraries
target_link_libraries(program PRIVATE ${ROOT_LIBRARIES} GSL::gsl m)
target_link_libraries(myroot PRIVATE ${ROOT_LIBRARIES} GSL::gsl m)
target_link_libraries(xsec_batch PRIVATE ${ROOT_LIBRARIES} GSL::gsl m)

# Custom clean target
add_custom_target(clean-all
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_clean.cmake
)

# Generate a cmake_clean.cmake file
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/cmake_clean.cmake
"
file(REMOVE_RECURSE
    \${CMAKE_BINARY_DIR}/CMakeFiles
    \${CMAKE_BINARY_DIR}/CMakeCache.txt
    \${CMAKE_BINARY_DIR}/cmake_install.cmake
    \${CMAKE_BINARY_DIR}/Makefile
    \${CMAKE_BINARY_DIR}/program
    \${CMAKE_BINARY_DIR}/myroot
    \${CMAKE_BINARY_DIR}/xsec_batch
    \${CMAKE_BINARY_DIR}/eventdict*
    \${CMAKE_BINARY_DIR}/event_3p3hdict*
)
")