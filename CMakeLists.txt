cmake_minimum_required(VERSION 3.10)

# set the project name
project(CHASE)

# bring headers into project
include_directories(include)
# set sources
file(GLOB SOURCES "src/*.cpp")

# add the executable
add_executable(CHASE ${SOURCES})
target_link_directories(CHASE PUBLIC lib)

# link library needed for the vertex parsing
target_link_libraries(CHASE pthread)
target_link_libraries(CHASE routingkit)
target_link_libraries(CHASE boost_filesystem)
target_link_libraries(CHASE boost_program_options)

# specify c++ standard
set_property(TARGET CHASE PROPERTY CXX_STANDARD 17)