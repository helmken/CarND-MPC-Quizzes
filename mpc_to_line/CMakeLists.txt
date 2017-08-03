project(mpc_to_line)

cmake_minimum_required (VERSION 3.5)

# -g allows for gdb debugging
# turn on -03 for best performance
add_definitions(-std=c++11 -DHAVE_CSTDDEF)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

set(sources
	src/FG_eval.h
	src/main.cpp
	src/matplotlibcpp.h
	src/MPC.cpp
	src/MPC.h
)

include_directories(src/Eigen-3.3)
include_directories(/usr/local/include)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
include_directories(../../../cppad-20170722.epl)
endif()

link_directories(/usr/local/lib)



add_executable(mpc_to_line 
	${sources}
)

message(${CMAKE_SYSTEM_NAME})

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
target_link_libraries(mpc_to_line 
    ipopt 
    python3.5m
)
else()
target_link_libraries(mpc_to_line 
    ipopt 
    python
)
endif()
