cmake_minimum_required(VERSION 3.28.3)

project(Filters)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/lib/DataHandler")

link_libraries(
	DataHandler
	Eigen3::Eigen
)

# EKF
add_executable(EKF 
	"${CMAKE_CURRENT_SOURCE_DIR}/src/ekf.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/filter.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp"
)

target_compile_definitions(EKF PUBLIC 
	EKF_TARGET
	DECOUPLED
)

# IEKF
add_executable(IEKF 
	"${CMAKE_CURRENT_SOURCE_DIR}/src/iekf.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/ekf.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/filter.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp"
)

target_compile_definitions(IEKF PUBLIC 
	IEKF_TARGET
	DECOUPLED
)

# Information Filter
add_executable(INFO 
	"${CMAKE_CURRENT_SOURCE_DIR}/src/information_filter.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/filter.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp"
)

target_compile_definitions(INFO PUBLIC 
	INFO_TARGET
	DECOUPLED
)



