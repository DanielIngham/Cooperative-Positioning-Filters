cmake_minimum_required(VERSION 3.28.3)

project(Filters)

if(NOT CMAKE_BUILD_TYPE)
	set(default_build_type "Release")
	message(STATUS "Set the build type to `${default_build_type}` as none was specified.")
	set(CMAKE_BUILD_TYPE ${default_build_type} CACHE STRING "Choose the build type." FORCE)
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
	"Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "${PROJECT_NAME} Build Type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)
find_package(yaml-cpp REQUIRED)
message(STATUS "Used config file: ${yaml-cpp_CONFIG}")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/lib/DataHandler")

option(RANGE_ONLY "Use bearing only update" OFF)
option(DECOUPLED "Use coupled correction step" OFF)
option(COUPLED "Use coupled correction step" ON)
option(ROBUST "Use Robust correction step" OFF)
option(MEASUREMENT_UPDATE "Include measurement update step when filtering." ON)
option(SAVE_INPUT "Saves the input data" OFF)

set (AGENT_SRC 
	"${CMAKE_CURRENT_SOURCE_DIR}/src/agent/agent.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/agent/adversarial_landmark.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/agent/landmark.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/agent/robot.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/agent/faulty_robot.cpp"
)

set (COMMON_SRC
	"${CMAKE_CURRENT_SOURCE_DIR}/src/common/config.cpp"
)


set (MODEL_SRC
	"${CMAKE_CURRENT_SOURCE_DIR}/src/models/process.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/models/measurement.cpp"
)

set (SENSOR_SRC 
    "${CMAKE_CURRENT_SOURCE_DIR}/src/sensors/odometry.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/sensors/odom_data.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/sensors/meas_data.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/sensors/meas_set.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/sensors/measurements.cpp"
)

set (UTILS_SRC
	"${CMAKE_CURRENT_SOURCE_DIR}/src/utils/matrix_operations.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/utils/performance_eval.cpp"
)

set (PROJ_SRC 
	${AGENT_SRC}
	${COMMON_SRC}
	${SENSOR_SRC}
	${MODEL_SRC}
	${UTILS_SRC}
)

function(add_filter_target TARGET_NAME)
    add_executable(${TARGET_NAME} 
	${ARGN}
	${PROJ_SRC}
    )

    target_compile_definitions(${TARGET_NAME} PRIVATE ${TARGET_NAME}_TARGET)

    configure_target_options(${TARGET_NAME})
endfunction()


function(configure_target_options TARGET_NAME)
    target_include_directories(${TARGET_NAME}
	PUBLIC
	"${CMAKE_CURRENT_SOURCE_DIR}/include"
    )

    target_link_libraries(${TARGET_NAME}
	PUBLIC
        DataHandler
        Eigen3::Eigen
	yaml-cpp::yaml-cpp
    )

    if (MEASUREMENT_UPDATE) 
        target_compile_definitions(${TARGET_NAME} PRIVATE MEASUREMENT_UPDATE)
    endif()

    if(ROBUST)
        target_compile_definitions(${TARGET_NAME} PRIVATE ROBUST)
    endif()

    if(COUPLED)
        target_compile_definitions(${TARGET_NAME} PRIVATE COUPLED)
    elseif(DECOUPLED)
        target_compile_definitions(${TARGET_NAME} PRIVATE DECOUPLED)
    else()
	target_compile_definitions(${TARGET_NAME} PRIVATE RANGE_ONLY)
    endif()

    if (SAVE_INPUT)
        target_compile_definitions(${TARGET_NAME} PRIVATE SAVE_INPUT)
    endif()
endfunction()


if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(STATUS "Configuring filters as an executable.")

    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	set(MAIN_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/src/debug.cpp")
    else()
	set(MAIN_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
    endif()

    # EKF
    add_filter_target(EKF 
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/filters/ekf.cpp"
	    ${MAIN_SOURCE}
    )

    # IEKF
    add_filter_target(IEKF 
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/filters/iekf.cpp"
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/filters/ekf.cpp"
	    ${MAIN_SOURCE}
    )

    # Converted measurement 
    add_filter_target(CMEKF 
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/filters/cmekf.cpp"
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/filters/ekf.cpp"
	    ${MAIN_SOURCE}
    )

    # Information Filter
    add_filter_target(INFO 
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/filters/information_filter.cpp"
	    ${MAIN_SOURCE}
    )

    # Information Filter
    add_filter_target(PARTICLE 
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/filters/particle.cpp"
	    ${MAIN_SOURCE}
    )

else()
    message(STATUS "Configuring filters as a library.")
    set(LIBRARY_NAME Filters)

    add_library(${LIBRARY_NAME} STATIC 
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/filters/iekf.cpp"
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/filters/ekf.cpp"
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/filters/cmekf.cpp"
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/filters/information_filter.cpp"
	    ${PROJ_SRC}
    )

    configure_target_options(${LIBRARY_NAME})

endif()
