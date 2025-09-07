cmake_minimum_required(VERSION 3.28.3)

project(Filters)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


function(add_filter_target TARGET_NAME)
    add_executable(${TARGET_NAME} ${ARGN})

    target_compile_definitions(${TARGET_NAME} PRIVATE ${TARGET_NAME}_TARGET)

    configure_target_options(${TARGET_NAME})
endfunction()


function(configure_target_options TARGET_NAME)
    target_link_libraries(${TARGET_NAME}
        PRIVATE
        DataHandler
        Eigen3::Eigen
    )

    if (MEASUREMENT_UPDATE) 
        target_compile_definitions(${TARGET_NAME} PRIVATE MEASUREMENT_UPDATE)
    endif()

    if(ROBUST)
        target_compile_definitions(${TARGET_NAME} PRIVATE ROBUST)
    endif()

    if(COUPLED)
        target_compile_definitions(${TARGET_NAME} PRIVATE COUPLED)
    else()
        target_compile_definitions(${TARGET_NAME} PRIVATE DECOUPLED)
    endif()

    if (SAVE_INPUT)
        target_compile_definitions(${TARGET_NAME} PRIVATE SAVE_INPUT)
    endif()
endfunction()

find_package(Eigen3 REQUIRED)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/lib/DataHandler")

option(COUPLED "Use coupled correction step" OFF)
option(ROBUST "Use Robust correction step" OFF)
option (MEASUREMENT_UPDATE "Include measurement update step when filtering." ON)
option(SAVE_INPUT "Saves the input data" OFF)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(STATUS "Configuring filters as an executable.")

    # EKF
    add_filter_target(EKF 
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/ekf.cpp"
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/filter.cpp"
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp"
    )

    # IEKF
    add_filter_target(IEKF 
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/iekf.cpp"
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/ekf.cpp"
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/filter.cpp"
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp"
    )

    # Information Filter
    add_filter_target(INFO 
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/information_filter.cpp"
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/filter.cpp"
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp"
    )
else()
    message(STATUS "Configuring filters as a library.")
    set(LIBRARY_NAME Filters)

    add_library(${LIBRARY_NAME} STATIC 
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/iekf.cpp"
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/ekf.cpp"
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/information_filter.cpp"
	    "${CMAKE_CURRENT_SOURCE_DIR}/src/filter.cpp"
    )

    configure_target_options(${LIBRARY_NAME})

endif()
