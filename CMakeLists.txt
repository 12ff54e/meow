cmake_minimum_required(VERSION 3.20)

project(meow VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_EXTENSIONS OFF)

option(MEOW_BUILD_TESTS "Build the meow test suite" ON)
option(MEOW_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" ON)
option(MEOW_BUILD_CUMES_INTEGRATION
       "Build optimizer targets and examples backed by cuMES" OFF)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_library(meow_trf src/trf.cpp)
add_library(meow::trf ALIAS meow_trf)

target_compile_features(meow_trf PUBLIC cxx_std_20)
target_link_libraries(meow_trf PUBLIC Eigen3::Eigen)
target_include_directories(
    meow_trf
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)

function(meow_enable_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4)
        if(MEOW_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE /WX)
        endif()
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
        if(MEOW_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE -Werror)
        endif()
    endif()
endfunction()

meow_enable_warnings(meow_trf)

if(MEOW_BUILD_CUMES_INTEGRATION)
    if(NOT TARGET cumes::solver)
        find_package(cuMES CONFIG REQUIRED)
    endif()

    add_library(meow_cumes_targets INTERFACE)
    add_library(meow::cumes_targets ALIAS meow_cumes_targets)
    target_compile_features(meow_cumes_targets INTERFACE cxx_std_20)
    target_include_directories(
        meow_cumes_targets
        INTERFACE
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>)
    target_link_libraries(meow_cumes_targets INTERFACE cumes::solver meow::trf)

    add_executable(cumes_meow_optimize examples/cumes_optimize.cpp)
    target_link_libraries(cumes_meow_optimize PRIVATE meow::cumes_targets)
    meow_enable_warnings(cumes_meow_optimize)

    add_executable(cumes_meow_qs_optimize examples/cumes_qs_optimize.cpp)
    target_link_libraries(cumes_meow_qs_optimize PRIVATE meow::cumes_targets)
    meow_enable_warnings(cumes_meow_qs_optimize)
endif()

if(MEOW_BUILD_TESTS)
    enable_testing()

    add_executable(test_trf tests/test_trf.cpp)
    target_link_libraries(test_trf PRIVATE meow::trf)
    meow_enable_warnings(test_trf)
    add_test(NAME trf_unit COMMAND test_trf)

    if(MEOW_BUILD_CUMES_INTEGRATION)
        foreach(target IN ITEMS plasma_size magnetic_gradient quasisymmetry)
            add_executable(test_${target}_target tests/test_${target}_target.cpp)
            target_include_directories(test_${target}_target PRIVATE
                                       ${CMAKE_CURRENT_SOURCE_DIR}/tests)
            target_link_libraries(test_${target}_target PRIVATE
                                  meow::cumes_targets)
            meow_enable_warnings(test_${target}_target)
            add_test(NAME ${target}_target_unit COMMAND test_${target}_target)
        endforeach()
    endif()

    add_executable(trf_reference_cases tests/reference_cases.cpp)
    target_link_libraries(trf_reference_cases PRIVATE meow::trf)
    meow_enable_warnings(trf_reference_cases)

    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        add_test(
            NAME vmec_namelist_converter_unit
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_vmec_namelist_converter.py)
        execute_process(
            COMMAND ${Python3_EXECUTABLE} -c "import scipy"
            RESULT_VARIABLE MEOW_SCIPY_RESULT
            OUTPUT_QUIET
            ERROR_QUIET)
        if(MEOW_SCIPY_RESULT EQUAL 0)
            add_test(
                NAME trf_scipy_comparison
                COMMAND ${Python3_EXECUTABLE}
                        ${CMAKE_CURRENT_SOURCE_DIR}/tests/verify_scipy.py
                        $<TARGET_FILE:trf_reference_cases>)
        endif()
    endif()
endif()

include(GNUInstallDirs)
install(TARGETS meow_trf EXPORT meowTargets)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(
    EXPORT meowTargets
    NAMESPACE meow::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/meow)
