#   SPDX-FileCopyrightText: (c) 2024 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
#   SPDX-License-Identifier: MIT License
#=========#=========#=========#=========#=========#=========#=========#=========#=========#=========


cmake_minimum_required(VERSION 3.0)

project(Test_Hyper_Sphere_Fitting_Methods)

set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../include)
set(SIGMA_3D_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../sigma_3D/include)
set(SIGMA_SRC_DIR ${SIGMA_3D_SRC_DIR}/../sigma/include)
set(EIGEN_SRC_DIR ${SIGMA_3D_SRC_DIR}/../eigen)

file(
	GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS
	${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/*.hpp
	${INCLUDE_DIR}/*.hpp
)

add_executable(
    Test_Hyper_Sphere_Fitting_Methods
    ${SRC_FILES}
)

target_include_directories(
    Test_Hyper_Sphere_Fitting_Methods
    PRIVATE
	${INCLUDE_DIR}
	${SIGMA_3D_SRC_DIR}
	${SIGMA_SRC_DIR}
	${EIGEN_SRC_DIR}
)

target_compile_definitions(
    Test_Hyper_Sphere_Fitting_Methods
    PRIVATE
    SGM_CXX_STANDARD=2017
)

target_compile_features(
	Test_Hyper_Sphere_Fitting_Methods PRIVATE 
	cxx_std_17
)