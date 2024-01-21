#   SPDX-FileCopyrightText: (c) 2024 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
#   SPDX-License-Identifier: MIT License
#=========#=========#=========#=========#=========#=========#=========#=========#=========#=========


cmake_minimum_required(VERSION 3.0)

project(Sphere_Fitting)

set(CMAKE_CXX_STANDARD 17)

add_executable(
    Sphere_Fitting
    main.cpp
)

set(EIGEN_LIBRARY_DIR $ENV{EIGEN_DIR})

target_include_directories(
    Sphere_Fitting
    PRIVATE
    ${EIGEN_LIBRARY_DIR}
    sigma/include
)