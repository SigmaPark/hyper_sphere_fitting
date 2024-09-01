::	SPDX-FileCopyrightText: (c) 2024 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
::	SPDX-License-Identifier: MIT License
::========::========::========::========::=======#::========::========::========::========::=======#

@echo off

powershell rm -recurse -force ./build

mkdir build
cd build
	
cmake .. 
cmake --build . --config Release

cd ..

pause