::	SPDX-FileCopyrightText: (c) 2024 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
::	SPDX-License-Identifier: MIT License
::========::========::========::========::=======#::========::========::========::========::=======#


powershell rm -recurse -force ./build

mkdir build
cd build
	
cmake .. 

cd ..

pause
