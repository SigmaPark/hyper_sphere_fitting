::	SPDX-FileCopyrightText: (c) 2024 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
::	SPDX-License-Identifier: MIT License
::========::========::========::========::=======#::========::========::========::========::=======#


powershell rm --recurse ./build

mkdir build
cd build
	
cmake .. -G "Visual Studio 16 2019"

cd ..

pause
