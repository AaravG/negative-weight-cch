@echo off
rem Build cch.exe with MSVC. Adjust the path if Visual Studio lives elsewhere.
call "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cl /nologo /O2 /std:c++20 /EHsc /openmp /Fe:cch.exe cch.cpp
