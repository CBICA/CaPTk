@echo off

REM batch file example of setting %PATH% to Matlab Runtime Components, then
REM calling LIBRA.exe
REM
REM The Windows variable "%~dp0" refers to the directory in which the batch
REM file itself is stored. The location of the Matlab Runtime Components
REM is known, relative to the %~dp0 directory.

set curdrive=%CD:~0,2%
set curdir=%~dp0
%curdrive%
cd %curdir%

PATH=%PATH%;%~dp0\Mathworks\MCR_R2013A_Win64_reduced\v81\runtime\win64
.\LIBRA.exe
