@echo off
setlocal enabledelayedexpansion

set RESULTENTRY=%1
set PYMOLSCRIPTLIBFILE=%2
set PYMOLINPUTFILE=%3
REM "pymol_input.py"
REM set RESULTENTRY=A0PJX2_Pos155

REM echo %RESULTENTRY%
REM goto :end

REM findstr /C:A0PJX2_Pos155 test_pymol.txt > %PYMOLSCRIPTFILE%
Rem Extract Line from the pymol script library file
for /f "delims=" %%i in ('findstr /C:%RESULTENTRY% %PYMOLSCRIPTLIBFILE%') do set t=%%i

REM set /p t=<%PYMOLSCRIPTFILE%
REM set t=%t:>=:%
REM echo %t%

Rem Insert comment
echo #LAUNCH PYMOL TO VIEW ALGORITHM RESULTS FOR ENTRY %RESULTENTRY% > %PYMOLINPUTFILE%

Rem Remove line identifier (resultentry)
for /f "tokens=1,* delims= " %%a in ("!t!") do ( set t=%%b )

Rem Split to lines ";"
:loop
for /f "tokens=1,* delims=;" %%a in ("!t!") do (
   REM if defined t echo %%a >> %PYMOLINPUTFILE%
   REM echo A='%%a'
   set line=%%a
   set t=%%b
   )

REM if defined t echo A=!line!
if defined t echo !line! >> %PYMOLINPUTFILE%
if defined t goto :loop

Rem Launch pymol
PymolWin.exe -d "@%PYMOLINPUTFILE%"
