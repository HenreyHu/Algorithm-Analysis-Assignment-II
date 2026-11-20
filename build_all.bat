@echo off
rem Windows shortcut for make. Double-click it, or run: build_all.bat [make target]
rem   no target, results/ complete : make           (build and run the 46 tests)
rem   no target, results/ incomplete: make results, then make mutation (recreates the whole results/ folder)
rem   a target (test, demo, bench, viz, results, mutation, ...): make <target>
rem All build logic lives in the Makefile.
setlocal
cd /d "%~dp0"

rem Opened by double-click (Explorer runs this file itself, outside CI): keep the window open at the end.
set KEEP_OPEN=0
if not defined CI echo %cmdcmdline% | find /i "%~f0" >nul && set KEEP_OPEN=1

where mingw32-make >nul 2>nul
if errorlevel 1 (
    echo mingw32-make not found. It ships with MinGW next to g++, so add the MinGW bin folder to PATH.
    set RC=1
    goto end
)

if not "%~1"=="" (
    call :run %1
    goto end
)

set MISSING=0
for %%F in (tests.csv bench_scaling.csv bench_carve.csv bench_flow.csv bench_improve.csv bench_sobel.csv seam_counts.csv mutation.csv) do if not exist "results\%%F" set MISSING=1
for %%I in (synthetic sailboat acacia broadway_tower starry_night) do if not exist "results\demo\%%I\summary.csv" set MISSING=1

if %MISSING%==0 (
    echo results/ is complete: running the tests only.
    call :run test
    goto end
)

echo results/ is missing or incomplete: recreating all of it. This takes several minutes.
call :run results
if not %RC%==0 goto end
call :run mutation
goto end

:run
echo.
echo ===== make %1 =====
mingw32-make %1
set RC=%errorlevel%
if %RC%==0 (echo DONE: make %1 succeeded) else (echo FAILED: make %1 returned %RC%)
exit /b %RC%

:end
if %KEEP_OPEN%==1 pause
exit /b %RC%
