REM https://codeforces.com/blog/entry/102287
REM WINDOWS
@echo off
for /l %%i in (1, 1, 100) do (
    echo %%i

    gen.exe > input.txt
    a.exe < input.txt > output.txt
    stl.exe < input.txt > answer.txt
    fc output.txt answer.txt || goto :out
) :out
echo on