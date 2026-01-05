:: https://codeforces.com/blog/entry/102287

@echo off

for /l %%i in (1, 1, 100) do (
    echo %%i

    gen.exe > input.txt
    a.exe < input.txt > output.txt
    stl.exe < input.txt > answer.txt

    fc output.txt answer.txt || goto :out
)

:out
echo on

:: linux
::for i in `seq 1 100`; do
::    echo $i 

::    ./ferris_gen $i 4 5 > input.txt
::    ./ferris < input.txt > output.txt
::    ./ferris_naive < input.txt > answer.txt

::    diff output.txt answer.txt || break
::done
