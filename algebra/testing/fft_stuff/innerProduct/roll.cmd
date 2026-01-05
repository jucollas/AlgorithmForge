@echo off


for /l %%i in (1,1,100) do (
	echo %%i
	gen.exe %1 > input.txt
	cntrl.exe < input.txt > rcntrl.txt
	fconv.exe < input.txt > rfconv.txt
	REM echo "xd"
	REM brute.exe < input.txt > rbrute.txt
	REM fc rcntrl.txt rbrute.txt || goto out
	fc rcntrl.txt rfconv.txt || goto out
)
:out