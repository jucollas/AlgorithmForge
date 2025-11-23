@echo off

FOR /R %%f IN (*.in) DO (
	echo %%f
	atcoder < %%f > atres.out
	montgom < %%f > montgom.out
	montatc < %%f > montatc.out

	mfc atres.out montgom.out
	mfc montatc.out atres.out
)

echo on