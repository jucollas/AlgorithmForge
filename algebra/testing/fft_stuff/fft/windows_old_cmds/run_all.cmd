@echo off

FOR /R %%f IN (*.in) DO (
	echo %%f
	atcoder < %%f > atres.out
	montgom < %%f > montgom.out
	montatc < %%f > montatc.out
	n_rad4  < %%f > n_rad4.out

	mfc atres.out montgom.out
	mfc montatc.out atres.out
	mfc n_rad4.out
)

echo on