mkdir -p tests
mkdir -p tests/out
for file in tests/*; do
    if [ -f "$file" ]; then
        echo "Processing: $file"
        # Add your commands here
        ./atcoder.out < $file > tests/out/atres.ans
        ./montgom.out < $file > tests/out/montgom.ans
        ./montatc.out < $file > tests/out/montatc.ans
        #./n_rad4.out < $file > atres.ans
        ./weird_rec.out < $file > tests/out/wrec.ans
        ./weird_rec_montg.out < $file > tests/out/wmrec.ans

        diff tests/out/atres.ans tests/out/montgom.ans
        diff tests/out/montatc.ans tests/out/atres.ans
        diff tests/out/wrec.ans tests/out/atres.ans
        diff tests/out/wrec.ans tests/out/wmrec.ans
    fi
done