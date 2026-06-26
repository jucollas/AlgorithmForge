# https://codeforces.com/blog/entry/102287
# linux
for i in {1..100}; do
    echo "on test $i"
    ./gen.out > input.txt
    ./a.out < input.txt > output.txt
    if [ $? -ne 0 ]; then
    	echo "Something wrong"
    	cat input.txt
    	cat output.txt
    	exit 1
    fi
    #./naive.out < input.txt > answer.txt
    #diff output.txt answer.txt || break
done
echo "All ok"