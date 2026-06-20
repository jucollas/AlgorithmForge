./a.out > inter.txt
echo "second" >> inter2.txt
cat inter.txt >> inter2.txt

echo "######################################### First run"
cat inter.txt
echo "######################################### Second run"

./a.out < inter2.txt
rm inter2.txt