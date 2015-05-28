
flexbar --reads test.csfastq --target result_right_nw --format csfastq --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-readlength 10 --adapter-threshold 1 --adapter-trim-end RIGHT > /dev/null

a=`diff correct_result_right_nw.csfastq result_right_nw.csfastq`

l1=`expr length "$a"`

if [ $l1 != 0 ]; then
echo "error testing mode csfastq, right, nw"
echo $a
exit -1
else
echo "Test 1 OK"
fi


flexbar --reads test.csfastq --target result_left_nw --format csfastq --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-readlength 10 --adapter-threshold 1 --adapter-trim-end LEFT > /dev/null

a=`diff correct_result_left_nw.csfastq result_left_nw.csfastq`

l1=`expr length "$a"`

if [ $l1 != 0 ]; then
echo "error testing mode cdfastq, left, nw"
echo $a
exit -1
else
echo "Test 2 OK"
fi


flexbar --reads test.csfastq --target result_any_nw --format csfastq --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-readlength 10 --adapter-threshold 1 --adapter-trim-end ANY > /dev/null

a=`diff correct_result_any_nw.csfastq result_any_nw.csfastq`

l1=`expr length "$a"`

if [ $l1 != 0 ]; then
echo "error testing mode any, left, nw"
echo $a
exit -1
else
echo "Test 3 OK"
fi


flexbar --reads test.csfastq --target result_left_tail_nw --format csfastq --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-readlength 10 --adapter-threshold 1 --adapter-trim-end LEFT_TAIL > /dev/null

a=`diff correct_result_left_tail_nw.csfastq result_left_tail_nw.csfastq`

l1=`expr length "$a"`

if [ $l1 != 0 ]; then
echo "error testing mode fasta, left_tail, nw"
echo $a
exit -1
else
echo "Test 4 OK"
fi


flexbar --reads test.csfastq --target result_right_tail_nw --format csfastq --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-readlength 10 --adapter-threshold 1 --adapter-trim-end RIGHT_TAIL > /dev/null

a=`diff correct_result_right_tail_nw.csfastq result_right_tail_nw.csfastq`

l1=`expr length "$a"`

if [ $l1 != 0 ]; then
echo "error testing mode fasta, right_tail, nw"
echo $a
exit -1
else
echo "Test 5 OK"
fi

echo ""

