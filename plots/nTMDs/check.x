baseout=`cat baseline.txt`

myjobs > out.txt

sed -n '2p' out.txt > newout.txt

jobsout=`cat newout.txt`

if [ "$baseout" = "$jobsout" ]; then
    echo "True"
else
    echo "False"
fi
