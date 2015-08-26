./src/kfilter.py -v -R -o ./$1.filt.csv ./$1
#./kfilter.py -v -k 21 -R -o $1.filt.csv $1 
./src/tester.py -i ./$1.filt.csv -o ./$1.test.csv
./src/parser.py -i ./$1.test.csv -o ./$1.$2.primer.csv -l $2 -c ./Ebola_Gire_All_my.deleme.fasta
for i in ./gire_fasta/*.fasta; do echo "hello $i"; done 
rm -r ./$1.test.csv


