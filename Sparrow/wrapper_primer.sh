./kfilter.py -v -R -o src/Ebola.filt.csv src/Ebola_reference.fasta
#./kfilter.py -v -k 21 -R -o $1.filt.csv $1 
./tester.py -i ./$1.csv -o ./$1.test.csv
./parser.py -i ./$1.test.csv -o ./$1.primer.csv -l $2
rm -r ./*test.csv


