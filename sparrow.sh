echo "Generating kmers of length $3..."
./src/kfilter.py -v -R -k $3 -o ./$1.filt.csv ./$1
#./kfilter.py -v -k 21 -R -o $1.filt.csv $1 

echo "Running primer quality tests..."
./src/tester.py -i ./$1.filt.csv -o ./$1.test.csv

echo "Scoring primers..."
./src/parser.py -i ./$1.test.csv -o ./$1.$2.primer.csv -l $2 -c ./Ebola_Gire_All_my.deleme.fasta

#for i in ./gire_fasta/*.fasta; do echo "hello $i"; done 
