#leave one out test

GL=1 #GA line

for gene in `cat gene.txt`
do

	num=1
	seq=`grep -c '>' ${gene}.faa`
	
	for name  in {1..'$seq'}
	do
		#extract hit fa
		grep '>' ${gene}.faa | awk '{print $1}' | sed 's/>//g'> tmp.ID
		cp tmp.ID tmp.1
		sed -i ${num}d tmp.1 #leave one
		seqkit grep -f tmp.1 ${gene}.faa > tmp.fa #extract seq according to ID

		#build db
		mafft --maxiterate 1000 --thread 10 --quiet tmp.fa > tmp_mafft.fa
		python covert.py #covert to sto format
		hmmbuild tmp.hmm tmp.stockholm

		#define GA TC NC
		GA="GA	`sed -n ''$GL'p' GA.txt`;"
		TC="TC	`sed -n ''$GL'p' GA.txt`;"
		NC="NC	`sed -n ''$GL'p' GA.txt`;"

		# add GA to the hmm profile
		sed "14i ${GA}" tmp.hmm > tmp.1.hmm
		sed -i "15i ${TC}" tmp.1.hmm
		sed -i "16i ${NC}" tmp.1.hmm

		#hmmsearch
		hmmsearch --cut_ga tmp.1.hmm test.fa > tmp.out

		#extract hit ID	
		sed -n "/full sequence/,/Domain annotation/p" tmp.out > tmp.2 #print hit result
		sed -i '1,3d' tmp.2
		sed -i '$d' tmp.2
		cat tmp.2| awk '{print $9}'|grep -v '^$' > tmp.3 #extract hit ID

		#Recall and Precission
		grep -f tmp.3 tmp.ID > tmp.4 #find true postive hit
		AP=`cat tmp.ID |wc -l` #All positive seq
		AH=`cat tmp.3 |wc -l` #All hit
		PH=`cat tmp.4 |wc -l` #Positive hit
		PHU=$[$PH*100]
		Recall=`expr $PHU / $AP`
		Precision=`expr $PHU / $AH`

		#output
		ID=`sed -n ${num}p tmp.ID`
		echo "${ID} ${Recall} ${Precision}" >> ${gene}.LOO.txt

		num=$[$num+1]
	done
	
	GL=$[$GL+1]
done

#rm tmp.*
