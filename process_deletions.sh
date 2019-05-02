#!/bin/bash

N=16
Nsq=256
counter=0
cutoff=0.80

echo -e "Chr\tStart\tEnd\tComparison\tn\tS\tcategory" > final_results.txt

for file in `ls *_deletions.bed`
do
	echo $file
	IFS=$'\n'
	for line in `cat $file`
	do
	    echo Now processing ${line}...
	    chr=`echo $line | awk '{print $1}'`
	    echo $line > tmp1.bed	
	    
	    bedtools intersect -a tmp1.bed -b tmpAnalyzed.bed -f 0.50 -r -wa -wb > tmpCheck.bed

	    if [[ ! -s tmpCheck.bed || "${counter}" == 0 ]]; then
		
		echo $line >> tmpAnalyzed.bed	
		counter=1

		grep -w "^${chr}" *_deletions.bed | awk -F':' '{print $2}' > tmp2.bed
				
		bedtools intersect -a tmp1.bed -b tmp2.bed -f 0.50 -r -wb > tmp3.bed

		n=`wc -l tmp3.bed | awk '{print $1}'`
		echo n $n
	
		x=`echo "scale=5; $n/$Nsq" | bc -l`		
#		echo x $x

		b=`echo "0.25 - $x" | bc -l`
#		echo b $b

		sqr=$(echo "scale=4; sqrt ( $b )" | bc -l)	
#		echo sqr $sqr
	
	        f=`echo "scale=4; 0.5 - $sqr" | bc -l`
#		echo f $f

		Nv=`echo "scale=4; $f*$N" | bc -l`
#		echo Nv $Nv

		for line2 in `cat tmp3.bed`
		do
			tumorSample=`echo ${line2} | awk '{print $NF}' | awk -F'___' '{print $1}'`
			echo ${tumorSample} >> tmpTumors.txt
		done

		maxNv=`sort tmpTumors.txt | uniq -c | sort -n | awk '{print $1}' | tail -n 1`
		echo maxNv $maxNv
		S=`echo "scale=4; $maxNv/$n" | bc -l`
		echo S $S

		if (( $(echo "$S >= $cutoff" | bc -l) )); then
			category="Somatic"
		else
			category="Germline"
		fi

		echo -e "${line}\t${n}\t${S}\t${category}" >> final_results.txt

		rm tmpTumors.txt
	    else
		echo "SV already analyzed.."
	    fi
	done
done
