find ./ -size 0 -print -delete

for i in *(guides|promoters).tsv; do LANG=en_EN sort -k 1,1 -o $i $i --parallel 12; done

for i in *promoters.tsv; do LANG=en_EN join -t " " -e"NA" -o auto -11 -21 -a1 -a2 ${i%%promoters.tsv}guides.tsv $i |  awk -v strain=${i%%.promoters.tsv} 'BEGIN{FS = " "; OFS = "\t";}{print strain, $2, $3}'> ${i%%promoters.tsv}context.tsv; done

for i in *context.tsv; do sort $i --parallel 12 | uniq -c | sed -e 's/^ *//;s/ /       /' | awk 'BEGIN{FS = " "; OFS = "\t"}{print $2, $3, "P1", $1}' > ${i%%context.tsv}counts.tsv; done

######
# Jasons Samples
for i in d*R1*gz; do
	fwd=$i;
	rev=$(echo $i|sed 's/\R1/R2/g')
	echo seal.sh overwrite=t in1=$fwd in2=$rev k=20 mm=f ambiguous=toss outm=stdout.fa hdist=1 trd=t rename fbm int=f ref=full_annotations.fasta \| awk \''$0~">"'\' \| sed \''s/>//g'\' \|  awk \''{winner = "Unknown"; score = 0; for (i=2; i<=NF; i++) {split ($i, parts, "="); if(parts[2]>score) {score=parts[2];winner=parts[1]}} print $1, winner; delete parts}'\' \> ${i%%_S*_R1_001.fastq.gz}.guides.tsv;
done > map_it_all2.sh;

awk 'BEGIN{FS = " "}{print FILENAME, $2}' dJMP*guides*tsv | sed 's/\.guides\.tsv//g' | sort --parallel 12 | uniq -c | sed -e 's/^ *//;s/ /       /' | awk 'BEGIN{FS = " "}{print $2, $3, "P1", $1/2}'> all_counts_seal2.tsv &
