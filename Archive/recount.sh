cd /mnt/bigdata/linuxhome/ryan.d.ward/Integrated_Mouse_Libraries/Round1
PATH=$PATH:/home/glbrc.org/ryan.d.ward/bin/bbmap/

for i in *R1*gz; do
	fwd=$i;
	rev=$(echo $i|sed 's/\R1/R2/g')
	echo seal.sh overwrite=t in=$fwd k=20 mm=f ambiguous=toss outm=stdout.fa hdist=1 trd=t rename fbm int=f ref=../GuidesOnly.fasta \| awk \''$0~">"'\' \| sed \''s/>//g'\' \|  awk \''{winner = "Unknown"; score = 0; for (i=2; i<=NF; i++) {split ($i, parts, "="); if(parts[2]>score) {score=parts[2];winner=parts[1]}} print $1, winner; delete parts}'\' \> ${i%%_S*_R1_001.fastq.gz}.guides.tsv;done > map_it_all.sh;

parallel -j8 'echo {} | sh' :::: map_it_all.sh  &

for i in *tsv; do
	awk 'BEGIN{FS = " "; OFS = " ";}{print FILENAME, $2}' $i |
	sort --parallel 12 -o $i.sorted;
done

cat *guides.tsv.sorted | sed 's/\.guides\.tsv//g'  | uniq -c | sed -e 's/^ *//;s/ /\t/g'  > all_counts_seal_round1_recount.tsv &
