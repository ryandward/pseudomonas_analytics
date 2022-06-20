export PATH=$PATH:/home/glbrc.org/ryan.d.ward/bin/bbmap/

for i in *R1*gz; do
	echo seal.sh overwrite=t in=$i k=20 mm=f ambiguous=toss outm=stdout.fa hdist=1 trd=t rename fbm int=f ref=full_annotations.fasta \| awk \''$0~">"'\' \| sed \''s/>//g'\' \|  awk \''{winner = "Unknown"; score = 0; for (i = 2; i <= NF; i++) {split ($i, parts, "="); if (parts[2] > score) {score = parts[2]; winner = parts[1]; }} print $1, winner; delete parts; }'\' \> ${i%%_S*_R1_001.fastq.gz}.guides.tsv
done > map_it_all.sh;

for i in *R2*gz; do
  echo seal.sh overwrite=t in=$i k=20 mm=f ambiguous=toss outm=stdout.fa hdist=1 trd=t rename fbm int=f ref=PromotersOnly.fasta \| awk \''$0~">"'\' \| sed \''s/>//g'\' \|  awk \''{winner = "Unknown";	score = 0; for (i = 2; i <= NF; i++) {split ($i, parts, "="); if (parts[2] > score) {score = parts[2]; winner = parts[1]; }} print $1, winner; delete parts; }'\' \> ${i%%_S*_R2_001.fastq.gz}.promoters.tsv
done >> map_it_all.sh;

parallel -j8 'echo {} | sh' :::: map_it_all.sh 

for i in *guides.tsv; do
  echo awk -vOFS='"\t"' \''FILENAME ~ "guides" {guides[$1] = $2} FILENAME ~ "promoters"  && $1 in guides{print $2, guides[$1]}'\' $i ${i%%guides.tsv}promoters.tsv \| sort \| uniq -c \| sed -e \''s/^ *//;s/ /\t/'\'  \| awk -vOFS=\'\\t\' \'{print \$2, \$3, \$1}\' \> ${i%%guides.tsv}context.tsv
done > contextify.sh

parallel -j8 'echo {} | sh' :::: contextify.sh

for i in *context*tsv; do awk '{name=FILENAME; gsub(".context.tsv", "", name); print $0, name }' $i; done > all_counts_seal.tsv
