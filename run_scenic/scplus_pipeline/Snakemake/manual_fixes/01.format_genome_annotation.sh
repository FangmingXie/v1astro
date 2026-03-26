
awk 'BEGIN{OFS="\t"} NR==1{print} NR>1{$1="chr"$1; print}' default_genome_annotation.tsv > genome_annotation.tsv
