
awk 'BEGIN{OFS="\t"} {$0=$1"\t0\t"substr($0,index($0,$2))} 1' mm10.chrom.sizes | sed '1i Chromosome\tStart\tEnd' - > chromsizes.tsv
