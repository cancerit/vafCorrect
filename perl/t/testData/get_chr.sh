samtools faidx chr1_input.fa 1:1-16900500| awk '($1!~/>/)'| perl -p -n -e 's/[aA-zZ]/N/g'  >chr1_out.fa
samtools faidx chr1_input.fa 1:16900501-16909560| awk '($1!~/>/)' >>chr1_out.fa
samtools faidx chr1_input.fa 1:16909561-2210520| awk '($1!~/>/)'| perl -p -n -e 's/[aA-zZ]/N/g'  >>chr1_out.fa
samtools faidx chr1_input.fa 1:2210521-2215500| awk '($1!~/>/)' >>chr1_out.fa

