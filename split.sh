#/bin/bash

for i in {1..22} X Y
do
  bcftools view $1##idx##$2 --regions ${i} -o $3_${i}.vcf.gz -Oz
done