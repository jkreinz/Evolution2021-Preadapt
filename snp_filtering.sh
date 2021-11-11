#merge individual regional files by chr
cat ../scafs | parallel -j 16 "mkdir temp_{} && cat ../{}:* | vcffirstheader | vcfstreamsort -w 10000 | vcfuniq | vcfallelicprimitives --keep-info | bcftools
sort -T temp_{} - > {}_herblastrerun.vcf"

#merge chr files into one
cat Scaffold_{1..16}_herblastrerun.vcf | vcffirstheader | vcfuniq | bcftools view -i 'F_MISSING<=0.3' - > allchroms_herblastrerun.vcf

#dustmask, hard filt for variants
cat <(head -n 1000 allchroms_commongarden.vcf | zgrep "#") <(bedtools subtract -a allchroms_commongarden.vcf -b /ohta/julia.kreiner/waterhemp/herbarium/femaleref/unscaled/dustmasked.bed) | vcfallelicprimitives --keep-info | vcffilter -s -f "QUAL > 10 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPL > 0 & RPR > 0"  > allchroms_commongarden_dustm_mallelesplt_hardfilt_noAB.vcf

#high depth based filters
vcffilter -f "QUAL / DP > 0.25" allchroms_commongarden.vcf | tee allchroms_commongarden.qualDP | cut -f8 | grep -oe "DP=[0-9]*" | sed -s 's/DP=//g' > allchroms_commongarden_qualDP.DEPTH
vcffilter -f "QUAL / DP > 0.25" allchroms_commongarden.vcf | grep -v "#" | cut -f1,2,6 > allchroms_commongarden_qualDP.QUAL

~/software/bin/bin/mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' allchroms_commongarden_qualDP.QUAL
python -c "print int(759.771+(3*(759.771**0.5)))"

vcffilter -f "QUAL / DP > 0.25" allchroms_commongarden.vcf | tee >(cut -f8 | grep -oe "DP=[0-9]*" | sed -s 's/DP=//g' > allchroms_commongarden_qualDP.DEPTH) >(grep -v "#" | cut -f1,2,6 > allchroms_commongarden_qualDP.QUAL)

#remove particularily high depth sites (depth is greater than the mean + 3(sqrt(mean)), that do not have a qual score < 2* depth
paste allchroms_commongarden_qualDP.QUAL allchroms_commongarden_qualDP.DEPTH | ~/software/bin/bin/mawk -v x=870 '$4 > x' | ~/software/bin/bin/mawk '$3 < 2 * $4' > allchroms_commongarden_qualDP.lowQD

vcftools --vcf allchroms_commongarden_dustm_mallelesplt_hardfilt_AC0_noAB.vcf --site-depth --exclude-positions allchroms_commongarden_qualDP.lowQD | cut -f3 > allchroms_commongarden_hardfilt_QD.justdepth

#plot this to also get max depth cuttoff
~/software/bin/bin/mawk '!/D/' allchroms_commongarden_hardfilt_QD.justdepth | ~/software/bin/bin/mawk -v x=188 '{print $1/x}' > meandepthpersite

vcftools --vcf allchroms_commongarden_dustm_mallelesplt_hardfilt_AC0_noAB.vcf --recode-INFO-all --max-meanDP XXXX --exclude-positions allchroms_commongarden_qualDP.lowQD --recode --min-alleles 2 --max-alleles 2 --remove-indels | vcffilter -s -f "AB < 0.01 | AB > 0.25 & AB < 0.75" > allchroms_commongarden_dustm_hardfilt_AC0_QualDP_biall_AB.vcf
