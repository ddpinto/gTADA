OUT_DIR=W_logLK_AUT
mkdir ${OUT_DIR}
##################
while read GENE_SET
do

rFile=testLogLKgeneSet.R

echo $rFile

MM=1000
iChain=1

bsub -P acc_psychgen -q premium -n $iChain -W 0:55 -R "rusage[mem=${MM}]" -R "span[hosts=1]" \
	-e ${OUT_DIR}/${rFile}.${MM}.err \
	-o ${OUT_DIR}/${rFile}.${MM}.out \
/hpc/users/nguyet26/R/bin/R CMD BATCH --no-save --no-restore \
    "--args ${GENE_SET} ${OUT_DIR}" ${rFile} \
${OUT_DIR}/${rFile}.${GENE_SET}.${MM}.Rout

done < listGeneSet.txt