rFile=increaseLLK1_removeLowLK1.R
rFile=increaseLLK1_removeLowLK1_useNLMINB.R
#OUTDIR=Test_LLK2
OUTDIR=Test_LLKnlminb

for ii in "ID"
#"EPI" "CHD" "DD" "ASD" "ID" "SCZ"
#"EPI" "CHD"
do
index=$(date|sed 's/ //g'|sed 's/:/_/g')
bsub -q premium -P acc_psychgen -M 20000 -R "rusage[mem=3000]" -W 22:50 -oo ${OUTDIR}/OUT \
  R CMD BATCH --no-save --no-restore "--args ${ii} ${index} ${OUTDIR}" ./${rFile} ${OUTDIR}/R.${ii}.RoutOUT

done
