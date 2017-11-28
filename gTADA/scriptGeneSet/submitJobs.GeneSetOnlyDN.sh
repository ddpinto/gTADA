export R_LIBS="/hpc/users/nguyet26/InstallSoftware/Rlibs"
module load R

OUT_DIR="TestGeneSetDNmcmc4"
OUT_DIR="TestGeneSetDNmcmc4_April3"
mkdir ${OUT_DIR}
RUN_TIME=$(date +"%H_%M_%h_%d_%Y")
rFile=gTADAsimulateData.R
rFile=gTADAsimulateDataSimulation.R

for RUN_TIME in $(seq 1 1 30)
do


for gDN1 in 30 40
#12 
#30 40
#12 20
do
for gDN2 in 2 
#5
do

for a0 in -2
#-2 -3 -4 -1 
do
for a1 in 5 4 3 2 1 0.5 0.2 10^-3
#5 0.5 1 2 3 4 5 10^-3 0 0.2 
do

for pi0 in 0.025 0.05 0.1 0.15 0.2 0.001
do
ntrio=1077
#/hpc/users/nguyet26/R/bin/
bsub -P acc_psychgen -q premium -n 1 -W 6:25 -R "span[hosts=1]" -R "rusage[mem=10000]" \
	-e ${OUT_DIR}/${rFile}.${ntrio}.${pi0}.${gDN1}.${gDN2}.${a0}.${a1}.${RUN_TIME}.err \
	-o ${OUT_DIR}/${rFile}.${ntrio}.${pi0}.${gDN1}.${gDN2}.${a0}.${a1}.${RUN_TIME}.out \
R CMD BATCH --no-save --no-restore \
    "--args ${ntrio} ${gDN1} ${gDN2} ${a0} ${a1} ${pi0} ${OUT_DIR} ${RUN_TIME}" ./${rFile} \
	${OUT_DIR}/${rFile}.${ntrio}.${pi0}.${gDN1}.${gDN2}.${a0}.${a1}.${RUN_TIME}.Rout

done
done
done
done
done
done

