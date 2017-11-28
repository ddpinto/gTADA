#module unload R 
#module unload R/3.2.3
#module unload R/3.3.0
#module purge
#module load R/3.3.1
#module load openssl
#outDir=$(pwd)/testRSTAN

export R_LIBS="/hpc/users/nguyet26/InstallSoftware/Rlibs"
module load R

outDir=$(pwd)/TestOutData
CURRENT_DIR=$(pwd)
echo $outDir

ii=500 #Burn-in
jj=10000 #Iteration
iThin=1

###Errors can happen here, if we set iChain > 1 (n cores > 1)
iChain=3

###
changeFile=101

for MM in 3000
do
##Errors can happen here if we set a large jj (e.g., jj > 30000)
for jj in 10000
#5000
#jj is running times
do

for sigmaPrior in $(seq 0.2 1 20)
#$(seq 0.22 0.4 5)
#$(seq 0.2 1 20)
#0.25 1 4 10 15
#6.5 8 10 11 15 20
#0.25 0.5 1 2 3 4 5 6
#$(seq 0.001 0.15 7)
do

for kF in $(seq 2 1 10)
#$(seq 1 1 10)
do

iThin=$(expr ${jj} / 500)
iThin=1
nThin=1
kk=$(date +"%H_%M_%h_%d_%Y")
#kk="REMOVEzero"

for lowerHyperGamma in 1
#0 0.5 0.75 1 2
do
while read line 
do

annotationType=$(echo $line|awk '{print $1}')
annotationType2=$(echo $line|awk '{print $2}')
annotationType3=$(echo $line|awk '{print $3}')


OUT_DATA_DIR=W_EPIgeneSetN2_BayesianPrior 
OUT_DATA_DIR=W_EPIgeneSetN2_BayesianPriorLaplace
OUT_DATA_DIR=W_AddGeneSetKF5
OUT_DATA_DIR=W_AddGeneSetKF5_LaplacePrior
OUT_DATA_DIR=W_AddGeneSetKF6_NormalPrior
#ElasticNcolinear
outDir=${OUT_DATA_DIR}
mkdir ${OUT_DATA_DIR}
mkdir ${outDir}

rFile=script_EPI_ID_denovo_Sep_2016_DN.R
#rFile=script_EPI_ID_denovo_Sep_2016_DN_oneClass.R
rFile=script_EPI_ID_denovo_Sep_2016_DN_useNew.R
rFile=script_EPI_ID_denovo_Sep_2016_DN_useNewGeneSets.R
rFile=script_EPI_ID_denovo_Sep_2016_DN_useNewGeneSetsCV_Kfold.R
#rFile=script_EPI_ID_denovo_Sep_2016_DN_useNewGeneSetsForEachGeneSets.R
#rFile=extTADAforAUT_CC.R
#rFile=extTADAforAUT_CC.R
#rFile=extTADAforAUT_DN.R
#################
nGroupDN=2  ##Number of categories
nGroupCC=2
adjustHyperBeta=1 #1 #1
swapData=2 #1 #2 for LoF
##################


echo $rFile

#/sc/orga/scratch/nguyet26/Re_annotate/script/TestOutData/outFile.out

bsub -P acc_psychgen -q premium -n $iChain -W 8:35 -R "rusage[mem=${MM}]" -R "span[hosts=1]" \
	-e ${outDir}/${rFile}.${iChain}.$jj.$ii.$kk.1.2.binomial.beta.cc.${changeFile}.${annotationType}.${annotationType2}.${annotationType3}.${MM}.err \
	-o ${outDir}/${rFile}.${iChain}.burnin.$jj.iteration.$ii.$kk.2.1.${changeFile}.${annotationType}.${annotationType2}.${annotationType3}.binomial.beta.cc.${MM}.out \
/hpc/users/nguyet26/R/bin/R CMD BATCH --no-save --no-restore \
    "--args $ii $jj $kk $iThin $iChain $annotationType $annotationType2 ${MM} ${CURRENT_DIR} ${lowerHyperGamma} ${annotationType3} ${OUT_DATA_DIR} ${nGroupDN} ${nGroupCC} ${adjustHyperBeta} ${swapData} ${sigmaPrior} ${kF}" ${CURRENT_DIR}/${rFile} \
$outDir/${rFile}.chain.${iChain}.${jj}.$kk.${changeFile}.${annotationType}.${annotationType2}.${annotationType3}.${MM}.Rout


done < listGroupEPIandID.txt2
#listGroup.txt1o
#list.3ClassesArray.txt1
#listDN.txt


done

done

done

done
done