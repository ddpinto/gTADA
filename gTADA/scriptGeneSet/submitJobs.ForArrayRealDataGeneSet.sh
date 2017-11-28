#module unload R 
#module unload R/3.2.3
#module unload R/3.3.0
#module purge
#module load R/3.3.1
#module load openssl
#outDir=$(pwd)/testRSTAN

#export R_LIBS="/hpc/users/nguyet26/InstallSoftware/Rlibs"

outDir=$(pwd)/TestOutData
CURRENT_DIR=$(pwd)
echo $outDir

ii=500 #Burn-in
jj=10000 #Iteration
iThin=10

###Errors can happen here, if we set iChain > 1 (n cores > 1)
iChain=3

###
changeFile=101


for MM in 5000
do
##Errors can happen here if we set a large jj (e.g., jj > 30000)
for jj in 1000 5000 10000
#500 1000 2000
#10000
#2000
do
nThin=200
iThin=$(expr $jj / 1000)
nThin=$(expr $jj / 1000)
#iThin=1
#nThin=1

kk=$(date +"%H_%M_%h_%d_%Y")
#kk="REMOVEzero"

for lowerHyperGamma in 1
#0 0.5 0.75 1 2
do
while read line 
do

simgaPrior=2

annotationType=$(echo $line|awk '{print $1}')
annotationType2=$(echo $line|awk '{print $2}')
annotationType3=$(echo $line|awk '{print $3}')

misType="missense"
misType="damaging"

#################
#ADD GENE SET
rFile=scriptSCZ2017usingGeneSet.R
#rFile=scriptAUTfromTADAusingGeneSet.R

OUT_DATA_DIR=W_SCZgeneSet3
#OUT_DATA_DIR=W_AUTgeneSet3
#OUT_DATA_DIR=W_SCZthesameHyperGamma
outDir=${OUT_DATA_DIR}
mkdir ${outDir}
mkdir ${OUT_DATA_DIR}

#########
nGroupDN=3  ##Number of categories
nGroupCC=1
adjustHyperBeta=1 #0 #1 #0 # 1 #0 # 1 #1
swapData=2 #1 #2 for LoF
singleClass="lof" 
#singleClass="missensedamaging"
#singleClass="cfpk"
singleClass="NO"

##################


echo $rFile

#/sc/orga/scratch/nguyet26/Re_annotate/script/TestOutData/outFile.out

bsub -P acc_psychgen -q premium -n $iChain -W 45:35 -M 25000 -R "rusage[mem=${MM}]" -R "span[hosts=1]" \
	-e ${outDir}/${rFile}.${iChain}.$jj.$ii.$kk.1.2.binomial.beta.cc.${changeFile}.${annotationType}.${annotationType2}.${annotationType3}.${MM}.err \
	-o ${outDir}/${rFile}.${iChain}.burnin.$jj.iteration.$ii.$kk.2.1.${changeFile}.${annotationType}.${annotationType2}.${annotationType3}.binomial.beta.cc.${MM}.out \
/hpc/users/nguyet26/R/bin/R CMD BATCH --no-save --no-restore \
    "--args $ii $jj $kk $iThin $iChain $annotationType $annotationType2 ${MM} ${CURRENT_DIR} ${lowerHyperGamma} ${annotationType3} ${OUT_DATA_DIR} ${nGroupDN} ${nGroupCC} ${adjustHyperBeta} ${misType} ${sigmaPrior} ${swapData} ${singleClass}" ${CURRENT_DIR}/${rFile} \
$outDir/${rFile}.chain.${iChain}.${jj}.$kk.${changeFile}.${annotationType}.${annotationType2}.${annotationType3}.${MM}.Rout

done < list2.txt
#listGroup.txt1
#list.3ClassesArray.txt1
#listDN.txt


done

done

done


