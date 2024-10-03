#!/bin/bash

#CTCF, DNase, H3K27ac, H3K4me3
mode=CTCF
genome=hg38

fileDir=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/MaxZ/Run2
scriptDir=/Users/ryanhagan/NoCoSMiCC/Scripts/ENCODE_cCRE_pipeline_replication

echo "Processing Files ..."
file=/Users/ryanhagan/NoCoSMiCC/Files/$mode-List_colon.txt
mkdir -p $fileDir/signal_output/$mode
q=$(wc -l $file | awk '{print $1}')
for j in `seq 1 1 $q`
do
    A=$(awk '{if (NR == '$j') print $1}' $file)
    B=$(awk '{if (NR == '$j') print $2}' $file)
    C=$( awk '{if (NR == '$j') print $3}' $file)
    mv $fileDir/signal_output/$A"-"$B.txt $fileDir/signal_output/$mode/
done

cd $fileDir/signal_output/$mode/

echo "Determining maxZ..."
python $scriptDir/max.zscore.array.py *.txt > $fileDir/$genome-$mode-maxZ.txt
mv *.txt ../



# cd to the script dir to execute the signal script
cd /Users/ryanhagan/NoCoSMiCC/Scripts/ENCODE_cCRE_pipeline_replication


## ATAC

./merged_CTCF_signal.sh /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/ATAC/merged_ATAC_signal.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/ATAC/ENCSR086OGH-ENCFF049WJI.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/ATAC/ENCSR355SGJ-ENCFF961XDO.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/ATAC/ENCSR386HAZ-ENCFF668GUI.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/ATAC/ENCSR404LLJ-ENCFF033RPN.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/ATAC/ENCSR548QCP-ENCFF796DRU.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/ATAC/ENCSR600ZHS-ENCFF057BIJ.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/ATAC/ENCSR668VCT-ENCFF509NRA.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/ATAC/ENCSR761TKU-ENCFF811ERB.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/ATAC/ENCSR846VLJ-ENCFF784HME.txt

mode=ATAC
genome=hg38

fileDir=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/$genome-rCARs/MaxZ
scriptDir=/Users/ryanhagan/NoCoSMiCC/Scripts/ENCODE_cCRE_pipeline_replication

cd /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/$mode/
python $scriptDir/max.zscore.array.py merged_ATAC_signal.txt > $fileDir/$genome-$mode-maxZ.txt

## single cell ATAC

# cd to the script dir to execute the signal script
cd /Users/ryanhagan/NoCoSMiCC/Scripts/ENCODE_cCRE_pipeline_replication

./merged_CTCF_signal.sh /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/scATAC/merged_scATAC_signal.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/scATAC/Transverse_ENCSR349XKD-Transverse_ENCSR349XKD.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/scATAC/Transverse_ENCSR434SXE-Transverse_ENCSR434SXE.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/scATAC/Transverse_ENCSR506YMX-Transverse_ENCSR506YMX.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/scATAC/Transverse_ENCSR997YNO-Transverse_ENCSR997YNO.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/scATAC/Left_ENCSR830FPR-Left_ENCSR830FPR.txt

mode=scATAC
genome=hg38

fileDir=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/$genome-rCARs/MaxZ
scriptDir=/Users/ryanhagan/NoCoSMiCC/Scripts/ENCODE_cCRE_pipeline_replication

cd /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/$mode/
python $scriptDir/max.zscore.array.py merged_scATAC_signal.txt > $fileDir/$genome-$mode-maxZ.txt


## CTCF

# cd to the script dir to execute the signal script
cd /Users/ryanhagan/NoCoSMiCC/Scripts/ENCODE_cCRE_pipeline_replication

./merged_CTCF_signal.sh /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/CTCF/merged_CTCF_signal.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/CTCF/ENCSR102CSD-ENCFF170KVV.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/CTCF/ENCSR236YGF-ENCFF686TFX.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/CTCF/ENCSR449SEF-ENCFF626YRZ.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/CTCF/ENCSR558HTE-ENCFF349LWW.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/CTCF/ENCSR608WPS-ENCFF646EZE.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/CTCF/ENCSR769WKR-ENCFF493XMW.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/CTCF/ENCSR833FWC-ENCFF517GQE.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/CTCF/ENCSR907BES-ENCFF435CDF.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/CTCF/ENCSR721AHD-ENCFF460ECX.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/CTCF/ENCSR857RJQ-ENCFF634JUC.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/CTCF/ENCSR925GDS-ENCFF154FOF.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/CTCF/ENCSR222SQE-ENCFF985EXZ.txt

mode=CTCF

cd /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/$mode/
python $scriptDir/max.zscore.array.py merged_CTCF_signal.txt > $fileDir/$genome-$mode-maxZ.txt


## H3K4me3
cd /Users/ryanhagan/NoCoSMiCC/Scripts/ENCODE_cCRE_pipeline_replication
./merged_CTCF_signal.sh /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K4me3/merged_H3K4me3_signal.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K4me3/ENCSR315EZG-ENCFF487CTD.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K4me3/ENCSR557OWY-ENCFF568IBR.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K4me3/ENCSR813ZEY-ENCFF252OBP.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K4me3/ENCSR933BVL-ENCFF339CRV.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K4me3/ENCSR172LVU-ENCFF886LUE.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K4me3/ENCSR421HUB-ENCFF423YBA.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K4me3/ENCSR793IKH-ENCFF800VAP.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K4me3/ENCSR900UIP-ENCFF237VMY.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K4me3/ENCSR960AAL-ENCFF122GSI.txt

mode=H3K4me3

cd /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/$mode/
python $scriptDir/max.zscore.array.py merged_H3K4me3_signal.txt > $fileDir/$genome-$mode-maxZ.txt


## H3K27ac
cd /Users/ryanhagan/NoCoSMiCC/Scripts/ENCODE_cCRE_pipeline_replication
./merged_CTCF_signal.sh /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K27ac/merged_H3K27ac_signal.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K27ac/ENCSR069EGE-ENCFF427MZX.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K27ac/ENCSR208QRN-ENCFF318ECM.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K27ac/ENCSR640XRV-ENCFF532ZGB.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K27ac/ENCSR792VLP-ENCFF741NZM.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K27ac/ENCSR268ZCF-ENCFF111DLN.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K27ac/ENCSR561YSH-ENCFF124AKT.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K27ac/ENCSR641SDI-ENCFF468UEP.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K27ac/ENCSR807XUB-ENCFF608MDC.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K27ac/ENCSR937EVN-ENCFF322NLT.txt

mode=H3K27ac

cd /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/$mode/
python $scriptDir/max.zscore.array.py merged_H3K27ac_signal.txt > $fileDir/$genome-$mode-maxZ.txt

## DNase
cd /Users/ryanhagan/NoCoSMiCC/Scripts/ENCODE_cCRE_pipeline_replication
./merged_CTCF_signal.sh /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/DNase/merged_DNase_signal.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/DNase/ENCSR340MRJ-ENCFF405NTZ.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/DNase/ENCSR504WYA-ENCFF291LZE.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/DNase/ENCSR923JYH-ENCFF753MXL.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/DNase/ENCSR276ITP-ENCFF299OOV.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/DNase/ENCSR279SXQ-ENCFF452PQB.txt \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/DNase/ENCSR867XFA-ENCFF613BDA.txt

mode=DNase

cd /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/$mode/
python $scriptDir/max.zscore.array.py merged_DNase_signal.txt > $fileDir/$genome-$mode-maxZ.txt




