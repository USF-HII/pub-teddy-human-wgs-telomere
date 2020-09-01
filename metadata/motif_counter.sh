#!/bin/bash

HELP="USAGE:\t `basename $0` -o <FILENAME> [OPTIONS] -h for help"

OPTIONS=":hi:o:pq:Q:rsuv"

HELP_FULL="\n$HELP\n\nmotif_counter 1.0\n\nFor further information, troubleshooting or to report bugs, contact:
\n\n\tMark Hills\tmhills@bccrc.ca
\n\t~~~~ ~~~~~\t~~~~~~~~~~~~~~~
\n\nSoftware to identify any user-supplied genic motif within a sequence library. Run the program, and input the sequence and number of repeats at the prompts.  The program will automatically convert any input sequence to its reverse complement and also scan for that sequence.  Program accepts nucleotides in the form G,C,A,T or g,c,a,t, and accepts the following base pair codes:\n\nR = G or A\t\tY = T or C\t\tM = A or C\nK = G or T\t\tS = G or C\t\tW = A or T\nH = A or C or T\t\tB = G or T or C\t\tV = G or C or A\nD = G or A or T\t\tN = G or A or T or C\n\n For cases in which the reverse complement is not needed (such as in directional libraries), use the -r option.  Analysis can be either for consecutive repeats within a read (e.g. scans for 3xTTAGGG) or non consecutive repeats (e.g. looks for 3 instances of TTAGGG.TTAGGG.TTAGGG where '.' represents any number of interspersed bases).  The former is standard analysis, while the latter must be specified with the -u option.  Note there are two quality score options.  The first (-q) is set in the counting of the total number of reads in a library.  This is set to a default of 20 to exculde poorly mapped reads.  Since many repeats will have a poor mapping quality due to multiple mapping locations, the second quality score (-Q) is set to 0.  If the region you are searching for is specific, this second variable can be changed to make analysis more stringent.

\n\nDEPENDENCIES:\n\t1. Requires samtools set in PATH\n

\n\nOPTIONS:-i PATH\t\tSpecifies location of bam files (default is current)\n
\t-o STR\t\tSpecifies output filename (default 'Pattern_Counter')\n
\t-p\t\tCreates file with corresponding pair end reads from motif read\n
\t-q INT\t\tQuality score for read counts (default is 20)\n
\t-Q INT\t\tQuality score for finding repeats (default is 0)\n
\t-r\t\tSuppresses reverse complement analysis\n
\t-s\t\tCreates a folder with motif-only bam files\n
\t-u\t\tUncoupled repeat analysis\n
\t-v\t\tVerbose output\n"

if ( ! getopts $OPTIONS opt); then
  echo -e $HELP;
#  exit $E_OPTERROR;
fi

# Check dependencies  #CAN I DO WHEREIS AND DIRECT TO PATH IF FOUND?
command -v samtools >/dev/null 2>&1 || { echo "Samtools is required to run this script but it's not installed. Aborting." >&2; exit 1; }

# RESET ALL VARIABLES
INPUT_BAM="."
OUTPUT_NAME="Pattern_Counter"
QUAL=20
QUAL2=0
RECOMP=0
SAVE_BAM=0
UNCOUPLED=0
PAIR_END=1

#extBAM=.bam
#export bam_dir=${bam_dir:-.bam}

let count=1
START_TIME=`date +%s`
#total_file="`ls $bam_dir/*.bam | wc -l`"

while getopts $OPTIONS opt; do
  case $opt in
    h)
     echo -e $HELP_FULL
     exit
    ;;
    i)
     INPUT_BAM=$OPTARG
    ;;
    o)
     OUTPUT_NAME=$OPTARG
    ;;
    p)
     PAIR_END=1
     SAVE_BAM=1
    ;;
    q)
     QUAL=$OPTARG
    ;;
    Q)
     QUAL2=$OPTARG
    ;;
    r)
     RECOMP=1
    ;;
    s)
     SAVE_BAM=1
    ;;
    u)
     UNCOUPLED=1
    ;;
    v)
     verbose=1
    ;;
    \?)
    echo -e "\n#############\nERROR! Invalid Parameter: -$OPTARG\nTry 'motif_counter.sh -h' for help\n#############" >&2
    exit 1
    ;;
  esac
done

if [ -f $OUTPUT_NAME.txt ]; then
    echo -en "$OUTPUT_NAME already exists!\nFile will be overwritten. Continue? [Y/N]  "
    read ans
        if [[ $ans == "Y" ]]; then
            echo -e "Overwriting file.." && rm $OUTPUT_NAME.txt
        elif [[ $ans == "y" ]]; then
            echo -e "Overwriting file.." && rm $OUTPUT_NAME.txt
        else
            echo -en "Please rename output file: "
            read ans2
            OUTPUT_NAME=$ans2
        fi
fi

echo -en "-> Please input sequence motif to search:  "
read PATTERN
echo -en "-> Please input number of repeats:  "
read NUMBER


PATTERN_FWD="`echo "$PATTERN" | awk '{print toupper($0)}' | sed 's/R/[GA]/g' | sed 's/Y/[TC]/g'| sed 's/M/[AC]/g'| sed 's/K/[GT]/g' | sed 's/S/[GC]/g' | sed 's/W/[AT]/g' | sed 's/H/[ACT]/g' | sed 's/B/[GTC]/g' | sed 's/V/[GCA]/g' | sed 's/D/[GAT]/g' | sed 's/N/[ACGT]/g'`"
PATTERN_FWD_FULL=`echo $PATTERN_FWD | printf "$PATTERN_FWD"%.s $(eval echo "{1..$NUMBER}")`

if [[ $RECOMP == 0 ]]; then
    PATTERN_REVCOMP="`echo "$PATTERN" | awk '{print toupper($0)}' | tr [GCTARYMKSWHDBV] [CGATYRKMWSDHVB] | rev`"
    PATTERN_REV="`echo "$PATTERN_REVCOMP" | sed 's/R/[GA]/g' | sed 's/Y/[TC]/g'| sed 's/M/[AC]/g'| sed 's/K/[GT]/g' | sed 's/S/[GC]/g' | sed 's/W/[AT]/g' | sed 's/H/[ACT]/g' | sed 's/B/[GTC]/g' | sed 's/V/[GCA]/g' | sed 's/D/[GAT]/g' | sed 's/N/[ACGT]/g'`"
    PATTERN_REV_FULL="`echo $PATTERN_REV | printf "$PATTERN_REV"%.s $(eval echo "{1..$NUMBER}")`"

    ((verbose)) && echo -e "\n-> -r not selected.  Scanning for ${NUMBER}x$PATTERN_FWD\n-> AND reverse complement ${NUMBER}x$PATTERN_REVCOMP"
else
    ((verbose)) && echo -e "\n-> -r selected.  Scanning for ${NUMBER}x$PATTERN\n-> WILL NOT ANALYZE REVERSE COMPLEMENT"
fi

if [[ $UNCOUPLED == 1 ]]; then
    echo -e "INDEX\tReads_(q>$QUAL)\t${NUMBER}x${PATTERN_FWD}_uncoupled\tPercent_fwd_uncoupled\t${NUMBER}x${PATTERN_REVCOMP}_uncoupled\tPercent_rev_uncoupled\t${NUMBER}x${PATTERN}_total_uncoupled\tPercent_total_uncoupled\t\t" > $OUTPUT_NAME.txt
    ((verbose)) && echo -e "-> -u selected! Uncoupling repeats!"
else
    echo -e "INDEX\tReads_(q>$QUAL)\t$PATTERN($NUMBER)_instances\tPercent_instances\t$PATTERN($NUMBER)_reads\tPercent.reads" > $OUTPUT_NAME.txt
fi

if [[ $SAVE_BAM = 1 ]]; then
    mkdir -p ${OUTPUT_NAME}_BAM_FILES
fi

#for f in "$INPUT_FOLDER"/*$extBAM
#do

((verbose)) && echo -e "\n-> Locating index..."
#IND="`echo "$f" | grep -o -E '[ACGT]{6}'`"

f=$INPUT_BAM
filename=$(basename $INPUT_BAM)
IND=${filename%%.*}

echo $filename
echo $IND

if [[ $IND == "" ]]; then
    ((verbose)) && echo -e "-> COULD NOT IDENTIFY INDEX! \n-> USING $f AS IDENTIFIER INSTEAD..."
    IND="`echo "$f" | sed 's/.*\///'`"
else
    ((verbose)) && echo -e "-> INDEX IDENTIFIED AS $IND"
fi

((verbose)) && echo -e "-> Analyzing $IND\n"

((verbose)) && echo -en "-> Counting reads q > $QUAL in $IND"

samtools view -q $QUAL -b $f > ${IND}_temp.bam
READ_COUNTS="`samtools view ${IND}_temp.bam | wc -l`"
((verbose)) && echo -e ": $READ_COUNTS"

if [[ $UNCOUPLED == 1 ]]; then
    if [[ $NUMBER > 1 ]]; then
        if [[ $RECOMP == 0 ]]; then

                AWK_REPEATER=`printf ".*"$PATTERN_FWD""%.s $(eval echo "{2..$NUMBER}")`

                if [[ $SAVE_BAM == 1 ]]; then
                    ((verbose)) && echo -e "-> Creating motif sam file for $IND with ${NUMBER}x$PATTERN_FWD (q > $QUAL2)"
                    samtools view -q $QUAL2 $f | awk '{if ($10 ~ /'"$PATTERN_FWD"''"$AWK_REPEATER"'/) print $0}' > ${OUTPUT_NAME}_BAM_FILES/$OUTPUT_NAME.fwd

                    ((verbose)) && echo -en "-> Counting uncoupled (non-consecutive) reads for ${NUMBER}x$PATTERN_FWD"
                    PAT_UNCOUP="`wc -l ${OUTPUT_NAME}_BAM_FILES/$OUTPUT_NAME.fwd | awk '{print $1}'`"
                    ((verbose)) && echo -e ": $PAT_UNCOUP"

                else
                    ((verbose)) && echo -en "-> Counting uncoupled (non-consecutive) reads with ${NUMBER}x$PATTERN_FWD"
                    PAT_UNCOUP="`samtools view -q $QUAL2 $f | awk '{if ($10 ~ /'"$PATTERN_FWD"''"$AWK_REPEATER"'/) print $0}' | wc -l | awk '{print $1}'`"
                    ((verbose)) && echo -e ": $PAT_UNCOUP"
                fi
                    UNCOUP_PER="`echo "scale=10; ($PAT_UNCOUP/$READ_COUNTS)*100" | bc`"

                AWK_REPEATER_REV=`printf ".*"$PATTERN_REV""%.s $(eval echo "{2..$NUMBER}")`

                if [[ $SAVE_BAM == 1 ]]; then
                    ((verbose)) && echo -e "-> Creating motif sam file for $IND with ${NUMBER}x$PATTERN_REV (q > $QUAL2) "
                    samtools view -q $QUAL2 $f | awk '{if ($10 ~ /'"$PATTERN_REV"''"$AWK_REPEATER_REV"'/) print $0}' > ${OUTPUT_NAME}_BAM_FILES/$OUTPUT_NAME.rev

                    ((verbose)) && echo -en "-> Counting uncoupled (non-consecutive) reads for ${NUMBER}x$PATTERN_REV"
                    PAT_UNCOUP_REV="`wc -l ${OUTPUT_NAME}_BAM_FILES/$OUTPUT_NAME.rev | awk '{print $1}'`"
                    ((verbose)) && echo -e ": $PAT_UNCOUP_REV"

                    ((verbose)) && echo -e "-> Creating motif bam file for $IND (q > $QUAL2)"
                    samtools view -H $f > ${OUTPUT_NAME}_BAM_FILES/$OUTPUT_NAME.header
                    cat ${OUTPUT_NAME}_BAM_FILES/$OUTPUT_NAME.header ${OUTPUT_NAME}_BAM_FILES/$OUTPUT_NAME.fwd ${OUTPUT_NAME}_BAM_FILES/$OUTPUT_NAME.rev > ${OUTPUT_NAME}_BAM_FILES/${OUTPUT_NAME}_${IND}_q$QUAL2.sam
                    samtools view -Sbh ${OUTPUT_NAME}_BAM_FILES/${OUTPUT_NAME}_${IND}_q$QUAL2.sam -o ${OUTPUT_NAME}_BAM_FILES/${OUTPUT_NAME}_${IND}_q$QUAL2.bam
                    rm ${OUTPUT_NAME}_BAM_FILES/${OUTPUT_NAME}_${IND}_q$QUAL2.sam
                else
                    ((verbose)) && echo -en "-> Counting uncoupled (non-consecutive) reads with ${NUMBER}x$PATTERN_REV"
                    PAT_UNCOUP_REV="`samtools view -q $QUAL2 $f | awk '{print $10}' |  awk '{if ($10 ~ /'"$PATTERN_REV"''"$AWK_REPEATER_REV"'/) print $0}' | wc -l | awk '{print $1}'`"
                    ((verbose)) && echo -e ": $PAT_UNCOUP_REV"
                fi
            UNCOUP_PER_REV="`echo "scale=10; ($PAT_UNCOUP_REV/$READ_COUNTS)*100" | bc`"

            ((verbose)) && echo -en "-> Combined total (non-consecutive) reads of ${NUMBER}x$PATTERN_FWD and ${NUMBER}x$PATTERN_REV"
            UNCOUP_COMBO="`echo "scale=10; ($PAT_UNCOUP + $PAT_UNCOUP_REV)" | bc`"
            UNCOUP_PER_COMBO="`echo "scale=10; ($UNCOUP_COMBO/$READ_COUNTS)*100" | bc`"
            ((verbose)) && echo -e ": $UNCOUP_COMBO"

            echo -e "$IND\t$READ_COUNTS\t$PAT_UNCOUP\t$UNCOUP_PER\t$PAT_UNCOUP_REV\t$UNCOUP_PER_REV\t$UNCOUP_COMBO\t$UNCOUP_PER_COMBO" >> $OUTPUT_NAME.txt
        else
            ((verbose)) && echo -en "-> Counting uncoupled (non-consecutive) reads of ${NUMBER}x$PATTERN_FWD (Not rev complement)"
            AWK_REPEATER=`printf ".*"$PATTERN_FWD""%.s $(eval echo "{2..$NUMBER}")`
            PAT_UNCOUP="`samtools view -q $QUAL2 $f | awk '{if ($10 ~ /'"$PATTERN_REV"''"$AWK_REPEATER_REV"'/) print $0}' | wc -l | awk '{print $1}'`"
            UNCOUP_PER="`echo "scale=10; ($PAT_UNCOUP/$READ_COUNTS)*100" | bc`"
            ((verbose)) && echo -e ": $PAT_UNCOUP"
            echo -e "$IND\t$READ_COUNTS\t$PAT_UNCOUP\t$UNCOUP_PER" >> $OUTPUT_NAME.txt
        fi
    else
        ((verbose)) && echo -e "-> -u selected, but option is redunant for one repeat!\n-> Not running uncoupled analysis"
    fi

elif [[ $RECOMP == 0 ]]; then
    ((verbose)) && echo -en "-> Counting number of instances of $PATTERN_FWD($NUMBER) and $PATTERN_REV($NUMBER)"
    PAT_STRING="`samtools view ${IND}_temp.bam | grep -o -E "$PATTERN_FWD_FULL|$PATTERN_REV_FULL" | wc -l`"
    ((verbose)) && echo -e ": $PAT_STRING"
    ((verbose)) && echo -en "-> Counting number of reads with $PATTERN_FWD($NUMBER) and $PATTERN_REV($NUMBER)"
    PAT_READS=`samtools view ${IND}_temp.bam | grep -o -E -c "$PATTERN_FWD_FULL|$PATTERN_REV_FULL"`
    ((verbose)) && echo -e ": $PAT_READS"
    PAT_STRING_PER="`echo "scale=10; ($PAT_STRING/$READ_COUNTS)*100" | bc`"
    PAT_READS_PER="`echo "scale=10; ($PAT_READS/$READ_COUNTS)*100" | bc`"
    echo -e "$IND\t$READ_COUNTS\t$PAT_STRING\t$PAT_STRING_PER\t$PAT_READS\t$PAT_READS_PER" >> $OUTPUT_NAME.txt

        #if [[ $SAVE_BAM == 1 ]]; then
        #    ((verbose)) && echo -e "-> Creating motif bam file for $IND"
        #    ${IND}_temp.bam |  awk '{if (($10 ~ /'"$PATTERN_FWD_FULL"'/) || ($10 ~ /'"$PATTERN_REV_FULL"'/)) print $0}' > ${OUTPUT_NAME}_BAM_FILES/$OUTPUT_NAME.fwd
        #    samtools view -H $f > ${OUTPUT_NAME}_BAM_FILES/$OUTPUT_NAME.header
        #    cat ${OUTPUT_NAME}_BAM_FILES/$OUTPUT_NAME.header ${OUTPUT_NAME}_BAM_FILES/$OUTPUT_NAME.fwd > ${OUTPUT_NAME}_BAM_FILES/${OUTPUT_NAME}.rev
        #    samtools view -Sbh ${OUTPUT_NAME}_BAM_FILES/${OUTPUT_NAME}.rev -o ${OUTPUT_NAME}_BAM_FILES/${OUTPUT_NAME}_${IND}_q$QUAL2.bam
        #fi
else
    ((verbose)) && echo -en "-> Counting number of instances of $PATTERN_FWD($NUMBER)"
    PAT_STRING="`samtools view -q $QUAL2 $f | grep -o -E "$PATTERN_FWD_FULL" | wc -l`"
    ((verbose)) && echo -e ": $PAT_STRING"
    ((verbose)) && echo -en "-> Counting number of reads with $PATTERN_FWD($NUMBER)"
    PAT_READS="`samtools view -q $QUAL2 $f | grep -o -E -c "$PATTERN_FWD_FULL"`"
    ((verbose)) && echo -e ": $PAT_READS"
    PAT_STRING_PER="`echo "scale=10; ($PAT_STRING/$READ_COUNTS)*100" | bc`"
    PAT_READS_PER="`echo "scale=10; ($PAT_READS/$READ_COUNTS)*100" | bc`"
    echo -e "$IND\t$READ_COUNTS\t$PAT_STRING\t$PAT_STRING_PER\t$PAT_READS\t$PAT_READS_PER" >> $OUTPUT_NAME.txt
        #if [[ $SAVE_BAM == 1 ]]; then
        #    ((verbose)) && echo -e "-> Creating motif bam file for $IND"
        #    samtools view -q $QUAL2 $f |  awk '{if ($10 ~ /'"$PATTERN_FWD_FULL"'/) print $0}' > ${OUTPUT_NAME}_BAM_FILES/$OUTPUT_NAME.fwd
        #    samtools view -H $f > ${OUTPUT_NAME}_BAM_FILES/$OUTPUT_NAME.header
        #    cat ${OUTPUT_NAME}_BAM_FILES/$OUTPUT_NAME.header ${OUTPUT_NAME}_BAM_FILES/$OUTPUT_NAME.fwd > ${OUTPUT_NAME}_BAM_FILES/${OUTPUT_NAME}_${IND}_q$QUAL2.rev
        #    samtools view -Sbh ${OUTPUT_NAME}_BAM_FILES/${OUTPUT_NAME}_${IND}_q$QUAL2.rev -o ${OUTPUT_NAME}_BAM_FILES/${OUTPUT_NAME}_${IND}_q$QUAL2.bam
        #fi
fi


rm -rf ${IND}_temp.bam
#if [[ $PAIR_END == 1 ]]; then
#
#    echo -e "\n-> Generating Index file from $f"
#    samtools view ${OUTPUT_NAME}_BAM_FILES/${OUTPUT_NAME}_${IND}_q$QUAL2.bam | awk '{k=$1; line[k]=$1; freq[k]++ } END { for (x in freq) if (freq[x]==1) print line[x]}' | awk '{if ($0 !~ /@/) print $0}' > unique_temp.sam
#
#    echo -e "\n-> Matching mate pair from $OUTPUT_NAME_$IND and $f"
#    samtools view $f | awk 'NR==FNR {a[$1] = $1; next}{if (a[$1]) print $0}' unique_temp.sam - > temp.sam
#
#    echo -e "\n-> Sorting mate pair sam file for $IND"
#    sort temp.sam > ${IND}_PET.sam
#    rm unique_temp.sam
#    rm temp.sam
#fi


MID_TIME=`date +%s`
dt=`expr $MID_TIME - $START_TIME`
ds=$((dt % 60))
dm=$(((dt / 60) % 60))
dh=$((dt / 3600))
time=`printf '%d:%02d:%02d' $dh $dm $ds`

((verbose)) &&  echo -e "\n-> Elapsed time for $IND samples: $time"
#remain=`expr $total_file - $count`
#frac=`echo "scale=2; ($count / $total_file) * 100" | bc`
#((verbose)) &&  echo -e "-> $remain samples remaining ($frac % complete)."
#count=`echo "$(($count + 1))"`


    #if [[ $SAVE_BAM == 1 ]]; then
    #        rm ${OUTPUT_NAME}_BAM_FILES/$OUTPUT_NAME.fwd
    #        rm ${OUTPUT_NAME}_BAM_FILES/$OUTPUT_NAME.rev
    #        rm ${OUTPUT_NAME}_BAM_FILES/$OUTPUT_NAME.header
    #fi

END_TIME=`date +%s`
dt=`expr $END_TIME - $START_TIME`
ds=$((dt % 60))
dm=$(((dt / 60) % 60))
dh=$((dt / 3600))
endtime=`printf '%d:%02d:%02d' $dh $dm $ds`
((verbose)) &&  echo -e "\n-> TOTAL ELAPSED TIME: $endtime\n"

