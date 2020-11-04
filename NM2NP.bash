#!/bin/bash
set -x
# Script to get matching NP protein accession given RefSeq accessions.
# Due to the occasional network errors (about 1 in 100), need to rerunning failed NM
# ... watch the script though
#
CAVA_DB=$1
# cava.ids.2NP.txt
OUTFILE=$2
PREFIX=$3
idsfile=${PREFIX}.transcripts.ids.txt
idstodo=${PREFIX}.transcripts.todo.ids.txt
gunzip -c ${CAVA_DB} | cut -f1 > ${idsfile}
TMPOUT=${PREFIX}.tr_pr.tmp.ids

NIDS=`wc -l ${idsfile} | cut -f1 -d ' '`
if [[ $NIDS -eq 0 ]]; then
    touch $OUTFILE
    rm $OUTFILE
    touch $OUTFILE
    quit
fi

touch ${OUTFILE}
NOUT=`wc -l ${OUTFILE} | cut -f1 -d ' '` 
echo "Creating Initial $idstodo file, allowing for restarts"
if [[ $NOUT -gt 0 ]] ; then
    gawk -F "\t" '(NF==1 || length($2)<1) {next}{print $0}' ${OUTFILE} > ${TMPOUT}
    cp ${TMPOUT} ${OUTFILE}
    echo "" | gawk -F "\t" -v OUTFILE=${OUTFILE} -v IDSFILE=${idsfile} 'BEGIN{while(getline<IDSFILE){NEED[$1]=1};while(getline < OUTFILE) {delete NEED[$1]}}{next}END{for( i in NEED) {print i}}' > ${idstodo}
else
    cp ${idsfile} ${idstodo}
fi

N=`wc -l ${idstodo} | cut -f1 -d ' '`

m=0
while [[ $N -gt 0 ]] ; do
    m=$(( m + 1 ))
    echo "Loop $m over $CAVA_DB to generate $N protein ids"
    IFS=$'\n'; for i in `cat ${idstodo} ` ; do NP=""; NP=$(IFS=$'\n'; esearch -db nuccore -query ${i} | efetch -format gpc | xtract -pattern INSDSeq -block INSDQualifier -if INSDQualifier_name -equals protein_id -element INSDQualifier_value); echo "${i} $NP" | sed -e 's/ /\t/' >> $OUTFILE; sleep 0.33;done
    NOUT=`wc -l ${OUTFILE} | cut -f1 -d ' '` 
    if [ "$NOUT" -gt 0 ] ; then
        gawk -F "\t" '(NF==1 || length($2)<1) {next}{print $0}' ${OUTFILE} > ${TMPOUT}
        cp ${TMPOUT} ${OUTFILE}
        echo "" | gawk -F "\t" -v OUTFILE=${OUTFILE} -v IDSFILE=${idsfile} 'BEGIN{while(getline<IDSFILE){NEED[$1]=1};while(getline < OUTFILE) {delete NEED[$1]}}{next}END{for( i in NEED) {print i}}' > ${idstodo}
    else
        cp ${idsfile} ${idstodo}
    fi

    N=`wc -l ${idstodo} | cut -f1 -d ' '`
done

# touch before deleting to avoid errors
touch ${TMPOUT}
rm ${TMPOUT}
touch ${idstodo}
rm ${idstodo}
