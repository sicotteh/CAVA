#IFS=$'\n'; for i in `cat cava.ids | gawk -F "\t" 'BEGIN{N=0}$1=="NM_182709.2"{N=1}N==1{print $1}' `; do NP=""; NP=$(IFS=$'\n'; esearch -db nuccore -query ${i} | efetch -format gpc | xtract -pattern INSDSeq -block INSDQualifier -if INSDQualifier_name -equals protein_id -element INSDQualifier_value); echo "${i} $NP" >> cava.ids.2NP.txt; sleep 0.33;done

# Loop over this section until the list of missing is empty
#gawk -F " " '(NF==1 || length($2)<1) {print $1}' cava.ids.2NP.txt > redos.ids.txt
# Add a step to clean up cava.ids.2NP.txt
IFS=$'\n'; for i in `cat redos2.ids.txt`; do NP=""; NP=$(IFS=$'\n'; esearch -db nuccore -query ${i} | efetch -format gpc | xtract -pattern INSDSeq -block INSDQualifier -if INSDQualifier_name -equals protein_id -element INSDQualifier_value); echo "${i} $NP" >> cava.ids.2NP.txt; sleep 0.33;done

# Add a final QC to make sure all cava.ids have an equivalent.
