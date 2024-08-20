#! /bin/sh

OPTSTRING=":i:o:"
pathReformat=/data/work/I2BC/hugo.pointier/msa_tools/tools/reformat.pl
pathReformat=~/programmation/stage/script_msa_tools/hhsuite/scripts/reformat.pl
while getopts ${OPTSTRING} opt; do
  case "${opt}" in
    i) 
      a3m=${OPTARG}
    ;;
    o) 
      fasta=${OPTARG}    
    ;;
    :)
      echo "Option -${OPTARG} requires an argument." >&2
      exit 1
    ;;
    ?) 
      echo "Invalid option: -${OPTARG}"
    ;;
  esac
done
  
firstline=$(head -n 1 ${a3m})
${pathReformat} -r a3m fas ${a3m} ${fasta}
if [[ ${firstline:0:1} == "#" ]] ; then sed  -i '1i '"$firstline" ${fasta} ; fi # add the header if there is one
