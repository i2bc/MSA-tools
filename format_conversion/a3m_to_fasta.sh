#! /bin/sh

# READ CONFIG FILE
script_dir=$(dirname $(dirname $(readlink -f "$0")))
pathReformat=$($script_dir/common.sh $script_dir/config.yaml pathReformat)

OPTSTRING=":i:o:"

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
