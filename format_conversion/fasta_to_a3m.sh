#! /bin/sh

# READ CONFIG FILE
script_dir=$(dirname $(dirname $(readlink -f "$0")))
pathReformat=$($script_dir/common.sh $script_dir/config.yaml pathReformat)

OPTSTRING=":i:o:"

while getopts ${OPTSTRING} opt; do
  case "${opt}" in
    i) 
      fasta=${OPTARG}
    ;;
    o) 
      a3m=${OPTARG}
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

${pathReformat} ${fasta} ${a3m}.tmp.a3m
echo "${a3m}"
# Then make sure you only have 1 line per sequence:
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ${a3m}.tmp.a3m > ${a3m}
rm ${a3m}.tmp.a3m
