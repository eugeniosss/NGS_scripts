#!/bin/bash

set -e

usage() {
  echo "Usage: $0 -g <GOALS_FILE> -c <CORRESPONDENCIES_FILE> -f <FINAL_FILE> -p <PATH_RAW>"
  echo "Example: $0 -g list.txt -c /media/jbod2/eugenio/correspondencies.csv -f dic.txt -p /media/jbod3/joel/ngs/raw_data/sequencing"
}

while getopts ":g:c:f:p:" opt; do
  case ${opt} in
    g )
      GOALS=$OPTARG
      ;;
    c )
      CORRESPONDENCIES=$OPTARG
      ;;
    f )
      FINAL=$OPTARG
      ;;
    p )
      PATH_RAW=$OPTARG
      ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      usage
      exit 1
      ;;
    : )
      echo "Invalid option: $OPTARG requires an argument" 1>&2
      usage
      exit 1
      ;;
  esac
done
shift $((OPTIND -1))

if [[ -z $GOALS || -z $CORRESPONDENCIES || -z $FINAL || -z $PATH_RAW ]]; then
  usage
  exit 1
fi

TMP=tmp.txt

touch $TMP

while IFS= read -r SAMPLE; do

        echo $SAMPLE,>> $TMP

        RAW=$(grep $SAMPLE $CORRESPONDENCIES | awk '{print $1}')

        FOLDER=$(grep $SAMPLE $CORRESPONDENCIES | awk '{print $5}')

        echo $RAW

        echo $FOLDER

        RAW_PATH_1=$(find ${PATH_RAW}/*${FOLDER}*/ -name ${RAW}* | grep -v "BKE" | grep -v "[R_]2\.")

        RAW_PATH_2=$(find ${PATH_RAW}/*${FOLDER}*/ -name ${RAW}* | grep -v "BKE" | grep "[R_]2\.")

        echo $RAW_PATH_1 | wc -w >> ver

        echo $RAW_PATH_2 | wc -w  >> ver

        echo $RAW_PATH_1 $RAW_PATH_2 >> $TMP

done < $GOALS

paste -d " " - - < $TMP > $FINAL

sed -i "s/, /,/g" $FINAL

rm $TMP
