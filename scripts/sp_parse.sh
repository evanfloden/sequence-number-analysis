LC_COLLATE=C sort --ignore-case -k1,1 data/spScores.*.csv > sp_sorted.csv
LC_COLLATE=C sort --ignore-case -k1,1 data/num_seqs.csv > seqs_sorted.csv
join sp_sorted.csv seqs_sorted.csv > sp_and_numseqs.csv
LC_COLLATE=C sort --ignore-case -k5,5 -k4,4 -k7,7 -k6,6n -k9,9n -k2,2n -k3,3n sp_and_numseqs.csv > spResults.csv

cat spResults.csv | cut -f 4,5,7,6 -d " " | LC_COLLATE=C sort --ignore-case -u > datasets.csv

while read line;
  do echo "";
  echo ""; 
  echo $line;
  grep "$line " spResults.csv > individual_datatsets.csv
  cat individual_datatsets.csv | cut -f 1,2 -d " " | LC_COLLATE=C sort --ignore-case -k1,1 -k2,2n -u > datapoints.csv
  cat individual_datatsets.csv | cut -f 1 -d " " | LC_COLLATE=C sort --ignore-case -k1,1 -u > ids.csv
  cat individual_datatsets.csv | cut -f 2 -d " " | LC_COLLATE=C sort --ignore-case -k1,1n -u > sizes.csv
  cat sizes.csv | tr '\n' ' '
  ID=''
  while read line2;
    do
    CURRENT_ID=$(echo $line2 | cut -d' ' -f1)
    if [ "$CURRENT_ID" != "$ID" ]; then
      echo " ";
      echo -ne "$CURRENT_ID ";
    fi
    SUM=$(grep "^$line2 " individual_datatsets.csv | cut -f 8 -d " " | paste -sd+ | bc);
    AVG=$( echo "scale=2; $SUM / 10"  | bc -l);
    echo -ne $AVG" ";
    ID=$CURRENT_ID;
  done < datapoints.csv
done <datasets.csv

rm datasets.csv datapoints.csv individual_datatsets.csv ids.csv sp_sorted.csv seqs_sorted.csv spResults.csv sp_and_numseqs.csv sizes.csv
 


