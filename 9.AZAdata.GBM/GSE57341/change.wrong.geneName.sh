command="sed"
for (( i = 15; i > 0; i-- )); do
	command="$command -e \"s/$i-Mar/MARCH$i/g\" -e \"s/$i-Sep/SEPT$i/g\" -e \"s/$i-Dec/DEC$i/g\" -e \"s/SEPT15/SEP15/g\""
done

eval "cat GSE57341_series_matrix.clean.probe.txt | $command > GSE57341_series_matrix.clean.probe.corrected.txt"
