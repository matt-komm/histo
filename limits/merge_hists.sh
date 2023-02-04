n_lines=$(wc -l < ../config/procs.txt)
#years=(2016 2017 2018)
years=(2016)
histPath="hists"
echo "number of lines" $n_lines
mkdir -p hists_merged

for ((i=1;i<=n_lines;i++)); do
    proc=$(awk "NR == $i" ../config/procs.txt)
    echo $i, $proc
    for year in "${years[@]}"; do
        echo $year
        files=($histPath/${proc}*_${year}*.root)
        if [ -e "${files[0]}" ]; then
            #rm hists_merged/${proc}_${year}.root
            hadd -f hists_merged/${proc}_${year}.root $histPath/${proc}*_${year}*.root
        else
            echo "skip"
        fi
    done
done

for year in "${years[@]}"; do
    #rm hists_merged2/${year}.root
    hadd -f hists_merged/${year}.root hists_merged/*_${year}.root
done

