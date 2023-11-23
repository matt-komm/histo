n_lines=$(wc -l < ../config/procs.txt)
#years=(2016 2017 2018)
years=($1)
histPath="hists_CWR3"
echo "number of lines" $n_lines
mkdir -p hists_merged3

for ((i=1;i<=n_lines;i++)); do
    proc=$(awk "NR == $i" ../config/procs.txt)
    for year in "${years[@]}"; do
        echo $i/$n_lines ${proc} $year
        files=($histPath/${proc}/*/*${year}*.root)
        if [ -e "${files[0]}" ]; then
            #rm hists_merged/${proc}_${year}.root
            if [ ! -f hists_merged3/${proc}_${year}.root ]; then
                hadd -f hists_merged3/${proc}_${year}.root $histPath/${proc}/*/*${year}*.root
            else
                echo "skip - output exists"
            fi
        else    
            echo "skip - input missing"
        fi
    done
done

for year in "${years[@]}"; do
    #rm hists_merged/${year}.root
    hadd -f hists_merged3/${year}.root hists_merged3/*_${year}.root
done

