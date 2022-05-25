n_lines=$(wc -l < ../config/procs.txt)
years=(2016 2017 2018)
#years=(2018)
echo "number of lines" $n_lines
mkdir -p hists_merged

for ((i=1;i<=n_lines;i++)); do
    proc=$(awk "NR == $i" ../config/procs.txt)
    echo $i, $proc
    for year in "${years[@]}"; do
    echo $year
        #rm hists_merged/${proc}_${year}.root
        hadd hists_merged/${proc}_${year}.root hists/${proc}*_${year}_chunk*.root
    done
done

#for year in "${years[@]}"; do
#    #rm hists_merged/${year}.root
#    #hadd hists_merged/${year}.root hists_merged2/muon_${year}.root hists_merged2/electron_${year}.root hists_merged2/dyjets_${year}.root hists_merged2/topbkg_${year}.root hists_merged2/qcd_${year}.root hists_merged2/wjets_${year}.root hists_merged2/vgamma_${year}.root
#done
