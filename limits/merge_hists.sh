n_lines=$(wc -l < config/procs.txt)
years=(2016)
echo "number of lines" $n_lines
mkdir -p hists_merged

for ((i=1;i<=n_lines;i++)); do
    proc=$(awk "NR == $i" config/procs.txt)
    echo $i, $proc
    for year in "${years[@]}"; do
    echo $year
        rm hists_merged/${proc}_${year}.root
        hadd hists_merged/${proc}_${year}.root hists/${proc}*_${year}.root
    done
done

for year in "${years[@]}"; do
    rm hists_merged/${year}.root
    hadd hists_merged/${year}.root hists_merged/muon_${year}.root hists_merged/electron_${year}.root hists_merged/dyjets_${year}.root hists_merged/topbkg_${year}.root hists_merged/qcd_${year}.root hists_merged/wjets_${year}.root hists_merged/vgamma_${year}.root
done
