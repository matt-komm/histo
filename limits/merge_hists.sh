n_lines=$(wc -l < config/procs.txt)
years=(2016 2017 2018)
echo "number of lines" $n_lines
mkdir -p limits/hists_merged

for ((i=1;i<=n_lines;i++)); do
    proc=$(awk "NR == $i" config/procs.txt)
    echo $i, $proc
    for year in "${years[@]}"; do
    echo $year
        rm limits/hists_merged/${proc}_${year}.root
        hadd limits/hists_merged/${proc}_${year}.root limits/hists/${proc}*_${year}.root
    done
done

for year in "${years[@]}"; do
    rm limits/hists_merged/${year}.root
    hadd limits/hists_merged/${year}.root limits/hists_merged/muon_${year}.root limits/hists_merged/electron_${year}.root limits/hists_merged/dyjets_${year}.root limits/hists_merged/topbkg_${year}.root limits/hists_merged/qcd_${year}.root limits/hists_merged/wjets_${year}.root limits/hists_merged/vgamma_${year}.root
done
