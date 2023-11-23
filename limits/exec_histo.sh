#!/bin/bash
hostname
date
cd /nfs/dust/cms/user/mkomm/HNL/histo/limits
export PATH=/nfs/dust/cms/user/mkomm/HNL/histo/env/bin:$PATH
source activate hnl
which python
which conda
function run_fct() {
    python -u make_hists.py "$@" || return 0
}
run_fct "$@"
date

