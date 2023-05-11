#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

        full_evap=${cwd}/${model}/evspsbl_${common}.nc
        full_prec=${cwd}/${model}/pr_${common}.nc
        full_qtend=${cwd}/${model}/qtend_${common}.nc
        full_qaht=${cwd}/${model}/qaht_${common}.nc

        python ${cwd}/make_qaht.py ${full_evap} ${full_prec} ${full_qtend} ${full_qaht} 
    fi

done
