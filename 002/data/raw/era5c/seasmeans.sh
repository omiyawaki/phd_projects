rea=("era5c")
typenames=("rad" "stf" "hydro")
yr_span="1979_2019"

# save current path
cwd=$(pwd)

for typename in ${typenames[@]}; do
    echo ${typename}

    cd ${cwd}/${typename}

    filename="${rea}_${typename}_${yr_span}"

    cdo -seasmean -selseas,djf ${filename}.nc ${filename}.djfmean.nc
    cdo -seasmean -selseas,mam ${filename}.nc ${filename}.mammean.nc
    cdo -seasmean -selseas,jja ${filename}.nc ${filename}.jjamean.nc
    cdo -seasmean -selseas,son ${filename}.nc ${filename}.sonmean.nc
done
