# This script can be used to download all data sets used for 
# this manuscript. The first few variables here can be toggled 
# to determine which data sets to download. Note that the 
# simulated DHFR data is NOT available publicly for download 
# and must be requested from the authors.

HEWL_DATA=true # 84 GB on SBGrid
PDZ2_DATA=true # 50 GB on SBGrid
DHFR_DATA=true # 5  GB on SBGrid

if [ "$HEWL_DATA" = true ]; then
    cd hewl; mkdir data; cd data
    rsync -av rsync://data.sbgrid.org/10.15785/SBGRID/1118 .
    cd 1118 ; shasum -c files.sha
    cd ../../../
fi

if [ "$PDZ2_DATA" = true ]; then
    cd pdz2; mkdir data; cd data
    rsync -av rsync://data.sbgrid.org/10.15785/SBGRID/1116 .
    cd 1116 ; shasum -c files.sha
    cd ../../../
fi

if [ "$DHFR_DATA" = true ]; then
    mkdir dhfr; cd dhfr; mkdir data; cd data
    rsync -av rsync://data.sbgrid.org/10.15785/SBGRID/1117 .
    cd 1117 ; shasum -c files.sha
    cd ../../../
fi
