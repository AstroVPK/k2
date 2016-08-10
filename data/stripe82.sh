#08/08/2016
#command to run script (no arguments):     ./stripe82.sh
#EDIT dir1 & dir2 to the paths to store lightcurves
#script to download vanderburg and Johnson calibrated lightcurves
#Parameters:
#-------------------------
#file: string
#         calls file with a list of EPIC IDS
#camp:  string
#        calls the campaign
#dir1 & dir2:  path 
#        directory path for downloaded fits lc file dump



file="./k08-stripe82quasars"
#file=$1
echo ${file}
#make directories
camp="8"
dir1="/Users/Jackster/Desktop/k2c0${camp}AGN/Vanderberg"
dir2="/Users/Jackster/Desktop/k2c0${camp}AGN/KTeam"
mkdir -p ${dir1}
mkdir -p ${dir2}

local_out=${dir1}
local_out2=${dir2}

while read line; do
	row=${line}
	fits_filename="hlsp_k2sff_k2_lightcurve_${row}-c0${camp}_kepler_v1_llc.fits"
	fits_address="http://archive.stsci.edu/missions/hlsp/k2sff/"
	fits_address="${fits_address}c0${camp}/${row:0:4}00000/${row:4:5}/${fits_filename}"
	#Andrew Vanderburg re-processed lightcurves
	wget -U 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:14.0) Gecko/20100101 Firefox/14.0.1'\
	--random-wait  -nd -q "$fits_address" -O "${local_out}/${fits_filename}" --content-disposition
	#Kepler Team Lightcurves for calculating errors
	fits_filename="ktwo${row}-c${camp}_llc.fits"
	fits_address="https://archive.stsci.edu/pub/k2/lightcurves/"
	fits_address="${fits_address}c${camp}/${row:0:4}00000/${row:4:5}/${fits_filename}"
	wget -U 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:14.0) Gecko/20100101 Firefox/14.0.1'\
	--random-wait  -nd -q "$fits_address" -O "${local_out2}/${fits_filename}" --content-disposition
	echo ${row}, ${line}
done < ${file}

#example wget to download fits file by parsing epic ID
#wget -q http://archive.stsci.edu/missions/hlsp/k2sff/c08/220100000/73631/hlsp_k2sff_k2_lightcurve_220173631-c08_kepler_v1_llc-default-aper.txt
#https://archive.stsci.edu/pub/k2/lightcurves/c8/220100000/33000/ktwo220133060-c08_llc.fits