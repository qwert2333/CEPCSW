#!/bin/bash
FILES=($(ls RecSample/Digi_nnHgg_*.root))

# Loop through each file with prefix 'Digi_nnHgg_'
for f in "${FILES[@]}"
do
    # Get the INDEX value from the file name
    INDEX=$(echo $f | sed 's/[^0-9]*//g')
    
    # Run the command with the corresponding INDEX value
    python Root2h5.py $f HDF5/Digi_nnHgg_$INDEX.h5
    #echo $f
    #echo $INDEX
done

