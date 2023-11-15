#!/bin/bash

rootOld=/vols/Scratch/abaram/shirley
rootNew=/vols/Scratch/abaram/shirley/alon

# copy preprocessed functional (smoothed and cleaned) data
cd ${rootOld}/data
counter=0
for s in sub-??; do
    let counter++
    # zero pad counter
    printf -v new_s "%02d" $counter
    for i in 1 2 3 4; do
        echo $s $counter $new_s
        $i
        mkdir ${rootNew}/preproc/sub-${new_s}/run-0${i}
        scp ${rootOld}/data/${s}/sess0${i}/r_cleaned_smoothed_0* ${rootNew}/preproc/sub-${new_s}/run-0${i}
    done
done

# copy strucural data
cd ${rootOld}/data
counter=0
for s in sub-??; do
    let counter++
    # zero pad counter
    printf -v new_s "%02d" $counter
    mkdir -p ${rootNew}/struct/sub-${new_s}
    scp ${rootOld}/data/${s}/structural/* ${rootNew}/struct/sub-${new_s}/
done

# copy behavioural data
cd $rootOld/beh
counter=0
for s_beh in s??Behave; do
    let counter++
    s=${s_beh:1:2}                  # e.g. 18
    printf -v new_s "%02d" $counter # e.g. new_s=01
    mkdir $rootNew/beh/sub-${new_s}
    for r in 1 2 3 4; do # runs/sessions
        scp $rootOld/beh/$s_beh/*_S${s}_s${r}.mat $rootNew/beh/sub-${new_s}/run-0${r}.mat
    done
done
