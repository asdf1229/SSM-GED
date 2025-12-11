#!/bin/bash

DATA_DIR=./test/dataset
EXEC_NAME=ssm_ged

#############################################
# BUILD VERSION 1
#############################################
echo "===== Building Version 1 (Approximate_Matching) ====="
rm -rf build_v1
mkdir build_v1
cd build_v1
cmake -DAPPROXIMATE_MATCHING_V2=OFF ..
make -j
cp ${EXEC_NAME} ../exec_v1
cd ..

#############################################
# BUILD VERSION 2
#############################################
echo "===== Building Version 2 (Approximate_Matching_V2) ====="
rm -rf build_v2
mkdir build_v2
cd build_v2
cmake -DAPPROXIMATE_MATCHING_V2=ON ..
make -j
cp ${EXEC_NAME} ../exec_v2
cd ..

#############################################
# RUN BOTH VERSIONS + SAVE OUTPUT + COMPARE COUNT
#############################################
echo "===== Running and Comparing count ====="

for folder in ${DATA_DIR}/G*/; do
    dirname=$(basename "$folder")   # G1, G2, ...
    gfile="${folder}/graph_g.txt"
    qfile="${folder}/graph_q.txt"

    echo "Dataset: $dirname"

    for t in {0..10}; do
        echo "  Running -t $t ..."

        mkdir -p ./result/${dirname}/v1
        mkdir -p ./result/${dirname}/v2

        ########################
        # Run V1
        ########################
        out1=./result/${dirname}/v1/t=${t}.txt
        ./exec_v1 -d "$gfile" -q "$qfile" -t $t > "$out1"
        count1=$(grep "count:" "$out1" | awk '{print $2}')

        ########################
        # Run V2
        ########################
        out2=./result/${dirname}/v2/t=${t}.txt
        ./exec_v2 -d "$gfile" -q "$qfile" -t $t > "$out2"
        count2=$(grep "count:" "$out2" | awk '{print $2}')

        ########################
        # Compare count
        ########################
        if [ "$count1" = "$count2" ]; then
            echo "    ✔ t=$t  SAME  (count=$count1)"
        else
            echo "    ✘ t=$t  DIFF  (v1=$count1 , v2=$count2)"
        fi

    done
done

echo "===== ALL DONE ====="
