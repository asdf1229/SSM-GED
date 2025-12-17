#!/bin/bash
set -e

DATA_DIR=./test/dataset
EXEC_NAME=ssm_ged

#############################################
# CREATE TIMESTAMPED RESULT DIRECTORY
#############################################
timestamp=$(date +"%Y%m%d_%H%M%S")
RESULT_DIR="./result/${timestamp}"
mkdir -p "${RESULT_DIR}"

echo "===== RESULTS WILL BE SAVED TO: ${RESULT_DIR} ====="


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

        mkdir -p "${RESULT_DIR}/${dirname}/v1"
        mkdir -p "${RESULT_DIR}/${dirname}/v2"

        ########################
        # Run V1 + time (seconds, 2 decimals)
        ########################
        out1="${RESULT_DIR}/${dirname}/v1/t=${t}.txt"
        start1=$(date +%s%N)
        ./exec_v1 -d "$gfile" -q "$qfile" -t $t > "$out1"
        end1=$(date +%s%N)
        raw1=$(echo "scale=4; ($end1 - $start1) / 1000000000" | bc)
        time1=$(printf "%.2f" "$raw1")   # <-- 自动补齐前导零

        count1=$(grep "count:" "$out1" | awk '{print $2}')

        ########################
        # Run V2 + time (seconds, 2 decimals)
        ########################
        out2="${RESULT_DIR}/${dirname}/v2/t=${t}.txt"
        start2=$(date +%s%N)
        ./exec_v2 -d "$gfile" -q "$qfile" -t $t > "$out2"
        end2=$(date +%s%N)
        raw2=$(echo "scale=4; ($end2 - $start2) / 1000000000" | bc)
        time2=$(printf "%.2f" "$raw2")   # <-- 自动补齐前导零

        count2=$(grep "count:" "$out2" | awk '{print $2}')

        ########################
        # Compare count + show time
        ########################
        if [ "$count1" = "$count2" ]; then
            echo "    ✔ t=$t  SAME  (count=$count1) (v1=${time1}s , v2=${time2}s)"
        else
            echo "    ✘ t=$t  DIFF  (v1=$count1 , v2=$count2) (v1=${time1}s , v2=${time2}s)"
        fi

    done
done

echo "===== ALL DONE ====="
echo "Results saved at: ${RESULT_DIR}"
