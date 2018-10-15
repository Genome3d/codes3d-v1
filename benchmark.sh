#!/usr/bin/env bash
NUM_RUNS=$1
FILENAME=$(basename $2)
CALL="${2%}"
RUNTIME=0
#python benchmark.py --preprocess codes3d/codes3d.py
for ARG in "${@:3}"; do
    CALL="$CALL $ARG";
done
TEMP_FILES=""
#mkdir error_logs
for RUN_NO in $(seq 1 $NUM_RUNS); do
    echo "BENCHMARK RUN #$RUN_NO/$NUM_RUNS"
    TEMP_FILES="$TEMP_FILES temp_$RUN_NO.log"
    START=$(date +%s)
    python /usr/local/lib/python2.7/dist-packages/kernprof.py -l $CALL #2> "error_logs/stderr_$RUN_NO"
    END=$(date +%s)
    rm -r codes3d_summary
    ((RUNTIME+= END-START))
    python /usr/local/lib/python2.7/dist-packages/line_profiler.py "${FILENAME}.lprof" > "temp_$RUN_NO.log";
    rm "${FILENAME%.py}.py.lprof"
    #if [ -s "error_logs/stderr_$RUN_NO" ]
    #    then
    #        continue
    #    else
    #        rm  "error_logs/stderr_$RUN_NO"
    #fi
done
echo "PROCESSING BENCHMARK RESULTS"
python benchmark.py --compile $TEMP_FILES 
rm $TEMP_FILES
#python benchmark.py --postprocess codes3d/codes3d.py
echo "PROCESSED BENCHMARK RESULTS"
AVG_RUN=$((RUNTIME/NUM_RUNS))
echo "AVERAGE RUNTIME: ${AVG_RUN}s" 
