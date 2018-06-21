#!/bin/bash -ex
echo "working directory =${PWD}"
    echo "home dir =${HOME}"

    # write the WDL script into a file
    cat >${DX_FS_ROOT}/source.wdl.uu64 <<'EOL'
dGFzayBCd2FNZW0gewogICAgRmlsZSBmYXN0cV8xIAogICAgRmlsZSBmYXN0cV8yIAogICAgU3RyaW5nPyBzYW1wbGVfaWQgID0gImVtcHR5aWQiCiAgICBJbnQ/IGluZGV4ICA9ICIwIgogICAgU3RyaW5nIGJhc2VfZmlsZV9uYW1lICA9ICJiYXNlIgogICAgU3RyaW5nIFJHX1BMICA9ICJJbGx1bWluYSIKICAgIEZpbGU/IHJlZl9mYXN0YSAgPSAiZHg6Ly9JbnRlbGxpc2VxIC0gUmVzb3VyY2VzOlJlZmVyZW5jZSBnZW5vbWVzL0dSQ2gzOC1uby1hbHQtYW5hbHlzaXMtc2V0L0dSQ2gzOC5ub19hbHRfYW5hbHlzaXNfc2V0LmZhIgogICAgRmlsZT8gcmVmX2Zhc3RhX2luZGV4ICA9ICJkeDovL0ludGVsbGlzZXEgLSBSZXNvdXJjZXM6UmVmZXJlbmNlIGdlbm9tZXMvR1JDaDM4LW5vLWFsdC1hbmFseXNpcy1zZXQvR1JDaDM4Lm5vX2FsdF9hbmFseXNpc19zZXQuZmEuZmFpIgogICAgRmlsZT8gcmVmX2RpY3QgID0gImR4Oi8vSW50ZWxsaXNlcSAtIFJlc291cmNlczpSZWZlcmVuY2UgZ2Vub21lcy9HUkNoMzgtbm8tYWx0LWFuYWx5c2lzLXNldC9HUkNoMzgubm9fYWx0X2FuYWx5c2lzX3NldC5kaWN0IgogICAgRmlsZT8gcmVmX2FsdCAKICAgIEZpbGU/IHJlZl9zYSAgPSAiZHg6Ly9JbnRlbGxpc2VxIC0gUmVzb3VyY2VzOlJlZmVyZW5jZSBnZW5vbWVzL0dSQ2gzOC1uby1hbHQtYW5hbHlzaXMtc2V0L0dSQ2gzOC5ub19hbHRfYW5hbHlzaXNfc2V0LmZhLnNhIgogICAgRmlsZT8gcmVmX2FtYiAgPSAiZHg6Ly9JbnRlbGxpc2VxIC0gUmVzb3VyY2VzOlJlZmVyZW5jZSBnZW5vbWVzL0dSQ2gzOC1uby1hbHQtYW5hbHlzaXMtc2V0L0dSQ2gzOC5ub19hbHRfYW5hbHlzaXNfc2V0LmZhLmFtYiIKICAgIEZpbGU/IHJlZl9id3QgID0gImR4Oi8vSW50ZWxsaXNlcSAtIFJlc291cmNlczpSZWZlcmVuY2UgZ2Vub21lcy9HUkNoMzgtbm8tYWx0LWFuYWx5c2lzLXNldC9HUkNoMzgubm9fYWx0X2FuYWx5c2lzX3NldC5mYS5id3QiCiAgICBGaWxlPyByZWZfYW5uICA9ICJkeDovL0ludGVsbGlzZXEgLSBSZXNvdXJjZXM6UmVmZXJlbmNlIGdlbm9tZXMvR1JDaDM4LW5vLWFsdC1hbmFseXNpcy1zZXQvR1JDaDM4Lm5vX2FsdF9hbmFseXNpc19zZXQuZmEuYW5uIgogICAgRmlsZT8gcmVmX3BhYyAgPSAiZHg6Ly9JbnRlbGxpc2VxIC0gUmVzb3VyY2VzOlJlZmVyZW5jZSBnZW5vbWVzL0dSQ2gzOC1uby1hbHQtYW5hbHlzaXMtc2V0L0dSQ2gzOC5ub19hbHRfYW5hbHlzaXNfc2V0LmZhLnBhYyIKICAgIFN0cmluZz8gYndhX2NvbW1hbmRsaW5lICA9ICItSyAxMDAwMDAwMDAgLXYgMyAtWSAkYmFzaF9yZWZfZmFzdGEiCiAgICBTdHJpbmcgZG9ja2VyX2ltYWdlIAogICAgSW50PyBudW1fY3B1ICA9ICIyIgogICAgY29tbWFuZCA8PDwKICAgIGVjaG8gIi0tLSBMT0c6IFNldHRpbmcgZW52aXJvbm1lbnRhbCB2YXJpYWJsZXMgLS0tIgogICAgUkdfSUQ9JChnemlwIC1jZCAke2Zhc3RxXzF9IHwgaGVhZCAtMSB8IGN1dCAtZCAnOicgLWYgMyw0IHwgc2VkICdzLzovXC4vZycpCiAgICBSR19QVT0iJFJHX0lEIiIuIiIke3NhbXBsZV9pZH0iCiAgICBSR19MQj0iJHtzYW1wbGVfaWR9IiIubGlicmFyeSIKICAgIFJHX1NNPSIke3NhbXBsZV9pZH0iCiAgICBlY2hvICItLS0gTE9HOiBHZXR0aW5nIGJ3YSB2ZXJzaW9uIC0tLSIKICAgIG1rZGlyIC90bXAvcmVzdWx0cwogICAgZHgtZG9ja2VyIHJ1biAtdiAvdG1wL3Jlc3VsdHM6L3RtcC9yZXN1bHRzIC0tZW50cnlwb2ludCAiL2Jpbi9zaCIgaW50ZWxsaXNlcS9id2E6bGF0ZXN0IC1jICJid2EgMj4mMSB8IGdyZXAgLWUgJ15WZXJzaW9uJyB8IHNlZCAncy9WZXJzaW9uOiAvLycgPiAvdG1wL3Jlc3VsdHMvdmVyc2lvbi50eHQiCiAgICBzZXQgLW8gLWUgLXggcGlwZWZhaWwKICAgIHNldCB0aGUgYmFzaCB2YXJpYWJsZSBuZWVkZWQgZm9yIHRoZSBjb21tYW5kLWxpbmUKICAgIGJhc2hfcmVmX2Zhc3RhPSR7cmVmX2Zhc3RhfQogICAgZHgtZG9ja2VyIHJ1biAtdiAvdG1wL3Jlc3VsdHM6L3RtcC9yZXN1bHRzIC0tZW50cnlwb2ludCAiL2Jpbi9zaCIgaW50ZWxsaXNlcS9id2E6bGF0ZXN0IC1jICJid2EgbWVtIC10ICR7bnVtX2NwdX0gLVIgIkBSR1x0SUQ6IiIkUkdfSUQiIlx0UFU6IiIkUkdfUFUiIlx0UEw6JHtSR19QTH1cdExCOiIiJFJHX0xCIiJcdFNNOiIiJFJHX1NNIiAke2J3YV9jb21tYW5kbGluZX0gJHtmYXN0cV8xfSAke2Zhc3RxXzJ9IDI+ICR7YmFzZV9maWxlX25hbWV9X3tpbmRleH0uYndhLnN0ZGVyci5sb2cgXAogICAgMT4gL3RtcC9yZXN1bHRzL215LmJhbSIKICAgICN8IHNhbWJsYXN0ZXIgMj4ge2Jhc2VfZmlsZV9uYW1lfV8ie2luZGV4fSIuc2FtYmxhc3Rlci5zdGRlcnIubG9nIFwKICAgICN8IHNhbXRvb2xzIHNvcnQgLUAgJHtudW1fY3B1fSAtID4ge2Jhc2VfZmlsZV9uYW1lfV8ie2luZGV4fSIuYWxuLmJhbSAyPiB7YmFzZV9maWxlX25hbWV9X3tpbmRleH0uc2FtdG9vbHMtc29ydC5zdGRlcnIubG9nCiAgICA+Pj4KICAgIG91dHB1dCB7CiAgICAgICAgRmlsZSBvdXRwdXRfYmFtID0gIi90bXAvcmVzdWx0cy9teS5iYW0iCiAgICAgICAgQXJyYXlbU3RyaW5nXSB2ZXJzaW9uID0gcmVhZF9saW5lcygiL3RtcC9yZXN1bHRzL3ZlcnNpb24udHh0IikKICAgIH0KfQ==
EOL
# decode the WDL script
base64 -d ${DX_FS_ROOT}/source.wdl.uu64 > ${DX_FS_ROOT}/source.wdl
# write the instance type DB
    cat >${DX_FS_ROOT}/instanceTypeDB.json <<'EOL'
{
  "instances": [{
    "diskGB": 28,
    "name": "mem1_ssd1_x2",
    "price": 1.0,
    "os": [["Ubuntu", "12.04"], ["Ubuntu", "14.04"], ["Ubuntu", "16.04"]],
    "memoryMB": 3766,
    "cpu": 2
  }, {
    "diskGB": 159,
    "name": "mem1_ssd2_x2",
    "price": 2.0,
    "os": [["Ubuntu", "12.04"], ["Ubuntu", "14.04"], ["Ubuntu", "16.04"]],
    "memoryMB": 3766,
    "cpu": 2
  }, {
    "diskGB": 27,
    "name": "mem2_ssd1_x2",
    "price": 3.0,
    "os": [["Ubuntu", "12.04"], ["Ubuntu", "14.04"], ["Ubuntu", "16.04"]],
    "memoryMB": 7225,
    "cpu": 2
  }, {
    "diskGB": 27,
    "name": "mem3_ssd1_x2",
    "price": 4.0,
    "os": [["Ubuntu", "12.04"], ["Ubuntu", "14.04"], ["Ubuntu", "16.04"]],
    "memoryMB": 15044,
    "cpu": 2
  }, {
    "diskGB": 77,
    "name": "mem1_ssd1_x4",
    "price": 5.0,
    "os": [["Ubuntu", "12.04"], ["Ubuntu", "14.04"], ["Ubuntu", "16.04"]],
    "memoryMB": 7225,
    "cpu": 4
  }, {
    "diskGB": 318,
    "name": "mem1_ssd2_x4",
    "price": 6.0,
    "os": [["Ubuntu", "12.04"], ["Ubuntu", "14.04"], ["Ubuntu", "16.04"]],
    "memoryMB": 7225,
    "cpu": 4
  }, {
    "diskGB": 72,
    "name": "mem2_ssd1_x4",
    "price": 7.0,
    "os": [["Ubuntu", "12.04"], ["Ubuntu", "14.04"], ["Ubuntu", "16.04"]],
    "memoryMB": 14785,
    "cpu": 4
  }, {
    "diskGB": 72,
    "name": "mem3_ssd1_x4",
    "price": 8.0,
    "os": [["Ubuntu", "12.04"], ["Ubuntu", "14.04"], ["Ubuntu", "16.04"]],
    "memoryMB": 30425,
    "cpu": 4
  }, {
    "diskGB": 157,
    "name": "mem1_ssd1_x8",
    "price": 9.0,
    "os": [["Ubuntu", "12.04"], ["Ubuntu", "14.04"], ["Ubuntu", "16.04"]],
    "memoryMB": 14785,
    "cpu": 8
  }, {
    "diskGB": 639,
    "name": "mem1_ssd2_x8",
    "price": 10.0,
    "os": [["Ubuntu", "12.04"], ["Ubuntu", "14.04"], ["Ubuntu", "16.04"]],
    "memoryMB": 14785,
    "cpu": 8
  }, {
    "diskGB": 147,
    "name": "mem2_ssd1_x8",
    "price": 11.0,
    "os": [["Ubuntu", "12.04"], ["Ubuntu", "14.04"], ["Ubuntu", "16.04"]],
    "memoryMB": 29905,
    "cpu": 8
  }, {
    "diskGB": 147,
    "name": "mem3_ssd1_x8",
    "price": 12.0,
    "os": [["Ubuntu", "12.04"], ["Ubuntu", "14.04"], ["Ubuntu", "16.04"]],
    "memoryMB": 61187,
    "cpu": 8
  }, {
    "diskGB": 302,
    "name": "mem1_ssd1_x16",
    "price": 13.0,
    "os": [["Ubuntu", "12.04"], ["Ubuntu", "14.04"], ["Ubuntu", "16.04"]],
    "memoryMB": 29900,
    "cpu": 16
  }, {
    "diskGB": 1278,
    "name": "mem1_ssd2_x16",
    "price": 14.0,
    "os": [["Ubuntu", "12.04"], ["Ubuntu", "14.04"], ["Ubuntu", "16.04"]],
    "memoryMB": 29900,
    "cpu": 16
  }, {
    "diskGB": 297,
    "name": "mem3_ssd1_x16",
    "price": 15.0,
    "os": [["Ubuntu", "12.04"], ["Ubuntu", "14.04"], ["Ubuntu", "16.04"]],
    "memoryMB": 122705,
    "cpu": 16
  }, {
    "diskGB": 637,
    "name": "mem1_ssd1_x32",
    "price": 16.0,
    "os": [["Ubuntu", "12.04"], ["Ubuntu", "14.04"], ["Ubuntu", "16.04"]],
    "memoryMB": 60139,
    "cpu": 32
  }, {
    "diskGB": 2877,
    "name": "mem1_ssd2_x36",
    "price": 17.0,
    "os": [["Ubuntu", "12.04"], ["Ubuntu", "14.04"], ["Ubuntu", "16.04"]],
    "memoryMB": 60139,
    "cpu": 36
  }, {
    "diskGB": 597,
    "name": "mem3_ssd1_x32",
    "price": 18.0,
    "os": [["Ubuntu", "12.04"], ["Ubuntu", "14.04"], ["Ubuntu", "16.04"]],
    "memoryMB": 245751,
    "cpu": 32
  }]
}
EOL

main() {
# Keep track of streaming files. Each such file
    # is converted into a fifo, and a 'dx cat' process
    # runs in the background.
    background_pids=()

    # evaluate input arguments, and download input files
    java -jar ${DX_FS_ROOT}/dxWDL.jar internal taskProlog ${DX_FS_ROOT}/source.wdl ${HOME}

    # uncomment to get more debugging outputs
    # ls -lR
    # cat ${HOME}/execution/meta/script

    # setup any file streams. Keep track of background
    # processes in the 'background_pids' array.
    # We 'source' the sub-script here, because we
    # need to wait for the pids. This can only be done
    # for child processes (not grand-children).
    if [[ -e ${HOME}/execution/meta/setup_streams ]]; then
       source ${HOME}/execution/meta/setup_streams > ${HOME}/execution/meta/background_pids.txt

       # reads the file line by line, and convert into a bash array
       mapfile -t background_pids < ${HOME}/execution/meta/background_pids.txt
       echo "Background processes ids: ${background_pids[@]}"
    fi

    # Run the shell script generated by the prolog.
    # Capture the stderr/stdout in files
    if [[ -e ${HOME}/execution/meta/script.submit ]]; then
        echo "docker submit script:"
        cat ${HOME}/execution/meta/script.submit
        ${HOME}/execution/meta/script.submit
    else
        /bin/bash ${HOME}/execution/meta/script
    fi

    # This section deals with streaming files.
    #
    # We cannot wait for all background processes to complete,
    # because the worker process may not read one of the fifo streams.
    # We want to make sure there were no abnormal terminations.
    #
    # Assumptions
    #  1) 'dx cat' returns zero status when a user reads only the beginning
    #  of a file
    for pid in ${background_pids[@]}; do
        p_status=0
        p_status=`ps --pid $pid --no-headers | wc -l`  || p_status=0

        if [[ $p_status == 0 ]]; then
            # the process is already dead, check correct exit status
            echo "wait $pid"
            rc=0
            wait $pid || rc=$?
            if [[ $rc != 0 ]]; then
                echo "Background download process $pid failed"
                exit $rc
            fi
        else
            echo "Warning: background download process $pid is still running."
            echo "Perhaps the worker process did not read it."
        fi
    done

    # Uncomment to see what the directory looks like after execution
    # ls -lR

    #  check return code of the script
    rc=`cat ${HOME}/execution/meta/rc`
    if [[ $rc != 0 ]]; then
        exit $rc
    fi

    # evaluate applet outputs, and upload result files
    java -jar ${DX_FS_ROOT}/dxWDL.jar internal taskEpilog ${DX_FS_ROOT}/source.wdl ${HOME}
}