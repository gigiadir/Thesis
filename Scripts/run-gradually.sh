#!/bin/bash

COMMANDS_FILE="/gpfs0/bgu-ofircohen/users/gigiadir/Scripts/submit_scDiffCom_jobs.sh"
MAX_JOBS=15
USER=$(whoami)

while IFS= read -r cmd; do
    # Skip shebang / blanks / comments from generated command files
    [[ -z "${cmd//[[:space:]]/}" ]] && continue
    [[ "$cmd" =~ ^[[:space:]]*# ]] && continue
    [[ ! "$cmd" =~ ^[[:space:]]*qsub ]] && continue
    # Wait until fewer than MAX_JOBS are running
    while true; do
        current_jobs=$(qstat -u "$USER" | tail -n +3 | wc -l)
        if [ "$current_jobs" -lt "$MAX_JOBS" ]; then
            break
        fi
        sleep 30
    done
    
    # Submit the job
    eval "$cmd"
    echo "Submitted: $cmd"
done < "$COMMANDS_FILE"

echo "All jobs submitted."