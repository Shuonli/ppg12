#!/bin/bash

echo "========================================"
echo "Starting the Golden Run List Generation"
echo "========================================"

# Create necessary directories
echo "Creating necessary directories..."
mkdir -p FileLists/
mkdir -p dst_list/
mkdir -p list/
echo "Directories created."
echo "----------------------------------------"

# Set the working directory
workplace=$(pwd)
echo "Working directory is set to $workplace"
echo "----------------------------------------"

# Step 1: Run the Python script to get the run numbers
echo "Step 1: Running Python script to generate initial run lists..."

python_output=$(python3 << EOF
import pyodbc

def get_all_run_numbers(cursor):
    query = """
    SELECT runnumber
    FROM datasets
    WHERE dsttype='DST_CALO_run2pp' and dataset = 'ana437_2024p007'
    GROUP BY runnumber
    HAVING SUM(events) >= 1000000 AND runnumber >= 47289;
    """
    cursor.execute(query)
    run_numbers = [row.runnumber for row in cursor.fetchall()]
    return run_numbers

def get_emcal_auto_golden_run_numbers(file_catalog_run_numbers, production_cursor):
    query = """
    SELECT runnumber
    FROM goodruns
    WHERE (emcal_auto).runclass = 'GOLDEN'
    """
    production_cursor.execute(query)
    emcal_auto_golden_run_numbers = {row.runnumber for row in production_cursor.fetchall()}
    return list(
        emcal_auto_golden_run_numbers.intersection(
            set(file_catalog_run_numbers)
        )
    )

def get_ihcal_auto_golden_run_numbers(file_catalog_run_numbers, production_cursor):
    query = """
    SELECT runnumber
    FROM goodruns
    WHERE (ihcal_auto).runclass = 'GOLDEN'
    """
    production_cursor.execute(query)
    ihcal_auto_golden_run_numbers = {row.runnumber for row in production_cursor.fetchall()}
    return list(
        ihcal_auto_golden_run_numbers.intersection(
            set(file_catalog_run_numbers)
        )
    )

def get_ohcal_auto_golden_run_numbers(file_catalog_run_numbers, production_cursor):
    query = """
    SELECT runnumber
    FROM goodruns
    WHERE (ohcal_auto).runclass = 'GOLDEN'
    """
    production_cursor.execute(query)
    ohcal_auto_golden_run_numbers = {row.runnumber for row in production_cursor.fetchall()}
    return list(
        ohcal_auto_golden_run_numbers.intersection(
            set(file_catalog_run_numbers)
        )
    )

def main():
    # Connect to the FileCatalog database
    file_catalog_conn = pyodbc.connect("DSN=FileCatalog;UID=phnxrc;READONLY=True")
    file_catalog_cursor = file_catalog_conn.cursor()

    # Get unique run numbers with at least 1 million total events
    file_catalog_run_numbers = get_all_run_numbers(file_catalog_cursor)
    file_catalog_run_numbers.sort()
    with open('list/list_runnumber_all.txt', 'w') as f:
        for run_number in file_catalog_run_numbers:
            f.write(f"{run_number}\n")
    print(f"TOTAL_RUNS:{len(file_catalog_run_numbers)}")

    # Close the FileCatalog database connection
    file_catalog_conn.close()

    # Connect to the Production database
    production_conn = pyodbc.connect("DSN=Production_write")
    production_cursor = production_conn.cursor()

    # Filter 'GOLDEN' run numbers
    emcal_auto_golden_run_numbers = get_emcal_auto_golden_run_numbers(file_catalog_run_numbers, production_cursor)
    ihcal_auto_golden_run_numbers = get_ihcal_auto_golden_run_numbers(file_catalog_run_numbers, production_cursor)
    ohcal_auto_golden_run_numbers = get_ohcal_auto_golden_run_numbers(file_catalog_run_numbers, production_cursor)

    # Get the intersection of all sets
    golden_run_numbers = list(
        set(emcal_auto_golden_run_numbers).intersection(
            set(ihcal_auto_golden_run_numbers),
            set(ohcal_auto_golden_run_numbers)
        )
    )
    golden_run_numbers.sort()
    with open('list/Full_ppGoldenRunList.txt', 'w') as f:
        for run_number in golden_run_numbers:
            f.write(f"{run_number}\n")
    print(f"COMBINED_GOLDEN_RUNS:{len(golden_run_numbers)}")

    # Close the Production database connection
    production_conn.close()

if __name__ == "__main__":
    main()
EOF
)

echo "Python script execution completed."
echo "----------------------------------------"

# Parse the output to get the counts
total_runs=$(echo "$python_output" | grep 'TOTAL_RUNS' | cut -d':' -f2)
combined_golden_runs=$(echo "$python_output" | grep 'COMBINED_GOLDEN_RUNS' | cut -d':' -f2)

echo "Summary of initial counts:"
echo "Total runs after firmware fix and >1M events: $total_runs"
echo "Number of runs passing Calo QA (all three GOLDEN statuses): $combined_golden_runs"
echo "----------------------------------------"

# Step 2: Check that the initial golden run list exists
echo "Step 2: Checking for the initial golden run list..."

golden_run_list="list/Full_ppGoldenRunList.txt"
if [[ ! -f "$golden_run_list" ]]; then
    echo "[ERROR] The file $golden_run_list does not exist. Please check the path and try again."
    exit 1
fi

echo "Initial golden run list found: $golden_run_list"
echo "----------------------------------------"

# Function to get total GL1 raw events for a list of runs
get_total_events() {
    input_file=$1
    total_events=0
    batch_size=100
    run_numbers=()
    while IFS= read -r runnumber; do
        if [[ -n "$runnumber" ]]; then
            run_numbers+=("$runnumber")
            if [[ ${#run_numbers[@]} -ge $batch_size ]]; then
                run_list=$(printf ",%s" "${run_numbers[@]}")
                run_list=${run_list:1}
                query="SELECT SUM(raw) FROM gl1_scalers WHERE runnumber IN (${run_list});"
                result=$(psql -h sphnxdaqdbreplica -d daq -t -c "$query")
                events=$(echo "$result" | xargs)
                if [[ "$events" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
                    total_events=$(echo "$total_events + $events" | bc)
                fi
                run_numbers=()
            fi
        fi
    done < "$input_file"
    # Process remaining runs
    if [[ ${#run_numbers[@]} -gt 0 ]]; then
        run_list=$(printf ",%s" "${run_numbers[@]}")
        run_list=${run_list:1}
        query="SELECT SUM(raw) FROM gl1_scalers WHERE runnumber IN (${run_list});"
        result=$(psql -h sphnxdaqdbreplica -d daq -t -c "$query")
        events=$(echo "$result" | xargs)
        if [[ "$events" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
            total_events=$(echo "$total_events + $events" | bc)
        fi
    fi
    echo "$total_events"
}

#################################
# Version 1: Calo QA + run time > 5 mins + livetime > 80%
#################################

echo "----------------------------------------"
echo "Processing Version 1: Calo QA + run time > 5 mins + livetime > 80%"
echo "----------------------------------------"

# Step 3a (Version 1): Apply run duration filter (>300 seconds)

input_file="list/Full_ppGoldenRunList.txt"  # Runs after Calo QA
output_file_duration_v1="list/list_runnumber_runtime_v1.txt"
> "$output_file_duration_v1"

total_runs_duration_v1=0
runs_dropped_runtime_v1=0

while IFS= read -r runnumber; do
    if [[ -z "$runnumber" ]]; then
        continue
    fi

    query="SELECT runnumber, EXTRACT(EPOCH FROM (ertimestamp - brtimestamp)) AS duration FROM run WHERE runnumber = ${runnumber};"

    result=$(psql -h sphnxdaqdbreplica -d daq -t -c "$query" | tr -d '[:space:]')

    duration=$(echo $result | awk -F '|' '{print $2}')

    if [[ $duration =~ ^[0-9]+([.][0-9]+)?$ ]]; then
        if (( $(echo "$duration > 300" | bc -l) )); then
            echo "$runnumber" >> "$output_file_duration_v1"
            total_runs_duration_v1=$((total_runs_duration_v1+1))
        else
            runs_dropped_runtime_v1=$((runs_dropped_runtime_v1+1))
        fi
    else
        runs_dropped_runtime_v1=$((runs_dropped_runtime_v1+1))
    fi

done < "$input_file"

echo "Total runs after run duration cut (>5 mins): $total_runs_duration_v1"
echo "Number of runs dropped due to run duration cut: $runs_dropped_runtime_v1"
echo "----------------------------------------"

# Step 4a (Version 1): Apply live/raw ratio filter for trigger index 10 (>80%)

input_file="$output_file_duration_v1"
output_file_livetime_v1="list/list_runnumber_livetime_v1.txt"
bad_file_livetime_v1="list/list_runnumber_bad_livetime_v1.txt"
> "$output_file_livetime_v1"
> "$bad_file_livetime_v1"

total_runs_livetime_v1=0
runs_dropped_livetime_v1=0

while IFS= read -r runnumber; do
    if [[ -z "$runnumber" ]]; then
        continue
    fi

    index_to_check=10

    query="SELECT index, raw, live FROM gl1_scalers WHERE runnumber = ${runnumber} AND index = ${index_to_check};"

    result=$(psql -h sphnxdaqdbreplica -d daq -t -c "$query")

    index_pass=false

    while IFS='|' read -r index raw live; do
        index=$(echo "$index" | xargs)
        raw=$(echo "$raw" | xargs)
        live=$(echo "$live" | xargs)

        if [[ "$raw" =~ ^[0-9]+$ ]] && [[ "$live" =~ ^[0-9]+$ ]] && [ "$raw" -ne 0 ]; then
            ratio=$(echo "scale=2; $live / $raw * 100" | bc -l)
            #echo "Run $runnumber, Index $index, Live = $live, Raw = $raw, Live/Raw Ratio = $ratio%"

            if [ "$index" -eq "$index_to_check" ]; then
                if (( $(echo "$ratio >= 80" | bc -l) )); then
                    index_pass=true
                fi
            fi
        fi
    done <<< "$result"

    if [[ "$index_pass" == true ]]; then
        echo "$runnumber" >> "$output_file_livetime_v1"
        total_runs_livetime_v1=$((total_runs_livetime_v1+1))
    else
        echo "$runnumber" >> "$bad_file_livetime_v1"
        runs_dropped_livetime_v1=$((runs_dropped_livetime_v1+1))
    fi

done < "$input_file"

echo "Total runs after livetime cut (>80% for trigger index 10): $total_runs_livetime_v1"
echo "Number of runs dropped due to livetime cut: $runs_dropped_livetime_v1"
echo "----------------------------------------"

# Step 5a (Version 1): Remove runs with missing bad tower maps

input_file="$output_file_livetime_v1"
output_file_final_v1="FileLists/Full_ppGoldenRunList_Version1.txt"

total_runs_input=$(wc -l < "$input_file")
echo "Total runs before removing runs with missing bad tower maps (Version 1): $total_runs_input"

bad_tower_runs=$(find /cvmfs/sphenix.sdcc.bnl.gov/calibrations/sphnxpro/cdb/CEMC_BadTowerMap -name "*p0*" | cut -d '-' -f2 | cut -d c -f1 | sort | uniq)
echo "$bad_tower_runs" > bad_tower_runs.txt

grep -Ff bad_tower_runs.txt "$input_file" > "$output_file_final_v1"

removed_runs=$(comm -23 <(sort "$input_file") <(sort "$output_file_final_v1"))
runs_dropped_bad_tower_v1=$(echo "$removed_runs" | wc -l)
total_runs_after_all_cuts_v1=$(wc -l < "$output_file_final_v1")

echo "Number of runs dropped due to missing bad tower maps (Version 1): $runs_dropped_bad_tower_v1"
echo "Total runs after removing runs with missing bad tower maps (Version 1): $total_runs_after_all_cuts_v1"
echo "----------------------------------------"

# Clean up temporary files
rm bad_tower_runs.txt

# Copy the final run numbers to dst_list folder
cp "$output_file_final_v1" dst_list/Final_RunNumbers_After_All_Cuts.txt
echo "Final run numbers after all cuts have been copied to dst_list/Final_RunNumbers_After_All_Cuts.txt"
echo "----------------------------------------"

# Step 6: Create .list file for CreateDstList.pl command
echo "Creating .list file for CreateDstList.pl..."
cp "$output_file_final_v1" Full_ppGoldenRunList_Version1.list
echo ".list file created."
echo "----------------------------------------"

# Step 7: Remove existing list files in dst_list directory
echo "Removing existing list files in dst_list directory..."
rm -f dst_list/*.list
echo "Existing list files removed."
echo "----------------------------------------"

# Step 8: Create DST lists
echo "Changing to dst_list directory to generate DST lists..."
cd dst_list

echo "Running CreateDstList.pl to generate the DST list..."

CreateDstList.pl --build ana437 --cdb 2024p007 DST_CALO_run2pp --list ../Full_ppGoldenRunList_Version1.list

echo "DST list generated and saved to dst_list."
echo "----------------------------------------"

cd ..

# Calculate total events at each stage
total_events_initial=$(get_total_events 'list/list_runnumber_all.txt')
total_events_calo_qa=$(get_total_events 'list/Full_ppGoldenRunList.txt')
total_events_after_runtime=$(get_total_events 'list/list_runnumber_runtime_v1.txt')
total_events_after_livetime=$(get_total_events 'list/list_runnumber_livetime_v1.txt')
total_events_after_all_cuts=$(get_total_events 'FileLists/Full_ppGoldenRunList_Version1.txt')

# Compute percentages relative to initial totals
percent_runs_calo_qa=$(echo "scale=2; $combined_golden_runs / $total_runs * 100" | bc)
percent_runs_after_runtime=$(echo "scale=2; $total_runs_duration_v1 / $total_runs * 100" | bc)
percent_runs_after_livetime=$(echo "scale=2; $total_runs_livetime_v1 / $total_runs * 100" | bc)
percent_runs_after_badtower=$(echo "scale=2; $total_runs_after_all_cuts_v1 / $total_runs * 100" | bc)

percent_events_calo_qa=$(echo "scale=2; $total_events_calo_qa / $total_events_initial * 100" | bc)
percent_events_after_runtime=$(echo "scale=2; $total_events_after_runtime / $total_events_initial * 100" | bc)
percent_events_after_livetime=$(echo "scale=2; $total_events_after_livetime / $total_events_initial * 100" | bc)
percent_events_after_badtower=$(echo "scale=2; $total_events_after_all_cuts / $total_events_initial * 100" | bc)

# Calculate percentages lost at each step
percent_runs_lost_calo_qa=$(echo "scale=2; 100 - $percent_runs_calo_qa" | bc)
percent_runs_lost_after_runtime=$(echo "scale=2; $percent_runs_calo_qa - $percent_runs_after_runtime" | bc)
percent_runs_lost_after_livetime=$(echo "scale=2; $percent_runs_after_runtime - $percent_runs_after_livetime" | bc)
percent_runs_lost_after_badtower=$(echo "scale=2; $percent_runs_after_livetime - $percent_runs_after_badtower" | bc)

percent_events_lost_calo_qa=$(echo "scale=2; 100 - $percent_events_calo_qa" | bc)
percent_events_lost_after_runtime=$(echo "scale=2; $percent_events_calo_qa - $percent_events_after_runtime" | bc)
percent_events_lost_after_livetime=$(echo "scale=2; $percent_events_after_runtime - $percent_events_after_livetime" | bc)
percent_events_lost_after_badtower=$(echo "scale=2; $percent_events_after_livetime - $percent_events_after_badtower" | bc)

# Final Summary
echo "========================================"
echo "Final Summary:"

printf "%-40s | %-9s | %-12s | %-15s | %-12s\n" "Stage" "Runs" "% of initial" "Events" "% of initial"
echo "------------------------------------------|-----------|--------------|-----------------|--------------"
printf "%-40s | %-9s | %-12s | %-15s | %-12s\n" "1) After firmware fix and >1M events" "$total_runs" "100%" "$total_events_initial" "100%"
printf "%-40s | %-9s | %-12s | %-15s | %-12s\n" "2) & pass Calo QA" "$combined_golden_runs" "${percent_runs_calo_qa}%" "$total_events_calo_qa" "${percent_events_calo_qa}%"
printf "%-40s | %-9s | %-12s | %-15s | %-12s\n" "3) & > 5 mins" "$total_runs_duration_v1" "${percent_runs_after_runtime}%" "$total_events_after_runtime" "${percent_events_after_runtime}%"
printf "%-40s | %-9s | %-12s | %-15s | %-12s\n" "4) & livetime > 80% of MB trigger" "$total_runs_livetime_v1" "${percent_runs_after_livetime}%" "$total_events_after_livetime" "${percent_events_after_livetime}%"
printf "%-40s | %-9s | %-12s | %-15s | %-12s\n" "5) & have bad tower maps" "$total_runs_after_all_cuts_v1" "${percent_runs_after_badtower}%" "$total_events_after_all_cuts" "${percent_events_after_badtower}%"

echo "========================================"
echo ""
echo "Percentage of runs lost at each step:"
echo "After Calo QA: ${percent_runs_lost_calo_qa}% of runs lost"
echo "After run time > 5 mins: ${percent_runs_lost_after_runtime}% of runs lost"
echo "After livetime >80%: ${percent_runs_lost_after_livetime}% of runs lost"
echo "After removing runs without bad tower maps: ${percent_runs_lost_after_badtower}% of runs lost"

echo ""
echo "Percentage of events lost at each step:"
echo "After Calo QA: ${percent_events_lost_calo_qa}% of events lost"
echo "After run time > 5 mins: ${percent_events_lost_after_runtime}% of events lost"
echo "After livetime >80%: ${percent_events_lost_after_livetime}% of events lost"
echo "After removing runs without bad tower maps: ${percent_events_lost_after_badtower}% of events lost"
echo "========================================"
