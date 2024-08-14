#!/bin/bash


# before running ./calculate_event_weights.sh ensure that the file is given permissions to run: chmod +x calculate_event_weights.sh

# This script calculates event weights for different datasets.
# It processes ROOT files containing event data and computes weights based on cross-section, 
# k-factor, filter efficiency, and luminosity. The results are saved to an output text file.
# The output file is to be used in run_analyses_sum_evntweight.sh in order to use the weights to analyse event sets (final paths).


# Step 1: Define the integrated luminosity in pb^-1.
# Luminosity is provided in fb^-1 (44.3 fb^-1) and is converted to pb^-1.
luminosity=$(echo "44.3 * 1000" | bc -l)

# Print the calculated luminosity for debugging purposes.
echo "Luminosity: $luminosity"

# Step 2: Define an associative array to store DSID, cross-section [pb], k-factor, and filter efficiency values.

data["700323"]="2221.3 1 0.0244"
data["700324"]="2221.3 1 0.13003"
data["700325"]="2221.3 1 0.846"

# Step 3: Define the base path where the ROOT files corresponding to these datasets are located.
# Adjust the path based on the dataset being analyzed (e.g., Zmumu, ttbar). This large file should contain background events for the DSIDs being analysed.
# The base path will have its whole tree explored, so ensure only relevant files are within the base path's tree
base_path="/path/to/background/Zmumu"

# Function: Convert a number from scientific notation to a decimal format.
# Required due to ROOT files contaniing scientific notation for weights.
sci_to_decimal() {
    local sci_number=$1
    local decimal_number=$(printf "%f\n" "$sci_number")
    echo "$decimal_number"
}

# Step 4: Function to extract the total weighted events from the sumWeights tree within a ROOT file.
# This function uses a ROOT macro to retrieve the value from the specified ROOT file.
get_total_weighted_events() {
    root_file=$1
    result=$(root -l -b -q "totWeight.C(\"$root_file\")" 2>/dev/null | tail -n 1)
    echo $result
}

# Step 5: Create a temporary ROOT macro to extract the "totalEventsWeighted" variable from the sumWeights tree.
# This macro will be used by the get_total_weighted_events function to retrieve the necessary data.
cat<<'EOF'> totWeight.C
#include <iostream>
#include <TFile.h>
#include <TTree.h>

void totWeight(const char* filename) {
    TFile file(filename);
    if (file.IsZombie()) {
        std::cerr<<"Cannot open file: "<<filename<<std::endl;
        return;
    }

    TTree* tree = (TTree*)file.Get("sumWeights");
    if (!tree) {
        std::cerr<<"Cannot find tree in file: "<<filename<<std::endl;
        return;
    }

    Float_t totalEventsWeighted;
    tree->SetBranchAddress("totalEventsWeighted", &totalEventsWeighted);
    tree->GetEntry(0); // Since there is only one entry

    std::cout<<totalEventsWeighted<<std::endl;
    file.Close();
}
EOF

# Step 6: Prepare the output file.
# If the output file already exists, clear its contents to ensure it only includes new data.
output_file="event_weights.txt"
> $output_file

# Step 7: Iterate over each DSID in the data array and compute the event weight for each corresponding ROOT file.
echo "DSID","ROOT FILE","EventWeight" >> $output_file
for dsid in "${!data[@]}"; do
    # Extract the cross-section, k-factor, and filter efficiency for the current DSID.
    IFS=' ' read -r crossSection k_factor filtEff <<< "${data[$dsid]}"
    echo "Processing DSID: $dsid with crossSection=$crossSection, k_factor=$k_factor, filtEff=$filtEff"

    # Loop through all ROOT files associated with this DSID and calculate the event weight for each file.
    for root_file in $(find $base_path/$dsid -name "*.root"); do
        echo "Processing ROOT file: $root_file"
        
        # Retrieve the total weighted events from the current ROOT file.
        total=$(get_total_weighted_events "$root_file")
        total_decimal=$(sci_to_decimal "$total")
        echo "Total from file $root_file: $total_decimal"
        
        echo "-----------------------------"
        echo "total Event Number for " $root_file ": " $total_decimal
        
        # Calculate the event weight using the provided formula.
        if [[ "$total_decimal" != "0" && -n "$total_decimal" ]]; then
            eventWeight=$(echo "scale=20; $crossSection * $filtEff * $k_factor * $luminosity / $total_decimal" | bc -l)
            echo "Calculated Event Weight for file $root_file: $eventWeight"
            
            # Append the event weight to the output file.
            echo "$dsid,$root_file, $eventWeight" >> "$output_file"
        else
            echo "File $root_file has zero total weighted events" >> "$output_file"
        fi
    done
done

# Step 8: Clean up the temporary ROOT macro created earlier.
rm totWeight.C
