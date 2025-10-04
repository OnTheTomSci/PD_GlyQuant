#!/bin/bash
pwd
# Create necessary directories
mkdir -p input_data output_data figures functions

# Print confirmation message
echo "Created directories:"
echo "- input_data"
echo "- output_data" 
echo "- figures"
echo "- functions"

cd output_data
mkdir -p peptidegroups_intensity
cd figures
mkdir -p peptidegroups_intensity


# Run peptidegroups_quant
Rscript peptidegroups_quant.R
echo "peptidegroups_quant.R completed successfully"


cd output_data
mkdir -p psm_counts
cd figures
mkdir -p psm_counts
# Run count_stats
Rscript /functions/count_stats.R
echo "count_stats.R completed successfully"