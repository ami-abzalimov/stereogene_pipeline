#!/bin/bash

# Define directories for input and output data
experiment_dir="./ex"
control_dir="./con"
p_value_dir="./p-value_peaks"
fold_dir="./fold_change_over_control"
alignments_dir="./alignments"

# Define genome sizes and blacklist files
genome_sizes="hg38.chrom.sizes"
blacklist_bed="ENCFF356LFX.bed"

# List of required commands for the script
required_commands=("samtools" "bamToBed" "genomeCoverageBed" "./bedGraphToBigWig" "./wigToBigWig" "bigwigCompare" "bedtools" "./bigWigToWig" "./bigBedToBed")

# Check if all required commands are available
for cmd in "${required_commands[@]}"; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "Error: Required command '$cmd' is not available."
        exit 1
    fi
done

# Function to process BAM files in a directory
process_bams() {
    # Input parameters: directory with BAM files, output BAM file, output BigWig file
    local bam_dir=$1
    local output_bam=$2
    local raw_experiment_bw=$3

    # Find all BAM files in the directory
    bam_files=($(find "$bam_dir" -name "*.bam"))
    
    # Check if any BAM files were found
    if [ ${#bam_files[@]} -eq 0 ]; then
        echo "Error: No BAM files found in directory '$bam_dir'."
        exit 1
    elif [ ${#bam_files[@]} -eq 1 ]; then
        # If only one BAM file, copy it to the output location
        cp "${bam_files[0]}" "$output_bam"
        raw_bam="${bam_files[0]}"
    else
        # If multiple BAM files, sort and merge them
        sorted_bams=()
        for bam in "${bam_files[@]}"; do
            sorted_bam="${bam%.bam}_sorted.bam"
            samtools sort -o "$sorted_bam" "$bam"
            sorted_bams+=("$sorted_bam")
        done

        # Merge sorted BAM files
        samtools merge -f "$output_bam" "${sorted_bams[@]}"
        raw_bam="$output_bam"

        # Remove intermediate sorted BAM files
        rm "${sorted_bams[@]}"
    fi
    
    # Convert BAM to BED and generate coverage files
    bamToBed -i "$raw_bam" > raw.bed
    genomeCoverageBed -bg -i raw.bed -g "$genome_sizes" > raw.bedGraph

    # Sort bedGraph and convert to BigWig format
    sort -k1,1 -k2,2n raw.bedGraph > raw_sorted.bedGraph
    mv raw_sorted.bedGraph raw.bedGraph

    ./bedGraphToBigWig raw.bedGraph "$genome_sizes" "$raw_experiment_bw"
    ./bigWigToWig "$raw_experiment_bw" raw_experiment.wig

    # Clean up temporary files
    rm raw.bed raw.bedGraph 
}

# Process experiment and control BAMs, saving results to alignments directory
mkdir -p alignments
process_bams "$experiment_dir" alignments/experiment_merged.bam alignments/raw_experiment.bw
process_bams "$control_dir" alignments/control_merged.bam alignments/raw_control.bw

# Index the merged BAM files
samtools index alignments/experiment_merged.bam
samtools index alignments/control_merged.bam

# Generate bedGraph and BigWig for experiment and control data
bamToBed -i alignments/experiment_merged.bam > alignments/experiment.bed
bamToBed -i alignments/control_merged.bam > alignments/control.bed

genomeCoverageBed -bg -i alignments/experiment.bed -g "$genome_sizes" > alignments/experiment.bedGraph
genomeCoverageBed -bg -i alignments/control.bed -g "$genome_sizes" > alignments/control.bedGraph

./bedGraphToBigWig alignments/experiment.bedGraph "$genome_sizes" alignments/filtered_experiment.bw
./bigWigToWig alignments/filtered_experiment.bw alignments/filtered_experiment.wig

./wigToBigWig alignments/experiment.bedGraph "$genome_sizes" alignments/experiment.bw
./wigToBigWig alignments/control.bedGraph "$genome_sizes" alignments/control.bw

# Normalize the experiment data using bigwigCompare
bigwigCompare -b1 alignments/experiment.bw -b2 alignments/control.bw --operation log2 -o alignments/normalized_experiment.bw
./bigWigToWig alignments/normalized_experiment.bw alignments/normalized_experiment.wig

# Clean up intermediate files for BAM processing
rm alignments/experiment_merged.bam alignments/experiment_merged.bam.bai alignments/control_merged.bam alignments/control_merged.bam.bai \
   alignments/experiment.bw alignments/control.bw alignments/normalized_experiment.wb alignments/experiment.bed \
   alignments/experiment.filtered.bedGraph alignments/experiment.bedGraph alignments/control.bed alignments/control.bedGraph

# Function to process .bw files in .wig format in a directory (for p-value_dir)
process_bw_to_wig() {
    local dir=$1
    bw_files=($(find "$dir" -name "*.bw"))
    
    if [ ${#bw_files[@]} -eq 0 ]; then
        echo "No .bw files found in $dir"
        return
    fi

    for bw_file in "${bw_files[@]}"; do
        wig_file="${bw_file%.bw}.wig"
        ./bigWigToWig "$bw_file" "$wig_file"
        echo "Converted $bw_file to $wig_file"
    done
}

# Function to process .bigBed files in .wig format in a directory (for fold_dir)
process_bigbed_to_wig() {
    local dir=$1
    bigbed_files=($(find "$dir" -name "*.bb"))

    if [ ${#bigbed_files[@]} -eq 0 ]; then
        echo "No .bigBed files found in $dir"
        return
    fi

    for bigbed_file in "${bigbed_files[@]}"; do
        bed_file="${bigbed_file%.bb}.bed"
        wig_file="${bigbed_file%.bb}.wig"

        # Convert bigBed -> bed -> wig (if required)
        ./bigBedToBed "$bigbed_file" "$bed_file"
        genomeCoverageBed -bg -i "$bed_file" -g "$genome_sizes" > temp.bedGraph
        
        # Sort bedGraph and convert to Wig
        sort -k1,1 -k2,2n temp.bedGraph > temp_sorted.bedGraph
        mv temp_sorted.bedGraph temp.bedGraph
        
        ./bedGraphToBigWig temp.bedGraph "$genome_sizes" temp.bw
        ./bigWigToWig temp.bw "$wig_file"

        # Remove temporary files and output result
        rm temp.bedGraph temp.bw "$bed_file"
        echo "Converted $bigbed_file to $wig_file"
    done
}

# Process p-value_dir (convert all .bw to .wig)
process_bw_to_wig "$p_value_dir"

# Process fold_dir (convert all .bb to .wig)
process_bigbed_to_wig "$fold_dir"

echo "Processing completed!"
