#!/bin/bash

# --------------------------------------------------------------------------- #
#
# ATAC-seq footprint pipeline                                         
# 
# --------------------------------------------------------------------------- #

#Flag-handling and input-handling:

#initialize variables
need_help=false
sam_file='unassigned'
sam_dir=''
out_dir=$( pwd)
starting_dir=$( pwd)
name='name'
genome="unassigned"
motifs="unassigned"

while getopts 'abho:n:s:g:m:v' flag; do
  case "${flag}" in
    s) sam_file="${OPTARG}" ;;
    o) out_dir="${OPTARG}" ;;
    g) genome="${OPTARG}" ;; 
    m) motifs="${OPTARG}" ;;
    n) name="${OPTARG}" ;;
    h) need_help=true ;;
    *) error "Unexpected option ${flag}" ;;
    esac
done

#check if all the helper scripts are in the right spot:
fimo_helper=false
filter_helper=false
centipede_helper=false
blacklist_check=false

#checking fimo helper
if [ -e ~/bin/ATAC_Pipeline/split_fimo_by_motif.py ]
then
    $filter_helper=true
fi

#checking filter helper
if [ -e ~/bin/ATAC_Pipeline/filter_by_postprobs.py ]
then
    $fimo_helper=true
fi

#checking centipede helper
if [ -e ~/bin/ATAC_Pipeline/centipede.R ]
then
    $centipede_helper=true
fi

#check blacklist file
if [ -e ~/bin/ATAC_Pipeline/ENCODE.hg19.blacklist.bed ]
then
    $blacklist_check=true
fi

#handle the help flag -h:
if [ "$need_help" = true ]
then
    echo "ATAC-seq footprint pipeline"
    echo ""
    echo "Usage: bash </path/to/atac_centipede.sh> [options]"
    echo ""
    echo "Tips:"
    echo "  Use absolute paths for all files and directories"
    echo "  Store helper scripts ~/bin/ATAC_Pipeline"
    echo "  Use the most recent versions of dependencies (bedtools, samtools, fimo...)"
    echo ""
    echo "  Options:"
    echo "      -s <sam file>"
    echo "      -o <path/to/output/dir>"
    echo "      -n <name> "
    echo "      -m <meme formatted motif file>"
    echo "      -g <referece_genome in fasta format> "
    echo "      -h <help>"
    echo ""

#throw an error message if necessary arguments are not present
elif [ $sam_file == 'unassigned' ] || [ $genome == 'unassigned' ] || [ $motifs == 'unassigned' ]
then
    echo "Error: not all required parameters have been initialized. -s, -g, and -m flags are required."
    echo "use -h option for help"

#are the helper scripts in the right spot?
elif [ $filter_helper = false ] || [ $fimo_helper = false ] || [ $centipede_helper = false ] || [ $blacklist_check = false ]
then
    echo "Error: not all helper scripts are in ~/bin/ATAC_Pipeline"
    echo "use -h option for help"

#processing
else 

    #make a new directory to keep all the mapped reads:
    mkdir "$out_dir/mapped_reads"
    
    #filter based on read quality, remove mitochondrial reads
    ext=".filtered.tmp.sam"
    filtered_sam="$name$ext"
    samtools view -h -f 2 -q 30 "$sam_file" | grep -v "chrM" > "$out_dir/mapped_reads/$filtered_sam"

    #convert from sam to bam:
    >&2 echo "convert from sam to bam"
    ext=".filtered.tmp.bam"
    filtered_bam="$name$ext"
    samtools view -bS "$out_dir/mapped_reads/$filtered_sam" > "$out_dir/mapped_reads/$filtered_bam"

    #sort bam file
    >&2 echo "sort bam file"
    ext=".sorted.tmp.bam"
    sorted_bam="$name$ext"
    samtools sort -@ 4 "$out_dir/mapped_reads/$filtered_bam" > "$out_dir/mapped_reads/$sorted_bam"

    #remove duplicates
    >&2 echo "remove duplicates"
    ext=".rmdup.bam"
    rmdup_bam="$name$ext"
    samtools rmdup "$out_dir/mapped_reads/$sorted_bam" "$out_dir/mapped_reads/$rmdup_bam"

    #index bam file:
    >&2 echo "index bam file"
    samtools index "$out_dir/mapped_reads/$rmdup_bam"

    #call peaks with macs2:
    >&2 echo "call peaks"
    mkdir "$out_dir/peaks"
    macs2 callpeak -t "$out_dir/mapped_reads/$rmdup_bam" --outdir "$out_dir/peaks" -n "$name" --nomodel --shift -100 --extsize 200 -B --keep-dup all

    #remove peaks that fall within ENCODE blacklist regions:
    >&2 echo "removing blacklisted peaks"
    ext="_peaks.narrowPeak"
    narrow_peak="$name$ext"
    black="$name_blacklisted.bed"
    bedtools intersect -v -a $out_dir/peaks/$narrow_peak -b ~/bin/ATAC_Pipeline/ENCODE.hg19.blacklist.bed -wa > $out_dir/peaks/$black

    #get fasta sequences from peaks:
    >&2 echo "get fasta sequences from peaks"
    ext=".peak_seqs.fa"
    peak_seqs="$name$ext"
    
    bedtools getfasta -fi "$genome" -bed "$out_dir/peaks/$black" > "$out_dir/peaks/$peak_seqs"

    #call fimo to get motifs
    >&2 echo "call fimo to get motifs"
    mkdir "$out_dir/motifs"
    ext=".fimo_motifs.txt"
    fimo_motifs="$name$ext"
    fimo --text "$motifs" "$out_dir/peaks/$peak_seqs" > "$out_dir/$fimo_motifs"

    #split big fimo motif into individual bed files for each motif:
    >&2 echo "split fimo file by motif"
    cd "$out_dir/motifs"
    python ~/bin/ATAC_Pipeline/split_fimo_by_motif.py "$out_dir/$fimo_motifs"

    #using atactk, create discrete matrices that will go into centipede:
    >&2 echo "make discrete matrices"
    ext="_matrices"
    mat_dir_name="$name$ext"
    cd "$out_dir"
    mkdir -p "$out_dir/matrices/$mat_dir_name"

    #loop through each motif to make matrices
    for motif in $out_dir/motifs/*
    do
        ext=".dmatrix.gz"
        mName=$(basename $motif .motif.bed)
        >&2 echo "creating $mName matrix"
        outName="$mName$ext"
        make_cut_matrix -d -b '(1-100000 1)' -p 4 "$out_dir/mapped_reads/$rmdup_bam" $motif | gzip -c > "$out_dir/matrices/$mat_dir_name/$outName"
    done

    #run centipede on each motif, output posterior probabilites
    >&2 echo "run centipede"
    mkdir "$out_dir/centipede_postprobs"
    cd "$out_dir/centipede_postprobs"
    for motif in "$out_dir/motifs"/*
    do 
        a=$(basename $motif .motif.bed)
        ext=".dmatrix.gz"
        matrix="$a$ext"
        >&2 echo "motif: $a"
        >&2 echo "matrix: $matrix"

        #run centipede:
        /usr/bin/Rscript --vanilla ~/bin/ATAC_Pipeline/centipede.R "$out_dir/matrices/$mat_dir_name/$matrix" "$motif"
    done
    cd "$out_dir"

    #filter out bed regions that satisfy posterior probability threshold (run for 0.99 and 0.95)
    mkdir "$out_dir/footprints_95"
    mkdir "$out_dir/footprints_99"
    for motif in $out_dir/motifs/*
    do
        foot=$( basename $motif .motif.bed)
        ext1="_footprints.bed"
        fName="$foot$ext1"
        ext2=".dmatrix.gz.txt"
        pName="$foot$ext2"

        #run python script to get footprint regions that meet threshold posterior probabilities
        >&2 echo "extracting footprints for $foot"
        python ~/bin/ATAC_Pipeline/filter_by_postprobs.py $motif $out_dir/centipede_postprobs/$pName 0.99 > $out_dir/footprints_99/$fName
        python ~/bin/ATAC_Pipeline/filter_by_postprobs.py $motif $out_dir/centipede_postprobs/$pName 0.95 > $out_dir/footprints_95/$fName 
        >&2 echo "done getting footprint regions for $foot"

    done

fi