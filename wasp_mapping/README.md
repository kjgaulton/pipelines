`python3 wasp-mapping.py --help`

    usage: wasp-mapping.py [-h] --anaconda_path ANACONDA_PATH [--bam_start]
                           --false_discovery_rate FALSE_DISCOVERY_RATE
                           --hdf5_directory HDF5_DIRECTORY --lane_names LANE_NAMES
                           --max_processes MAX_PROCESSES --memory_limit
                           MEMORY_LIMIT --operating_directory OPERATING_DIRECTORY
                           [--paired_end] --reference_genome REFERENCE_GENOME
                           [--rna_seq] --rsid_directory RSID_DIRECTORY
                           [--STAR_genomeDir STAR_GENOMEDIR]
                           [--STAR_sjdbGTFfile STAR_SJDBGTFFILE]
                           [--text_based_snp_files] --variant_positions_file
                           VARIANT_POSITIONS_FILE --wasp_directory WASP_DIRECTORY

    Implementation of the re-mapping procedure to eliminate referencebias detailed
    in WASP.

    optional arguments:
      -h, --help            show this help message and exit
      --anaconda_path ANACONDA_PATH
                            Path to Anaconda distribution
      --bam_start           Input files are BAM files, skip initial alignment
      --false_discovery_rate FALSE_DISCOVERY_RATE
                            Set false discovery rate for binomial tests
      --hdf5_directory HDF5_DIRECTORY
                            Path to HDF5 snp files directory
      --lane_names LANE_NAMES
                            Path to text file listing lane_names to be considered
      --max_processes MAX_PROCESSES
                            Maximum number of processes allowed
      --memory_limit MEMORY_LIMIT
                            Approximate memory limit in gigabytes
      --operating_directory OPERATING_DIRECTORY
                            Path to operating directory
      --paired_end          Input data is paired-end
      --reference_genome REFERENCE_GENOME
                            Path to reference genome
      --rna_seq             Use splice-aware alignment for RNA-seq reads
      --rsid_directory RSID_DIRECTORY
                            Path to RSID directory
      --STAR_genomeDir STAR_GENOMEDIR
                            Path to STAR genome directory
      --STAR_sjdbGTFfile STAR_SJDBGTFFILE
                            Path to GTF file
      --text_based_snp_files
                            use text-based snp files instead
      --variant_positions_file VARIANT_POSITIONS_FILE
                            Path to variant positions file
      --wasp_directory WASP_DIRECTORY
                            Path to WASP directory
