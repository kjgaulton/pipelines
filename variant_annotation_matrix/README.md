`python3 vam.py --help`


    usage: vam.py [-h] -v VARIANTS -a ANNOTATIONS -o OUTPUT
                  [-p {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25}]
                  [-d DELIMITER] [--chromosome CHROMOSOME] [--position POSITION]
                  [--add-header ADD_HEADER | --replace-header REPLACE_HEADER] [-z]
                  [-f] [--params PARAMS]

    Construct a binary variant-annotation matrix from a table of input variants
    and a .bed file defining annotations.

    optional arguments:
      -h, --help            show this help message and exit

    I/O arguments:
      -v VARIANTS, --variants VARIANTS
                            path to variants file
      -a ANNOTATIONS, --annotations ANNOTATIONS
                            path to .bed file giving annotation data
      -o OUTPUT, --output OUTPUT
                            path to output file

    multiprocessing arguments:
      -p {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25}, --processes {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25}
                            number of processes to launch

    formatting arguments:
      -d DELIMITER, --delimiter DELIMITER
                            delimiter for input file, if not tab
      --chromosome CHROMOSOME
                            provide the (1-based) index of the column giving the
                            chromosome of each variant in the variants file, or a
                            string giving its header
      --position POSITION   provide the (1-based) index of the column giving the
                            position of each variant in the variants file, or a
                            string giving its header
      --add-header ADD_HEADER
                            Provide a string to be used as the header for a
                            variants file that does not have one
      --replace-header REPLACE_HEADER
                            Provide a string to be used as a replacement for the
                            variants file's header

    output compression arguments:
      -z, --compress-output
                            compress output file with gzip

    optional FGWAS-specific arguments:
      -f, --fine-mapping    produce output sorted by SEGNUMBER for fine mapping
      --params PARAMS       path to .params file contatining a subset of
                            annotations to use

    The ( --chromosome / --position ) argument is not required if the
    corresponding column is identified in the variants file header by a (case-
    insensitive) prefix of ("chromosome" / "position") at least three characters
    long. Examples: "chr" "POS" "ChRoM" "positio"
