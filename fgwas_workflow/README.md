`fgwas-workflow.py --help`

    usage: fgwas-workflow.py [-h] -i INPUT -o OUTPUT -t THRESHOLD -p PROCESSES
                             [-b BED] [-f] [-k WINDOW] [-w ANNOTATIONS] [-l LLK]
                             [--restrict-to-defined-ci-annotations]
                             [--positive-estimates-only]

    Apply the workflow suggested by the fgwas manual

    optional arguments:
      -h, --help            show this help message and exit

    required arguments:
      -i INPUT, --input INPUT
                            Path to fgwas input file
      -o OUTPUT, --output OUTPUT
                            Prefix for output files
      -t THRESHOLD, --threshold THRESHOLD
                            Likelihood threshold for model improvement
      -p PROCESSES, --processes PROCESSES
                            Maximum number of fgwas processes to launch

    optional arguments:
      -b BED, --bed BED     Path to .bed file containing regional definitions
      -f, --fine            Use for fine mapping input
      -k WINDOW, --window WINDOW
                            Window size
      -w ANNOTATIONS, --annotations ANNOTATIONS
                            "+"-separated list of annotations to use as a model
                            starting point.
      -l LLK, --llk LLK     Required if -w is used. Likelihood of model with the
                            annotation list given by -w.

    alternate workflow arguments:
      --restrict-to-defined-ci-annotations
                            Exclude from analysis any annotations whose individual
                            results have undefined confidence intervals.
      --positive-estimates-only
                            Add an annotation to the joint model only if it has a
                            positive parameter estimate.
