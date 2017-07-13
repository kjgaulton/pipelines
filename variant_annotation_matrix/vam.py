#!/usr/bin/env python3
#------------------------------------------------------------------------------#
#                                  vam.py                                      #
#------------------------------------------------------------------------------#

# Construct a matrix from a variants file and a .bed annotation file




#-------------------------------- Imports -------------------------------------#

import argparse
from operator import attrgetter
from multiprocessing import Pool
import gzip




#----------------------- Class Bootstrapping Function -------------------------#

def bootstrap(cls):
    cls.bootstrap()
    return(cls)




#----------------------------- Class Definitions ------------------------------#

@bootstrap
class Chromosome():
    '''A chromosome'''
    
    @classmethod
    def bootstrap(cls):
        for comparison_string in (
            '__eq__',
            '__ne__',
            '__lt__',
            '__le__',
            '__gt__',
            '__ge__'
        ):
            integer_comparison = getattr(int, comparison_string)
            def chromosome_comparison(
                self,
                chr,
                integer_comparison=integer_comparison
            ):
                return integer_comparison(self.int, chr.int)
            setattr(cls, comparison_string, chromosome_comparison)
    
    def __init__(self, chromosome_name):
        self.char = (
            str(chromosome_name)
            .casefold()
            .replace('chr', '')
            .replace('23', 'x')
            .replace('24', 'y')
            .replace('25', 'm')
        )
        try:
            self.int = int(self.char)
            if self.int not in range(1,26):
                self.raise_chromosome_name_error(chromosome_name)
        except ValueError:
            try:
                self.int = {'x': 23, 'y': 24, 'm': 25}[self.char]
            except KeyError:
                self.raise_chromosome_name_error(chromosome_name)
        self.variants = []
        self.annotations = {}
        self.variants_have_been_sorted = False
        self.annotations_have_been_sorted = False
    
    def __int__(self):
        return self.int
    
    def sort_variants(self):
        self.variants.sort(key=attrgetter('position'))
        self.variants_have_been_sorted = True
    
    def sort_annotations(self):
        for name, interval_list in self.annotations.items():
            interval_list.sort()
        self.annotations_have_been_sorted = True
    
    def annotate_variants(self, annotation_subset):
        if (
            self.variants_have_been_sorted
            and
            self.annotations_have_been_sorted
        ):
            for variant in self.variants:
                for annotation_name, interval_list in self.annotations.items():
                    empty_interval_count = 0
                    for interval in interval_list:
                        if interval[1] < variant.position:
                            empty_interval_count += 1
                        elif interval[0] > variant.position:
                            break
                        else:
                            if annotation_subset:
                                if annotation_name in annotation_subset:
                                    variant.annotations.add(annotation_name)
                            else:
                                variant.annotations.add(annotation_name)
                            break
                    for x in range(empty_interval_count):
                        interval_list.pop(0)
        else:
            raise Exception(
                'Variants and annotations must be sorted before variants can '
                'be annotated'
            )
                
    def raise_chromosome_name_error(self, chromosome_name):
        raise SystemExit(
            'Error: {} is not a valid chromosome name. Valid names include '
            'integers 1-25 or characters XxYyMm and may or not be prefixed '
            'with chr, CHR, Chr, etc.'.format(str(chromosome_name))
        )




class HeaderErrors():
    '''Errors for header parsing'''
    
    def raise_no_segnumber(self):
        raise SystemExit(
            'Error: The -f flag was used, but no SEGNUMBER column was '
            'found in the variants file header.'
        )
    
    def raise_column_index_out_of_range(self, column_string, index, ncol):
        raise SystemExit(
            (
                'Error: The index provided to --{} was {}, but '
                'there are only {} columns in the variants file.'
            ).format(
                column_string,
                index + 1,
                ncol
            )
        )
    
    def raise_ambiguity(self, column_string):
        raise SystemExit(
             (
                'Error: No {} column was specified and the variants file '
                'header is ambiguous.'
            ).format(
                column_string
            )
        )
        
    def raise_bad_column_header(self, column_string):
        raise SystemExit(
            (
                'Error: The string provided to --{} is not in the variants '
                'file header.'
            ).format(
                column_string
            )
        )
        
    def raise_possible_delimiter_error(self):
        raise SystemExit(
            (
                'Error: There doesn\'t seem to be enough information in the '
                'variants file. Did you forget to set the delimiter with -d ?'
            )
        ) 




class VariantsHeader():
    '''The header of an input variants file'''
    
    errors = HeaderErrors()
    
    def __init__(
        self,
        header_string,
        delimiter,
        chromosome_arg,
        position_arg,
        fine_mapping_arg
    ):
        self.tuple = tuple(
            header_string.replace(
                '\n', ''
            ).replace(
                '\\t', '\t'
            ).split(
                delimiter
            )
        )
        if len(self.tuple) < 2:
            self.errors.raise_possible_delimiter_error()
        if fine_mapping_arg:
            try:
                self.segnumber_index = self.tuple.index('SEGNUMBER')
            except ValueError:
                self.errors.raise_no_segnumber()
        for column_string, arg in (
            ('chromosome', chromosome_arg),
            ('position', position_arg)
        ):
            self.get_index(column_string, arg)
        self.sorting_attrs = tuple(
            attr
            for
            attr
            in
            ('segnumber', 'chromosome', 'position')
            if
            hasattr(self, '{}_index'.format(attr))
        )
    
    def get_index(self, column_string, arg):
        attr = '{}_index'.format(column_string)
        try:
            index = int(arg) - 1
            if index > len(self.tuple):
                self.errors.raise_column_index_out_of_range(
                    column_string,
                    index,
                    len(self.tuple)
                )
            else:
               setattr(self, attr, index)
        except ValueError:
            try:
                setattr(self, attr, self.tuple.index(arg) - 1)
            except ValueError:
                self.errors.raise_bad_column_header(column_string)
        except TypeError:
            for index in range(len(self.tuple)):
                prefix = self.tuple[index].casefold()
                if (
                    (
                        column_string[0:3] == prefix[0:3]
                    )
                    and
                    (
                        column_string[3:].startswith(prefix[3:])
                    )
                ):
                    try:
                        getattr(self, '{}_index'.format(column_string))
                        self.errors.raise_ambiguity(column_string)
                    except AttributeError:
                        setattr(self, attr, index)
            try:
                getattr(self, '{}_index'.format(column_string))
            except AttributeError:
                self.errors.raise_ambiguity(column_string)




class Variant():
    '''A variant'''
    
    raise_possible_delimiter_error = (
        HeaderErrors().raise_possible_delimiter_error
    )
    
    def __init__(self, tup, header):
        self.tuple = tup
        if len(self.tuple) < 2:
            self.raise_possible_delimiter_error()
        self.chromosome = Chromosome(
            self.tuple[getattr(header, 'chromosome_index')]
        ).int
        for attr in (
            attr for attr in header.sorting_attrs if attr != 'chromosome'
        ):
            try:
                setattr(
                    self,
                    attr,
                    int(
                        self.tuple[getattr(header, '{}_index'.format(attr))]
                    )
                )
            except ValueError:
                raise SystemExit(
                    (
                        'Error: The {} of a variant with the following tuple '
                        'representation could not be parsed as an integer:\n{}'
                    ).format(attr, repr(self.tuple))
                )
        self.annotations = set()




class VariantDistributor():
    '''
    A container for variants that can distribute them over the chromosomes.
    '''

    def __init__(self, variant_generator):
        self.generator = variant_generator
    
    def distribute(self, genome):
        for variant in self.generator:
            try:
                genome[variant.chromosome].variants.append(variant)
                genome[variant.chromosome].variants_have_been_sorted = False
            except KeyError:
                genome[variant.chromosome] = Chromosome(variant.chromosome)
                genome[variant.chromosome].variants.append(variant)
                genome[variant.chromosome].variants_have_been_sorted = False




class AnnotationDistributor():
    '''
    A container for annotations that can distribute them over the chromosomes.
    '''
    
    def __init__(self, annotation_generator):
       self.generator = annotation_generator
    
    def distribute(self, genome):
        for annotation in self.generator:
            try:
                (
                    genome[Chromosome(annotation[0]).int]
                    .annotations[annotation[3]]
                    .append((int(annotation[1]), int(annotation[2])))
                )
                (
                    genome[Chromosome(annotation[0]).int]
                    .annotations_have_been_sorted
                ) = (
                    False
                )
            except KeyError:
                try:
                    (
                        genome[Chromosome(annotation[0]).int]
                        .annotations[annotation[3]]
                    ) = (
                        [(int(annotation[1]), int(annotation[2]))]
                    )
                    (
                        genome[Chromosome(annotation[0]).int]
                        .annotations_have_been_sorted
                    ) = (
                        False
                    )
                except KeyError:
                    pass
            except (ValueError, IndexError):
                raise SystemExit(
                    (
                        'Error: The annotation file does not appear to be in '
                        'headerless BED format at a line with the following '
                        'tuple representation:\n{}'
                    ).format(repr(annotation))
                )




class MultiprocessingHandler():
    '''
    A device to perform multiprocessing tasks
    '''
    
    def __init__(self, genome):
        self.genome = genome
    
    def sort_variants(self):
        keys = []
        enumerated_variant_positions = []
        for key, chromosome in self.genome.items():
            keys.append(key)
            enumerated_variant_positions.append(
                tuple(
                    (elem, n)
                    for
                    (n, elem)
                    in
                    enumerate(
                        variant.position
                        for
                        variant
                        in
                        chromosome.variants
                    )
                )
            )
        with Pool(processes=args.processes) as pool:
            enumerated_variant_positions = pool.map(
                sorted,
                enumerated_variant_positions
            )
        for index in range(len(keys)):
            self.genome[keys[index]].variants = tuple(
                self.genome[keys[index]].variants[i]
                for
                i
                in
                (
                    n
                    for
                    (elem, n)
                    in
                    enumerated_variant_positions[index]
                )
            )
            self.genome[keys[index]].variants_have_been_sorted = True
    
    def sort_annotations(self):
        keys = []
        interval_lists = []
        for key, chromosome in self.genome.items():
            for name, interval_list in chromosome.annotations.items():
                keys.append((key, name))
                interval_lists.append(list(interval_list))
        with Pool(processes=args.processes) as pool:
            interval_lists = pool.map(
                sorted,
                interval_lists,
                chunksize = min(1, int(len(interval_lists) / args.processes))
            )
        for index in range(len(keys)):
            self.genome[keys[index][0]].annotations[keys[index][1]] = (
                list(interval_lists[index])
            )
            self.genome[keys[index][0]].annotations_have_been_sorted = True
    
    def annotate_variants(self, annotation_subset):
        keys = []
        argument_list = []
        for key, chromosome in self.genome.items():
            if (
                chromosome.variants_have_been_sorted
                and
                chromosome.annotations_have_been_sorted
            ):
                keys.append(key)
                argument_list.append(
                    (
                        tuple(
                            variant.position for variant in chromosome.variants
                        ),
                        dict(chromosome.annotations),
                        annotation_subset
                    )
                )
            else:
                raise Exception(
                    'Variants and annotations must be sorted before variants '
                    'can be annotated'
                )
                
        with Pool(processes=args.processes) as pool:
            annotations_per_variant = pool.starmap(
                annotate_variants,
                argument_list
            )
        for index in range(len(keys)):
            for i in range(len(self.genome[keys[index]].variants)):
                (
                    self
                    .genome[keys[index]]
                    .variants[i]
                    .annotations
                ) = annotations_per_variant[index][i]





#--------------------------- Function Definitions -----------------------------#

### Convenience functions for parsing arguments

# Construct annotation susbset
def get_annotation_subset(params_filename):
    try:
        with open(params_filename, 'r') as params_file:
            f.readline()
            return {
                line.split(' ')[0].replace('_ln', '')
                for
                line
                in
                params_file
                if
                float(line.split(' ')[2]) > 0
            }
    except TypeError:
        print('Using all annotations since no --params argument was given.')
        return None




# Handle a .gz extension on the output path
def output_handle_gz_extension(output_file_path):
    return ''.join(
        (
            output_file_path[0:-3],
            (output_file_path[-3:] != '.gz') * output_file_path[-3:]
        )
    )




### Loading data from the input files

# Load variants
def load_variants(args):
    with open(args.variants, 'r') as variants_file:
        if args.add_header:
            variants_header = VariantsHeader(
                args.add_header,
                args.delimiter,
                args.chromosome,
                args.position,
                args.fine_mapping
            )
        elif args.replace_header:
            variants_header = VariantsHeader(
                args.replace_header,
                args.delimiter,
                args.chromosome,
                args.position,
                args.fine_mapping
            )
            variants_file.readline()
        else:
            variants_header = VariantsHeader(
                variants_file.readline(),
                args.delimiter,
                args.chromosome,
                args.position,
                args.fine_mapping
            )
        variant_distributor = VariantDistributor(
            Variant(
                tup,
                variants_header
            )
            for
            tup
            in
            {
                tuple(
                    line.replace('\n', '').split(args.delimiter)
                )
                for
                line
                in
                variants_file
            }
        )
    return variants_header, variant_distributor




# Load annotations
def load_annotations(args):
    with open(args.annotations, 'r') as annotations_file:
        annotation_distributor = AnnotationDistributor(
            annotation
            for
            annotation
            in
            {
                tuple(
                    line.replace('\n', '').split()
                )
                for
                line
                in
                annotations_file
            }
        )
    return annotation_distributor




### Multiprocessing function

# Annotate variants
def annotate_variants(variant_positions, annotations, annotation_subset):
    annotations_per_variant = []
    for variant_position in variant_positions:
        annotations_per_variant.append(set())
        for annotation_name, interval_list in annotations.items():
            empty_interval_count = 0
            for interval in interval_list:
                if interval[1] < variant_position:
                    empty_interval_count += 1
                elif interval[0] > variant_position:
                    break
                else:
                    if annotation_subset:
                        if annotation_name in annotation_subset:
                            annotations_per_variant[-1].add(annotation_name)
                    else:
                        annotations_per_variant[-1].add(annotation_name)
                    break
            for x in range(empty_interval_count):
                interval_list.pop(0)
    return annotations_per_variant




### Generating and writing the matrix

# Generate lines of the matrix
def generate_matrix(variants_header, genome):
    chromosomes = tuple(
        sorted(
            chromosome
            for
            key, chromosome
            in
            genome.items()
        )
    )
    annotations = tuple(
        sorted(
            {
                annotation_name
                for
                chromosome
                in
                chromosomes
                for
                annotation_name
                in
                chromosome.annotations.keys()
            }
        )
    )
    yield variants_header.tuple + annotations
    for chromosome in chromosomes:
        for variant in chromosome.variants:
            yield (
                variant.tuple
                +
                tuple(
                    str(
                        int(
                            annotation
                            in
                            variant.annotations
                        )
                    )
                    for
                    annotation
                    in
                    annotations
                )
            )




# Write to output
def write_output(args, variants_header, matrix_generator):
    if args.fine_mapping:
        header = next(matrix_generator)
        matrix_generator = (
            row
            for
            row
            in
            (
                (
                    header,
                )
                + 
                tuple(
                    sorted(
                        matrix_generator,
                        key = lambda row: int(
                            row[variants_header.segnumber_index]
                        )
                    )
                )
            )
        )
    output_data = (
        '\n'.join(
            (
                args.delimiter.join(
                    row
                )
                for
                row
                in
                matrix_generator
            )
        )
        +
        '\n'
    )
    if args.compress_output:
        with gzip.open('{}.gz'.format(args.output), 'wb') as f:
           f.write(bytes(output_data, 'utf8'))
    else:
        with open(args.output, 'w') as f:
            f.write(output_data)

    


# Main
def main(args):
    annotation_subset = get_annotation_subset(args.params)
    genome = {}
    print('Loading variants')
    variants_header, variant_distributor = load_variants(args)
    print('Distributing variants')
    variant_distributor.distribute(genome)
    print('Loading annotations')
    annotation_distributor = load_annotations(args)
    print('Distributing annotations')
    annotation_distributor.distribute(genome)
    print('Annotating variants')
    if args.processes > 1:
        multiprocessing_handler = MultiprocessingHandler(genome)
        multiprocessing_handler.sort_variants()
        multiprocessing_handler.sort_annotations()
        multiprocessing_handler.annotate_variants(annotation_subset)
    else:
        for key, chromosome in genome.items():
            chromosome.sort_variants()
            chromosome.sort_annotations()
            chromosome.annotate_variants(annotation_subset)
    print('Preparing output matrix')
    matrix_generator = generate_matrix(variants_header, genome)
    print('Writing to file')
    write_output(args, variants_header, matrix_generator)
    print('Done')




# Parse arguments
def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            'Construct a binary variant-annotation matrix from a table of '
            'input variants and a .bed file defining annotations.'
        ),
        epilog=(
            'The  ( --chromosome / --position ) argument is not required if '
            'the corresponding column is identified in the variants file '
            'header by a (case-insensitive) prefix of ("chromosome" / '
            '"position") at least three characters long. Examples: "chr" "POS" '
            '"ChRoM" "positio"'
        )
    )
    
    io_group = parser.add_argument_group('I/O arguments')
    io_group.add_argument(
        '-v',
        '--variants',
        required=True,
        help='path to variants file'
    )
    io_group.add_argument(
        '-a',
        '--annotations',
        required=True,
        help='path to .bed file giving annotation data'
    )
    io_group.add_argument(
        '-o',
        '--output',
        required=True,
        help='path to output file'
    )
    
    multiprocessing_group = parser.add_argument_group(
        'multiprocessing arguments'
    )
    multiprocessing_group.add_argument(
        '-p',
        '--processes',
        type=int,
        default=1,
        choices=range(1, 26),
        help='number of processes to launch'
    )
    
    formatting_group = parser.add_argument_group('formatting arguments')
    formatting_group.add_argument(
        '-d',
        '--delimiter',
        default='\t',
        help='delimiter for input file, if not tab'
    )
    formatting_group.add_argument(
        '--chromosome',
        help=(
            'provide the (1-based) index of the column giving the chromosome '
            'of each variant in the variants file, or a string giving its '
            'header'
        )
    )
    formatting_group.add_argument(
        '--position',
        help=(
            'provide the (1-based) index of the column giving the position '
            'of each variant in the variants file, or a string giving its '
            'header'
        )
    )
    provide_header_group = formatting_group.add_mutually_exclusive_group()
    provide_header_group.add_argument(
        '--add-header',
        help=(
            'Provide a string to be used as the header for a variants file '
            'that does not have one'
        )
    )
    provide_header_group.add_argument(
        '--replace-header',
        help=(
            'Provide a string to be used as a replacement for the variants '
            'file\'s header'
        )
    )
    
    output_compression_group = parser.add_argument_group(
        'output compression arguments'
    )
    output_compression_group.add_argument(
        '-z',
        '--compress-output',
        action='store_true',
        help='compress output file with gzip'
    )
    
    fgwas_group = parser.add_argument_group(
        'optional FGWAS-specific arguments'
    )
    fgwas_group.add_argument(
        '-f',
        '--fine-mapping',
        action='store_true',
        help='produce output sorted by SEGNUMBER for fine mapping'
    )
    fgwas_group.add_argument(
        '--params',
        help='path to .params file contatining a subset of annotations to use'
    )
    
    args = parser.parse_args()
    if args.compress_output:
        args.output = output_handle_gz_extension(args.output)
    return args




#-------------------------------- Execute -------------------------------------#

if __name__ == '__main__':
    args = parse_arguments()
    main(args)
