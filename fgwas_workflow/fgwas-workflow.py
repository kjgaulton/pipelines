#!/usr/bin/env python3
#------------------------------------------------------------------------------#
#                             fgwas-workflow.py                                #
#------------------------------------------------------------------------------#

# Apply the workflow suggested by the fgwas manual




#-------------------------------- Imports -------------------------------------#

import argparse
import gzip
import os
import subprocess
import operator
import itertools
from multiprocessing import Pool




#---------------------------- Class definitions -------------------------------#

class InputFileHeader():
    '''
    The header of the input file
    '''

    non_annotation_column_names = {
            'SNPID',
            'CHR',
            'POS',
            'Z',
            'F',
            'N',
            'NCASE',
            'NCONTROL',
            'LNBF',
            'SEGNUMBER'
        }
    
    def __init__(self, args):
        self.args = args
        with gzip.open(args.input, 'rt') as f:
            self.tuple = tuple(f.readline().replace('\n', '').split(' '))
        self.annotations = {
            column_name
            for
            column_name
            in
            self.tuple
            if
            column_name
            not
            in
            self.non_annotation_column_names
        }
    
    def validate(self):
        '''
        Validate the input file header
        '''
        for column_name in self.non_annotation_column_names:
            if self.tuple.count(column_name) > 1:
                raise SystemExit(
                    (
                        'ERROR: Invalid header on input file: Either {0} is '
                        'duplicated, or there is an annotation called {0}.'
                    ).format(column_name)
                )
            elif (
                (
                    column_name
                    not
                    in
                    {'N', 'NCASE', 'NCONTROL', 'LNBF', 'SEGNUMBER'}
                )
                and
                (
                    column_name
                    not
                    in
                    self.tuple
                )
            ):
                raise SystemExit(
                    (
                        'ERROR: Invalid header on input file: {} is missing.'
                    ).format(column_name)
                )
        if {'N', 'NCASE', 'NCONTROL'} <= set(self.tuple):
            raise SystemExit(
                'ERROR: N, NCASE, and NCONTROL were all found in the input '
                'file header. Please use EITHER quantitative OR case-control '
                'format.'
            )




class IndividualAnnotationResults():
    '''
    Individual annotation results
    '''
    
    def __init__(self, args, header):
        self.args = args
        self.annotations = header.annotations
        self.data = None
        self.defined_ci_annotations = set()
        self.header = header
    
    def collect(self):
        '''
        Collect individual annotation results
        '''
        with Pool(processes=min(args.processes, len(self.annotations))) as pool:
            self.data = (
                (
                    ('parameter', 'ln(lk)', 'CI_lo', 'estimate', 'CI_hi'),
                )
                +
                tuple(
                    sorted(
                        pool.starmap(
                            collect_individual_annotation_result,
                            (
                                (args, annotation, self.header)
                                for
                                annotation
                                in
                                self.annotations
                            )
                        ),
                        key=operator.itemgetter(1),
                        reverse=True
                    )
                )
            )
    
    def export(self):
        '''
        Export individual annotation results
        '''
        if self.data:
             with open(
                '{}-individual-annotation-results.txt'.format(args.output),
                'w'
            ) as f:
                f.write(
                    '\n'.join(
                        '\t'.join((str(entry) for entry in line))
                        for
                        line
                        in
                        self.data
                    )
                    +
                    '\n'
                )
        else:
            raise Exception(
                'Cannot export individual annotation results with no data '
                'present.'
            )
    
    def identify_defined_ci_annotations(self):
        if self.data:
            for parameter, llk, ci_lo, estimate, ci_hi in self.data[1:]:
                try:
                    float(ci_lo)
                    float(ci_hi)
                    self.defined_ci_annotations.add(
                        parameter.replace('_ln', '')
                    )
                except ValueError:
                    pass
        if not self.defined_ci_annotations:
            raise SystemExit(
                'None of the provided annotations were viable (all '
                'individual results had undefined confidence intervals). '
                'Terminating analysis.'
            )
    
    def identify_best_starting_state(self):
        if self.data:
            best_individual_llk = max(
                float(llk)
                for
                parameter, llk, ci_lo, estimate, ci_hi
                in
                self.data[1:]
            )
            best_individual_annotations = []
            for parameter, llk, ci_lo, estimate, ci_hi in self.data[1:]:
                if (
                    float(llk) == best_individual_llk
                    and
                    (
                        (
                            parameter.replace('_ln', '')
                            in
                            self.defined_ci_annotations
                        )
                        if
                        self.defined_ci_annotations
                        else
                        True
                    )
                ):
                    best_individual_annotations.append(
                        (parameter, llk, ci_lo, estimate, ci_hi)
                    )
            if len(best_individual_annotations) == 1:
                best_individual_annotation = best_individual_annotations[0]
                self.best_starting_state = {
                    'annotations': [best_individual_annotation[0]],
                    'llk': best_individual_annotation[1],
                    'estimates': {
                        best_individual_annotation[0]:
                            best_individual_annotation[2:]
                    },
                    'xvl': float('-inf'),
                    'xv_penalty': 0,
                    'output_files': {}
                }
            elif len(best_individual_annotations) > 1:
                annotation_combinations = tuple(
                    itertools.chain.from_iterable(
                        itertools.combinations(
                            (
                                parameter
                                for
                                parameter, llk, ci_lo, estimate, ci_hi
                                in
                                best_individual_annotations
                            ),
                            length
                        )
                        for
                        length
                        in
                        range(2, len(best_individual_annotations) + 1)
                    )
                )
                with Pool(
                    processes=min(
                        self.args.processes,
                        len(annotation_combinations)
                    )
                ) as pool:
                    pool.map(
                        call_fgwas,
                        (
                            {
                                'args': self.args,
                                'annotation':
                                    '+'.join(
                                        self.annotations
                                        +
                                        list(annotation_combination)
                                    ),
                                'tag': 
                                    '+'.join(
                                        annotation_combination
                                    ),
                                'header': individual_results.header
                            }
                            for
                            annotation_combination
                            in
                            annotation_combinations
                        )
                    )
                llk_list = []
                for annotation_combination in annotation_combinations:
                    with open(
                        '{}-{}.llk'
                        .format(
                            args.output,
                            '+'.join(annotation_combination)
                        )
                    ) as f:
                        llk_list.append(
                            (
                                float(
                                    f
                                    .readline()
                                    .replace('ln(lk): ', '')
                                    .replace('\n', '')
                                ),
                                annotation_combination
                            )
                        )
                best_combo_llk = max(
                    llk
                    for
                    llk, annotation_combination
                    in
                    llk_list
                )
                best_combos = tuple(
                    annotation_combination
                    for
                    llk, annotation_combination
                    in
                    llk_list
                    if
                    llk == best_combo_llk
                )
                if (
                    (best_combo_llk > best_individual_llk)
                    and
                    len(best_combos) == 1
                ):
                    best_combo = list(best_combos[0])
                    estimates = {}
                    with open(
                        '{}-{}.params'
                        .format(args.output, '+'.join(best_combo))
                    ) as f:
                        for line in f:
                            parsed_line = tuple(
                                line.replace('\n', '').split(' ')
                            )
                            parameter = parsed_line[0][:-3]
                            for annotation in self.annotations:
                                if annotation == parameter:
                                    estimates[annotation] = tuple(
                                        parsed_line[1:]
                                    )
                                    break
                    self.best_starting_state = {
                        'annotations': best_combo,
                        'llk': best_combo_llk,
                        'estimates': estimates,
                        'xvl': float('-inf'),
                        'xv_penalty': 0,
                        'output_files': {}
                    }
                else:
                    print(
                        'Next best annotation is ambiguous, taking a null '
                        'step and proceeding to cross-validation phase'
                    )
                for annotation in set(
                    '+'.join(combo)
                    for
                    combo
                    in
                    annotation_combinations
                ):
                    clean_up_intermediate_files(
                        args,
                        '{}-{}.',
                        annotation
                    )
        else:
            raise Exception(
                'Cannot identify annotations with with well-defined confidence '
                'intervals with no data present.'
            )
        




class FgwasModel():
    '''
    An FGWAS model
    '''
    
    def __init__(self, args):
        self.args = args
        if self.args.initialize:
            self.annotations = self.args.annotations.split('+')
            self.llk = args.llk
        else:
            self.annotations = []
            self.llk = float('-inf')
        self.estimates = {}
        self.xvl = float('-inf')
        self.cross_validation_penalty = 0
        self.output_files = {}
        self.cache()
    
    def cache(self):
        '''
        Cache the model state
        '''
        self.annotations_cache = list(self.annotations)
        self.llk_cache = float(self.llk)
        self.estimates_cache = dict(self.estimates)
        self.xvl_cache = float(self.xvl)
        self.cross_validation_penalty_cache = float(
            self.cross_validation_penalty
        )
        self.output_files_cache = dict(self.output_files)
    
    def revert(self):
        '''
        Revert to the cached model state
        '''
        self.annotations = self.annotations_cache
        self.llk = self.llk_cache
        self.estimates = self.estimates_cache
        self.xvl = self.xvl_cache
        self.cross_validation_penalty = self.cross_validation_penalty_cache
        self.output_files = self.output_files_cache
        print('Reverted last step')
    
    def clear(self):
        '''
        Clear the model to an empty state
        '''
        self.cache()
        self.annotations = []
        self.llk = float('-inf')
        self.estimates = {}
        self.xvl = float('-inf')
        self.cross_validation_penalty = 0
        self.output_files = {}
    
    def set_state(self, state):
        '''
        Set the model state according to a dictionary definition
        '''
        self.cache()
        self.annotations = state['annotations']
        self.llk = state['llk']
        self.estimates = state['estimates']
        self.xvl = state['xvl']
        self.cross_validation_penalty = state['xv_penalty']
        self.output_files = state['output_files']
    
    def collect_output_files(self, annotation):
        '''
        Collect the output files tagged with a specific annotation
        '''
        for extension in ('llk', 'params', 'ridgeparams'):
            with open(
                '{}-{}.{}'.format(self.args.output, annotation, extension),
                'r'
            ) as f:
                self.output_files[extension] = f.read()
        for extension in ('bfs', 'segbfs'):
            with gzip.open(
                '{}-{}.{}.gz'.format(self.args.output, annotation, extension),
                'rb'
            ) as f:
                self.output_files[extension] = f.read()
    
    def export(self, tag):
        '''
        Export the output file state
        '''
        for extension in ('llk', 'params', 'ridgeparams'):
            with open(
                '{}-{}.{}'.format(self.args.output, tag, extension),
                'w'
            ) as f:
                f.write(self.output_files[extension])
        for extension in ('bfs', 'segbfs'):
            with gzip.open(
                '{}-{}.{}.gz'.format(self.args.output, tag, extension),
                'wb'
            ) as f:
                f.write(self.output_files[extension])
    
    def append(self, annotation, header):
        '''
        Append an annotation to the model and determine the resulting joint
        model
        '''
        self.cache()
        self.annotations.append(annotation)
        call_fgwas(
            {
                'args': self.args,
                'annotation': '+'.join(self.annotations),
                'header': header
            }
        )
        with open(
            '{}-{}.llk'.format(args.output, annotation.split('+')[-1])
        ) as f:
            self.llk = float(
                f.readline().replace('ln(lk): ', '').replace('\n', '')
            )
        
        with open(
            '{}-{}.params'.format(args.output, annotation.split('+')[-1])
        ) as f:
            for line in f:
                parsed_line = tuple(line.replace('\n', '').split(' '))
                parameter = parsed_line[0][:-3]
                for annotation in self.annotations:
                    if annotation == parameter:
                        self.estimates[annotation] = tuple(parsed_line[1:])
                        break
        clean_up_intermediate_files(args, '{}-{}.', annotation)
        
    def append_best_annotation(self, individual_results):
        '''
        Identify the next best annotation, append it to the model, and
        determine the resulting joint model
        '''
        self.cache()
        remaining_annotations = {
            annotation
            for
            annotation
            in
            (
                individual_results.defined_ci_annotations
                if
                individual_results.defined_ci_annotations
                else
                individual_results.annotations
            )
            if
            annotation
            not
            in
            self.annotations
        }
        with Pool(
            processes=min(
                self.args.processes,
                len(remaining_annotations)
            )
        ) as pool:
            pool.map(
                call_fgwas,
                (
                    {
                        'args': self.args,
                        'annotation': '+'.join(self.annotations + [annotation]),
                        'header': individual_results.header
                    }
                    for
                    annotation
                    in
                    remaining_annotations
                )
            )
        llk_list = []
        for annotation in remaining_annotations:
            positive_estimate = False
            with open(
                '{}-{}.params'.format(args.output, annotation)
            ) as f:
                for line in f:
                    parameter, ci_lo, estimate, ci_hi = (
                        line
                        .replace('\n', '')
                        .split(' ')
                    )
                    if (
                        (parameter.replace('_ln', '') == annotation)
                        and 
                        (float(estimate) > 0)
                    ):
                        positive_estimate = True
            if (
                positive_estimate
                if
                individual_results.args.positive_estimates_only
                else
                True
            ):
                with open(
                    '{}-{}.llk'.format(args.output, annotation)
                ) as f:
                    llk_list.append(
                        (
                            float(
                                f
                                .readline()
                                .replace('ln(lk): ', '')
                                .replace('\n', '')
                            ),
                            annotation
                        )
                    )
        if individual_results.args.positive_estimates_only:
            print(
                '{} annotations can be added with positive parameter estimates'
                .format(len(llk_list))
            )
        if len(llk_list) == 0:
            print('No more annotations can be added')
        else:
            best_llk = max(llk for llk, annotation in llk_list)
            best_annotations = tuple(
                annotation
                for
                llk, annotation
                in
                llk_list
                if
                llk == best_llk
            )
            if len(best_annotations) == 1:
                best_annotation = best_annotations[0]
                self.llk = best_llk
                self.annotations.append(best_annotation)
                with open(
                    '{}-{}.params'.format(args.output, best_annotation)
                ) as f:
                    for line in f:
                        parsed_line = tuple(line.replace('\n', '').split(' '))
                        parameter = parsed_line[0][:-3]
                        for annotation in self.annotations:
                            if annotation == parameter:
                                self.estimates[annotation] = tuple(
                                    parsed_line[1:]
                                )
                                break
                self.collect_output_files(best_annotation)
                for annotation in remaining_annotations:
                    clean_up_intermediate_files(args, '{}-{}.', annotation)
                print(
                    'Added {} to joint model (llk: {})'
                    .format(best_annotation, str(best_llk))
                )
            elif len(best_annotations) > 1:
                annotation_combinations = tuple(
                    itertools.chain.from_iterable(
                        itertools.combinations(best_annotations, length)
                        for
                        length
                        in
                        range(2, len(best_annotations) + 1)
                    )
                )
                with Pool(
                    processes=min(
                        self.args.processes,
                        len(annotation_combinations)
                    )
                ) as pool:
                    pool.map(
                        call_fgwas,
                        (
                            {
                                'args': self.args,
                                'annotation':
                                    '+'.join(
                                        self.annotations
                                        +
                                        list(annotation_combination)
                                    ),
                                'tag': 
                                    '+'.join(
                                        annotation_combination
                                    ),
                                'header': individual_results.header
                            }
                            for
                            annotation_combination
                            in
                            annotation_combinations
                        )
                    )
                llk_list = []
                for annotation_combination in annotation_combinations:
                    with open(
                        '{}-{}.llk'
                        .format(
                            args.output,
                            '+'.join(annotation_combination)
                        )
                    ) as f:
                        llk_list.append(
                            (
                                float(
                                    f
                                    .readline()
                                    .replace('ln(lk): ', '')
                                    .replace('\n', '')
                                ),
                                annotation_combination
                            )
                        )
                best_combo_llk = max(
                    llk
                    for
                    llk, annotation_combination
                    in
                    llk_list
                )
                best_combos = tuple(
                    annotation_combination
                    for
                    llk, annotation_combination
                    in
                    llk_list
                    if
                    llk == best_combo_llk
                )
                if (best_combo_llk > best_llk) and len(best_combos) == 1:
                    best_combo = list(best_combos[0])
                    self.llk = best_combo_llk
                    self.annotations.extend(best_combo)
                    with open(
                        '{}-{}.params'
                        .format(args.output, '+'.join(best_combo))
                    ) as f:
                        for line in f:
                            parsed_line = tuple(
                                line.replace('\n', '').split(' ')
                            )
                            parameter = parsed_line[0][:-3]
                            for annotation in self.annotations:
                                if annotation == parameter:
                                    self.estimates[annotation] = tuple(
                                        parsed_line[1:]
                                    )
                                    break
                    self.collect_output_files('+'.join(best_combo))
                    print(
                        'Added {} to joint model (llk: {})'
                        .format('+'.join(best_combo), str(best_llk))
                    )
                else:
                    print(
                        'Next best annotation is ambiguous, taking a null '
                        'step and proceeding to cross-validation phase'
                    )
                for annotation in remaining_annotations.union(
                    set(
                        '+'.join(combo)
                        for
                        combo
                        in
                        annotation_combinations
                    )
                ):
                    clean_up_intermediate_files(
                        args,
                        '{}-{}.',
                        annotation
                    )

    def calibrate_cross_validation_penalty(self, header):
        '''
        Calibrate the cross validation penalty
        '''
        xv_penalties = tuple(
            0.5 * x / max(args.processes, 10)
            for
            x
            in
            range(1, max(args.processes, 10) + 1)
        )
        with Pool(processes=args.processes) as pool:
            pool.map(
                call_fgwas,
                tuple(
                    {
                        'args': self.args,
                        'annotation': '+'.join(self.annotations),
                        'header': header,
                        'xv_penalty': xv_penalty
                    }
                    for
                    xv_penalty
                    in
                    xv_penalties
                )
            )
        xv_likelihoods = []
        for xv_penalty in xv_penalties:
            with open(
                '{}-p{}.ridgeparams'.format(args.output, str(xv_penalty)),
                'r'
            ) as f:
                xv_likelihoods.append(
                    (
                        float(
                            f
                            .read()
                            .split('\n')[-2]
                            .replace('X-validation penalize ln(lk): ', '')
                            .replace('\n', '')
                        ),
                        xv_penalty
                    )
                )
            clean_up_intermediate_files(args, '{}-p{}.', str(xv_penalty))
        best_xvl_penalty = (
            sorted(xv_likelihoods, key=operator.itemgetter(0), reverse=True)
            [0]
        )
        self.xvl = best_xvl_penalty[0]
        self.cross_validation_penalty = best_xvl_penalty[1]
        
    def remove_worst_annotation(self, header):
        '''
        Identify the annotation contributing least to the model and remove it
        '''
        if len(self.annotations) < 2:
            raise Exception(
                'Can\'t remove annotations from a model with only one '
                'annotation.'
            )
        self.cache()
        with Pool(
            processes=min(
                self.args.processes,
                len(self.annotations)
            )
        ) as pool:
            pool.map(
                call_fgwas,
                (
                    {
                        'args': self.args,
                        'annotation': annotation_0,
                        'header': header,
                        'xv_penalty': self.cross_validation_penalty,
                        'drop': '+'.join(
                            annotation_1
                            for
                            annotation_1
                            in
                            self.annotations
                            if
                            (annotation_0 != annotation_1)
                        )
                    }
                    for
                    annotation_0
                    in
                    self.annotations
                )
            )
        xvl_list = []
        for annotation in self.annotations:
            with open(
                '{}-{}.ridgeparams'.format(args.output, annotation),
                'r'
            ) as f:
                xvl_list.append(
                    (
                        float(
                            f
                            .read()
                            .split('\n')[-2]
                            .replace('X-validation penalize ln(lk): ', '')
                            .replace('\n', '')
                        ),
                        annotation
                    )
                )
        best_xvl = max(xvl for xvl, annotation in xvl_list)
        worst_annotations = tuple(
            annotation
            for
            xvl, annotation
            in
            xvl_list
            if
            xvl == best_xvl
        )
        if len(worst_annotations) == 1:
            worst_annotation = worst_annotations[0]
            self.xvl = best_xvl
            self.annotations.remove(worst_annotation)
            with open(
                '{}-{}.params'.format(args.output, worst_annotation)
            ) as f:
                for line in f:
                    parsed_line = tuple(line.replace('\n', '').split(' '))
                    parameter = parsed_line[0][:-3]
                    for annotation in self.annotations:
                        if annotation == parameter:
                            self.estimates[annotation] = tuple(
                                parsed_line[1:]
                            )
                            break
            self.collect_output_files(worst_annotation)
            for annotation in self.annotations_cache:
                clean_up_intermediate_files(args, '{}-{}.', annotation)
            print(
                'Dropped {} from joint model (xvl: {})'
                .format(worst_annotation, str(best_xvl))
            )
        elif len(worst_annotations) > 1:
            annotation_combinations = tuple(
                itertools.chain.from_iterable(
                    itertools.combinations(worst_annotations, length)
                    for
                    length
                    in
                    range(2, len(worst_annotations) + 1)
                )
            )
            with Pool(
                processes=min(
                    self.args.processes,
                    len(annotation_combinations)
                )
            ) as pool:
                pool.map(
                    call_fgwas,
                    (
                        {
                            'args': self.args,
                            'annotation':
                                '+'.join(
                                    self.annotations
                                    +
                                    list(annotation_combination)
                                ),
                            'tag': 
                                '+'.join(
                                    sorted(annotation_combination)
                                ),
                            'header': header,
                            'xv_penalty': self.cross_validation_penalty,
                            'drop': '+'.join(
                                annotation
                                for
                                annotation
                                in
                                self.annotations
                                if
                                (annotation not in annotation_combination)
                            )
                        }
                        for
                        annotation_combination
                        in
                        annotation_combinations
                    )
                )
            xvl_list = []
            for annotation_combination in annotation_combinations:
                with open(
                    '{}-{}.ridgeparams'
                    .format(
                        args.output,
                        '+'.join(annotation_combination)
                    )
                ) as f:
                    xvl_list.append(
                        (
                            float(
                                f
                                .read()
                                .split('\n')[-2]
                                .replace(
                                    'X-validation penalize ln(lk): ',
                                    ''
                                )
                                .replace('\n', '')
                            ),
                            annotation_combination
                        )
                    )
            best_combo_xvl = max(
                xvl
                for
                xvl, annotation_combination
                in
                xvl_list
            )
            best_combos = [
                annotation_combination
                for
                xvl, annotation_combination
                in
                xvl_list
                if
                xvl == best_combo_xvl
            ]
            if best_combo_xvl >= best_xvl:
                best_combo = sorted(
                    list(
                        {
                            annotation
                            for
                            combo
                            in
                            best_combos
                            for
                            annotation
                            in
                            combo
                        }
                    )
                )
                self.xvl = best_combo_xvl
                for annotation in best_combo:
                    self.annotations.remove(annotation)
                with open(
                    '{}-{}.params'
                    .format(args.output, '+'.join(best_combo))
                ) as f:
                    for line in f:
                        parsed_line = tuple(
                            line.replace('\n', '').split(' ')
                        )
                        parameter = parsed_line[0][:-3]
                        for annotation in self.annotations:
                            if annotation == parameter:
                                self.estimates[annotation] = tuple(
                                    parsed_line[1:]
                                )
                                break
                self.collect_output_files('+'.join(best_combo))
                print(
                    'Dropped {} from joint model (xvl: {})'
                    .format('+'.join(best_combo), str(best_xvl))
                )
            else:
                print(
                    'Next annotation to drop is ambiguous, taking a null '
                    'step and terminating analysis'
                )
            for annotation in set(self.annotations_cache).union(
                set(
                    '+'.join(combo)
                    for
                    combo
                    in
                    annotation_combinations
                )
            ):
                clean_up_intermediate_files(
                    args,
                    '{}-{}.',
                    annotation
                )
        
        
        

#--------------------------- Function definitions -----------------------------#

# Call fgwas
def call_fgwas(args_dict):
    '''
    A function for calling FGWAS. Its argument is a dictionary, which should
    at least include keys "args", "annotation" and "header", and may also
    include keys "tag", "xv_penalty" and/or "drop." "xv_penalty" is required if
    "drop" is included.
    '''
    args = args_dict['args']
    annotation = args_dict['annotation']
    header = args_dict['header']
    tag = (
        args_dict['tag']
        if
        'tag'
        in
        args_dict.keys()
        else
        annotation.split('+')[-1]
    )
    xv_penalty = (
        args_dict['xv_penalty']
        if
        'xv_penalty'
        in
        args_dict.keys()
        else
        None
    )
    drop = (
        args_dict['drop']
        if
        'drop'
        in
        args_dict.keys()
        else
        ''
    )
    fgwas_command_line = (
        (
            'fgwas',
            '-i', args.input,
            '-w', annotation * (not bool(drop)) + drop * bool(drop),
            '-o', (
                (
                    '{}-{}'.format(args.output, tag)
                    *
                    (
                        (not bool(xv_penalty))
                        or
                        bool(drop)
                    )
                )
                +
                (
                    '{}-p{}'.format(args.output, str(xv_penalty))
                    *
                    (
                        bool(xv_penalty)
                        and
                        (not bool(drop))
                    )
                )
            )
        )
        +
        (('-xv', '-p', str(xv_penalty)) * bool(xv_penalty))
        +
        (('-bed', args.bed) * bool(args.bed))
        +
        (('-cc',) * ({'NCASE', 'NCONTROL'} <= set(header.tuple)))
        +
        (('-fine',) * ('SEGNUMBER' in header.tuple))
        +
        (('-k', args.window) * bool(args.window))
        +
        ('-print',)
    )
    subprocess.call(
        fgwas_command_line,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )




# Clean up intermediate files
def clean_up_intermediate_files(args, format_string, format_arg):
    '''
    Remove the intermediate output files created by the fgwas run
    '''
    for intermediate_file in tuple(
        (format_string + extension).format(args.output, format_arg)
        for
        extension
        in
        ('llk', 'params', 'ridgeparams', 'segbfs.gz', 'bfs.gz')
    ):
        os.remove(intermediate_file)




# Collect an individual annotation result
def collect_individual_annotation_result(args, annotation, header):
    '''
    Collect the result for an individual annotation
    '''
    model = FgwasModel(args)
    model.clear()
    model.append(annotation, header)
    return (annotation, model.llk) + model.estimates[annotation]




# Main
def main(args):
    print('Reading input file header')
    header = InputFileHeader(args)
    print('Validating input file header')
    header.validate()
    model = FgwasModel(args)
    individual_results = IndividualAnnotationResults(args, header)
    print('Collecting individual annotation results')
    individual_results.collect()
    print('Exporting individual annotation results')
    individual_results.export()
    if args.restrict_to_defined_ci_annotations:
        print(
            'Identifying annotations with well-defined confidence intervals'
        )
        individual_results.identify_defined_ci_annotations()
        print(
            '{} annotations with well-defined confidence intervals.'
            .format(len(individual_results.defined_ci_annotations))
        )
    if not args.initialize:
        individual_results.identify_best_starting_state()
        model.set_state(individual_results.best_starting_state)
    print(
        'Constructing joint model, beginning with: {} (llk: {})'
        .format(', '.join(model.annotations), model.llk)
    )
    for iteration in range(
        len(
            individual_results.defined_ci_annotations
            if
            individual_results.defined_ci_annotations
            else
            individual_results.annotations
        )
        -
        len(
            model.annotations
        )
    ):
        model.append_best_annotation(individual_results)
        if model.llk < (model.llk_cache):
            model.revert()
            break
    print('Exporting pre-cross-validation results')
    model.export('pre-xv')
    print('Calibrating cross-validation penalty')
    model.calibrate_cross_validation_penalty(header)
    print('Beginning cross-validation phase')
    number_of_annotations = len(model.annotations)
    for iteration in range(number_of_annotations - 1):
        model.remove_worst_annotation(header)
        if model.xvl <= model.xvl_cache:
            model.revert()
            break
    print('Exporting post-cross-validation results')
    model.export('post-xv')
    print('Workflow complete.')




# Parse arguments
def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            'Apply the workflow suggested by the fgwas manual'
        )
    )
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument(
        '-i',
        '--input',
        required=True,
        help='Path to fgwas input file'
    )
    required_group.add_argument(
        '-o',
        '--output',
        required=True,
        help='Prefix for output files'
    )
    required_group.add_argument(
        '-t',
        '--threshold',
        required=True,
        type=float,
        help='Likelihood threshold for model improvement'
    )
    required_group.add_argument(
        '-p',
        '--processes',
        required=True,
        type=int,
        help='Maximum number of fgwas processes to launch'
    )
    optional_group = parser.add_argument_group('optional arguments')
    optional_group.add_argument(
        '-b',
        '--bed',
        help='Path to .bed file containing regional definitions'
    )
    optional_group.add_argument(
        '-k',
        '--window',
        help='Window size'
    )
    optional_group.add_argument(
        '-w',
        '--annotations',
        help=(
            '"+"-separated list of annotations to use as a model starting '
            'point.'
        )
    )
    optional_group.add_argument(
        '-l',
        '--llk',
        type=float,
        help=(
            'Required if -w is used. Likelihood of model with the annotation '
            'list given by -w.'
        )
    )
    alternate_workflow_group = parser.add_argument_group(
        'alternate workflow arguments'
    )
    alternate_workflow_group.add_argument(
        '--restrict-to-defined-ci-annotations',
        action='store_true',
        help=(
            'Exclude from analysis any annotations whose individual results '
            'have undefined confidence intervals.'
        )
    )
    alternate_workflow_group.add_argument(
        '--positive-estimates-only',
        action='store_true',
        help=(
            'Add an annotation to the joint model only if it has a positive '
            'parameter estimate.'
        )
    )
    args = parser.parse_args()
    if args.annotations and (not args.llk):
        raise SystemExit(
            'ERROR: When a starting annotation list is provided, a starting '
            'model likelihood must also be provided.'
        )
    elif args.llk and (not args.annotations):
        raise SystemExit(
            'ERROR: A starting likelihood was provided without a starting '
            'annotaiton list.'
        )
    elif args.llk and args.annotations:
        args.initialize = True
    else:
        args.initialize = False
    if args.threshold < 0:
        raise SystemExit(
            'ERROR: The likelihood threshold must be nonnegative.'
        )
    return args




#-------------------------------- Execute -------------------------------------#

if __name__ == '__main__':
    args = parse_arguments()
    main(args)
