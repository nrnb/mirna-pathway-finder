#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import inspect
from mirna_pathway_finder import MirnaPathwayFinder

arg_spec = inspect.getargspec(MirnaPathwayFinder)
arg_defaults = dict(zip(reversed(arg_spec.args), reversed(arg_spec.defaults)))
default_column_number = arg_defaults['query_value_list_column_index'] + 1


def main():
    parser = argparse.ArgumentParser(
        description='''Given a list of mirbase names,
        return an HTML page with a list of matching hits from WikiPathways.''')
    parser.add_argument('ids',
                        type=str,
                        help='identifier or file path to identifier list')
    parser.add_argument('-c', '--column',
                        default=default_column_number,
                        type=int,
                        help='''column number for identifiers in identifier list file
                        (default = {0})'''.format(default_column_number))
    parser.add_argument('-s', '--source',
                        help='''source file path with mappings from
                        WikiPathways identifiers to contained miRNAs.
                        (default = {0})
                        '''.format(arg_defaults['mappings_path']))
    parser.add_argument('-t', '--type',
                        default='rna',
                        help='''input type (rna or protein; default = {0})
                        '''.format(arg_defaults['node_type']))
    parser.add_argument('-o', '--output',
                        default='.',
                        help='''output directory path
                        (default = current working directory)''')
    parser.add_argument('--cache',
                        default=arg_defaults['cache'],
                        type=bool,
                        help='''Cache pathway to miRNA mappings and use in
                        subsequent runs to reduce parse time
                        (default = {0})'''.format(arg_defaults['cache']))
    parser.add_argument('-d', '--debug',
                        default=arg_defaults['debug'],
                        type=bool,
                        help='''Show debug messages (default = {0})
                        '''.format(arg_defaults['debug']))
    args = parser.parse_args()

    query_values = args.ids
    mappings_path = args.source
    query_value_list_column_index = args.column - 1
    node_type = args.type
    output_dir = args.output
    cache = args.cache
    debug = args.debug

    return MirnaPathwayFinder(
        mappings_path=mappings_path,
        query_values=query_values,
        query_value_list_column_index=query_value_list_column_index,
        node_type=node_type,
        output_dir=output_dir,
        cache=cache,
        debug=debug)

if __name__ == '__main__':
    main()
