#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from mirna_pathway_finder import MirnaPathwayFinder


def main():

    # TODO update the input args to whatever we actually need

    parser = argparse.ArgumentParser(
        description='''Given a list of mirbase names, return an HTML page with a list of matching hits from WikiPathways.''')
    parser.add_argument('ids',
                        type=str,
                        help='identifier or file path to identifier list')
    parser.add_argument('-c', '--column',
                        default=1,
                        type=int,
                        help='''column number for identifiers in identifier list file
                        (default = 1)''')
    parser.add_argument('-s', '--source',
                        default='./pathway-to-mirna-mappings.json',
                        help='''source file path with mappings from
                        WikiPathways identifiers to contained miRNAs.
                        (default = file named "pathway-to-mirna-mappings.json"
                            in current working directory)''')
    parser.add_argument('-t', '--type',
                        default='rna',
                        help='input type (rna or protein; default = rna)')
    parser.add_argument('-o', '--output',
                        default='.',
                        help='''output directory path
                        (default = current working directory)''')
    parser.add_argument('--cache',
                        default=True,
                        type=bool,
                        help='''Cache pathway to miRNA mappings and use in
                        subsequent runs to reduce parse time
                        (default = True)''')
    parser.add_argument('-d', '--debug',
                        default=False,
                        type=bool,
                        help='Show debug messages (default = False)')
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
