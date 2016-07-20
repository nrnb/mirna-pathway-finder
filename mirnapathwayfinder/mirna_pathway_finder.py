#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
# import numpy as np
import os
import pandas as pd
import pystache
import re
import urllib

current_path = os.path.dirname(os.path.abspath(__file__))


def MirnaPathwayFinder(
        mappings_path=None,
        query_values=None,
        query_value_list_column_index=0,
        node_type='rna',
        output_dir='.',
        cache=True,
        debug=False):

    def print_debug(message):
        if debug:
            print message

    def generate_widget_uri(mapping, highlight_string):
        return ''.join([
            'http://www.wikipathways.org/wpi/PathwayWidget.php?id=',
            mapping['identifier'],
            '&',
            highlight_string,
        ])

    def generate_pathways_table(html_template_input, query_value_list):
        # TODO handle the case where the query values are NOT display names
        highlight_values = map(
            lambda query_value: 'label[]=' + urllib.quote(
                query_value
            ), query_value_list
        )
        highlight_string = str.join('&', highlight_values) + '&colors=red'

        f = open(current_path + '/table-template.html', 'r')
        table_template = f.read()
        widget_uri = generate_widget_uri(
            html_template_input[0], highlight_string)
        initial_html_string = pystache.render(
            table_template, html_template_input)
        html_string_with_widget_url = initial_html_string.replace(
            'widget_uri',
            widget_uri
        )

        update_widget_path = os.path.join(current_path, 'update-widget.js')
        with open(update_widget_path, 'r') as update_widget:
            update_widget_string = ''.join([
                'var highlightString = \'',
                highlight_string,
                '\';\n',
                update_widget.read()
            ])

        html_string_with_update_widget = html_string_with_widget_url.replace(
            'update_widget_string',
            update_widget_string
        )
        f = open(os.path.join(output_dir, 'pathways.html'), 'w')
        f.write(html_string_with_update_widget)
        return html_string_with_update_widget

    def has_targeter(row):
        columns_to_check = (
            'stem_loop_name',
            'mature_name',
            'mirbase',
            'mirbase.mature',
            'ncbigene',
            'hgnc',
            'targeter_stem_loop_name',
            'targeter_mature_name',
            'targeter_mirbase',
            'targeter_mirbase.mature',
            'targeter_ncbigene',
            'targeter_hgnc',
            )

        possible_targeters = set(filter(lambda y: isinstance(y, str),
                                 map(lambda x: row[x], columns_to_check)))
        return len(possible_targeters.intersection(query_value_list)) > 0

    def has_target(row):
        columns_to_check = (
            'targeter_stem_loop_name',
            'targeter_mature_name',
            'targeter_mirbase',
            'targeter_mirbase.mature',
            'targeter_ncbigene',
            'targeter_hgnc',
            )

        possible_targeters = set(filter(lambda y: isinstance(y, str),
                                 map(lambda x: row[x], columns_to_check)))
        return len(possible_targeters.intersection(query_value_list)) > 0

    results_limit = 20

    query_value_list = set()
    if os.path.isfile(query_values):
        with open(query_values, 'rb') as csvfile:
            query_value_list_reader = csv.reader(
                csvfile, delimiter='\t', quotechar='|')
            for row in query_value_list_reader:
                query_value_list.add(row[query_value_list_column_index])
    else:
        if hasattr(query_values, '__iter__'):
            query_value_list = set(query_values)
        else:
            query_value_list.add(query_values)

    query_value_list = map(lambda x: re.sub(
        '^http:\/\/identifiers.org\/(hgnc|ncbigene|mirbase|mirbase\.mature)\/',
        '',
        x
        ), query_value_list)

    if mappings_path is None:
        # TODO remove the date part of the file name
        mappings_path = os.path.join(
            current_path, '..', 'wp-mir-table-hs-20160715.tsv')

    # TODO integrate this old code into the current code. Specifically,
    #      handle if the mapping data is provided as a Python object.
    #    # parse wp-mir-table-hs.csv (or other file, if specified)
    #    # to get mappings between pathways and mirnas,
    #    # including for each pathway:
    #    # * genes: all gene products in the pathway, annotated as genes
    #    # * mirna_hits_as_gene_specified: miRNAs actually existing in pathway,
    #    #                                 annotated as genes
    #    # * mirna_hits_as_mirna_specified: miRNAs actually shown on pathway,
    #    #                                  annotated as miRNAs
    #    # * mirna_hits_as_mirna_inferred: miRNAs NOT actually specified on the
    #    #       pathway but
    #    #       inferred to exist on the pathway because they target genes or
    #    #       proteins that DO actually exist on the pathway
    #    pathway_to_mirna_mappings = mappings_path
    #    pathway_to_mirna_mappings_list = []
    #    if os.path.isfile(pathway_to_mirna_mappings):
    #        with open(pathway_to_mirna_mappings, 'rb') as csvfile:
    #            pathway_to_mirna_mappings_reader = csv.DictReader(csvfile)
    #            for row in pathway_to_mirna_mappings_reader:
    #                genes = parse_hits_field(row['genes'])
    #                mirna_hits_as_gene_specified = parse_hits_field(
    #                                                               row['ghits'])
    #                mirna_hits_as_mirna_specified = parse_hits_field(
    #                                                                row['mhits'])
    #                mirna_hits_as_mirna_inferred = parse_hits_field(
    #                                                               row['mthits'])
    #
    #                wp_identifier = re.search('WP\d+', row['link']).group(0)
    #                parsed_row = {
    #                    'name': row['name'],
    #                    'identifier': wp_identifier,
    #                    'id': row['link'],
    #                    'genes': genes,
    #                    'mirna_hits_as_gene_specified':
    #                        mirna_hits_as_gene_specified,
    #                    'mirna_hits_as_mirna_specified':
    #                        mirna_hits_as_mirna_specified,
    #                    'mirna_hits_as_mirna_inferred':
    #                        mirna_hits_as_mirna_inferred,
    #                }
    #                pathway_to_mirna_mappings_list.append(parsed_row)
    #    else:
    #        if hasattr(pathway_to_mirna_mappings, '__iter__'):
    #            pathway_to_mirna_mappings_list += pathway_to_mirna_mappings
    #        else:
    #            pathway_to_mirna_mappings_list.append(pathway_to_mirna_mappings)
    wp_mirna = pd.read_csv(mappings_path,
                           sep='\t',
                           dtype=str)

    # TODO remove this. It's just for dev.
    # wp_mirna = wp_mirna.head(1000)

    with_targeter = wp_mirna[wp_mirna.apply(
        lambda d: has_targeter(d), axis=1)]

    with_targeter_by_pwy = with_targeter.groupby(['wikipathways'])

    # get targeter count by pathway
    targeter_n_by_pwy = with_targeter_by_pwy[
        'shown_or_inferred_mature_name'].nunique()

    # get shown_targeter count by pathway
    shown_targeter_n_by_pwy = with_targeter_by_pwy[
        'mature_name'].nunique()

    # get target count by pathway
    with_target = with_targeter[with_targeter.apply(
        lambda d: has_target(d), axis=1)]
    with_target_by_pwy = with_target.groupby(['wikipathways'])
    target_n_by_pwy = with_target_by_pwy['ncbigene'].nunique()

    pathways = with_targeter['wikipathways'].unique()
    d = {
        'wikipathways': pathways,
        'targeter_count': targeter_n_by_pwy,
        'shown_targeter_count': shown_targeter_n_by_pwy,
        'target_count': target_n_by_pwy
    }

    wp_counts = pd.DataFrame(
        data=d
        ).fillna(
            value=0
        ).nlargest(
            results_limit,
            ['shown_targeter_count', 'target_count', 'targeter_count'])

    results = with_targeter.join(wp_counts,
                                 on='wikipathways',
                                 how='left',
                                 lsuffix='',
                                 rsuffix='_r',
                                 sort=False)

    pathway_level_columns = [
        'shown_targeter_count',
        'target_count',
        'targeter_count',
        'wikipathways',
        'link',
        'pathway_name'
    ]
    results_sorted = results.sort_values(by=pathway_level_columns,
                                         axis=0,
                                         ascending=False,
                                         inplace=False,
                                         kind='quicksort',
                                         na_position='last')

    results_by_pwy = results_sorted.groupby(pathway_level_columns, sort=False)

    html_template_input = []
    for name_by_pwy, group_by_pwy in results_by_pwy:
        targets_by_targeters = []
        shown_targeters = []
        mature_names = filter(lambda x: isinstance(x, str),
                              group_by_pwy['mature_name'])
        if len(mature_names) > 0:
            shown_targeters = shown_targeters + map(
                lambda x: {'name': x}, mature_names)
        result = {
            'id': name_by_pwy[pathway_level_columns.index('link')],
            'identifier': name_by_pwy[
                pathway_level_columns.index('wikipathways')],
            'name': name_by_pwy[pathway_level_columns.index('pathway_name')],
            'targets_by_targeters': targets_by_targeters,
            'shown_targeters': shown_targeters,
            'shown_targeter_count': name_by_pwy[
                pathway_level_columns.index('shown_targeter_count')],
            'target_count': name_by_pwy[
                pathway_level_columns.index('target_count')],
            'targeter_count': name_by_pwy[
                pathway_level_columns.index('targeter_count')]
        }
        by_targeter = group_by_pwy.groupby('targeter_mature_name')
        for name_by_targeter, group_by_targeter in by_targeter:
            targets_by_targeters.append({
                'targeter': name_by_targeter,
                'targets': ', '.join(group_by_targeter[
                    'ncbigene'].unique().tolist())
            })
        html_template_input.append(result)

    generate_pathways_table(html_template_input, query_value_list)
