#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
from functional import seq
from math import fsum
import os
import pystache
import re
from rx import Observable, Observer
import urllib


def MirnaPathwayFinder(
        mappings_path=None,
        query_values=None,
        query_value_list_column_index=0,
        node_type='rna',
        output_dir='.',
        cache=True,
        debug=False):

    def parse_hits_field(field):
        return filter(lambda x: x != '', field.split(','))

    def print_debug(message):
        if debug:
            print message

    current_path = os.path.dirname(os.path.abspath(__file__))
    if mappings_path is None:
        mappings_path = os.path.join(current_path, '..', 'wp-mir-table-hs.csv')

    # Get mappings from gene to miRNA
    gene_mir_map_path = os.path.join(
        current_path, '..', 'wp-mir-table-builder', 'gene-mir-map.txt'
    )
    gene_to_mirna_mappings = {}
    with open(gene_mir_map_path, 'rU') as csvfile:
        gene_to_mirna_mappings_reader = csv.DictReader(
            csvfile, delimiter='\t')
        for row in gene_to_mirna_mappings_reader:
            mirna = row['miRBase ID(s)']
            gene = row['EntrezGene ID']
            if (mirna is not None) and (gene is not None):
                gene_to_mirna_mappings[gene] = mirna

    # Get gene targets for each miRNA
    mir_gene_targets_path = os.path.join(
        current_path, '..', 'wp-mir-table-builder', 'mir-gene-targets.txt'
    )
    targeted_genes_by_mirna = {}
    targeted_genes_by_mirna_handled_mirnas = set()
    targeted_genes_by_mirna_entries = set()
    with open(mir_gene_targets_path, 'rU') as csvfile:
        targeted_genes_by_mirna_reader = csv.DictReader(
            csvfile,
            delimiter='\t'
        )
        for row in targeted_genes_by_mirna_reader:
            mirna = row['miRNA']
            gene = row['target']
            if mirna not in targeted_genes_by_mirna_handled_mirnas:
                targeted_genes_by_mirna_handled_mirnas.add(mirna)
                targeted_genes_by_mirna[mirna] = set()

            key = (mirna + gene)
            if key not in targeted_genes_by_mirna_entries:
                targeted_genes_by_mirna_entries.add(key)
                targeted_genes_by_mirna[mirna].add(gene)

    targeted_genes_by_mirna_keys = targeted_genes_by_mirna.keys()

    # parse wp-mir-table-hs.csv (or other file, if specified)
    # to get mappings between pathways and mirnas,
    # including for each pathway:
    # * genes: all gene products in the pathway, annotated as genes
    # * mirna_hits_as_gene_specified: miRNAs actually existing in pathway,
    #   annotated as genes
    # * mirna_hits_as_mirna_specified: miRNAs actually existing in pathway,
    #   annotated as miRNAs
    # * mirna_hits_as_mirna_inferred: miRNAs NOT actually specified on the
    #       pathway but
    #       inferred to exist on the pathway because they target genes or
    #       proteins that DO actually exist on the pathway
    pathway_to_mirna_mappings = mappings_path
    pathway_to_mirna_mappings_list = []
    if os.path.isfile(pathway_to_mirna_mappings):
        with open(pathway_to_mirna_mappings, 'rb') as csvfile:
            pathway_to_mirna_mappings_reader = csv.DictReader(csvfile)
            for row in pathway_to_mirna_mappings_reader:
                genes = parse_hits_field(row['genes'])
                mirna_hits_as_gene_specified = parse_hits_field(row['ghits'])
                mirna_hits_as_mirna_specified = parse_hits_field(row['mhits'])
                mirna_hits_as_mirna_inferred = parse_hits_field(row['mthits'])

                wp_identifier = re.search('WP\d+', row['link']).group(0)
                parsed_row = {
                    'name': row['name'],
                    'identifier': wp_identifier,
                    'id': row['link'],
                    'genes': genes,
                    'mirna_hits_as_gene_specified':
                        mirna_hits_as_gene_specified,
                    'mirna_hits_as_mirna_specified':
                        mirna_hits_as_mirna_specified,
                    'mirna_hits_as_mirna_inferred':
                        mirna_hits_as_mirna_inferred,
                }
                pathway_to_mirna_mappings_list.append(parsed_row)
    else:
        if hasattr(pathway_to_mirna_mappings, '__iter__'):
            pathway_to_mirna_mappings_list += pathway_to_mirna_mappings
        else:
            pathway_to_mirna_mappings_list.append(pathway_to_mirna_mappings)

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

    # TODO handle the case where the query values are NOT display names
    highlight_values = map(
        lambda query_value: 'label[]=' + urllib.quote(
            query_value
        ), query_value_list
    )
    highlight_string = str.join('&', highlight_values) + '&colors=red'

    class MyObserver(Observer):
        def on_next(self, x):
            print_debug("Got: %s" % x)

        def on_error(self, e):
            print_debug("Got error: %s" % e)

        def on_completed(self):
            print_debug("Sequence completed")

    def generate_widget_uri(mapping):
        return ''.join([
            'http://www.wikipathways.org/wpi/PathwayWidget.php?id=',
            mapping['identifier'],
            '&',
            highlight_string,
        ])

    def get_matches_and_counts(mapping):
        # Find matching miRNAs in pathways, annotated as genes or miRNAs
        mirna_hits_as_gene_specified = set(
            mapping['mirna_hits_as_gene_specified'])
        mirna_hits_as_mirna_specified = set(
            mapping['mirna_hits_as_mirna_specified'])

        mirna_matches_as_gene_specified = query_value_list.intersection(
            mirna_hits_as_gene_specified
        )
        mirna_matches_as_mirna_specified = query_value_list.intersection(
            mirna_hits_as_mirna_specified
        )

        mirna_matches_specified = mirna_matches_as_gene_specified.union(
            mirna_matches_as_mirna_specified
        )

        mapping['mirna_matches_specified'] = map(
            lambda mirna: {'name': mirna},
            mirna_matches_specified
        )

        mapping['mirna_match_specified_count'] = len(
            mapping['mirna_matches_specified']
        )

        # Find all miRNAs in pathways, including both those actually in
        # the pathway and those that are inferred based on the
        # presence of an miRNA target (a gene or protein) in the pathway
        mirna_matches_as_mirna_inferred = query_value_list.intersection(
            mirna_hits_as_mirna_inferred
        )

        mirna_matches_all = mirna_matches_as_mirna_specified.union(
            mirna_matches_as_mirna_inferred
        )

        if len(mirna_matches_as_gene_specified) > 0:
            # user specified desired miRNA(s) as gene(s), so
            # we need to specify it as an miRNA so that we can
            # look at the targeting
            for gene in mirna_matches_as_gene_specified:
                mirna_matches_all.add(gene_to_mirna_mappings[gene])

        mapping['mirna_matches_all'] = map(
            lambda mirna: {'name': mirna},
            mirna_matches_all
        )

        mapping['mirna_match_all_count'] = len(
            mapping['mirna_matches_all']
        )

        # Find miRNA targets (genes and proteins) in pathways
        # TODO from mir-gene-targets.txt, get the gene targets.
        # Then get the intersection between that result and genes

        genes = set(mapping['genes'])

        def create_target_title_string(mirna):
            if mirna not in targeted_genes_by_mirna_keys:
                return 'none'

            targeted_genes = targeted_genes_by_mirna[mirna]
            if targeted_genes is None:
                return 'none'

            targeted_genes_on_pathway = targeted_genes.intersection(genes)
            if len(targeted_genes_on_pathway) == 0:
                return 'none'

            return ', '.join(
                targeted_genes_on_pathway
            )

        if len(mirna_matches_all) > 0:
            mapping['target_matches_by_mirna'] = (
                seq(mirna_matches_all)
                .map(lambda mirna: {
                    'mirna': mirna,
                    'targets': create_target_title_string(mirna)
                })
                )

            target_matches_all = (
                seq(mirna_matches_all)
                .filter(lambda mirna: mirna in targeted_genes_by_mirna_keys)
                .map(lambda mirna: targeted_genes_by_mirna[mirna])
                .filter(lambda targets: hasattr(targets, '__iter__'))
                .aggregate(set(), lambda acc, targets: acc.union(targets))
                ).intersection(genes)
        else:
            target_matches_all = []

        if len(target_matches_all) > 0:
            mapping['target_matches'] = map(
                lambda target_match: {'name': target_match},
                target_matches_all
            )
        else:
            mapping['target_matches'] = []

        mapping['target_match_count'] = len(mapping['target_matches'])

        return mapping

    def has_match(mapping):
        match_count = fsum([
            mapping['mirna_match_all_count'],
            mapping['target_match_count'],
        ])
        return match_count > 0

    def sort_by_match_counts(mappings):
        return sorted(
            mappings,
            key=lambda mapping: (
                mapping['mirna_match_specified_count'],
                mapping['mirna_match_all_count'],
                mapping['target_match_count']
            ),
            reverse=True
        )

    def take_top_matches(mappings):
        return mappings[0:20]

    def generate_pathways_table(mappings):
        f = open(current_path + '/table-template.html', 'r')
        table_template = f.read()
        widget_uri = generate_widget_uri(mappings[0])
        initial_html_string = pystache.render(table_template, mappings)
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

    # Run and generate HTML output
    mappings_source = Observable.from_(pathway_to_mirna_mappings_list)

    mappings = mappings_source.map(
        get_matches_and_counts
    ).filter(
        has_match
    ).to_list(
    ).filter(
        lambda x: len(x) > 0
    ).map(
        sort_by_match_counts
    ).map(
        take_top_matches
    ).map(generate_pathways_table)

    mappings.subscribe(MyObserver())
