#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import json
import os
import operator
from os import listdir
from os.path import isfile, join
import pystache
import re
import rx
from rx import Observable, Observer
import urllib

def MirnaPathwayFinder(
        mappings_path='./pathway-to-mirna-mappings.json',
        query_values=None,
        query_value_list_column_index=0,
        node_type='rna',
        output_dir='.',
        cache=True,
        debug=False):

    debug = True

    def print_debug(message):
        if debug:
            print message

    def has_file_extension(filename, expected_file_extension):
        expected_index = len(filename) - len(expected_file_extension)
        actual_index = filename.find(
            expected_file_extension, expected_index, len(filename))
        return actual_index == expected_index

    def has_matching_node(current_mapping_graph, query_value):
        verified_query_value = None
        if current_mapping_graph.has_node(query_value):
            verified_query_value = query_value
        else:
            for onenode in current_mapping_graph.nodes():
                current_node = current_mapping_graph.node[onenode]
                if ((('identifiers' in current_node)
                    and (query_value in current_node['identifiers']))
                    or (('@label' in current_node)
                        and (query_value == current_node['@label']))
                        or (('label' in current_node)
                            and (query_value == current_node['label']))
                        or (('mimat id' in current_node)
                            and (query_value == current_node['mimat id']))
                        or (('name' in current_node)
                            and (query_value == current_node['name']))):
                    if (('mimat id' in current_node)
                            and (query_value == current_node['mimat id'])):
                        verified_query_value = onenode
                        break
                    elif (('identifiers' in current_node)
                            and (query_value in current_node['identifiers'])):
                        verified_query_value = onenode
                    elif not verified_query_value:
                        verified_query_value = onenode
        return verified_query_value

    log_result = dict()
    log_result['total_result_count'] = 0
    log_result['skipped_count'] = 0
    log_result['results_by_source'] = {}

    pathway_to_mirna_mappings = './wp-mir-table-hs.csv'
    pathway_to_mirna_mappings_list = []
    if os.path.isfile(pathway_to_mirna_mappings):
        with open(pathway_to_mirna_mappings, 'rb') as csvfile:
            pathway_to_mirna_mappings_reader = csv.DictReader(csvfile)
            for row in pathway_to_mirna_mappings_reader:
                gene_hits = row['ghits'].split(',')
                mirna_hits = row['mhits'].split(',')
                mirna_target_hits = row['mthits'].split(',')
                wp_identifier = re.search('WP\d+', row['link']).group(0)
                parsed_row = dict([('name', row['name']), ('identifier', wp_identifier), ('id', row['link']), ('gene_hits', gene_hits), ('mirna_hits', mirna_hits), ('mirna_target_hits', mirna_target_hits)])
                pathway_to_mirna_mappings_list.append(parsed_row)
    else:
        if hasattr(pathway_to_mirna_mappings, '__iter__'):
            pathway_to_mirna_mappings_list += pathway_to_mirna_mappings
        else:
            pathway_to_mirna_mappings_list.append(pathway_to_mirna_mappings)

    query_value_list = []
    if os.path.isfile(query_values):
        with open(query_values, 'rb') as csvfile:
            query_value_list_reader = csv.reader(
                csvfile, delimiter='\t', quotechar='|')
            for row in query_value_list_reader:
                query_value_list.append(row[query_value_list_column_index])
    else:
        if hasattr(query_values, '__iter__'):
            query_value_list += query_values
        else:
            query_value_list.append(query_values)

    # TODO handle the case where the query values are NOT display names
    highlight_values = map(lambda query_value: 'label[]=' + urllib.quote(query_value), query_value_list)
    highlight_string = str.join('&', highlight_values) + '&colors=red'

    log_result['query_count'] = len(query_value_list)

    class MyObserver(Observer):
        def on_next(self, x):
            print("Got: %s" % x)

        def on_error(self, e):
            print("Got error: %s" % e)

        def on_completed(self):
            print("Sequence completed")

    def generate_widget_uri(mapping):
        return 'http://www.wikipathways.org/wpi/PathwayWidget.php?id=' + mapping['identifier'] + '&' + highlight_string

    def get_hit_counts(mapping):
        mapping['matching_gene_hits'] = map(lambda matching_gene_hit: dict([('name', matching_gene_hit)]), set(mapping['gene_hits']).intersection(query_value_list))
        mapping['matching_mirna_hits'] = map(lambda matching_gene_hit: dict([('name', matching_gene_hit)]), set(mapping['mirna_hits']).intersection(query_value_list))
        mapping['matching_mirna_target_hits'] = map(lambda matching_gene_hit: dict([('name', matching_gene_hit)]), set(mapping['mirna_target_hits']).intersection(query_value_list))
        return mapping

    def has_hit(mapping):
        hit_count = len(mapping['matching_gene_hits']) + len(mapping['matching_mirna_hits']) + len(mapping['matching_mirna_target_hits'])
        return hit_count > 0

    def uri_encode_mirna_label(mirna):
        mirna['labelUriEncoded'] = urllib.quote(mirna['label'])
        return mirna

    def uri_encode_mirna_label(mirna):
        mirna['labelUriEncoded'] = urllib.quote(mirna['label'])
        return mirna

    def uri_encode_mirna_labels(mapping):
        mirnas_uri_encoded = map(uri_encode_mirna_label, mapping['mirnas'])
        mapping['mirnas_uri_encoded'] = mirnas_uri_encoded
        return mapping

    table_template = '''<!DOCTYPE html>
        <html>
            <body>
                <table>
                    <tr>
                        <th>Name</th>
                        <th>Identifier</th>
                        <th>Matching Genes</th>
                        <th>Matching miRNAs</th>
                        <th>Matching miRNA Targets</th>
                    </tr>
                    {{#.}}<tr>
                        <td id="identifier"><a href="#wikipathways-widget-anchor">{{identifier}}</a></td>
                        <td id="matching-gene-hits">{{#matching_gene_hits}}{{name}},{{/matching_gene_hits}}</td>
                        <td id="matching-mirna-hits">{{#matching_mirna_hits}}{{name}},{{/matching_mirna_hits}}</td>
                        <td id="matching-mirna-target-hits">{{#matching_mirna_target_hits}}{{name}},{{/matching_mirna_target_hits}}</td>
                    </tr>{{/.}}
                </table>
                <a name="wikipathways-widget-anchor"></a>
                <div id="wikipathways-widget">
                    <iframe
                        src="widget_uri"
                        width="600px"
                        height="300px"
                        style="overflow:hidden;">
                    </iframe>
                </div>
                <script src="https://cdnjs.cloudflare.com/ajax/libs/d3/3.5.6/d3.min.js"></script>
                <script src="../mirnapathwayfinder/update-widget.js"></script>
            </body>
        </html>'''

    def generate_pathways_table(mappings):
        widget_uri = generate_widget_uri(mappings[0])
        html_string = pystache.render(table_template, mappings)
        full_html_string = html_string.replace('widget_uri', widget_uri)
        f = open('./demos/index.html', 'w')
        f.write(full_html_string)
        return full_html_string

    mappings_source = Observable.from_(pathway_to_mirna_mappings_list)
    matching_mappings = mappings_source.map(get_hit_counts).filter(has_hit).to_list().map(generate_pathways_table)
    matching_mappings.subscribe(MyObserver())
