#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
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

    def print_debug(message):
        if debug:
            print message

    current_path = os.path.dirname(os.path.abspath(__file__))
    if mappings_path is None:
        mappings_path = os.path.join(current_path, '..', 'wp-mir-table-hs.csv')

    # parse wp-mir-table-hs.csv (or other file, if specified)
    # to get mappings between pathways and mirnas,
    # including for each pathway:
    # * gene_hits: targeted genes
    # * mirna_hits: miRNAs actually existing in pathway
    # * targeting_mirna_hits: miRNAs actually existing in the pathway or
    #       inferred to exist in the pathway because they target genes or
    #       proteins that actually exist in the pathway
    pathway_to_mirna_mappings = mappings_path
    pathway_to_mirna_mappings_list = []
    if os.path.isfile(pathway_to_mirna_mappings):
        with open(pathway_to_mirna_mappings, 'rb') as csvfile:
            pathway_to_mirna_mappings_reader = csv.DictReader(csvfile)
            for row in pathway_to_mirna_mappings_reader:
                gene_hits = row['ghits'].split(',')
                mirna_hits = row['mhits'].split(',')
                targeting_mirna_hits = row['mthits'].split(',')
                wp_identifier = re.search('WP\d+', row['link']).group(0)
                parsed_row = {
                    'name': row['name'],
                    'identifier': wp_identifier,
                    'id': row['link'],
                    'gene_hits': gene_hits,
                    'mirna_hits': mirna_hits,
                    'targeting_mirna_hits':
                        targeting_mirna_hits,
                }
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
    # TODO use the gene_hits, mirna hits and mirna targets hits instead of
    # just the query values to create this string
    highlight_values = map(lambda query_value: 'label[]=' + urllib.quote(query_value), query_value_list)
    highlight_string = str.join('&', highlight_values) + '&colors=red'

    class MyObserver(Observer):
        def on_next(self, x):
            print_debug("Got: %s" % x)

        def on_error(self, e):
            print_debug("Got error: %s" % e)

        def on_completed(self):
            print_debug("Sequence completed")

    def generate_widget_uri(mapping):
        return 'http://www.wikipathways.org/wpi/PathwayWidget.php?id=' + mapping['identifier'] + '&' + highlight_string

    def get_hits_and_counts(mapping):
        # NOTE: we're not currently using the matching gene hits.
        # They are intended to
        # be used for specifying which node(s) to highlight,
        # in the case that miRNAs
        # are annotated with gene ids.
        gene_query_value_list = set(mapping['gene_hits']).intersection(query_value_list)
        mapping['gene_hits'] = map(lambda gene_hit: {'name': gene_hit}, gene_query_value_list)
        mapping['gene_hit_count'] = len(mapping['gene_hits'])

        mirna_query_value_list = set(mapping['mirna_hits']).intersection(query_value_list)
        mapping['mirna_hits'] = map(lambda mirna_hit: {'name': mirna_hit}, mirna_query_value_list)
        mapping['mirna_hit_count'] = len(mapping['mirna_hits'])

        targeting_mirna_query_value_list = set(mapping['targeting_mirna_hits']).intersection(query_value_list)
        mapping['targeting_mirna_hits'] = map(lambda targeting_mirna_hit: {'name': targeting_mirna_hit}, targeting_mirna_query_value_list)
        mapping['targeting_mirna_hit_count'] = len(mapping['targeting_mirna_hits'])

        return mapping

    def has_hit(mapping):
        hit_count = mapping['gene_hit_count'] + mapping['mirna_hit_count'] + mapping['targeting_mirna_hit_count']
        return hit_count > 0

    table_template = '''<!DOCTYPE html>
        <html>
            <style>
                body {
                    font: 11px tahoma,arial,helvetica,sans-serif;
                }
                #pathway-to-mirna {
                        font-family: "tahoma,arial,helvetica,sans-serif"
                        font-size: 11px;
                        background: none repeat scroll 0% 0% #FFF;
                        width: 600px;
                        border-collapse: collapse;
                        text-align: left;
                        margin: 20px;
                }
                #pathway-to-mirna th {
                        border-bottom: 1px solid #AAA;
                        padding: 0px 0px 3px;
                }
                #pathway-to-mirna tr:hover {
                        background-color:#DDD;
                }
                a:link {
                    text-decoration: none;
                }

                a:visited {
                    text-decoration: none;
                }

                a:hover {
                    text-decoration: underline;
                }

                a:active {
                    text-decoration: underline;
                }
                p {
                        max-width:600px;
                }
            </style>
            <body>
                <h2>Pathway Finder</h2>
                <p>The table lists pathways containing your miRNAs of interest and/or gene product targets of your miRNAs. The first column lists a clickable pathway title that updates the Interactive Pathway Viewer below. The second column lists pathway identifiers that link to WikiPathways.org. The third column shows the miRNA count in each pathway. The last column shows the count of the gene products targeted by miRNAs and, in parentheses, the count of the miRNAs doing the targeting.</p>
                <p>The list is sorted by "miRNAs" (primary) and by "miRNA Targets" (secondary) found on each pathway. The results are limited to the top 20.</p>
                <h3>Table of Pathway Results</h3>
                <table id="pathway-to-mirna">
                    <tr>
                        <th>Pathway Title <i>(click to view pathway)</i></th>
                        <th>Linkout</th>
                        <th>miRNAs on Pathway</th>
                        <th>Targets on Pathway (Targeting miRNAs)</th>
                    </tr>
                    {{#.}}<tr>
                        <td id="name" title="view pathway"><a href="#wikipathways-widget-anchor">{{name}}</a></td>
                        <td id="identifier" title="view pathway at WikiPathways.org"><a href="{{id}}" target="_blank">{{identifier}}</a></td>
                        <td align="center" id="matching-mirna-hits" title="{{#mirna_hits}} {{name}}\n{{/mirna_hits}}">{{mirna_hit_count}}</td>
                        <td align="center" id="matching-mirna-target-hits" title="{{#targeting_mirna_hits}} {{name}}\n{{/targeting_mirna_hits}}">0 ({{targeting_mirna_hit_count}})</td>
                    </tr>{{/.}}
                </table>
                <a name="wikipathways-widget-anchor"></a>
                <div id="wikipathways-widget">
                  <h3>Interactive Pathway Viewer</h3>
                    <iframe
                        src="widget_uri"
                        width="600px"
                        height="300px"
                        style="overflow:hidden;">
                    </iframe>
                </div>
                <script src="https://cdnjs.cloudflare.com/ajax/libs/d3/3.5.6/d3.min.js"></script>
                <script type="text/javascript">
                    update_widget_string
                </script>
            </body>
        </html>'''

    def sort_by_hit_counts(mappings):
        return sorted(mappings, key=lambda mapping: (mapping['mirna_hit_count'], mapping['targeting_mirna_hit_count']), reverse=True)

    def take_top_hits(mappings):
        return mappings[0:20]

    def generate_pathways_table(mappings):
        widget_uri = generate_widget_uri(mappings[0])
        initial_html_string = pystache.render(table_template, mappings)
        html_string_with_widget_url = initial_html_string.replace('widget_uri', widget_uri)

        update_widget_path = os.path.join(current_path, 'update-widget.js')
        with open(update_widget_path, 'r') as update_widget:
            update_widget_string = 'var highlightString = \'' + highlight_string + '\';\n' + update_widget.read()
        html_string_with_update_widget = html_string_with_widget_url.replace('update_widget_string', update_widget_string)
        f = open(os.path.join(output_dir, 'pathways.html'), 'w')
        f.write(html_string_with_update_widget)
        return html_string_with_update_widget

    # Run to return result
    mappings_source = Observable.from_(pathway_to_mirna_mappings_list)
    mappings = mappings_source.map(get_hits_and_counts).filter(has_hit).to_list().map(sort_by_hit_counts).map(take_top_hits).map(generate_pathways_table)
    mappings.subscribe(MyObserver())
