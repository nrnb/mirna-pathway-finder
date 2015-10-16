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
        mappings_path='./wp-mir-table-hs.csv',
        query_values=None,
        query_value_list_column_index=0,
        node_type='rna',
        output_dir='.',
        cache=True,
        debug=False):

    def print_debug(message):
        if debug:
            print message

    pathway_to_mirna_mappings = mappings_path
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
    # TODO use the gene_hits, mirna hits and mirna targets hits instead of just the query values to create this string
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
        # NOTE: we're not currently using the matching gene hits. They are intended to
        # be used for specifying which node(s) to highlight, in the case that miRNAs
        # are annotated with gene ids.
        mapping['matching_gene_hits'] = map(lambda matching_gene_hit: dict([('name', matching_gene_hit)]), set(mapping['gene_hits']).intersection(query_value_list))
        mapping['matching_gene_hit_count'] = len(mapping['matching_gene_hits'])
        mapping['matching_mirna_hits'] = map(lambda matching_mrna_hit: dict([('name', matching_mrna_hit)]), set(mapping['mirna_hits']).intersection(query_value_list))
        mapping['matching_mirna_hit_count'] = len(mapping['matching_mirna_hits'])
        mapping['matching_mirna_target_hits'] = map(lambda mirna_target_hit: dict([('name', mirna_target_hit)]), set(mapping['mirna_target_hits']).intersection(query_value_list))
        mapping['matching_mirna_target_hit_count'] = len(mapping['matching_mirna_target_hits'])
        return mapping

    def has_hit(mapping):
        hit_count = mapping['matching_gene_hit_count'] + mapping['matching_mirna_hit_count'] + mapping['matching_mirna_target_hit_count']
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
		<p>The table lists pathways containing your miRNAs of interest and/or protein targets of your miRNAs. The first column lists a clickable pathway title that updates the Interactive Pathway Viewer below. The second column lists pathway identifers that link to WikiPathways.org. The list is sorted by "miRNAs" (primary) and by "miRNA Targets" (secondary) found on each pathway. The top 20 results are listed.</p>
		<h3>Table of Pathway Results</h3>
                <table id="pathway-to-mirna">
                    <tr>
                        <th>Pathway Title <i>(click to view pathway)</i></th>
                        <th>Linkout</th>
                        <th>miRNAs</th>
                        <th>miRNA Targets</th>
                    </tr>
                    {{#.}}<tr>
                        <td id="name" title="view pathway"><a href="#wikipathways-widget-anchor">{{name}}</a></td>
                        <td id="identifier" title="view pathway at WikiPathways.org"><a href="{{id}}" target="_blank">{{identifier}}</a></td>
                        <td align="center" id="matching-mirna-hits" title="{{matching_mirna_hits}}">{{matching_mirna_hit_count}}</td>
                        <td align="center" id="matching-mirna-target-hits" title="{{matching_mirna_target_hits}}">{{matching_mirna_target_hit_count}}</td>
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
        return sorted(mappings, key=lambda mapping: (mapping['matching_mirna_hit_count'], mapping['matching_mirna_target_hit_count']), reverse=True)

    def take_top_hits(mappings):
        return mappings[0:19]

    def generate_pathways_table(mappings):
        widget_uri = generate_widget_uri(mappings[0])
        initial_html_string = pystache.render(table_template, mappings)
        html_string_with_widget_url = initial_html_string.replace('widget_uri', widget_uri)
        with open('./mirnapathwayfinder/update-widget.js', 'r') as update_widget:
            update_widget_string = 'var highlightString = \'' + highlight_string + '\';\n' + update_widget.read()
        html_string_with_update_widget = html_string_with_widget_url.replace('update_widget_string', update_widget_string)
        f = open(output_dir + '/pathways.html', 'w')
        f.write(html_string_with_update_widget)
        return html_string_with_update_widget

    mappings_source = Observable.from_(pathway_to_mirna_mappings_list)
    matching_mappings = mappings_source.map(get_hits_and_counts).filter(has_hit).to_list().map(sort_by_hit_counts).map(take_top_hits).map(generate_pathways_table)
    matching_mappings.subscribe(MyObserver())
