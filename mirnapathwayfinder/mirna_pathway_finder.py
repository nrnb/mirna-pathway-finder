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
    #highlight_values = map(lambda query_value: 'display-name:' + urllib.quote(query_value), query_value_list)
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

#import csv
#import os
#from os import listdir
#from os.path import isfile, join
#import networkx as nx
#import networkxgmml
#import re
#import regexes
#
#def TargetInteractionFinder(
#        source_xgmml='./source_xgmml/',
#        node_ids=None,
#        node_id_list_column_index=0,
#        node_type='rna',
#        output_dir='.',
#        cache=True,
#        debug=False):
#
#    def print_debug(message):
#        if debug:
#            print message
#
#    def get_file_paths(filelike_object, file_extension):
#        if hasattr(filelike_object, 'name'):
#            return filelike_object.name
#        else:
#            source_xgmml_file_paths = []
#            if os.path.isdir(filelike_object):
#                source_xgmml_file_paths += [filelike_object + f for f in
#                                            listdir(filelike_object)
#                                            if isfile(join(filelike_object, f))
#                                            and has_file_extension(
#                                                f, file_extension)]
#            elif os.path.isfile(filelike_object):
#                source_xgmml_file_paths += [filelike_object]
#            else:
#                print filelike_object
#                raise ValueError('File or directory input error.')
#
#            return source_xgmml_file_paths
#
#    def has_file_extension(filename, expected_file_extension):
#        expected_index = len(filename) - len(expected_file_extension)
#        actual_index = filename.find(
#            expected_file_extension, expected_index, len(filename))
#        return actual_index == expected_index
#
#    def has_matching_node(current_mapping_graph, node_id):
#        verified_node_id = None
#        if current_mapping_graph.has_node(node_id):
#            verified_node_id = node_id
#        else:
#            for onenode in current_mapping_graph.nodes():
#                current_node = current_mapping_graph.node[onenode]
#                if ((('identifiers' in current_node)
#                    and (node_id in current_node['identifiers']))
#                    or (('@label' in current_node)
#                        and (node_id == current_node['@label']))
#                        or (('label' in current_node)
#                            and (node_id == current_node['label']))
#                        or (('mimat id' in current_node)
#                            and (node_id == current_node['mimat id']))
#                        or (('name' in current_node)
#                            and (node_id == current_node['name']))):
#                    if (('mimat id' in current_node)
#                            and (node_id == current_node['mimat id'])):
#                        verified_node_id = onenode
#                        break
#                    elif (('identifiers' in current_node)
#                            and (node_id in current_node['identifiers'])):
#                        verified_node_id = onenode
#                    elif not verified_node_id:
#                        verified_node_id = onenode
#        return verified_node_id
#
#    results_file_path = output_dir + 'interactions.csv'
#    results_file = open(results_file_path, 'w')
#    results_file.write('queryid,targetid,score,pvalue,pmid,datasource\n')
#    results_file.close()
#
#    logfile_path = output_dir + 'summary.log'
#    open(logfile_path, 'w').close()
#
#    log_result = dict()
#    log_result['total_result_count'] = 0
#    log_result['skipped_count'] = 0
#    log_result['results_by_source'] = {}
#
#    file_extension = '.xgmml'
#    source_xgmml_file_paths = []
#    if hasattr(source_xgmml, '__iter__'):
#        for source_xgmml_filelike_object in source_xgmml:
#            source_xgmml_file_paths += get_file_paths(
#                source_xgmml_filelike_object, file_extension)
#    else:
#        source_xgmml_file_paths += get_file_paths(source_xgmml, file_extension)
#
#    print_debug('source_xgmml_file_paths')
#    print_debug(source_xgmml_file_paths)
#
#    versions = []
#
#    node_id_list = []
#    if os.path.isfile(node_ids):
#        with open(node_ids, 'rb') as csvfile:
#            node_id_list_reader = csv.reader(
#                csvfile, delimiter='\t', quotechar='|')
#            for row in node_id_list_reader:
#                node_id_list.append(row[node_id_list_column_index])
#    else:
#        if hasattr(node_ids, '__iter__'):
#            node_id_list += node_ids
#        else:
#            node_id_list.append(node_ids)
#
#    log_result['query_count'] = len(node_id_list)
#
#    for source_xgmml_file_path in source_xgmml_file_paths:
#        source_xgmml_dirname = os.path.dirname(source_xgmml_file_path)
#        source_xgmml_basename = os.path.basename(source_xgmml_file_path)
#
#        version = source_xgmml_basename.replace('.xgmml', '')
#        versions.append(version)
#
#        cached_graph_filename = '.' + source_xgmml_basename.replace(
#            '.xgmml', '.p')
#        cached_graph_filepath = os.path.join(
#            source_xgmml_dirname, cached_graph_filename)
#        if cache and os.path.isfile(cached_graph_filepath):
#            print_debug('Using cached version of source_xgmml.')
#            current_mapping_graph = nx.read_gpickle(cached_graph_filepath)
#        else:
#            current_xgmml_file = open(source_xgmml_file_path)
#            current_mapping_graph = networkxgmml.XGMMLReader(
#                current_xgmml_file)
#            nx.write_gpickle(current_mapping_graph, cached_graph_filepath)
#
#        log_result[version + '_result_count'] = 0
#
#        datasource_subgraph_node_ids = []
#
#        id_mappings = {}
#
#        for node_id in node_id_list:
#            verified_node_id = None
#            name_or_identifier = None
#            if not node_id:
#                log_result['skipped_count'] += 1
#                if debug:
#                    print '''Warning: No node_id in line. Please provide a
#                    name (e.g., hsa-miR-542-3p),
#                    identifier (MIMAT0003389) or
#                    IRI (http://identifiers.org/mirbase.mature/MIMAT0003389)'''
#            else:
#                if regexes.identifiers_org.match(node_id):
#                    name_or_identifier = re.sub(
#                        regexes.identifiers_org, '', node_id)
#
#                verified_node_id = has_matching_node(
#                    current_mapping_graph, node_id)
#                # If provided node_id is an identifiers.org IRI,
#                # but it doesn't match any nodes, we can
#                # try to find nodes by matching the identifier
#                # alone (stripping off the namespace).
#                if name_or_identifier and (not verified_node_id):
#                    verified_node_id = has_matching_node(
#                        current_mapping_graph, name_or_identifier)
#
#                    if not verified_node_id:
#                        print 'Warning: No node found as a match for '
#                        print '\t\tnode_id "' + node_id + '"'
#                        print '\t\tin datasource "' + version + '"'
#
#                # if node is verified to exist in this graph
#                if verified_node_id:
#                    if (verified_node_id != node_id):
#                        print_debug('adding to mappings')
#                        print_debug('verified_node_id')
#                        print_debug(verified_node_id)
#                        id_mappings[verified_node_id] = node_id
#
#                    # if node is an RNA, we want to look for proteins,
#                    # which are targets (sucessors)
#                    if node_type == 'rna':
#                        neighbor_verified_node_ids = current_mapping_graph.successors(
#                            verified_node_id)
#                    # otherwise, node is a protein, so we'll look for RNAs,
#                    # which are sources (predecessors)
#                    else:
#                        neighbor_verified_node_ids = current_mapping_graph.predecessors(
#                            verified_node_id)
#
#                    current_source_log_result = log_result['results_by_source']
#                    if node_id not in current_source_log_result:
#                        current_source_log_result[node_id] = {}
#                    current_source_log_result[node_id][version] = len(
#                        neighbor_verified_node_ids)
#
#                    log_result['total_result_count'] += len(
#                        neighbor_verified_node_ids)
#
#                    log_result[version + '_result_count'] += len(
#                        neighbor_verified_node_ids)
#
#                    item_subgraph_node_ids = neighbor_verified_node_ids
#                    item_subgraph_node_ids.append(verified_node_id)
#                    datasource_subgraph_node_ids += item_subgraph_node_ids
#
#        if (len(id_mappings) > 0):
#            print_debug('id_mappings')
#            print_debug(id_mappings)
#            tmp_subgraph = current_mapping_graph.subgraph(
#                datasource_subgraph_node_ids)
#            datasource_subgraph = nx.relabel_nodes(
#                tmp_subgraph, id_mappings, copy=False)
#            print_debug('relabeled')
#        else:
#            datasource_subgraph = current_mapping_graph.subgraph(
#                datasource_subgraph_node_ids)
#
#        output_graph_name = version + '-subgraph.xgmml'
#        output_xgmml_file_path = output_dir + output_graph_name
#        output_xgmml = open(output_xgmml_file_path, 'w')
#        networkxgmml.XGMMLWriter(
#            output_xgmml, datasource_subgraph, output_graph_name)
#
#        with open(results_file_path, 'a') as csvfile:
#            csvwriter = csv.writer(csvfile, delimiter=',',
#                                   quotechar='|', quoting=csv.QUOTE_MINIMAL)
#            for oneedge in datasource_subgraph.edges(data=True):
#                print_debug('oneedge')
#                print_debug(oneedge)
#                edge_attributes = oneedge[2]
#                if oneedge[0] in node_id_list:
#                    row_value = [oneedge[0], oneedge[1]]
#                else:
#                    row_value = [oneedge[1], oneedge[0]]
#
#                if 'score' in edge_attributes:
#                    row_value.append(edge_attributes['score'])
#                else:
#                    row_value.append('')
#
#                if 'pValue' in edge_attributes:
#                    row_value.append(edge_attributes['pValue'])
#                else:
#                    row_value.append('')
#
#                if 'References (PMID)' in edge_attributes:
#                    row_value.append(edge_attributes['References (PMID)'])
#                else:
#                    row_value.append('')
#
#                if 'datasource' in edge_attributes:
#                    row_value.append(edge_attributes['datasource'])
#                else:
#                    row_value.append('')
#
#                print_debug('row_value')
#                print_debug(row_value)
#                csvwriter.writerow(row_value)
#
#    print_debug('log_result')
#    print_debug(log_result)
#    with file(logfile_path, 'w') as f:
#        print >>f, '# Queried {} line(s).'.format(
#            log_result['query_count'])
#        print >>f, '# Skipped {} line(s).'.format(
#            log_result['skipped_count'])
#        print >>f, '# Found {} results.'.format(
#            log_result['total_result_count'])
#        header = ['queryid'] + versions
#        print >>f, '\t'.join(header)
#        for one_node_id in log_result['results_by_source']:
#            row = [one_node_id]
#            for oneversion in versions:
#                logs_by_version = log_result['results_by_source'][one_node_id]
#                if oneversion in logs_by_version:
#                    row.append(str(logs_by_version[oneversion]))
#                else:
#                    row.append('0')
#            print >>f, '\t'.join(row)
#
#    return log_result
