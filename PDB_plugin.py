# Metadata for Plugin Manager
#
# Authors: PDBe, Artur.Klauser@computer.org
# License: BSD-2-Clause
# Version: 2.0.0a
# Citation-Required: No
"""
Described at PyMOL wiki:
https://pymolwiki.org/index.php/PDB_plugin

--- original plugin ---
Authors : PDBe
Date    : 2016

--- enhancements ---
Authors : Artur.Klauser@computer.org
Date    : 2019

DESCRIPTION

    Fetches the specified data set from PDB and breaks it into multiple
    objects based on the criteria named by the function.

USAGE

    PDB_Analysis_Assemblies pdb_id
    PDB_Analysis_Molecules pdb_id
    PDB_Analysis_Domains pdb_id
    PDB_Analysis_Validation pdb_id
    count_chains selection

ARGUMENTS

    pdb_id = string: 4-character PDB entry ID
    selection = string: selection-expression or name-pattern {default: (all)}.

EXAMPLES

    PDB_Analysis_Molecules 3mzw

NOTES

    For detailed descriptions see help on individual commands.
"""

from __future__ import print_function

# parsing the input
import json
import datetime
import socket
import re, os
import time
import sys
import random
import logging

import pymol
from pymol import cmd
from pymol import stored
from chempy.cif import *

# API calls
server_root = 'https://www.ebi.ac.uk/pdbe/api/'
summaryURL = server_root + 'pdb/entry/summary/'
poly_seq_schemeURL = server_root + 'pdb/entry/residue_listing/'
protein_mappingURL = server_root + 'mappings/'
nucleic_mappingURL = server_root + 'nucleic_mappings/'
# seqURL = server_root + 'mappings/sequence_domains/entry/'
moleculesURL = server_root + 'pdb/entry/molecules/'
validationURL = server_root + 'validation/global-percentiles/entry/'
validation_residueURL = server_root + 'validation/protein-RNA-DNA-geometry-outlier-residues/entry/'
validation_ramaURL = server_root + 'validation/protein-ramachandran-sidechain-outliers/entry/'

# ftp site
ebi_ftp = 'ftp://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/divided/mmCIF/%s/%s.cif.gz'
# clean mmcif
updated_ftp = 'https://www.ebi.ac.uk/pdbe/static/entry/%s_updated.cif.gz'

# dictionaries
stored.domain_dict = {}
stored.seq_scheme = {}
stored.molecule_dict = {}
stored.poly_count = 0
stored.residue_dict = {}
stored.ca_p_only = []
stored.entities = {}

# cmd.set('connect_mode', 4)

empty_cif_list = ['', '.', '?', None]


class ApiCall:
    """functions specific to API calls"""

    def __init__(self, url):
        self.url = url
        self.data = {}

    def return_data(self):
        try:
            import requests
        except Exception as e:
            logging.exception(e)
            return self.data
        self.encode_url()
        r = requests.get(url=self.url, timeout=60)
        if r.status_code == 200:
            json_data = r.json()
            self.data = json_data
        elif r.status_code == 404:
            self.data = {}

        else:
            logging.debug(r.status_code, r.reason)
            self.data = {}

        return self.data

    def encode_url(self):
        self.url = re.sub(' ', '%20', self.url)


# url handling
def url_response_urllib2(url, description):
    """try several times to get a response from the API"""
    data = None
    try:
        import urllib2
    except Exception as e:
        logging.exception(e)
        return data
    data_response = False
    tries = 1
    limit = 5
    logging.debug(url)
    date = datetime.datetime.now().strftime('%Y,%m,%d  %H:%M')
    while tries < limit:
        try:
            response = urllib2.urlopen(url, None, 60)
        except urllib2.HTTPError as e:
            logging.debug('##############')
            logging.debug('%s HTTP API error - %s, error code - %s, try %s' % (date, description, e.code, str(tries)))
            # logging.debug(e.code)
            if e.code == 404:
                data_response = True
                data = {}
                break
            else:
                # logging.debug(entry_url)
                time.sleep(random.randint(1, 20))  # sleep 1-20 seconds
                tries += 1
                logging.debug(e)
                data = {}
        except urllib2.URLError as e:
            logging.debug('##############')
            logging.debug('%s URL API error - %s, try %s' % (date, e, str(tries)))
            # logging.debug(entry_url)
            time.sleep(random.randint(1, 20))  # sleep 1-20 seconds
            tries += 1
            logging.debug(e)
            data = {}
        except socket.timeout as e:
            logging.debug('##############')
            logging.debug('%s Timeout API error - %s, try %s' % (date, e, str(tries)))
            # logging.debug(entry_url)
            time.sleep(random.randint(1, 20))  # sleep 1-20 seconds
            tries += 1
            logging.debug(e)
            data = {}
        else:
            logging.debug('received a response from the ' + description + ' API after ' + str(tries) + ' tries.')
            data = json.load(response)
            data_response = True
            break

    if not data_response:
        logging.debug('No response from the %s API' % description)
        # cmd.quit(1)

    return data


# url handling
def url_response(url, description):
    logging.debug(description)
    try:
        import requests
        data = ApiCall(url=url).return_data()
    except:
        data = url_response_urllib2(url=url, description=description)
    return data


class Color(object):
    """Manage use of colors."""

    _object_colors = [
        'red', 'green', 'blue', 'yellow', 'magenta', 'cyan', 'tv_red',
        'tv_green', 'tv_blue', 'tv_yellow', 'lightmagenta', 'palecyan',
        'raspberry', 'chartreuse', 'marine', 'paleyellow', 'hotpink',
        'aquamarine', 'darksalmon', 'splitpea', 'slate', 'yelloworange', 'pink',
        'greencyan', 'salmon', 'smudge', 'lightblue', 'limon', 'lightpink',
        'teal', 'deepsalmon', 'palegreen', 'skyblue', 'wheat', 'dirtyviolet',
        'deepteal', 'warmpink', 'limegreen', 'purpleblue', 'sand', 'violet',
        'lightteal', 'firebrick', 'lime', 'deepblue', 'violetpurple', 'ruby',
        'limon', 'density', 'purple', 'chocolate', 'forest', 'deeppurple',
        'brown',
    ]
    _validation_colors = ['yellow', 'orange', 'red']

    @classmethod
    def set_object_color(cls, color_num, obj):
        """Colors object from predetermined set, reusing as necessary."""
        color = cls._object_colors[color_num % len(cls._object_colors)]
        cmd.color(color, obj)

    @classmethod
    def set_validation_color(cls, color_num, obj):
        """Colors object from predetermined set, reusing as necessary."""
        color = cls._validation_colors[color_num % len(cls._validation_colors)]
        cmd.color(color, obj)

    @classmethod
    def set_validation_background_color(cls, obj):
        """Colors object with validation background color."""
        color = 'validation_background_color'
        # Define this as a new color name.
        cmd.set_color(color, [0.4, 1.0, 0.4])
        cmd.color(color, obj)


def poly_display_type(asym, mol_type, length):
    """return either ribbon or cartoon depending on number of polymer chains"""
    ##start with it being cartoon - then change as needed.
    display_type = 'cartoon'
    if stored.poly_count > 50:
        display_type = 'ribbon'
    elif asym in stored.ca_p_only:
        logging.debug('set ribbon trace on')
        cmd.set('ribbon_trace_atoms', 1)
        display_type = 'ribbon'
    else:
        if re.search('polypeptide', mol_type):
            if length < 20:
                display_type = 'sticks'
        if re.search('nucleotide', mol_type):
            if length < 3:
                display_type = 'sticks'
        if re.search('saccharide', mol_type):
            display_type = 'sticks'

    # logging.debug('asym: %s, mol_type: %s, length: %s, display_type: %s' %(asym, mol_type, length, display_type))
    return display_type


class worker_functions():
    def trans_setting(self, selection, trans_value):
        cmd.set('cartoon_transparency', trans_value, selection)
        cmd.set('ribbon_transparency', trans_value, selection)
        cmd.set('stick_transparency', trans_value, selection)
        cmd.set('sphere_transparency', trans_value, selection)

    def count_poly(self, pdbid):
        url = moleculesURL + pdbid
        # global molecule_dict
        if stored.molecule_dict == {}:
            stored.molecule_dict = url_response(url, 'molecules')
        if pdbid in stored.molecule_dict:
            for entity in stored.molecule_dict[pdbid]:
                # add ca only list
                if entity['ca_p_only'] != False:
                    for a in entity['in_struct_asyms']:
                        if a not in stored.ca_p_only:
                            stored.ca_p_only.append(a)
                if entity['molecule_type'] not in ['water', 'bound']:
                    stored.poly_count += 1


def poly_seq_scheme(pdbid):
    """build a dictionary like poly_seq_scheme"""
    url = poly_seq_schemeURL + pdbid
    data = url_response(url, 'seq_scheme')
    if pdbid in data:
        for entity in data[pdbid]['molecules']:
            for chain in entity['chains']:
                chain_id = chain['chain_id']
                asym_id = chain['struct_asym_id']
                for residue in chain['residues']:
                    CIFnum = residue['residue_number']
                    PDBnum = residue['author_residue_number']
                    PDBinsCode = residue['author_insertion_code']
                    residueName = residue['residue_name']
                    if residue['observed_ratio'] != 0:
                        PDBobs = True
                    else:
                        PDBobs = False

                    try:
                        stored.seq_scheme[asym_id].update({CIFnum: {'PDBnum': PDBnum, 'PDBinsCode': PDBinsCode,
                                                                    'observed': PDBobs, 'residueName': residueName,
                                                                    'chainID': chain_id}})
                    except:
                        stored.seq_scheme.update({asym_id: {}})
                        stored.seq_scheme[asym_id].update({CIFnum: {'PDBnum': PDBnum, 'PDBinsCode': PDBinsCode,
                                                                    'observed': PDBobs, 'residueName': residueName,
                                                                    'chainID': chain_id}})
                        # logging.debug(seq_scheme)


def check_range_observed(asym_id, start, end, unobs):
    range_dict = []
    firstResidue = 1
    lastResidue = 1
    for r in stored.seq_scheme[asym_id]:
        if r <= firstResidue:
            firstResidue = r
        if r >= lastResidue:
            lastResidue = r

    if end > lastResidue:
        # logging.debug('ERROR: residue %s from SIFTS is greater than the lastSeqRes number %s' %(end, lastSeqRes))
        end = lastResidue

    # logging.debug('first residue: %s, last residue %s' %(firstResidue, lastResidue))
    # logging.debug(seq_scheme[asym_id][start])
    if start in stored.seq_scheme[asym_id]:
        chain = stored.seq_scheme[asym_id][start]['chainID']

        while stored.seq_scheme[asym_id][start]['observed'] == False:
            # logging.debug('PDB start is not observed, adding 1 to residue number %s' %(start))
            if start == lastResidue:
                # logging.debug('domain not observed')
                unobs = True
                break
            else:
                start += 1
                if start == lastResidue:
                    # logging.debug('domain not observed')
                    unobs = True
                    break
                elif start == end:
                    # logging.debug('domain not observed')
                    unobs = True
                    break

    if end in stored.seq_scheme[asym_id]:
        while stored.seq_scheme[asym_id][end]['observed'] == False:
            # logging.debug('PDB end is not observed, minusing 1 from residue number %s' %(end))
            if end == firstResidue:
                # logging.debug('domain not observed')
                unobs = True
                break
            else:
                end -= 1
                if end == firstResidue:
                    # logging.debug('domain not observed')
                    unobs = True
                    break
                elif start == end:
                    # logging.debug('domain not observed')
                    unobs = True
                    break

    if unobs == False:
        # logging.debug('CIF start: %s, CIF end %s' %(start, end))
        # logging.debug('PDB start: %s, PDB end %s' %(PDBstart, PDBend))

        # need to do something about non continuous ranges which pymol doesn't cope with.
        """find all sections where the residue number increases by one. If it jumps then store this as a separate
        residue block while number is less than endbegin with number = start plus 1. see how the numbering differs
        from startif its 1 then its continuous - move to the next residue
        if its zero then there must be insert codes - this is ok - move to the next residue
        if its greater than 1 then store the previous residues as a block.
        """
        currentCIF = start
        blockStartCIF = start
        nextCIF = start + 1
        while nextCIF < end:
            currentPDB = stored.seq_scheme[asym_id][currentCIF]['PDBnum']
            nextPDB = stored.seq_scheme[asym_id][nextCIF]['PDBnum']
            if currentPDB not in empty_cif_list:
                if nextPDB in empty_cif_list:
                    startPDB = stored.seq_scheme[asym_id][blockStartCIF]['PDBnum']
                    swap = checkOrder(startPDB, currentPDB)
                    startPDBins = insert_code(asym_id, blockStartCIF)
                    currentPDB = insert_code(asym_id, currentCIF)
                    if swap == False:
                        range_dict.append({'PDBstart': startPDBins, 'PDBend': currentPDB, 'chainID': chain})
                    else:
                        range_dict.append({'PDBstart': currentPDB, 'PDBend': startPDBins, 'chainID': chain})
                    # move start to next block
                    blockStartCIF = nextCIF

                else:
                    currentPDB = int(stored.seq_scheme[asym_id][currentCIF]['PDBnum'])
                    nextPDB = int(stored.seq_scheme[asym_id][nextCIF]['PDBnum'])

                    if nextPDB - currentPDB > 1:
                        # logging.debug('numbering not continues, positive jump - store as section')
                        startPDB = stored.seq_scheme[asym_id][blockStartCIF]['PDBnum']
                        swap = checkOrder(startPDB, currentPDB)
                        startPDBins = insert_code(asym_id, blockStartCIF)
                        currentPDB = insert_code(asym_id, currentCIF)
                        if swap == False:
                            range_dict.append({'PDBstart': startPDBins, 'PDBend': currentPDB, 'chainID': chain})
                        else:
                            range_dict.append({'PDBstart': currentPDB, 'PDBend': startPDBins, 'chainID': chain})
                            # move start to next block
                        blockStartCIF = nextCIF
                    elif currentPDB - nextPDB > 1:
                        # logging.debug('numbering not continues, negative jump - store as section')
                        startPDB = stored.seq_scheme[asym_id][blockStartCIF]['PDBnum']
                        swap = checkOrder(startPDB, currentPDB)
                        startPDBins = insert_code(asym_id, blockStartCIF)
                        currentPDB = insert_code(asym_id, currentCIF)
                        if swap == False:
                            range_dict.append({'PDBstart': startPDBins, 'PDBend': currentPDB, 'chainID': chain})
                        else:
                            range_dict.append({'PDBstart': currentPDB, 'PDBend': startPDBins, 'chainID': chain})
                        # move start to next block
                        blockStartCIF = nextCIF
            else:
                blockStartCIF = nextCIF
            nextCIF += 1
            currentCIF += 1

        startPDB = stored.seq_scheme[asym_id][blockStartCIF]['PDBnum']
        swap = checkOrder(blockStartCIF, nextCIF)
        startPDBins = insert_code(asym_id, blockStartCIF)
        nextPDB = insert_code(asym_id, nextCIF)
        if swap == False:
            range_dict.append({'PDBstart': startPDBins, 'PDBend': nextPDB, 'chainID': chain})
        else:
            range_dict.append({'PDBstart': nextPDB, 'PDBend': startPDBins, 'chainID': chain})

            # logging.debug(range_dict)
    else:
        logging.debug('domain unobserved')

    return range_dict


def checkOrder(start, end):
    swap = False
    # logging.debug('PRE: start: %s, end %s, swap %s' %(start, end, swap))
    if type(start) is str:
        start = start.split('\\')[-1]
    if type(end) is str:
        end = end.split('\\')[-1]

    if int(start) > int(end):
        swap = True
    # logging.debug('POST: start: %s, end: %s, swap %s' %(start, end, swap))
    return swap


def insert_code(asym_id, CIFnum):
    if stored.seq_scheme[asym_id][CIFnum]['PDBinsCode'] == None:
        PDBnum = stored.seq_scheme[asym_id][CIFnum]['PDBnum']
    else:
        PDBnum = '%s%s' % (
            stored.seq_scheme[asym_id][CIFnum]['PDBnum'], stored.seq_scheme[asym_id][CIFnum]['PDBinsCode'])

    if re.search('-', str(PDBnum)):
        PDBnum = re.sub('-', '\-', str(PDBnum))

    return PDBnum



class Validation:
    def launch_validation(self, pdbid):

        url = validationURL + pdbid
        val_data = url_response(url, 'validation')

        if val_data:
            logging.debug('There is validation for this entry')

            url = validation_residueURL + pdbid
            res_data = url_response(url, 'residue validation')

            url = validation_ramaURL + pdbid
            rama_data = url_response(url, 'rama validation')

            self.per_residue_validation(pdbid, res_data, rama_data)
        else:
            logging.debug('No validation for this entry')

            ###this takes too long for really large entries.
            # Validation().per_chain_per_residue_validation(pdbid, res_data, rama_data)

    """validation of all polymeric entries"""

    def validation_selection(self, selection, display_id):
        # logging.debug(selection)

        # uses a list so 0 means that there is one outlier for this residue.

        if selection in stored.residue_dict:
            color_num = stored.residue_dict[selection] + 1
            stored.residue_dict[selection] = color_num
        else:
            stored.residue_dict.update({selection: 0})
            color_num = 0

        if color_num > 2:
            color_num = 2
        Color.set_validation_color(color_num, selection)

    def geometric_validation(self, pdbid, res_data, chain=False, model=1):
        """check for geometric validation outliers in res_data """

        if res_data:
            if pdbid in res_data:
                if 'molecules' in res_data[pdbid]:
                    if res_data[pdbid]['molecules']:
                        for entity in res_data[pdbid]['molecules']:
                            # logging.debug(entity)
                            entity_id = entity['entity_id']
                            for chain in entity['chains']:
                                chain_id = chain['chain_id']
                                # logging.debug(chain_id)
                                for model in chain['models']:
                                    model_id = int(model['model_id'])
                                    for outlier_type in model['outlier_types']:
                                        outliers = model['outlier_types'][outlier_type]
                                        # logging.debug(outlier_type)
                                        # logging.debug(outliers)
                                        for outlier in outliers:
                                            PDB_res_num = outlier['author_residue_number']
                                            PDB_ins_code = outlier['author_insertion_code']
                                            alt_code = outlier['alt_code']
                                            if PDB_ins_code not in [None, ' ']:
                                                PDB_res_num = '%s%s' % (PDB_res_num, PDB_ins_code)

                                            # logging.debug(PDB_res_num)
                                            selection = 'chain %s and resi %s' % (chain_id, PDB_res_num)
                                            if model == 1 or model_id:
                                                if chain == False or chain_id:
                                                    self.validation_selection(selection, pdbid)

                    else:
                        logging.debug('no residue validation for this entry')

    def ramachandran_validation(self, pdbid, rama_data, chain=False, model=1):
        """display ramachandran outliers"""

        if rama_data:
            if pdbid in rama_data:
                for k in rama_data[pdbid]:
                    v = rama_data[pdbid][k]
                    if v:
                        logging.debug('ramachandran %s for this entry' % k)
                        for outlier in rama_data[pdbid][k]:
                            # logging.debug(outlier)
                            entity_id = outlier['entity_id']
                            model_id = int(outlier['model_id'])
                            chain_id = outlier['chain_id']
                            PDB_res_num = outlier['author_residue_number']
                            PDB_ins_code = outlier['author_insertion_code']
                            alt_code = outlier['alt_code']

                            if PDB_ins_code not in [None, ' ']:
                                PDB_res_num = '%s%s' % (PDB_res_num, PDB_ins_code)

                            selection = 'chain %s and resi %s' % (chain_id, PDB_res_num)
                            if model == 1 or model_id:
                                if chain == False or chain_id:
                                    self.validation_selection(selection, pdbid)
                            else:
                                logging.debug('is multimodel!!!, outlier not in model 1, not shown.')

                    else:
                        logging.debug('no %s' % (k))

    def per_residue_validation(self, pdbid, res_data, rama_data):
        """validation of all outliers, colored by number of outliers"""

        # display only polymers
        if stored.molecule_dict == {}:
            url = moleculesURL + pdbid
            stored.molecule_dict = url_response(url, 'molecules')
        if pdbid in stored.molecule_dict:
            for entity in stored.molecule_dict[pdbid]:
                # logging.debug(entity)
                mol_type = entity['molecule_type']
                ca_only_list = []
                if entity['ca_p_only'] != False:
                    ca_only_list = entity['ca_p_only']
                if mol_type not in ['Water', 'Bound']:
                    for a in entity['in_struct_asyms']:
                        if stored.seq_scheme == {}:
                            poly_seq_scheme(pdbid)
                        unobs = False
                        start = 1
                        end = entity['length']
                        length = entity['length']
                        asym_id = a
                        display_type = poly_display_type(asym_id, 'polypeptide', length)
                        # logging.debug(asym_id)
                        x = check_range_observed(asym_id, start, end, unobs)
                        if x:
                            for y in x:
                                start = y['PDBstart']
                                end = y['PDBend']
                                chain = y['chainID']
                                selection = 'chain %s and resi %s-%s and %s' % (chain, start, end, pdbid)
                                # logging.debug(selection)
                                cmd.show(display_type, selection)

        Color.set_validation_background_color(pdbid)
        # display(pdbid, image_type)
        cmd.enable(pdbid)

        self.geometric_validation(pdbid, res_data)
        self.ramachandran_validation(pdbid, rama_data)

        # logging.debug(stored.residue_dict)


def show_entities(pdbid):
    logging.debug('Display entities')
    cmd.set('cartoon_transparency', 0.3, pdbid)
    cmd.set('ribbon_transparency', 0.3, pdbid)
    url = moleculesURL + pdbid
    num = 1
    if stored.molecule_dict == {}:
        stored.molecule_dict = url_response(url, 'molecules')
    if pdbid in stored.molecule_dict:
        for entity in stored.molecule_dict[pdbid]:
            entity_name = ''
            if entity['molecule_name']:
                if len(entity['molecule_name']) > 1:
                    for instance, mol in enumerate(entity['molecule_name']):
                        if mol != entity['molecule_name'][-1]:
                            mol = mol + '-'
                            entity_name += mol
                        else:
                            entity_name += mol
                    entity_name = entity_name + '_chimera'
                else:
                    entity_name = entity['molecule_name'][0]
            else:
                logging.debug('No name from the API')
                entity_name = 'entity_%s' % (entity['entity_id'])
            # logging.debug('entity %s: %s' %(entity['entity_id'], entity_name))

            # logging.debug(entity['entity_id'])
            ####see if its polymeric
            mol_type = entity['molecule_type']
            # entity_name = entity['molecule_name']
            if entity_name:
                objectName = re.sub(' ', '_', entity_name)
                objectName = re.sub(r'\W+', '', objectName)
            else:
                objectName = 'entity_%s' % (entity['entity_id'])
            # logging.debug(objectName)
            display_type = ''
            object_selection = []
            pymol_selection = ''
            if mol_type != 'Water':
                # use asym to find residue number and chain ID
                for a in entity['in_struct_asyms']:
                    if stored.seq_scheme == {}:
                        poly_seq_scheme(pdbid)
                    if mol_type == 'Bound':
                        ##bound molecules have no CIFresidue number and this is defaulted to 1 in the API
                        # logging.debug(entity_name)
                        # logging.debug(stored.seq_scheme[a])
                        for res in stored.seq_scheme[a]:
                            # logging.debug(stored.seq_scheme[a][res])
                            short = stored.seq_scheme[a][res]
                            chain = short['chainID']
                            res = ''
                            display_type = 'spheres'
                            # logging.debug(short)
                            if short['PDBinsCode'] == None:
                                res = str(short['PDBnum'])
                            else:
                                res = str(short['PDBnum']) + short['PDBinsCode']
                            selection = 'chain %s and resi %s and %s' % (chain, res, pdbid)
                            object_selection.append(selection)
                            # logging.debug(selection)
                    else:
                        ###find the first and last residues
                        unobs = False
                        start = 1
                        end = entity['length']
                        asym_id = a
                        ###need to work out if cartoon is the right thing to display
                        length = entity['length']
                        display_type = poly_display_type(asym_id, mol_type, length)
                        # logging.debug(asym_id)
                        x = check_range_observed(asym_id, start, end, unobs)
                        # logging.debug(x)
                        if x:
                            for y in x:
                                start = y['PDBstart']
                                end = y['PDBend']
                                chain = y['chainID']
                                selection = 'chain %s and resi %s-%s and %s' % (chain, start, end, pdbid)
                                # logging.debug(selection)
                                object_selection.append(selection)

                for o in object_selection:
                    if len(object_selection) > 1:
                        if o == object_selection[-1]:
                            pymol_selection += '(%s)' % (o)
                        else:
                            pymol_selection += '(%s) or ' % (o)
                    else:
                        pymol_selection += '%s' % (o)

                try:
                    stored.entities[pdbid]['entity'][entity['entity_id']]['selection']['selection_list'].append(
                        pymol_selection)
                    logging.debug('appended %s' % (pymol_selection))
                except:
                    stored.entities.setdefault(pdbid, {}) \
                        .setdefault('entity', {}) \
                        .setdefault(entity['entity_id'], {}) \
                        .setdefault('selection', {}) \
                        .update({'display_type': display_type, 'selection_list': [pymol_selection]})

                # logging.debug(pymol_selection)
                if len(objectName) > 250:
                    objectName = objectName[:249]
                cmd.select('test_select', pymol_selection)
                cmd.create(objectName, 'test_select')
                # logging.debug(display_type)
                cmd.show(display_type, objectName)

                # color by entity.
                Color.set_object_color(int(entity['entity_id']), objectName)

    cmd.delete('test_select')
    # cmd.show('cartoon', pdbid)


def show_assemblies(pdbid, mmCIF_file):
    """iterate through the assemblies and output images"""
    logging.info('Generating assemblies')
    # cmd.hide('everything')
    worker_functions().trans_setting('all', 0.0)

    try:
        assemblies = cmd.get_assembly_ids(pdbid)  # list or None
        logging.debug(assemblies)
        if assemblies:
            for a_id in assemblies:
                logging.debug('Assembly: %s' % (a_id))
                assembly_name = pdbid + '_assem_' + a_id
                logging.debug(assembly_name)
                cmd.set('assembly', a_id)
                cmd.load(mmCIF_file, assembly_name, format='cif')
                logging.debug('finished Assembly: %s' % (a_id))
    except:
        logging.debug('pymol version does not support assemblies')


def mapping(pdbid):
    """make domain objects"""
    if not stored.seq_scheme:
        poly_seq_scheme(pdbid)
    domain_to_make = ['CATH', 'SCOP', 'Pfam', 'Rfam']
    segment_identifier = {'CATH': 'domain', 'SCOP': 'scop_id', 'Pfam': '', 'Rfam': ''}
    protein_url = protein_mappingURL + pdbid
    nucleic_url = nucleic_mappingURL + pdbid
    protein_domains = url_response(protein_url, 'protein_domains')
    nucleic_domains = url_response(nucleic_url, 'nucleic_domains')
    for data in [protein_domains, nucleic_domains]:
        if data == {}:
            logging.debug('no domain information for this entry')
        else:
            for domain_type in data[pdbid]:
                if domain_type in domain_to_make:
                    for domain in data[pdbid][domain_type]:
                        if data[pdbid][domain_type][domain]['mappings']:
                            identifier = data[pdbid][domain_type][domain]['identifier']
                            # logging.debug(domain_type)
                            # logging.debug(domain)
                            for mapping in data[pdbid][domain_type][domain]['mappings']:
                                unobs = False
                                PDBstart = ''
                                PDBend = ''
                                domain_segment_id = ''
                                domain_name = ''
                                # check segment id
                                if 'segment_id' in mapping:
                                    domain_segment_id = mapping['segment_id']
                                else:
                                    # force it to be one if its missing
                                    domain_segment_id = 1
                                if segment_identifier[domain_type] in mapping:
                                    domain_name = str(mapping[segment_identifier[domain_type]])
                                # logging.debug('%s: %s' %(domain_name, domain_segment_id))

                                # then check it exits
                                start = mapping['start']['residue_number']
                                end = mapping['end']['residue_number']
                                chain = mapping['chain_id']
                                entity_id = mapping['entity_id']
                                asym_id = mapping['struct_asym_id']
                                if asym_id in stored.seq_scheme:
                                    x = check_range_observed(asym_id, start, end, unobs)
                                    if x:
                                        for y in x:
                                            start = y['PDBstart']
                                            end = y['PDBend']
                                            try:
                                                stored.domain_dict[domain_type][domain][domain_name].append(
                                                    {'asym_id': asym_id, 'chain': chain, 'start': start, 'end': end,
                                                     'entity_id': entity_id})
                                            except:
                                                stored.domain_dict.setdefault(domain_type, {}) \
                                                    .setdefault(domain, {}) \
                                                    .update({domain_name: [
                                                    {'asym_id': asym_id, 'chain': chain, 'start': start, 'end': end,
                                                     'entity_id': entity_id}]})


def show_domains(pdbid):
    mapping(pdbid)
    obj_dict = {}
    chain_dict = {}
    if stored.domain_dict:
        # logging.debug(domain_dict)
        for domain_type in stored.domain_dict:
            for domain in stored.domain_dict[domain_type]:
                # logging.debug(domain)
                for instance in stored.domain_dict[domain_type][domain]:
                    asym_list = []
                    entity_list = []
                    # logging.debug(instance)
                    object_selection = []
                    pymol_selection = ''
                    for segment in stored.domain_dict[domain_type][domain][instance]:
                        PDBstart = segment['start']
                        PDBend = segment['end']
                        chain = segment['chain']
                        asym_id = segment['asym_id']
                        entity_id = segment['entity_id']

                        asym_list.append(asym_id)
                        entity_list.append(entity_id)
                        chain_dict.update({asym_id: {'chain': chain, 'entity_id': entity_id}})

                        selection = 'chain %s and resi %s-%s and %s' % (chain, PDBstart, PDBend, pdbid)
                        object_selection.append(selection)
                    for o in object_selection:
                        if len(object_selection) > 1:
                            if o == object_selection[-1]:
                                pymol_selection += '(%s)' % (o)
                            else:
                                pymol_selection += '(%s) or ' % (o)
                        else:
                            pymol_selection += '%s' % (o)
                    # logging.debug(pymol_selection)
                    objectName = '%s_%s_%s' % (domain_type, domain, instance)
                    obj_dict.update({objectName: {'asym_list': asym_list, 'entity_list': entity_list}})
                    cmd.select('test_select', pymol_selection)
                    cmd.create(objectName, 'test_select')

            for chain in chain_dict:
                logging.debug(chain)
                c_select = 'chain %s and %s' % (chain_dict[chain]['chain'], pdbid)

                entity_id = chain_dict[chain]['entity_id']
                length = None
                for entity in stored.molecule_dict[pdbid]:
                    if entity_id == entity['entity_id']:
                        length = entity['length']
                display_type = poly_display_type(chain, 'polypeptide', length)
                cmd.show(display_type, c_select)
                cmd.color('grey', c_select)
            num = 1
            # logging.debug(obj_dict)
            for obj in obj_dict:
                # logging.debug(obj)
                cmd.enable(obj)
                # domain can span multiple asym_id's, could be different type
                for i, a in enumerate(obj_dict[obj]['asym_list']):
                    entity_id = obj_dict[obj]['entity_list'][i]
                    length = None
                    for entity in stored.molecule_dict[pdbid]:
                        if entity_id == entity['entity_id']:
                            length = entity['length']
                    display_type = poly_display_type(a, 'polypeptide', length)
                    cmd.show(display_type, obj)
                    Color.set_object_color(num, obj)
                    num += 1

            cmd.delete('test_select')


def PDBe_startup(pdbid, method, mmCIF_file=None, file_path=None):
    if pdbid:
        pdbid = pdbid.lower()
    if mmCIF_file:
        filename = os.path.basename(mmCIF_file)
        if re.search('_', filename):
            pdbid = re.split('_', filename)[0].lower()
        else:
            pdbid = re.split('\.', filename)[0].lower()

    try:
        cmd.set('cif_keepinmemory')
    except:
        logging.exception('version of pymol does not support keeping all cif items')

    # check the PDB code actually exists.
    url = summaryURL + pdbid
    summary = url_response(url, 'summary')

    if summary:

        # clear the dictionaries
        stored.domain_dict = {}
        stored.seq_scheme = {}
        stored.molecule_dict = {}
        stored.residue_dict = {}

        logging.debug('pdbid: %s' % (pdbid))
        mid_pdb = pdbid[1:3]

        obj_list = cmd.get_object_list('all')
        if pdbid not in obj_list:
            if mmCIF_file:
                file_path = mmCIF_file
            else:
                # if local file exists then use it
                if os.path.exists('%s.cif' % (pdbid)):
                    logging.debug('CIF found in current directory')
                    file_path = '%s.cif' % (pdbid)
                else:
                    logging.debug('fetching %s' % (pdbid))
                    try:
                        # connect mode 4 works only with version 1.7 of pymol
                        cmd.set('assembly', '')
                        file_path = updated_ftp % pdbid
                        logging.debug('setting connect mode to mode 4')
                        cmd.set('connect_mode', 4)
                    except:
                        logging.exception('pymol version does not support assemblies or connect mode 4')
                        file_path = ebi_ftp % (mid_pdb, pdbid)

            logging.debug('File to load: %s' % file_path)
            cmd.load(file_path, pdbid, format='cif')
        cmd.hide('everything', pdbid)
        worker_functions().count_poly(pdbid)

        if method == 'molecules':
            show_entities(pdbid)
        elif method == 'domains':
            show_entities(pdbid)
            show_domains(pdbid)
        elif method == 'validation':
            Validation().launch_validation(pdbid)
        elif method == 'assemblies':
            show_assemblies(pdbid, file_path)
        elif method == 'all':
            show_entities(pdbid)
            show_domains(pdbid)
            show_assemblies(pdbid, file_path)
            Validation().launch_validation(pdbid)
        else:
            logging.warning('provide a method')
        cmd.zoom(pdbid, complete=1)

    elif mmCIF_file:
        logging.warning('no PDB ID, show assemblies from mmCIF file')
        cmd.load(mmCIF_file, pdbid, format='cif')
        # worker_functions().count_poly(pdbid)
        show_assemblies(pdbid, mmCIF_file)

    else:
        logging.error('please provide a 4 letter PDB code')


def GetPdbId(label):
    from pymol.Qt import QtWidgets

    pdbid, ok_pressed = QtWidgets.QInputDialog.getText(
        None, 'PDB Loader Service',
        label + '\nPlease enter a 4-digit PDB ID:',
        QtWidgets.QLineEdit.Normal, '')

    if ok_pressed:
        return pdbid
    else:
        return None


def PDBeLoaderDialog():
    pdbid = GetPdbId(
        'Highlight chemically distinct molecules, domains and assemblies in a'
        ' PDB entry.')
    if pdbid: PDBe_startup(pdbid, 'all')


def PDBeMoleculeDialog():
    pdbid = GetPdbId('Highlight chemically distinct molecules in a PDB entry.')
    if pdbid: PDBe_startup(pdbid, 'molecules')


def PDBeDomainDialog():
    pdbid = GetPdbId(
        'Display Pfam, SCOP, CATH and Rfam domains on a PDB entry.')
    if pdbid: PDBe_startup(pdbid, 'domains')


def PDBeValidationDialog():
    pdbid = GetPdbId('Display geometric outliers on a PDB entry.')
    if pdbid: PDBe_startup(pdbid, 'validation')


def PDBeAssemblyDialog():
    pdbid = GetPdbId('Display assemblies for a PDB entry.')
    if pdbid: PDBe_startup(pdbid, 'assemblies')


class PdbIdShortcut(cmd.Shortcut):
    """
    PDB entry ID autocomplete helper: we don't have the whole set of valid keys
    in the PDB, but we try to guide the user to filling in a valid key and
    verify full-length keys with the PDB, showing the entry title as guidance.
    """
    def interpret(self, key, mode=0):
        # The filler is a 'space' but ordered after visible ASCII chars. We use
        # it to convince autocomplete that there is still more than 1 option
        # available, even though we only want to show one, which is just a sort
        # of 'regular expression' to guide the user.
        filler = '\x7f'
        if len(key) == 0:
            # Returning a list indicates that there are still >1 options.
            return ['[0-9]' + '[alphanumeric]' * 3, filler]
        elif len(key) < 4:
            pdb_id_pattern = re.compile(r'\d[a-zA-Z0-9]{0,3}$')
            if re.match(pdb_id_pattern, key) is None:
                return None  # key doesn't match PDB ID format
            # Returning a list indicates that there are still >1 options.
            return [key + '[alphanumeric]' * (4 - len(key)), filler]
        else:
            url = summaryURL + key
            summary = url_response(url, 'summary')
            # If the key doesn't return a valid description return None to
            # indicate that the key is invalid.
            if not summary:
                return None
            key_list = summary.get(key.lower(), None)
            if not key_list:
                return None
            if len(key_list) != 1:
                return None
            title = key_list[0].get('title', None)
            if not title:
                return None
            # We've got a valid key. Show the title of that PDB entry to the
            # user and make the key the only possible autocomplete match.
            print(key + ':', title)
            return key


@cmd.extendaa([lambda: PdbIdShortcut(), 'PDB Entry Id', ''])
def PDB_Analysis_Molecules(pdbid):
    """
DESCRIPTION

    Fetches the specified entry from PDB and highlights the chemically
    distinct molecules, including proteins, nucleic acids and ligands. Each
    molecule is given a unique color and selected as a different object.
    Polymers (protein, DNA and RNA chains) are shown as cartoon or ribbon
    representation and ligands are shown as spheres.

USAGE

    PDB_Analysis_Molecules pdb_id

ARGUMENTS

    pdb_id = string: 4-character PDB entry ID

EXAMPLES

    PDB_Analysis_Molecules 3l2p
    """
    PDBe_startup(pdbid, 'molecules')


@cmd.extendaa([lambda: PdbIdShortcut(), 'PDB Entry Id', ''])
def PDB_Analysis_Domains(pdbid):
    """
DESCRIPTION

    Fetches the specified entry from PDB and highlights the Pfam, Rfam, SCOP
    and CATH domains. Each domain is colored in a different color and selected
    as a different object. This enables each domain to be turned on and off with
    PyMOL. The Domains Plugin also highlights chemically distinct molecules and
    domains are overlaid on these molecules.

USAGE

    PDB_Analysis_Domains pdb_id

ARGUMENTS

    pdb_id = string: 4-character PDB entry ID

EXAMPLES

    PDB_Analysis_Domains 3b43
    """
    PDBe_startup(pdbid, 'domains')


@cmd.extendaa([lambda: PdbIdShortcut(), 'PDB Entry Id', ''])
def PDB_Analysis_Validation(pdbid):
    """
    '''
DESCRIPTION

    Fetches the specified entry from PDB and overlays geometric validation.
    Geometry outliers are colored as per the PDB validation report:
    Residues with no geometric outliers are shown in green. Residues with one
    geometric outliers are shown in yellow. Residues with two geometric
    outliers are shown in orange. Residues with three or more geometric
    outliers are shown in red.

USAGE

    PDB_Analysis_Validation pdb_id

ARGUMENTS

    pdb_id = string: 4-character PDB entry ID

EXAMPLES

    PDB_Analysis_Validation 2gc2
    """
    PDBe_startup(pdbid, 'validation')


@cmd.extendaa([lambda: PdbIdShortcut(), 'PDB Entry Id', ''])
def PDB_Analysis_Assemblies(pdbid):
    """
DESCRIPTION

    Fetches the specified entry from PDB and highlights the assemblies.

USAGE

    PDB_Analysis_Assemblies pdb_id

ARGUMENTS

    pdb_id = string: 4-character PDB entry ID

EXAMPLES

    PDB_Analysis_Assemblies 5j96
    """
    PDBe_startup(pdbid, 'assemblies')


@cmd.extendaa(cmd.auto_arg[0]['count_atoms'])
def count_chains(selection='(all)', state=0):
    """
DESCRIPTION

    "count_chains" returns the number of chains in the selection.

USAGE

    count_chains [ selection [, state ]]


ARGUMENTS

    selection = string: selection-expression or name-pattern {default: (all)}.
    state = int: object state, -1 for current state, 0 for all states
    {default: 0 (all states)}

EXAMPLES

    count_chains visible

SEE ALSO

    get_chains
    """
    n = len(cmd.get_chains(selection, state))
    print('count_chains: %s chains' % n)
    return n


def Initialize():
    # get preferences
    pref_loglevel = 'PDB_PLUGIN_LOGLEVEL'
    loglevel = pymol.plugins.pref_get(pref_loglevel, None)
    if loglevel is None:
        # supported log levels: DEBUG, INFO, WARNING, ERROR, CRITICAL
        loglevel = 'WARNING'
        pymol.plugins.pref_set(pref_loglevel, loglevel)
        pymol.plugins.pref_save()

    # set up logging of 'root' logger
    logger = logging.getLogger()
    numeric_loglevel = getattr(logging, loglevel.upper(), None)
    if numeric_loglevel is None:
      logger.setLevel(logging.WARNING)
      logging.error('Preference ' + pref_loglevel + ' is invalid;' +
                    ' using WARNING instead.')
    else:
      logger.setLevel(numeric_loglevel)


# Run when used as a plugin.
def __init_plugin__(app=None):
    Initialize()
    try:
        # Simply add the menu entry and callback
        from pymol.plugins import addmenuitemqt

        addmenuitemqt('PDB Analysis|All', PDBeLoaderDialog)
        addmenuitemqt('PDB Analysis|Molecules', PDBeMoleculeDialog)
        addmenuitemqt('PDB Analysis|Domains', PDBeDomainDialog)
        addmenuitemqt('PDB Analysis|Validation', PDBeValidationDialog)
        addmenuitemqt('PDB Analysis|Assemblies', PDBeAssemblyDialog)

    except Exception as e:
        logging.error('unable to make menu items')
        logging.exception(e)


def usage():
    usage = """
        Usage
            pymol PDB_plugin.py mmCIF_file_name
        """
    print(usage)
    # pymol.cmd.quit()


# Run when used from the command line.
if __name__ == '__main__':
    Initialize()
    mmCIF_file = None
    pdbid = None
    if len(sys.argv) > 1:
        logging.debug(sys.argv)
        for arr in sys.argv:
            pdb_id_re = re.compile(r'\d[A-z0-9]{3}')
            mmCIF_file_RE = re.compile(r'mmCIF_file=(.*)')

            if pdb_id_re.match(arr):
                pdbid = pdb_id_re.match(arr).group(0)
            if mmCIF_file_RE.match(arr):
                mmCIF_file = mmCIF_file_RE.match(arr).group(1)

        if pdbid or mmCIF_file:
            logging.debug('pdbid: %s' % pdbid)
            logging.debug('mmCIF file: %s' % mmCIF_file)
            PDBe_startup(pdbid, 'all', mmCIF_file=mmCIF_file)
        else:
            logging.error('Please provide a pdbid or mmCIF_file=FILE after -- ')

# vi:expandtab:smarttab:sw=4:tw=80
