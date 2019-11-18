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

import datetime
import importlib  # needs at least python 2.7
import json  # parsing the input
import logging
import os
import random
import re
import socket
import sys
import time

import pymol
import pymol.plugins
from pymol import cmd
from pymol import stored
from chempy.cif import *  # noqa: F401,F403

# ftp site
_EBI_FTP = ('ftp://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/divided/'
            'mmCIF/%s/%s.cif.gz')
# clean mmcif
_UPDATED_FTP = 'https://www.ebi.ac.uk/pdbe/static/entry/%s_updated.cif.gz'
_EMPTY_CIF = frozenset(['', '.', '?', None])

# dictionaries
stored.domains = {}
stored.seq_scheme = {}
stored.molecules = {}
stored.poly_count = 0
stored.residues = {}
stored.ca_p_only = set()


class PdbFetcher(object):
    """Downloads PDB data from URLs.

    This class tries to use one of several libraries to get the data, depending
    on which are installed on the system.
    """

    def __init__(self):
        self._modules = {}  # Modules loaded by this class - like sys.modules.
        self._fetcher = None  # Fetcher method to use for get_data.

        # Figure out which library we can use to fetch data.
        import_errors = []
        for lib in ('requests', 'urllib2', 'urllib.request', 'urllib.error'):
            try:
                self._modules[lib] = importlib.import_module(lib)
            except Exception as e:
                import_errors.append(e)
                continue

            if 'requests' in self._modules:
                # Module 'requests' works on python 2.x and 3.x.
                logging.debug('using requests module')
                self._fetcher = self._get_data_with_requests
                break
            elif 'urllib2' in self._modules:
                # Module 'urllib2' works on python 2.x.
                # However, the fetcher code uses the python 3.x split up naming
                # conventions, so point those names to the right place.
                logging.debug('using urllib2 module in python 2.x')
                self._modules['urllib.request'] = self._modules['urllib2']
                self._modules['urllib.error'] = self._modules['urllib2']
                self._fetcher = self._get_data_with_urllib
                break
            elif ('urllib.request' in self._modules and
                  'urllib.error' in self._modules):
                # Module 'urllib.[request,errors]' works on python 3.x.
                logging.debug('using urllib module in python 3.x')
                self._fetcher = self._get_data_with_urllib
                break
        if not self._fetcher:
            raise Exception(
                __name__.split('.')[-1] +
                ": PDB can't be accessed due to missing python libraries:\n" +
                '\n'.join(['    ' + str(error) for error in import_errors]) +
                '\n==> Install one of them!')

    def get_data(self, url, description):
        """Returns PDB data from the given URL."""
        logging.debug(description)
        url = self._quote(url)
        return self._fetcher(url, description)

    @staticmethod
    def _quote(url):
        """Returns URL with space escaped."""
        return re.sub(' ', '%20', url)

    def _get_data_with_requests(self, url, description):
        """Uses requests module to fetch and return data from the PDB URL."""
        requests = self._modules['requests']
        response = requests.get(url=url, timeout=60)
        if response.status_code == 200:
            data = response.json()
        elif response.status_code == 404:
            data = {}
        else:
            logging.debug(response.status_code, response.reason)
            data = {}

        return data

    def _get_data_with_urllib(self, url, description):
        """Uses urllib module to fetch and return data from the PDB URL."""
        urllib2_request = self._modules['urllib.request']
        urllib2_error = self._modules['urllib.error']
        data = {}
        data_response = False
        limit = 5
        logging.debug(url)
        date = datetime.datetime.now().strftime('%Y,%m,%d  %H:%M')
        for tries in range(1, limit + 1):
            try:
                response = urllib2_request.urlopen(url, None, 60)
            except urllib2_error.HTTPError as e:
                logging.debug(
                    '%s HTTP API error - %s, error code - %s, try %d' %
                    (date, description, e.code, tries))
                # logging.debug(e.code)
                if e.code == 404:
                    data_response = True
                    break
                else:
                    # logging.debug(entry_url)
                    time.sleep(random.randint(1, 20))  # sleep 1-20 seconds
                    logging.debug(e)
            except urllib2_error.URLError as e:
                logging.debug('%s URL API error - %s, try %d' %
                              (date, e, tries))
                # logging.debug(entry_url)
                time.sleep(random.randint(1, 20))  # sleep 1-20 seconds
                logging.debug(e)
            except socket.timeout as e:
                logging.debug('%s Timeout API error - %s, try %d' %
                              (date, e, tries))
                # logging.debug(entry_url)
                time.sleep(random.randint(1, 20))  # sleep 1-20 seconds
                logging.debug(e)
            else:
                logging.debug(
                    'received a response from the %s API after %d tries' %
                    (description, tries))
                data = json.load(response)
                data_response = True
                break

        if not data_response:
            logging.error('No response from the %s API' % description)

        return data


class PdbApi(object):
    """Handles getting data from the PDB API service."""

    def __init__(self, server_root='https://www.ebi.ac.uk/pdbe/api'):
        self._server_root = server_root.rstrip('/')
        self._fetcher = PdbFetcher()

    def get_summary(self, pdbid):
        url = self._get_url('pdb/entry/summary', pdbid)
        return self._fetcher.get_data(url, 'summary')

    def get_molecules(self, pdbid):
        url = self._get_url('pdb/entry/molecules', pdbid)
        return self._fetcher.get_data(url, 'molecules')

    def get_seq_scheme(self, pdbid):
        url = self._get_url('pdb/entry/residue_listing', pdbid)
        return self._fetcher.get_data(url, 'seq_scheme')

    def get_protein_domains(self, pdbid):
        url = self._get_url('mappings', pdbid)
        return self._fetcher.get_data(url, 'protein_domains')

    def get_nucleic_domains(self, pdbid):
        url = self._get_url('nucleic_mappings', pdbid)
        return self._fetcher.get_data(url, 'nucleic_domains')

    def get_validation(self, pdbid):
        url = self._get_url('validation/global-percentiles/entry', pdbid)
        return self._fetcher.get_data(url, 'validation')

    def get_residue_validation(self, pdbid):
        url = self._get_url(
            'validation/protein-RNA-DNA-geometry-outlier-residues/entry', pdbid)
        return self._fetcher.get_data(url, 'residue validation')

    def get_rama_validation(self, pdbid):
        url = self._get_url(
            'validation/protein-ramachandran-sidechain-outliers/entry', pdbid)
        return self._fetcher.get_data(url, 'rama validation')

    def _get_url(self, api_url, pdbid):
        url = '/'.join((self._server_root, api_url, pdbid))
        logging.debug('url:', url)
        return url


# FIXME(r2r): don't leave this as global variable
pdb = PdbApi()


class Color(object):
    """Manage use of colors."""

    _object_colors = [
        'red',
        'green',
        'blue',
        'yellow',
        'magenta',
        'cyan',
        'tv_red',
        'tv_green',
        'tv_blue',
        'tv_yellow',
        'lightmagenta',
        'palecyan',
        'raspberry',
        'chartreuse',
        'marine',
        'paleyellow',
        'hotpink',
        'aquamarine',
        'darksalmon',
        'splitpea',
        'slate',
        'yelloworange',
        'pink',
        'greencyan',
        'salmon',
        'smudge',
        'lightblue',
        'limon',
        'lightpink',
        'teal',
        'deepsalmon',
        'palegreen',
        'skyblue',
        'wheat',
        'dirtyviolet',
        'deepteal',
        'warmpink',
        'limegreen',
        'purpleblue',
        'sand',
        'violet',
        'lightteal',
        'firebrick',
        'lime',
        'deepblue',
        'violetpurple',
        'ruby',
        'limon',
        'density',
        'purple',
        'chocolate',
        'forest',
        'deeppurple',
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


def poly_display_type(asym, molecule_type, length):
    """return either ribbon or cartoon depending on number of polymer chains"""
    # start with it being cartoon - then change as needed.
    display_type = 'cartoon'
    if stored.poly_count > 50:
        display_type = 'ribbon'
    elif asym in stored.ca_p_only:
        logging.debug('set ribbon trace on')
        cmd.set('ribbon_trace_atoms', 1)
        display_type = 'ribbon'
    else:
        if re.search('polypeptide', molecule_type):
            if length < 20:
                display_type = 'sticks'
        if re.search('nucleotide', molecule_type):
            if length < 3:
                display_type = 'sticks'
        if re.search('saccharide', molecule_type):
            display_type = 'sticks'

    # logging.debug(
    #     'asym: %s, molecule_type: %s, length: %s, display_type: %s' %
    #     (asym, molecule_type, length, display_type))
    return display_type


class WorkerFunctions(object):

    @staticmethod
    def set_transparency(selection, transparency):
        cmd.set('cartoon_transparency', transparency, selection)
        cmd.set('ribbon_transparency', transparency, selection)
        cmd.set('stick_transparency', transparency, selection)
        cmd.set('sphere_transparency', transparency, selection)

    @staticmethod
    def count_poly(pdbid):
        # global molecules
        if not stored.molecules:
            stored.molecules = pdb.get_molecules(pdbid)
        for molecule in stored.molecules.get(pdbid, []):
            # add ca only list
            if molecule['ca_p_only']:
                for a in molecule['in_struct_asyms']:
                    stored.ca_p_only.add(a)
            if molecule['molecule_type'] not in ['water', 'bound']:
                stored.poly_count += 1


def poly_seq_scheme(pdbid):
    """build a dictionary like poly_seq_scheme"""
    data = pdb.get_seq_scheme(pdbid)
    if pdbid in data:
        for molecule in data[pdbid]['molecules']:
            for chain in molecule['chains']:
                chain_id = chain['chain_id']
                asym_id = chain['struct_asym_id']
                for residue in chain['residues']:
                    cif_num = residue['residue_number']
                    pdb_num = residue['author_residue_number']
                    pdb_ins_code = residue['author_insertion_code']
                    residue_name = residue['residue_name']
                    if residue['observed_ratio'] != 0:
                        is_observed = True
                    else:
                        is_observed = False

                    stored.seq_scheme.setdefault(asym_id, {})[cif_num] = {
                        'PDBnum': pdb_num,
                        'PDBinsCode': pdb_ins_code,
                        'observed': is_observed,
                        'residueName': residue_name,
                        'chainID': chain_id
                    }
                    # logging.debug(seq_scheme)


def get_ranges(asym_id, start, end):  # noqa: C901 too complex

    def order(start, end, start_code, end_code):
        """Returns (start_code, end_code) ordered such that start < end."""
        # logging.debug('PRE: start: %s, end %s' % (start, end))
        if type(start) is str:
            start = start.split('\\')[-1]
        if type(end) is str:
            end = end.split('\\')[-1]

        if int(start) > int(end):
            result = (end_code, start_code)
        else:
            result = (start_code, end_code)
        # logging.debug('POST: start: %s, end: %s, %s' % (start, end, result))
        return result

    def insert_code(asym_id, cif_num):
        if not stored.seq_scheme[asym_id][cif_num]['PDBinsCode']:
            pdb_num = stored.seq_scheme[asym_id][cif_num]['PDBnum']
        else:
            pdb_num = '%s%s' % (
                stored.seq_scheme[asym_id][cif_num]['PDBnum'],
                stored.seq_scheme[asym_id][cif_num]['PDBinsCode'])

        if re.search('-', str(pdb_num)):
            pdb_num = re.sub('-', '\\-', str(pdb_num))

        return pdb_num

    first_residue = 1
    last_residue = 1
    for r in stored.seq_scheme[asym_id]:
        if r <= first_residue:
            first_residue = r
        if r >= last_residue:
            last_residue = r

    if end > last_residue:
        # logging.debug('ERROR: residue %s from SIFTS is greater than the '
        #               'lastSeqRes number %s' % (end, lastSeqRes))
        end = last_residue

    # logging.debug('first residue: %s, last residue %s' %
    #               (first_residue, last_residue))
    # logging.debug(seq_scheme[asym_id][start])
    is_observed = True
    if start in stored.seq_scheme[asym_id]:
        chain = stored.seq_scheme[asym_id][start]['chainID']

        while not stored.seq_scheme[asym_id][start]['observed']:
            # logging.debug('PDB start is not observed, adding 1 to residue '
            #               'number %s' % start)
            if start == last_residue:
                # logging.debug('domain not observed')
                is_observed = False
                break
            else:
                start += 1
                if start == last_residue:
                    # logging.debug('domain not observed')
                    is_observed = False
                    break
                elif start == end:
                    # logging.debug('domain not observed')
                    is_observed = False
                    break

    if end in stored.seq_scheme[asym_id]:
        while not stored.seq_scheme[asym_id][end]['observed']:
            # logging.debug('PDB end is not observed, minusing 1 from residue '
            #               'number %s' % end)
            if end == first_residue:
                # logging.debug('domain not observed')
                is_observed = False
                break
            else:
                end -= 1
                if end == first_residue:
                    # logging.debug('domain not observed')
                    is_observed = False
                    break
                elif start == end:
                    # logging.debug('domain not observed')
                    is_observed = False
                    break

    ranges = []
    if not is_observed:
        logging.debug('domain unobserved')
        return ranges

    # logging.debug('CIF start: %s, CIF end %s' % (start, end))
    # logging.debug('PDB start: %s, PDB end %s' % (PDBstart, PDBend))

    # Need to do something about non continuous ranges which pymol doesn't
    # cope with.
    #
    # Find all sections where the residue number increases by one. If it
    # jumps then store this as a separate residue block while number is less
    # than end begin with number = start plus 1. See how the numbering
    # differs from start if its 1 then its continuous - move to the next
    # residue if its zero then there must be insert codes - this is ok -
    # move to the next residue if its greater than 1 then store the previous
    # residues as a block.
    block_start_cif = start
    for current_cif in range(start, end - 1):
        next_cif = current_cif + 1
        current_pdb = stored.seq_scheme[asym_id][current_cif]['PDBnum']
        next_pdb = stored.seq_scheme[asym_id][next_cif]['PDBnum']
        if current_pdb in _EMPTY_CIF:
            block_start_cif = next_cif
            continue

        if next_pdb in _EMPTY_CIF:
            start_pdb = stored.seq_scheme[asym_id][block_start_cif]['PDBnum']
            start_code, end_code = order(start_pdb, current_pdb,
                                         insert_code(asym_id, block_start_cif),
                                         insert_code(asym_id, current_cif))
            ranges.append({
                'PDBstart': start_code,
                'PDBend': end_code,
                'chainID': chain
            })
            # move start to next block
            block_start_cif = next_cif
            continue

        current_pdb = int(stored.seq_scheme[asym_id][current_cif]['PDBnum'])
        next_pdb = int(stored.seq_scheme[asym_id][next_cif]['PDBnum'])

        if next_pdb - current_pdb > 1:
            # logging.debug('numbering not continues, positive '
            #               'jump - store as section')
            start_pdb = stored.seq_scheme[asym_id][block_start_cif]['PDBnum']
            start_code, end_code = order(start_pdb, current_pdb,
                                         insert_code(asym_id, block_start_cif),
                                         insert_code(asym_id, current_cif))
            ranges.append({
                'PDBstart': start_code,
                'PDBend': end_code,
                'chainID': chain
            })
            # move start to next block
            block_start_cif = next_cif
        elif current_pdb - next_pdb > 1:
            # logging.debug('numbering not continues, negative '
            #               'jump - store as section')
            start_pdb = stored.seq_scheme[asym_id][block_start_cif]['PDBnum']
            start_code, end_code = order(start_pdb, current_pdb,
                                         insert_code(asym_id, block_start_cif),
                                         insert_code(asym_id, current_cif))
            ranges.append({
                'PDBstart': start_code,
                'PDBend': end_code,
                'chainID': chain
            })
            # move start to next block
            block_start_cif = next_cif

    start_pdb = stored.seq_scheme[asym_id][block_start_cif]['PDBnum']
    start_code, end_code = order(block_start_cif, end,
                                 insert_code(asym_id, block_start_cif),
                                 insert_code(asym_id, end))
    ranges.append({
        'PDBstart': start_code,
        'PDBend': end_code,
        'chainID': chain
    })
    # logging.debug(ranges)

    return ranges


class Validation(object):
    """Validation methods."""

    @classmethod
    def launch_validation(cls, pdbid):
        val_data = pdb.get_validation(pdbid)

        if val_data:
            logging.debug('There is validation for this entry')

            res_data = pdb.get_residue_validation(pdbid)
            rama_data = pdb.get_rama_validation(pdbid)

            cls.per_residue_validation(pdbid, res_data, rama_data)
        else:
            logging.debug('No validation for this entry')

            # This takes too long for really large entries.
            # Validation.per_chain_per_residue_validation(
            #   pdbid, res_data, rama_data)

    """validation of all polymeric entries"""

    @staticmethod
    def validation_selection(selection, display_id):
        # logging.debug(selection)

        # uses a list so 0 means that there is one outlier for this residue.

        if selection in stored.residues:
            color_num = stored.residues[selection] + 1
            stored.residues[selection] = color_num
        else:
            stored.residues.update({selection: 0})
            color_num = 0

        if color_num > 2:
            color_num = 2
        Color.set_validation_color(color_num, selection)

    @classmethod
    def geometric_validation(cls, pdbid, res_data):
        """check for geometric validation outliers in res_data """
        try:
            molecules = res_data[pdbid]['molecules']
        except Exception:
            molecules = []
            logging.debug('no residue validation for this entry')

        for molecule in molecules:
            # logging.debug(molecule)
            for chain in molecule['chains']:
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
                            if PDB_ins_code not in [None, ' ']:
                                PDB_res_num = '%s%s' % (PDB_res_num,
                                                        PDB_ins_code)

                            # logging.debug(PDB_res_num)
                            selection = 'chain %s and resi %s' % (chain_id,
                                                                  PDB_res_num)
                            if (model == 1 or model_id) and (not chain or
                                                             chain_id):
                                cls.validation_selection(selection, pdbid)

    @classmethod
    def ramachandran_validation(cls, pdbid, rama_data, chain=False, model=1):
        """display ramachandran outliers"""
        try:
            rama = rama_data[pdbid]
        except Exception:
            rama = []

        for k in rama:
            v = rama_data[pdbid][k]
            if v:
                logging.debug('ramachandran %s for this entry' % k)
                for outlier in rama_data[pdbid][k]:
                    # logging.debug(outlier)
                    model_id = int(outlier['model_id'])
                    chain_id = outlier['chain_id']
                    PDB_res_num = outlier['author_residue_number']
                    PDB_ins_code = outlier['author_insertion_code']

                    if PDB_ins_code not in [None, ' ']:
                        PDB_res_num = '%s%s' % (PDB_res_num, PDB_ins_code)

                    selection = 'chain %s and resi %s' % (chain_id, PDB_res_num)
                    if model == 1 or model_id:
                        if not chain or chain_id:
                            cls.validation_selection(selection, pdbid)
                    else:
                        logging.debug(
                            'is multimodel!!!, outlier not in model 1,'
                            ' not shown.')

            else:
                logging.debug('no %s' % (k))

    @classmethod
    def per_residue_validation(cls, pdbid, res_data, rama_data):
        """validation of all outliers, colored by number of outliers"""

        # display only polymers
        if not stored.molecules:
            stored.molecules = pdb.get_molecules(pdbid)
        for molecule in stored.molecules.get(pdbid, []):
            # logging.debug(molecule)
            molecule_type = molecule['molecule_type']
            if molecule_type not in ['Water', 'Bound']:
                for a in molecule['in_struct_asyms']:
                    if not stored.seq_scheme:
                        poly_seq_scheme(pdbid)
                    start = 1
                    end = molecule['length']
                    length = molecule['length']
                    asym_id = a
                    display_type = poly_display_type(asym_id, 'polypeptide',
                                                     length)
                    # logging.debug(asym_id)
                    ranges = get_ranges(asym_id, start, end)
                    for r in ranges:
                        start = r['PDBstart']
                        end = r['PDBend']
                        chain = r['chainID']
                        selection = 'chain %s and resi %s-%s and %s' % (
                            chain, start, end, pdbid)
                        # logging.debug(selection)
                        cmd.show(display_type, selection)

        Color.set_validation_background_color(pdbid)
        # display(pdbid, image_type)
        cmd.enable(pdbid)

        cls.geometric_validation(pdbid, res_data)
        cls.ramachandran_validation(pdbid, rama_data)

        # logging.debug(stored.residues)


def show_molecules(pdbid):  # noqa: C901 too complex
    logging.debug('Display entities')
    cmd.set('cartoon_transparency', 0.3, pdbid)
    cmd.set('ribbon_transparency', 0.3, pdbid)
    if not stored.molecules:
        stored.molecules = pdb.get_molecules(pdbid)
    for molecule in stored.molecules.get(pdbid, []):
        if molecule['molecule_name']:
            entity_name = '-'.join(molecule['molecule_name'])
            if len(molecule['molecule_name']) > 1:
                entity_name += '_chimera'
        else:
            logging.debug('No name from the API')
            entity_name = 'entity_%s' % (molecule['entity_id'])
        # logging.debug('molecule %s: %s' %(molecule['entity_id'], entity_name))

        # logging.debug(molecule['entity_id'])
        # see if its polymeric
        molecule_type = molecule['molecule_type']
        # entity_name = molecule['molecule_name']
        if entity_name:
            object_name = re.sub(' ', '_', entity_name)
            object_name = re.sub(r'\W+', '', object_name)
        else:
            object_name = 'entity_%s' % (molecule['entity_id'])
        # logging.debug(object_name)
        display_type = ''
        object_selection = []
        pymol_selection = ''
        if molecule_type != 'Water':
            # use asym to find residue number and chain ID
            for a in molecule['in_struct_asyms']:
                if not stored.seq_scheme:
                    poly_seq_scheme(pdbid)
                if molecule_type == 'Bound':
                    # bound molecules have no CIFresidue number and this is
                    # defaulted to 1 in the API
                    # logging.debug(entity_name)
                    # logging.debug(stored.seq_scheme[a])
                    for res in stored.seq_scheme[a]:
                        # logging.debug(stored.seq_scheme[a][res])
                        short = stored.seq_scheme[a][res]
                        chain = short['chainID']
                        res = ''
                        display_type = 'spheres'
                        # logging.debug(short)
                        if not short['PDBinsCode']:
                            res = str(short['PDBnum'])
                        else:
                            res = str(short['PDBnum']) + short['PDBinsCode']
                        selection = 'chain %s and resi %s and %s' % (chain, res,
                                                                     pdbid)
                        object_selection.append(selection)
                        # logging.debug(selection)
                else:
                    # find the first and last residues
                    start = 1
                    end = molecule['length']
                    asym_id = a
                    # need to work out if cartoon is the right thing to
                    # display
                    length = molecule['length']
                    display_type = poly_display_type(asym_id, molecule_type,
                                                     length)
                    # logging.debug(asym_id)
                    ranges = get_ranges(asym_id, start, end)
                    # logging.debug(ranges)
                    for r in ranges:
                        start = r['PDBstart']
                        end = r['PDBend']
                        chain = r['chainID']
                        selection = 'chain %s and resi %s-%s and %s' % (
                            chain, start, end, pdbid)
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

            logging.debug(pymol_selection)
            if len(object_name) > 250:
                object_name = object_name[:249]
            cmd.select('test_select', pymol_selection)
            cmd.create(object_name, 'test_select')
            # logging.debug(display_type)
            cmd.show(display_type, object_name)

            # color by molecule.
            Color.set_object_color(int(molecule['entity_id']), object_name)

    cmd.delete('test_select')
    # cmd.show('cartoon', pdbid)


def show_assemblies(pdbid, mm_cif_file):
    """iterate through the assemblies and output images"""
    logging.info('Generating assemblies')
    # cmd.hide('everything')
    WorkerFunctions.set_transparency('all', 0.0)

    try:
        assemblies = cmd.get_assembly_ids(pdbid)  # list or None
        logging.debug(assemblies)
        if assemblies:
            for assembly_id in assemblies:
                logging.debug('Assembly: %s' % (assembly_id))
                assembly_name = pdbid + '_assem_' + assembly_id
                logging.debug(assembly_name)
                cmd.set('assembly', assembly_id)
                cmd.load(mm_cif_file, assembly_name, format='cif')
                logging.debug('finished Assembly: %s' % (assembly_id))
    except Exception:
        logging.debug('pymol version does not support assemblies')


def map_domains(pdbid):  # noqa: C901 too complex
    """make domain objects"""
    if not stored.seq_scheme:
        poly_seq_scheme(pdbid)
    segment_identifier = {
        'CATH': 'domain',
        'SCOP': 'scop_id',
        'Pfam': '',
        'Rfam': ''
    }
    domain_to_make = frozenset(segment_identifier.keys())

    def store_domain(asym_id, chain, start, end, entity_id):
        stored.domains.setdefault(domain_type,
                                  {}).setdefault(domain, {}).setdefault(
                                      domain_name, []).append({
                                          'asym_id': asym_id,
                                          'chain': chain,
                                          'start': start,
                                          'end': end,
                                          'entity_id': entity_id
                                      })

    protein_domains = pdb.get_protein_domains(pdbid)
    nucleic_domains = pdb.get_nucleic_domains(pdbid)
    for data in [protein_domains, nucleic_domains]:
        if not data:
            logging.debug('no domain information for this entry')
            continue
        for domain_type in data[pdbid]:
            if domain_type not in domain_to_make:
                continue
            for domain in data[pdbid][domain_type]:
                if not data[pdbid][domain_type][domain]['mappings']:
                    continue
                # logging.debug(domain_type)
                # logging.debug(domain)
                for mapping in data[pdbid][domain_type][domain]['mappings']:
                    # check segment id
                    # domain_segment_id = mapping.get('segment_id', 1)
                    domain_name = str(
                        mapping.get(segment_identifier[domain_type], ''))
                    # logging.debug('%s: %s' % (domain_name, domain_segment_id))

                    # then check it exits
                    start = mapping['start']['residue_number']
                    end = mapping['end']['residue_number']
                    chain = mapping['chain_id']
                    entity_id = mapping['entity_id']
                    asym_id = mapping['struct_asym_id']
                    if asym_id in stored.seq_scheme:
                        ranges = get_ranges(asym_id, start, end)
                        for r in ranges:
                            start = r['PDBstart']
                            end = r['PDBend']
                            store_domain(asym_id, chain, start, end, entity_id)


def show_domains(pdbid):  # noqa: C901 too complex
    map_domains(pdbid)
    if not stored.domains:
        return

    objects = {}
    chains = {}
    # logging.debug(stored.domains)
    for domain_type in stored.domains:
        for domain in stored.domains[domain_type]:
            # logging.debug(domain)
            for instance in stored.domains[domain_type][domain]:
                asym_list = []
                entity_list = []
                # logging.debug(instance)
                object_selection = []
                pymol_selection = ''
                for segment in stored.domains[domain_type][domain][instance]:
                    pdb_start = segment['start']
                    pdb_end = segment['end']
                    chain = segment['chain']
                    asym_id = segment['asym_id']
                    entity_id = segment['entity_id']

                    asym_list.append(asym_id)
                    entity_list.append(entity_id)
                    chains[asym_id] = {'chain': chain, 'entity_id': entity_id}

                    selection = 'chain %s and resi %s-%s and %s' % (
                        chain, pdb_start, pdb_end, pdbid)
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
                object_name = '%s_%s_%s' % (domain_type, domain, instance)
                objects[object_name] = {
                    'asym_list': asym_list,
                    'entity_list': entity_list
                }
                cmd.select('test_select', pymol_selection)
                cmd.create(object_name, 'test_select')

        for chain in chains:
            logging.debug(chain)
            c_select = 'chain %s and %s' % (chains[chain]['chain'], pdbid)

            entity_id = chains[chain]['entity_id']
            length = None
            for molecule in stored.molecules[pdbid]:
                if entity_id == molecule['entity_id']:
                    length = molecule['length']
            display_type = poly_display_type(chain, 'polypeptide', length)
            cmd.show(display_type, c_select)
            cmd.color('grey', c_select)
        num = 1
        # logging.debug(objects)
        for object_name in objects:
            # logging.debug(object_name)
            cmd.enable(object_name)
            # domain can span multiple asym_id's, could be different type
            for i, a in enumerate(objects[object_name]['asym_list']):
                entity_id = objects[object_name]['entity_list'][i]
                length = None
                for molecule in stored.molecules[pdbid]:
                    if entity_id == molecule['entity_id']:
                        length = molecule['length']
                display_type = poly_display_type(a, 'polypeptide', length)
                cmd.show(display_type, object_name)
                Color.set_object_color(num, object_name)
                num += 1

        cmd.delete('test_select')


def PDBe_startup(  # noqa: 901 too complex
    pdbid, method, mm_cif_file=None, file_path=None):
    if pdbid:
        pdbid = pdbid.lower()
    if mm_cif_file:
        filename = os.path.basename(mm_cif_file)
        if re.search('_', filename):
            pdbid = re.split('_', filename)[0].lower()
        else:
            pdbid = re.split('[.]', filename)[0].lower()

    try:
        cmd.set('cif_keepinmemory')
    except Exception:
        logging.exception(
            'version of pymol does not support keeping all cif items')

    # check the PDB code actually exists.
    summary = pdb.get_summary(pdbid)

    if summary:
        # clear the dictionaries
        stored.domains = {}
        stored.seq_scheme = {}
        stored.molecules = {}
        stored.residues = {}

        logging.debug('pdbid: %s' % (pdbid))
        mid_pdb = pdbid[1:3]

        obj_list = cmd.get_object_list('all')
        if pdbid not in obj_list:
            if mm_cif_file:
                file_path = mm_cif_file
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
                        file_path = _UPDATED_FTP % pdbid
                        logging.debug('setting connect mode to mode 4')
                        cmd.set('connect_mode', 4)
                    except Exception:
                        logging.exception(
                            'pymol version does not support assemblies or '
                            'connect mode 4')
                        file_path = _EBI_FTP % (mid_pdb, pdbid)

            logging.debug('File to load: %s' % file_path)
            cmd.load(file_path, pdbid, format='cif')
        cmd.hide('everything', pdbid)
        WorkerFunctions.count_poly(pdbid)

        if method == 'molecules':
            show_molecules(pdbid)
        elif method == 'domains':
            show_molecules(pdbid)
            show_domains(pdbid)
        elif method == 'validation':
            Validation.launch_validation(pdbid)
        elif method == 'assemblies':
            show_assemblies(pdbid, file_path)
        elif method == 'all':
            show_molecules(pdbid)
            show_domains(pdbid)
            show_assemblies(pdbid, file_path)
            Validation.launch_validation(pdbid)
        else:
            logging.warning('provide a method')
        cmd.zoom(pdbid, complete=1)

    elif mm_cif_file:
        logging.warning('no PDB ID, show assemblies from mmCIF file')
        cmd.load(mm_cif_file, pdbid, format='cif')
        # WorkerFunctions.count_poly(pdbid)
        show_assemblies(pdbid, mm_cif_file)

    else:
        logging.error('please provide a 4 letter PDB code')


class PdbeGui(object):
    """Handles the GUI aspects of this pluing."""

    def __init__(self):
        # Make sure we can load the GUI library at startup.
        try:
            self._qt = importlib.import_module('pymol.Qt')
        except Exception as e:
            logging.exception(e)
            raise Exception(
                __name__.split('.')[-1] +
                ": Can't start GUI due to missing python libraries:\n" +
                '    ' + str(e))

    def _get_pdbid(self, label):
        """Gets a PDB entry ID from a dialog window and returns it."""
        pdbid, ok_pressed = self._qt.QtWidgets.QInputDialog.getText(
            None, 'PDB Loader Service',
            label + '\nPlease enter a 4-digit PDB ID:',
            self._qt.QtWidgets.QLineEdit.Normal, '')

        if ok_pressed:
            return pdbid
        else:
            return None

    def analyze_all(self):
        pdbid = self._get_pdbid(
            'Highlight chemically distinct molecules, domains and assemblies'
            ' in a PDB entry.')
        if pdbid:
            PDBe_startup(pdbid, 'all')

    def analyze_molecules(self):
        pdbid = self._get_pdbid(
            'Highlight chemically distinct molecules in a PDB entry.')
        if pdbid:
            PDBe_startup(pdbid, 'molecules')

    def analyze_domains(self):
        pdbid = self._get_pdbid(
            'Display Pfam, SCOP, CATH and Rfam domains on a PDB entry.')
        if pdbid:
            PDBe_startup(pdbid, 'domains')

    def analyze_validation(self):
        pdbid = self._get_pdbid('Display geometric outliers on a PDB entry.')
        if pdbid:
            PDBe_startup(pdbid, 'validation')

    def analyze_assemblies(self):
        pdbid = self._get_pdbid('Display assemblies for a PDB entry.')
        if pdbid:
            PDBe_startup(pdbid, 'assemblies')


class PdbIdAutocomplete(cmd.Shortcut):
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
        key = key.lower()
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
            summary = pdb.get_summary(key)
            # If the key doesn't return a valid description return None to
            # indicate that the key is invalid.
            try:
                title = summary[key][0]['title']
            except Exception:
                return None
            # We've got a valid key. Show the title of that PDB entry to the
            # user and make the key the only possible autocomplete match.
            print(key + ':', title)
            return key


@cmd.extendaa([lambda: PdbIdAutocomplete(), 'PDB Entry Id', ''])
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


@cmd.extendaa([lambda: PdbIdAutocomplete(), 'PDB Entry Id', ''])
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


@cmd.extendaa([lambda: PdbIdAutocomplete(), 'PDB Entry Id', ''])
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


@cmd.extendaa([lambda: PdbIdAutocomplete(), 'PDB Entry Id', ''])
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
    chains = cmd.get_chains(selection, state)
    # PyMOL 2.0 returns -1 instead of [] if no chains are found.
    if isinstance(chains, list):
        n = len(chains)
    else:
        logging.warning("get_chains didn't return a list: %s" % chains)
        n = 0
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
    numeric_loglevel = None
    if isinstance(loglevel, str):
        numeric_loglevel = getattr(logging, loglevel.upper(), None)
    if numeric_loglevel is None:
        logger.setLevel(logging.WARNING)
        logging.error('Invalid preference %s = "%s"; using WARNING instead.' %
                      (pref_loglevel, loglevel))
    else:
        logger.setLevel(numeric_loglevel)


# Run when used as a plugin.
def __init_plugin__(app=None):
    Initialize()
    gui = PdbeGui()
    try:
        # Simply add the menu entry and callback
        from pymol.plugins import addmenuitemqt

        submenu = 'PDB Analysis|'
        addmenuitemqt(submenu + 'All', lambda: gui.analyze_all())
        addmenuitemqt(submenu + 'Molecules', lambda: gui.analyze_molecules())
        addmenuitemqt(submenu + 'Domains', lambda: gui.analyze_domains())
        addmenuitemqt(submenu + 'Validation', lambda: gui.analyze_validation())
        addmenuitemqt(submenu + 'Assemblies', lambda: gui.analyze_assemblies())

    except Exception as e:
        logging.error('unable to make menu items')
        logging.exception(e)


def usage():
    usage = """
Usage
    pymol PDB_plugin.py [pdbid | mmCIF_file=<file>]

        pdbid ... PDB entry ID; CIF data will be downloaded
        file .... CIF file for a PDB entry.
                  The file name should start with PDB ID.
"""
    print(usage)
    # pymol.cmd.quit()


# Run when used from the command line.
def main(argv=sys.argv):
    Initialize()
    mm_cif_file = None
    pdbid = None
    logging.debug(argv)
    for arg in argv:
        pdb_id_re = re.compile(r'\d[A-z0-9]{3}')
        mm_cif_file_re = re.compile(r'mmCIF_file=(.*)')

        match = pdb_id_re.match(arg)
        if match:
            pdbid = match.group(0)
        match = mm_cif_file_re.match(arg)
        if match:
            mm_cif_file = match.group(1)

    if pdbid or mm_cif_file:
        logging.debug('pdbid: %s' % pdbid)
        logging.debug('mmCIF file: %s' % mm_cif_file)
        PDBe_startup(pdbid, 'all', mm_cif_file=mm_cif_file)
    else:
        logging.error('Please provide a pdbid or mmCIF_file=FILE after -- ')
        usage()


if __name__ == '__main__':
    main()

# vi:expandtab:smarttab:sw=4:tw=80
