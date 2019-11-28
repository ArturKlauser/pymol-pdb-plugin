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

# Documentation resources:
#   mmCIF: (macromolecular Crystallographic Information Framework)
#          http://mmcif.wwpdb.org/
# Naming Cheat Sheet:
#   entity: synonymous with molecule or molecular entity
#     within molecule:
#       in_chains: corresponds to PyMOL chain name
#       in_struct_asyms: corresponds to PyMOL segment name
#       top level molecule name/key: corresponds to PyMOL object name
#   asym_id: unique identifier within the crystallographic asymmetric unit

from __future__ import print_function

from collections import namedtuple
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

# ftp site
_EBI_FTP = ('ftp://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/divided/'
            'mmCIF/%s/%s.cif.gz')
# clean mmcif
_UPDATED_FTP = 'https://www.ebi.ac.uk/pdbe/static/entry/%s_updated.cif.gz'
_EMPTY_PDB_NUM = frozenset(['', '.', '?', None])


def clear_analysis_data():
    # Note that pymol.stored is just a convenient global variable where this
    # data is attached to. It is not necessary to have it there for interaction
    # with pymol.

    # --- inputs from PDB API
    stored.molecules = {}  # get_molecules()
    Sequences.clear()

    # --- results of analysis
    stored.residues = {}  # validation_selection()
    stored.polymer_count = 0  # Molecules()._process_molecules()
    stored.ca_p_only_segments = set()  # Molecules()._process_molecules()


def extendaa(*arg, **kw):
    """Wrapper across cmd.extendaa that also adds new command to pymol.cmd."""

    def wrapper(func):
        _self = kw.get('_self', cmd)
        name = func.__name__
        setattr(_self, name, func)
        return _self.extendaa(*arg, **kw)(func)

    return wrapper


class PdbFetcher(object):
    """Downloads PDB json data from URLs.

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

    def get_data(self, url, description, **kw):
        """Returns PDB data from the given URL."""
        logging.debug(description)
        url = self._quote(url)
        return self._fetcher(url, description, **kw)

    @staticmethod
    def _quote(url):
        """Returns URL with space escaped."""
        return url.replace(' ', '%20')

    def _get_data_with_requests(self, url, description):
        """Uses requests module to fetch and return data from the PDB URL."""
        requests = self._modules['requests']
        response = requests.get(url=url, timeout=60)
        if response.status_code == 200:
            data = response.json()
        elif response.status_code == 404:
            data = {}
        else:
            logging.debug('%d %s' % (response.status_code, response.reason))
            data = {}

        return data

    def _get_data_with_urllib(self,
                              url,
                              description,
                              sleep_min=1,
                              sleep_max=20):
        """Uses urllib module to fetch and return data from the PDB URL.

        Argument sleep_(min,max)=0 is mostly interesting for testing.
        """
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
                    time.sleep(random.randint(sleep_min, sleep_max))
                    logging.debug(e)
            except urllib2_error.URLError as e:
                logging.debug('%s URL API error - %s, try %d' %
                              (date, e, tries))
                # logging.debug(entry_url)
                time.sleep(random.randint(sleep_min, sleep_max))
                logging.debug(e)
            except socket.timeout as e:
                logging.debug('%s Timeout API error - %s, try %d' %
                              (date, e, tries))
                # logging.debug(entry_url)
                time.sleep(random.randint(sleep_min, sleep_max))
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
        """Returns a dictionary of molecule data.

        The contents is a conglomerate of data found in the mmCIF entity group
        http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Groups/entity_group.html
        (Better reference would be good)

        Format:
          <molecules> = dict(<PDB_ID>: list(<entity>))
            <entity> = dict(<key>: <value>)
            some keys:
              entity_id: unique numeric ID (key)
              molecule_name: list(name of this molecule)
              in_chains: list(<chain_name>)
              in_struct_asyms: list(<segment_name>)
              number_of_copies: number of duplicate entities
              molecule_type: enum(polypeptide|Bound|Water|...)
                for polypeptide type:
                  sequence: sequence string
                  length: sequence length
              ca_p_only: true/false; only CA and/or P atom positions available
        """
        url = self._get_url('pdb/entry/molecules', pdbid)
        return self._fetcher.get_data(url, 'molecules')

    def get_sequences(self, pdbid):
        """Returns a dictionary of sequence data.

        The contents is similar to mmCIF poly_seq_scheme.
        http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/pdbx_poly_seq_scheme.html

        Format:
          <sequences> = dict(<PDB_ID>: dict("molecules": list(<molecule>)))
            <molecule> = dict(<key>: <value)
              available keys:
                entity_id: molecule number (see get_molecules above)
                chains: list(<chain>)
              <chain> = dict(<key>: <value>)
                available keys:
                  chain_id: PyMOL's chain name (e.g in selection)
                  struct_asym_id: PyMOL's segment name (e.g. in selection)
                  residues: list(<residue>)
                      <residue>: dict(<key>: <value>)
                        available keys:
                          residue_number: int; sequential within this dict
                          author_residue_number: int; PyMOL residue number
                                                 within chain
                          observer_ration: float; ?
                          residue_name: char; 3-letter residue code
        """
        url = self._get_url('pdb/entry/residue_listing', pdbid)
        return self._fetcher.get_data(url, 'sequences')

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

    def get_ramachandran_validation(self, pdbid):
        url = self._get_url(
            'validation/protein-ramachandran-sidechain-outliers/entry', pdbid)
        return self._fetcher.get_data(url, 'ramachandran validation')

    def _get_url(self, api_url, pdbid):
        url = '/'.join((self._server_root, api_url, pdbid))
        logging.debug('url: %s' % url)
        return url


# FIXME(r2r): don't leave this as global variable
pdb = PdbApi()


class Presentation(object):
    """Manage PyMOL presentation properties."""

    _OBJECT_COLORS = [
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
    _VALIDATION_COLORS = ['yellow', 'orange', 'red']

    @classmethod
    def set_object_color(cls, color_num, obj):
        """Colors object from predetermined set, reusing as necessary."""
        color = cls._OBJECT_COLORS[color_num % len(cls._OBJECT_COLORS)]
        cmd.color(color, obj)

    @classmethod
    def set_validation_color(cls, color_num, obj):
        """Colors object from predetermined set, clamping to max."""
        color = cls._VALIDATION_COLORS[min(color_num,
                                           len(cls._VALIDATION_COLORS) - 1)]
        cmd.color(color, obj)

    @classmethod
    def set_validation_background_color(cls, obj):
        """Colors object with validation background color."""
        color = 'validation_background_color'
        # Define this as a new color name.
        cmd.set_color(color, [0.4, 1.0, 0.4])
        cmd.color(color, obj)

    @staticmethod
    def set_transparency(selection, transparency):
        cmd.set('cartoon_transparency', transparency, selection)
        cmd.set('ribbon_transparency', transparency, selection)
        cmd.set('stick_transparency', transparency, selection)
        cmd.set('sphere_transparency', transparency, selection)


def get_polymer_display_type(segment_id, molecule_type, length):
    """Returns display type depending on molecule complexity."""
    # Start with it being cartoon - then change as needed.
    display_type = 'cartoon'
    if stored.polymer_count > 50:
        display_type = 'ribbon'
    elif segment_id in stored.ca_p_only_segments:
        logging.debug('set ribbon trace on')
        cmd.set('ribbon_trace_atoms', 1)
        display_type = 'ribbon'
    elif 'polypeptide' in molecule_type and length < 20:
        display_type = 'sticks'
    elif 'nucleotide' in molecule_type and length < 3:
        display_type = 'sticks'
    elif 'saccharide' in molecule_type:
        display_type = 'sticks'

    # logging.debug(
    #     'segment_id: %s, molecule_type: %s, length: %s, display_type: %s' %
    #     (segment_id, molecule_type, length, display_type))
    return display_type


class Sequences(object):
    """Polymer sequences."""

    Residue = namedtuple('Residue',
                         'chain_id pdb_num pdb_residue_num is_observed')
    Range = namedtuple('Range', 'chain_id start_residue_num end_residue_num')

    # TODO(r2r): For refactoring purposes temporarily moved stored.sequences
    # into a Sequences class-level object. This should really be instance-level.
    _sequences = {}  # get_sequences() reformatted by Sequences()._build()
    """Stores polymer sequence data.

    Format:
      <sequences> = dict(<segment_id>: <residues>)
        <residues> = dict(<residue_num>: <residue>)
          <residue_num> = sequential numeric residue id within this sequence
          <residue> = namedtuple
            chain_id: PyMOL chain name
            pdb_num: PDB residue number (excluding insertion codes)
            pdb_residue_num: PyMOL residue number (including insertion codes),
                             i.e. 'resi' in selections
            is_observed: bool TODO(r2r): exact meaning unclear
    """

    def __init__(self, pdbid):
        self._pdbid = pdbid
        if not self._sequences:
            self._build()

    @classmethod
    def clear(cls):
        cls._sequences = {}

    @staticmethod
    def get_pdb_residue_num(pdb_num, pdb_insertion_code):
        """Returns a residue number as used in PDB data."""
        if not pdb_insertion_code or pdb_insertion_code == ' ':
            pdb_residue_num = str(pdb_num)
        else:
            pdb_residue_num = '%s%s' % (pdb_num, pdb_insertion_code)
        # TODO(r2r): not clear why this substitution would be necessary.
        pdb_residue_num = pdb_residue_num.replace('-', '\\-')
        return pdb_residue_num

    def _build(self):
        """Builds a dictionary of sequence residues."""
        data = pdb.get_sequences(self._pdbid)
        if self._pdbid not in data:
            return
        for molecule in data[self._pdbid]['molecules']:
            for chain in molecule['chains']:
                chain_id = chain['chain_id']
                segment_id = chain['struct_asym_id']
                for residue in chain['residues']:
                    residue_num = residue['residue_number']
                    pdb_num = residue['author_residue_number']
                    pdb_ins_code = residue['author_insertion_code']
                    pdb_residue_num = self.get_pdb_residue_num(
                        pdb_num, pdb_ins_code)
                    if residue['observed_ratio'] != 0:
                        is_observed = True
                    else:
                        is_observed = False

                    self._sequences.setdefault(segment_id,
                                               {})[residue_num] = self.Residue(
                                                   chain_id, pdb_num,
                                                   pdb_residue_num, is_observed)
                    # logging.debug(self._sequences)

    @staticmethod
    def _order_range(start, end, start_pdb_residue_num, end_pdb_residue_num):
        """Returns (start_pdb_residue_num, end_pdb_residue_num) ordered such
        that start < end.
        """
        # logging.debug('PRE: start: %s, end %s' % (start, end))
        if type(start) is str:
            start = start.split('\\')[-1]
        if type(end) is str:
            end = end.split('\\')[-1]

        if int(start) > int(end):
            result = (end_pdb_residue_num, start_pdb_residue_num)
        else:
            result = (start_pdb_residue_num, end_pdb_residue_num)
        # logging.debug('POST: start: %s, end: %s, %s' % (start, end, result))
        return result

    @classmethod
    def _append_range(cls, start_residue, end_residue, ranges):
        if not start_residue:
            return  # Range was not started yet - ignore.

        start_pdb_residue_num, end_pdb_residue_num = cls._order_range(
            start_residue.pdb_num, end_residue.pdb_num,
            start_residue.pdb_residue_num, end_residue.pdb_residue_num)
        ranges.append(
            cls.Range(start_residue.chain_id, start_pdb_residue_num,
                      end_pdb_residue_num))

    @staticmethod
    def _get_trimmed_range(sequence, start_residue_num, end_residue_num):
        """Trims the sequence range of unobserved residues."""
        # Get the range of existing residue numbers.
        residue_numbers = sequence.keys()
        first_residue_num = min(residue_numbers)
        last_residue_num = max(residue_numbers)
        # Bound start/end residue number range to existing residue numbers.
        start_residue_num = max(start_residue_num, first_residue_num)
        end_residue_num = min(end_residue_num, last_residue_num)
        # logging.debug('first_residue_num: %s, last_residue_num: %s' %
        #               (first_residue_num, last_residue_num))
        # logging.debug('start_residue_num: %s, end_residue_num: %s' %
        #               (start_residue_num, end_residue_num))

        # Trim unobserved ends of sequence.
        while not sequence[start_residue_num].is_observed:
            start_residue_num += 1
            if start_residue_num > end_residue_num:
                logging.debug('domain unobserved')
                return (None, None)
        while not sequence[end_residue_num].is_observed:
            end_residue_num -= 1
            if start_residue_num > end_residue_num:
                logging.debug('domain unobserved')
                return (None, None)
        return (start_residue_num, end_residue_num)

    @classmethod
    def get_ranges(cls, segment_id, start_residue_num, end_residue_num):
        """Returns contiguous residue ranges present in the sequence.

        Pymol doesn't cope with non-contiguous ranges.

        This function finds all ranges where the residue number increases by 1.
        If it jumps then a separate residue range is generated.
        """
        ranges = []
        if segment_id not in cls._sequences:
            return ranges

        sequence = cls._sequences[segment_id]
        start_residue_num, end_residue_num = cls._get_trimmed_range(
            sequence, start_residue_num, end_residue_num)
        # logging.debug('start_residue_num: %s, end_residue_num: %s' %
        #               (start_residue_num, end_residue_num))
        if not start_residue_num:
            return ranges

        # Look at difference of pdb_num between neighboring residues in
        # start-end range:
        # - If its 1 then its contiguous; continue increasing current range.
        # - If its 0 then there must be insert codes; that's ok, continue
        #   increasing current range.
        # - if its > 1 or < 0 then there is a discontinuity; store the previous
        #   residues as a range and start a new range on the next residue.
        # Also, if we run into a residue that doesn't have a valid pdb_num then
        # close off the previous range.
        range_start_residue = sequence[start_residue_num]
        for current_residue_num in range(start_residue_num, end_residue_num):
            current_residue = sequence[current_residue_num]
            current_pdb_num = current_residue.pdb_num
            next_pdb_num = sequence[current_residue_num + 1].pdb_num

            # Skip over residues that don't have valid PDB numbers.
            if current_pdb_num in _EMPTY_PDB_NUM:
                continue
            # Start a new range if none is started yet.
            if not range_start_residue:
                range_start_residue = current_residue
            if next_pdb_num in _EMPTY_PDB_NUM:
                cls._append_range(range_start_residue, current_residue, ranges)
                range_start_residue = None
                continue

            pdb_num_jump = int(next_pdb_num) - int(current_pdb_num)
            if pdb_num_jump > 1 or pdb_num_jump < 0:
                # logging.debug('numbering not contiguous, jump %d - '
                #               'store as range' % pdb_num_jump)
                cls._append_range(range_start_residue, current_residue, ranges)
                range_start_residue = None

        # Append the last open range (if any) until end of residues of interst.
        cls._append_range(range_start_residue, sequence[end_residue_num],
                          ranges)
        # logging.debug(ranges)

        return ranges

    def get_range_selection(self, rng):
        """Returns a PyMOL selection for the range."""
        selection = 'chain %s and resi %s-%s and %s' % (
            rng.chain_id, rng.start_residue_num, rng.end_residue_num,
            self._pdbid)
        # logging.debug(selection)
        return selection

    def append_residue_selections(self, segment_id, selections):
        """Adds PyMOL selections strings for the segment to selections list."""
        for residue in self._sequences[segment_id].values():
            # logging.debug(residue)
            selection = 'chain %s and resi %s and %s' % (
                residue.chain_id, residue.pdb_residue_num, self._pdbid)
            selections.append(selection)
            # logging.debug(selection)


class Validation(object):
    """Validation methods."""

    @classmethod
    def launch_validation(cls, pdbid):
        val_data = pdb.get_validation(pdbid)

        if val_data:
            logging.debug('There is validation for this entry')

            residue_data = pdb.get_residue_validation(pdbid)
            ramachandran_data = pdb.get_ramachandran_validation(pdbid)

            cls.per_residue_validation(pdbid, residue_data, ramachandran_data)
        else:
            logging.debug('No validation for this entry')

            # This takes too long for really large entries.
            # Validation.per_chain_per_residue_validation(
            #   pdbid, residue_data, ramachandran_data)

    """validation of all polymeric entries"""

    @staticmethod
    def validation_selection(selection, display_id):
        # logging.debug(selection)

        # uses a list so 0 means that there is one outlier for this residue.
        if selection in stored.residues:
            color_num = stored.residues[selection] + 1
        else:
            color_num = 0
        stored.residues[selection] = color_num
        Presentation.set_validation_color(color_num, selection)

    @classmethod
    def geometric_validation(cls, pdbid, residue_data):
        """check for geometric validation outliers in residue_data """
        try:
            molecules = residue_data[pdbid]['molecules']
        except Exception:
            logging.debug('no residue validation for this entry')
            return

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
                            pdb_residue_num = Sequences.get_pdb_residue_num(
                                outlier['author_residue_number'],
                                outlier['author_insertion_code'])
                            # logging.debug(pdb_residue_num)
                            selection = 'chain %s and resi %s' % (
                                chain_id, pdb_residue_num)
                            if (model == 1 or model_id) and (not chain or
                                                             chain_id):
                                cls.validation_selection(selection, pdbid)

    @classmethod
    def ramachandran_validation(cls,
                                pdbid,
                                ramachandran_data,
                                chain=False,
                                model=1):
        """display ramachandran outliers"""

        if pdbid not in ramachandran_data:
            return

        for key in ramachandran_data[pdbid]:
            outliers = ramachandran_data[pdbid][key]
            if not outliers:
                logging.debug('no %s' % key)
                continue

            logging.debug('ramachandran %s for this entry' % key)
            for outlier in outliers:
                # logging.debug(outlier)
                model_id = int(outlier['model_id'])
                chain_id = outlier['chain_id']
                pdb_residue_num = Sequences.get_pdb_residue_num(
                    outlier['author_residue_number'],
                    outlier['author_insertion_code'])
                selection = 'chain %s and resi %s' % (chain_id, pdb_residue_num)
                if model == 1 or model_id:
                    if not chain or chain_id:
                        cls.validation_selection(selection, pdbid)
                else:
                    logging.debug('is multimodel!!!, outlier not in model 1,'
                                  ' not shown.')

    @classmethod
    def per_residue_validation(cls, pdbid, residue_data, ramachandran_data):
        """validation of all outliers, colored by number of outliers"""

        sequences = Sequences(pdbid)
        # display only polymers
        if not stored.molecules:
            stored.molecules = pdb.get_molecules(pdbid)
        for molecule in stored.molecules.get(pdbid, []):
            # logging.debug(molecule)
            molecule_type = molecule['molecule_type']
            if molecule_type in ['Water', 'Bound']:
                continue
            for segment_id in molecule['in_struct_asyms']:
                length = molecule['length']
                display_type = get_polymer_display_type(segment_id,
                                                        'polypeptide', length)
                # logging.debug(segment_id)
                ranges = sequences.get_ranges(segment_id, 1, length)
                for rng in ranges:
                    selection = sequences.get_range_selection(rng)
                    cmd.show(display_type, selection)

        Presentation.set_validation_background_color(pdbid)
        # display(pdbid, image_type)
        cmd.enable(pdbid)

        cls.geometric_validation(pdbid, residue_data)
        cls.ramachandran_validation(pdbid, ramachandran_data)

        # logging.debug(stored.residues)


class Molecules(object):
    """Analyze and visualize molecules."""

    def __init__(self, pdbid):
        self._pdbid = pdbid
        # Analysis depends on some globally available data; load it.
        if not stored.molecules:
            stored.molecules = pdb.get_molecules(pdbid)
        self._process_molecules()
        self._sequences = Sequences(pdbid)

    def _process_molecules(self):
        """Generate derivative data from raw molecules.

        * Creates ca_p_only_segments set which contains all segment_ids that
          only have CA and/or P atom positions available in location data.
        * Creates a global count of polymers in the data.
        """
        for molecule in stored.molecules.get(self._pdbid, []):
            # add ca only list
            if molecule['ca_p_only']:
                for segment_id in molecule['in_struct_asyms']:
                    stored.ca_p_only_segments.add(segment_id)
            if molecule['molecule_type'] not in ['Water', 'Bound']:
                stored.polymer_count += 1

    def _process_molecule(self, molecule):
        """Returns the display type and selection criteria for the molecule."""
        display_type = ''
        selections = []
        molecule_type = molecule['molecule_type']
        for segment_id in molecule['in_struct_asyms']:
            if molecule_type == 'Bound':
                # Use segment_id to find residue number and chain ID.
                # logging.debug(entity_name)
                display_type = 'spheres'
                self._sequences.append_residue_selections(
                    segment_id, selections)
            else:
                # Need to work out if cartoon is the right thing to display.
                length = molecule['length']
                display_type = get_polymer_display_type(segment_id,
                                                        molecule_type, length)
                # logging.debug(segment_id)
                ranges = self._sequences.get_ranges(segment_id, 1, length)
                for rng in ranges:
                    selection = self._sequences.get_range_selection(rng)
                    selections.append(selection)

        return display_type, selections

    def show(self):
        logging.debug('Display molecules')
        cmd.set('cartoon_transparency', 0.3, self._pdbid)
        cmd.set('ribbon_transparency', 0.3, self._pdbid)
        for molecule in stored.molecules.get(self._pdbid, []):
            if molecule['molecule_type'] == 'Water':
                continue  # Don't show water.

            molecule_name = molecule['molecule_name']
            if molecule_name:
                entity_name = '-'.join(molecule_name)
                if len(molecule_name) > 1:
                    entity_name += '_chimera'
            else:
                logging.debug('No name from the API')
                entity_name = 'entity_%s' % molecule['entity_id']
            # logging.debug('molecule %s: %s' %
            #               (molecule['entity_id'], entity_name))

            object_name = entity_name.replace(' ', '_')
            object_name = re.sub(r'\W+', '', object_name)
            object_name = object_name[:250]
            # logging.debug(object_name)

            display_type, selections = self._process_molecule(molecule)

            pymol_selection = ' or '.join(['(%s)' % x for x in selections])
            logging.debug(pymol_selection)
            cmd.select('temp_select', pymol_selection)
            cmd.create(object_name, 'temp_select')
            # logging.debug(display_type)
            cmd.show(display_type, object_name)

            # Color by molecule.
            Presentation.set_object_color(int(molecule['entity_id']),
                                          object_name)

        cmd.delete('temp_select')


def show_assemblies(pdbid, mm_cif_file):
    """iterate through the assemblies and output images"""
    logging.info('Generating assemblies')
    # cmd.hide('everything')
    Presentation.set_transparency('all', 0.0)

    try:
        assemblies = cmd.get_assembly_ids(pdbid)  # list or None
        assemblies = assemblies if assemblies else []
        logging.debug(assemblies)
        for assembly_id in assemblies:
            logging.debug('Assembly: %s' % assembly_id)
            assembly_name = pdbid + '_assem_' + assembly_id
            logging.debug(assembly_name)
            cmd.set('assembly', assembly_id)
            cmd.load(mm_cif_file, assembly_name, format='cif')
            logging.debug('finished Assembly: %s' % assembly_id)
    except Exception:
        logging.debug('pymol version does not support assemblies')


class Domains(object):
    """Analyze and visualize domains."""
    _SEGMENT_ID_NAME = {
        # domain_type -> segment_id key string
        'CATH': 'domain',
        'SCOP': 'scop_id',
        'Pfam': '',
        'Rfam': ''
    }
    Segment = namedtuple('Segment', 'entity_id chain_id segment_id start end')

    def __init__(self, pdbid):
        self._pdbid = pdbid
        # Analysis depends on some globally available data; load it.
        if not stored.molecules:
            stored.molecules = pdb.get_molecules(pdbid)
        self._sequences = Sequences(pdbid)

    def _map_all(self):
        """Make all domains."""
        # mapped_domains = dict(<domain_type>: <typed_domain>)
        # <typed_domain> = dict(<domain_id>: <named_domain>)
        # <named_domain> = dict(<domain_name>: list(Segment)
        mapped_domains = {}

        self._map_domains(pdb.get_protein_domains(self._pdbid), mapped_domains)
        self._map_domains(pdb.get_nucleic_domains(self._pdbid), mapped_domains)
        return mapped_domains

    def _map_domains(self, domains, mapped_domains):
        if not domains:
            logging.debug('no domain information for this entry')
            return
        for domain_type in domains[self._pdbid]:
            try:
                segment_id_name = self._SEGMENT_ID_NAME[domain_type]
            except Exception:
                # Ignore all domain types not in _SEGMENT_ID_NAME.
                continue
            for domain_id, domain in domains[self._pdbid][domain_type].items():
                # logging.debug(domain_type)
                # logging.debug(domain_id)
                for mapping in domain.get('mappings', []):
                    domain_name = str(mapping.get(segment_id_name, ''))
                    start_residue_num = mapping['start']['residue_number']
                    end_residue_num = mapping['end']['residue_number']
                    chain_id = mapping['chain_id']
                    entity_id = mapping['entity_id']
                    segment_id = mapping['struct_asym_id']
                    ranges = self._sequences.get_ranges(segment_id,
                                                        start_residue_num,
                                                        end_residue_num)
                    for rng in ranges:
                        mapped_domains.setdefault(domain_type, {}).setdefault(
                            domain_id, {}).setdefault(domain_name, []).append(
                                self.Segment(entity_id, chain_id, segment_id,
                                             rng.start_residue_num,
                                             rng.end_residue_num))

    def show(self):
        mapped_domains = self._map_all()
        if not mapped_domains:
            return

        # Precompute molecule lengths indexed by entity_id.
        molecule_length = {}  # entity_id -> length
        for molecule in stored.molecules[self._pdbid]:
            entity_id = molecule['entity_id']
            molecule_length[entity_id] = molecule.get('length', None)

        Object = namedtuple('Object', 'name entity_ids segment_ids')
        Chain = namedtuple('Chain', 'entity_id chain_id segment_id')
        # logging.debug(mapped_domains)
        for domain_type, typed_domain in mapped_domains.items():
            objects = []  # list of Object
            chains = []  # list of Chain
            for domain_id, named_domain in typed_domain.items():
                # logging.debug(domain_id)
                for domain_name, segments in named_domain.items():
                    # logging.debug(domain_name)
                    segment_ids = []
                    entity_ids = []
                    object_selections = []
                    for segment in segments:
                        segment_ids.append(segment.segment_id)
                        entity_ids.append(segment.entity_id)
                        chains.append(
                            Chain(segment.entity_id, segment.chain_id,
                                  segment.segment_id))

                        selection = 'chain %s and resi %s-%s and %s' % (
                            segment.chain_id, segment.start, segment.end,
                            self._pdbid)
                        object_selections.append(selection)

                    # Create an object containing all segments.
                    pymol_selection = ' or '.join(
                        ['(%s)' % x for x in object_selections])
                    logging.debug(pymol_selection)
                    object_name = '%s_%s_%s' % (domain_type, domain_id,
                                                domain_name)
                    objects.append(Object(object_name, entity_ids, segment_ids))
                    cmd.select('temp_select', pymol_selection)
                    cmd.create(object_name, 'temp_select')

            # Show all original chains in grey as default background.
            for chain in chains:
                logging.debug(chain.segment_id)
                c_select = 'chain %s and %s' % (chain.chain_id, self._pdbid)
                entity_id = chain.entity_id
                length = molecule_length[entity_id]
                display_type = get_polymer_display_type(chain.segment_id,
                                                        'polypeptide', length)
                cmd.show(display_type, c_select)
                cmd.color('grey', c_select)

            # Show each mapped object in a different color.
            num = 1
            for obj in objects:
                # logging.debug(obj.name)
                cmd.enable(obj.name)
                # A domain can span multiple segment_ids which could be
                # of different display type.
                # TODO(r2r): Seems broken. It's just overriding the object's
                # display type and color multiple times - last one wins.
                for i, segment_id in enumerate(obj.segment_ids):
                    entity_id = obj.entity_ids[i]
                    length = molecule_length[entity_id]
                    display_type = get_polymer_display_type(
                        segment_id, 'polypeptide', length)
                    cmd.show(display_type, obj.name)
                    Presentation.set_object_color(num, obj.name)
                    num += 1

            cmd.delete('temp_select')


def PDBe_startup(  # noqa: 901 too complex
    pdbid, method, mm_cif_file=None, file_path=None):
    if pdbid:
        pdbid = pdbid.lower()
    if mm_cif_file:
        filename = os.path.basename(mm_cif_file)
        if '_' in filename:
            pdbid = filename.split('_')[0].lower()
        else:
            pdbid = filename.split('.')[0].lower()

    try:
        cmd.set('cif_keepinmemory')
    except Exception:
        logging.exception(
            'version of pymol does not support keeping all cif items')

    # check the PDB code actually exists.
    summary = pdb.get_summary(pdbid)

    if summary:
        # start over fresh
        clear_analysis_data()

        logging.debug('pdbid: %s' % pdbid)
        mid_pdb = pdbid[1:3]

        obj_list = cmd.get_object_list('all')
        if pdbid not in obj_list:
            if mm_cif_file:
                file_path = mm_cif_file
            else:
                # if local file exists then use it
                if os.path.exists('%s.cif' % pdbid):
                    logging.debug('CIF found in current directory')
                    file_path = '%s.cif' % pdbid
                else:
                    logging.debug('fetching %s' % pdbid)
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
        molecules = Molecules(pdbid)

        if method == 'molecules':
            molecules.show()
        elif method == 'domains':
            molecules.show()
            Domains(pdbid).show()
        elif method == 'validation':
            Validation.launch_validation(pdbid)
        elif method == 'assemblies':
            show_assemblies(pdbid, file_path)
        elif method == 'all':
            molecules.show()
            Domains(pdbid).show()
            show_assemblies(pdbid, file_path)
            Validation.launch_validation(pdbid)
        else:
            logging.warning('provide a method')
        cmd.zoom(pdbid, complete=1)

    elif mm_cif_file:
        logging.warning('no PDB ID, show assemblies from mmCIF file')
        cmd.load(mm_cif_file, pdbid, format='cif')
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
        return pdbid if ok_pressed else None

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

    MAX_CACHE_SIZE = 100

    def __init__(self, get_summary=None):
        """Constructor.
        Args:
            get_summary: Function that takes the pdbid as its only argument
                and returns a PDB summary or None if the pdbid is not valid.
        """

        self._pdb_id_re = re.compile(r'\d[a-zA-Z0-9]{0,3}$')
        # for dependency injection
        self._get_summary = (get_summary if get_summary else
                             lambda key: pdb.get_summary(key))
        self._title_cache = {}

    def _get_title_from_pdb_with_caching(self, key):
        # Try to get it from the cache.
        if key in self._title_cache:
            return self._title_cache[key]

        summary = self._get_summary(key)
        # If the key doesn't return a valid description return None to
        # indicate that the key is invalid.
        try:
            title = summary[key][0]['title']
        except Exception:
            title = None
        # Purge the cache if it gets too large.
        if len(self._title_cache) >= self.MAX_CACHE_SIZE:
            self._title_cache = {}
        # Put it into the cache, even for negative results.
        self._title_cache[key] = title
        return title

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
            if re.match(self._pdb_id_re, key) is None:
                return None  # key doesn't match PDB ID format
            # Returning a list indicates that there are still >1 options.
            return [key + '[alphanumeric]' * (4 - len(key)), filler]
        else:
            title = self._get_title_from_pdb_with_caching(key)
            # If we don't have a valid title return None to indicate that the
            # key is invalid.
            if title is None:
                return None
            # We've got a valid key. Show the title of that PDB entry to the
            # user and make the key the only possible autocomplete match.
            print(key + ':', title)
            return key


PDB_ID_AUTOCOMPLETE = [
    lambda: PdbIdAutocomplete(get_summary=lambda key: pdb.get_summary(key)),
    'PDB Entry Id', ''
]


@extendaa(PDB_ID_AUTOCOMPLETE)
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


@extendaa(PDB_ID_AUTOCOMPLETE)
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


@extendaa(PDB_ID_AUTOCOMPLETE)
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


@extendaa(PDB_ID_AUTOCOMPLETE)
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


@extendaa(cmd.auto_arg[0]['count_atoms'])
def count_chains(selection='(all)', state=0, quiet=False):
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
    n = len(chains)
    if not quiet:
        print('count_chains: %s chains' % n)
    return n


def initialize():
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
    initialize()
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
    initialize()
    mm_cif_file = None
    pdbid = None
    logging.debug(argv)
    pdb_id_re = re.compile(r'\d[A-z0-9]{3}')
    mm_cif_file_re = re.compile(r'mmCIF_file=(.*)')

    for arg in argv:
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
