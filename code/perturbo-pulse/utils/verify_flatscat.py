#!/bin/python3

import argparse
import pprint

from typing import Any, Dict, List, Tuple


def idx2band(index: int, num_bands: int) -> Tuple[int, int, int]:
    '''
    See boltz_utils.f90, function idx2band() for details.
    Note that Fortran indexes are 1-based, but Python indexes are 0-based.
    '''
    # mu <= nmodes, nk <= nbnd, mkq <= nbnd
    # index = (mu-1)*nbnd*nbnd + (nk-1)*nbnd + (mkq-1)
    # band(1) = mod( idx, nbnd ) + 1
    # band(2) = mod( idx/nbnd, nbnd ) + 1
    # band(3) = idx / (nbnd*nbnd) + 1
    return (index % num_bands + 1,
            (index // num_bands) % num_bands + 1,
            (index // (num_bands * num_bands)) + 1)


def band2idx(bands: Tuple[int, int, int], num_bands: int) -> int:
    '''
    See boltz_utils.f90, function band2idx() for details.
    Note that Fortran indexes are 1-based, but Python indexes are 0-based.
    '''
    # index = (mu-1)*nbnd*nbnd + (nk-1)*nbnd + (mkq-1)
    # idx = ( band(3)-1 )*nbnd*nbnd + ( band(2)-1 )*nbnd + ( band(1)-1 )
    return ((bands[0] - 1) +
            (bands[1] - 1) * num_bands +
            (bands[2] - 1) * num_bands * num_bands)


def make_rowdata(keynames: List[str], values: List[Any]) -> Dict[str, Any]:
    rowdata = {}
    for (k, v) in zip(keynames, values):
        rowdata[k] = v

    return rowdata


def load_scat(path):
    data = []
    with open(path) as f:
        line = f.readline()
        col_names = line.strip().split()

        while True:
            # These are the scalar values in the scat array entry
            line = f.readline()
            if not line:
                break # All done!

            values = line.strip().split()
            values = [int(v) for v in values]

            rowdata = make_rowdata(col_names, values)

            # Band indexes
            line = f.readline()
            if not line:
                raise IOError('Expected row of data containing band indexes')
            
            values = line.strip().split()
            if values[0] != 'bands_index':
                raise IOError('Expected row to start with "bands_index" string.')
            
            rowdata['bands_index'] = [int(v) for v in values[1:]]

            # eph_g2 values
            line = f.readline()
            if not line:
                raise IOError('Expected row of data containing eph_g2 values')
            
            values = line.strip().split()
            if values[0] != 'eph_g2':
                raise IOError('Expected row to start with "eph_g2" string.')
            
            rowdata['eph_g2'] = [float(v) for v in values[1:]]

            # Done with this row
            data.append(rowdata)
    
    return data


def load_flatscat(path):
    data = []
    with open(path) as f:
        line = f.readline()
        col_names = line.strip().split()

        while True:
            line = f.readline()
            if not line:
                break # All done!

            values = line.strip().split()
            values = [int(v) for v in values[:-1]] + [float(values[-1])]
            data.append(make_rowdata(col_names, values))

    return data


def load_scat_targets(path):
    data = []
    with open(path) as f:
        line = f.readline()
        col_names = line.strip().split()

        while True:
            # These are the scalar values in the scat array entry
            line = f.readline()
            if not line:
                break # All done!

            values = line.strip().split()
            values = [int(v) for v in values]

            rowdata = make_rowdata(col_names, values)

            # Source indexes (indexed against flatscat)
            line = f.readline()
            if not line:
                raise IOError('Expected row of data containing source indexes')
            
            values = line.strip().split()
            if values[0] != 'sources':
                raise IOError('Expected row to start with "sources" string.')
            
            rowdata['sources'] = [int(v) for v in values[1:]]

            # Done with this row
            data.append(rowdata)
    
    return data


def verify_scat_internal(scat_data):
    issues: List[str] = []
    i_row = 1
    for rowdata in scat_data:
        if len(rowdata['bands_index']) != rowdata['nchl']:
            issues.append(f'Row {i_row}:  nchl ({rowdata["nchl"]}) doesn\'t match number of elements in bands_index ({len(rowdata["bands_index"])})')
        
        if len(rowdata['eph_g2']) != rowdata['nchl']:
            issues.append(f'Row {i_row}:  nchl ({rowdata["nchl"]}) doesn\'t match number of elements in eph_g2 ({len(rowdata["eph_g2"])})')
        
        # Note that each scat entry occupies 3 rows.
        i_row += 3

    return issues


def verify_scat_targets_internal(scat_tgts_data):
    issues: List[str] = []
    i_row = 1
    for rowdata in scat_tgts_data:
        if len(rowdata['sources']) != rowdata['count']:
            issues.append(f'Row {i_row}:  count ({rowdata["count"]}) doesn\'t match number of elements in sources ({len(rowdata["sources"])})')

        # Note that each scat_targets entry occupies 2 rows.
        i_row += 2

    return issues


def verify_scat_flatscat(scat_data, flatscat_data, max_issues=20):
    issues: List[str] = []

    # We need to assume the scat and flatscat datasets are in the same order;
    # otherwise it's just prohibitavely expensive to compare the two datasets

    i_flatscat = 0
    for i_scat in range(len(scat_data)):
        scat_row = scat_data[i_scat]
        for i_chl in range(scat_row['nchl']):
            if i_flatscat >= len(flatscat_data):
                issues.append(f'Reached end of flatscat (row {i_flatscat+1}) before reaching end of scat (row {i_scat+1}, channel {i_chl+1}); ending test.')
                return issues

            flatscat_row = flatscat_data[i_flatscat]

            (scat_m, scat_n, scat_im) = idx2band(scat_row['bands_index'][i_chl], scat_row['kg_numb'])
            scat_eph_g2 = scat_row['eph_g2'][i_chl]

            if (scat_row['ik']  != flatscat_row['ik'] or
                scat_row['ikq'] != flatscat_row['ikq'] or
                scat_row['iq']  != flatscat_row['iq'] or
                scat_m          != flatscat_row['m'] or
                scat_n          != flatscat_row['n'] or
                scat_im         != flatscat_row['im'] or
                scat_eph_g2     != flatscat_row['eph_g2']):
                issues.append(f'Mismatch between scat row {i_scat+1} channel {i_chl+1} and flatscat row {i_flatscat+1}\n' +
                    'scat=' + pprint.pformat(scat_row) + '\nflatscat=' + pprint.pformat(flatscat_row)
                )
            
            if len(issues) >= max_issues:
                issues.append(f'Reached max-issue limit of {max_issues}; ending test.')
                return issues
            
            i_flatscat += 1

    if i_flatscat < len(flatscat_data):
        issues.append(f'Found {len(flatscat_data) - i_flatscat} extra rows at end of flatscat at end of scat dataset.')

    return issues


def verify_flatscat_scattgts(flatscat_data, scat_tgts_data):
    issues: List[str] = []

    # We need to assume the scat and flatscat datasets are in the same order;
    # otherwise it's just prohibitavely expensive to compare the two datasets

    src_indexes = set()
    num_flatscat = len(flatscat_data)

    for i_tgt in range(len(scat_tgts_data)):
        tgt_row = scat_tgts_data[i_tgt]
        tgt_i1 = tgt_row['i1']
        tgt_i2 = tgt_row['i2']

        for i_src in tgt_row['sources']:
            if i_src in src_indexes:
                issues.append(f'Source-index {i_src} appears multiple times')
            
            # Note that Fortran indexing starts at 1, but Python indexing starts at 0.
            src_indexes.add(i_src)

            if i_src < 0:
                if -i_src > num_flatscat:
                    issues.append(f'Source-index {i_src} is out of range for {num_flatscat} flatscat rows')

                if (flatscat_data[-i_src - 1]['m'] != tgt_i1 or flatscat_data[-i_src - 1]['ikq'] != tgt_i2):
                    issues.append(f'flatscat row {-i_src} doesn\'t specify target-element ({tgt_i1},{tgt_i2}) for (m,ikq)')
            
            elif i_src > 0:
                if i_src > num_flatscat:
                    issues.append(f'Source-index {i_src} is out of range for {num_flatscat} flatscat rows')

                if (flatscat_data[i_src - 1]['n'] != tgt_i1 or flatscat_data[i_src - 1]['ik'] != tgt_i2):
                    issues.append(f'flatscat row {i_src} doesn\'t specify target-element ({tgt_i1},{tgt_i2}) for (n,ik)')
            
            else:
                issues.append('Source-index of {i_src} encountered!')
    
    for i in range(num_flatscat):
        if (i + 1) not in src_indexes:
            issues.append(f'scat_targets data didn\'t specify an index {i+1}')

        if -(i + 1) not in src_indexes:
            issues.append(f'scat_targets data didn\'t specify an index {-(i+1)}')

    return issues


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Verify data-dumps of internal Perturbo data structures.')

    parser.add_argument('--scat', dest='scat_filepath', nargs='?',
        help='Specify path and filename of scat data-dump file')

    parser.add_argument('--flatscat', dest='flatscat_filepath', nargs='?',
        help='Specify path and filename of scat data-dump file')

    parser.add_argument('--scat_tgts', dest='scat_tgts_filepath', nargs='?',
        help='Specify path and filename of scat_targets data-dump file')

    args = parser.parse_args()

    scat_data = None    
    if args.scat_filepath:
        print(f'Loading scat data-structure contents from path:  {args.scat_filepath}')
        scat_data = load_scat(args.scat_filepath)
    else:
        print(f'No scat data file specified.')

    flatscat_data = None
    if args.flatscat_filepath:
        print(f'Loading flatscat data-structure contents from path:  {args.flatscat_filepath}')
        flatscat_data = load_flatscat(args.flatscat_filepath)
    else:
        print(f'No flatscat data file specified.')

    scat_tgts_data = None    
    if args.scat_tgts_filepath:
        print(f'Loading scat_targets data-structure contents from path:  {args.scat_tgts_filepath}')
        scat_tgts_data = load_scat_targets(args.scat_tgts_filepath)
    else:
        print(f'No scat data file specified.')

    issues_found = False

    if scat_data is not None:
        print('Verifying internal consistency of scat data.')
        issues = verify_scat_internal(scat_data)
        if issues:
            issues_found = True
            print('Issues found:')
            for iss in issues:
                print(f' * {iss}')

    if scat_tgts_data is not None:
        print('Verifying internal consistency of scat_targets data.')
        issues = verify_scat_targets_internal(scat_tgts_data)
        if issues:
            issues_found = True
            print('Issues found:')
            for iss in issues:
                print(f' * {iss}')

    if scat_data is not None and flatscat_data is not None:
        print('Verifying consistency between scat and flatscat data.')
        issues = verify_scat_flatscat(scat_data, flatscat_data)
        if issues:
            issues_found = True
            print('Issues found:')
            for iss in issues:
                print(f' * {iss}')

    if flatscat_data is not None and scat_tgts_data is not None:
        print('Verifying consistency between flatscat and scat_targets data.')
        issues = verify_flatscat_scattgts(flatscat_data, scat_tgts_data)
        if issues:
            issues_found = True
            print('Issues found:')
            for iss in issues:
                print(f' * {iss}')

    if not issues_found:
        print('No issues found, hooray!')
