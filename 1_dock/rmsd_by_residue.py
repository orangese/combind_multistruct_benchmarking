__doc__ = """
Calculates the RMSD for each residue of two or more identical or similar
proteins.

The number of residues in each structure must be the same length if -d is not
set. If structures are non-conformers, specify -d to identify similar residues
by distances. First structure is the reference. Multiple structures may be
specified as separate files or in one file. Output is a CSV file with the
RMSD values.

Copyright Schrodinger, LLC. All rights reserved.
"""

#- Imports -------------------------------------------------------------------

import os
import sys
from schrodinger import structure
from schrodinger.structutils import analyze, measure, rmsd
from schrodinger.protein import analysis
import schrodinger.utils.cmdline as cmdline
import csv
import math
import collections


#- Functions -----------------------------------------------------------------
def get_res_info(res):
    """
    Given a _Residue object, return a string representation of the residue.
    """
    return '%s(%s)' % (str(res), res.pdbres.strip())


def get_parser():
    """
    Returns schrodinger.utils.cmdline SingleDashOptionParser object with
    our application's flags and defaults defined.

    """
    usage = "$SCHRODINGER/run %prog -o <output.csv> <input1.mae> <input2.mae>"

    parser = cmdline.SingleDashOptionParser(
        usage=usage,
        description=__doc__,)
    parser.add_option(
        "-o",
        dest="output",
        default='rmsd_by_residue-out.csv',
        help="Output csv file.")
    parser.add_option(
        "-asl",
        dest="asl",
        default='all',
        help='ASL used to limit RMSD calculations. This ASL expression will '
        'limit the atoms in each residue that are used to calculate '
        'RMSD. For example if you want to calculate RMSD for just the '
        'backbone atoms the command line would be "-asl backbone". '
        'If you want the RMSD for sidechain atoms the commandline '
        'would be "-asl sidechain". Default: "all".')
    parser.add_option(
        "-d",
        dest="dist",
        type='float',
        help='Distance cutoff for C-alpha matching. If this is set, only '
        'the residues from <input2.mae> that are closest to each residue '
        'in <input1.mae> will be used (as long as the C-alphas are within '
        'the cutoff distance). If this is not set each residue will be '
        'looped over by residue order.')
    parser.add_option(
        "-t",
        dest="torsional",
        action='store_true',
        help='This will compute per residue torsional rmsd between the two '
        'structures provided.  For example, if a GLY is compared to a CYS then'
        ' only phi and psi angles will be used in the rmsd calculation.  For '
        ' example, if ALA is compared to CYS then only phi, psi, and chi1 will'
        ' be used in the rmsd calculation.  For example. if PHE is compared to'
        ' CYS then phi, psi, chi1, and chi2 will be used in the rmsd '
        'calculation.')
    parser.add_option(
        "-w",
        dest="write_bfactor_struct",
        type='string',
        help='Provide a filename to write a structure which will contain the '
        'resulting per residue torsional rmsd in the bfactor column scaled '
        'from the minimum 0 to the maximum 180. This can be only be used when '
        'the -t option is invoked.')
    parser.add_option(
        "-m",
        dest="write_matrix",
        action='store_true',
        help='Write the output CSV file as a matrix where columns are residues '
        'as named in the reference structure, and rows are structures (by '
        'default, there is a row for each residue pair interaction). In '
        'the matrix mode, only RMSD values are reported.')

    return parser


def calculate_in_place_torsion_rmsd(torsion_list1, torsion_list2):
    """
    This will calculate the torsional rmsd between the torsions in list1 and
    torsions in list2.  If torsions in list1 are missing then we skip that
    torsion and will not consider it.  For example torsion_list1 for a glycine
    would be [ '100', '23', '-', '-'] and if torsion list2 was for a lysine it
    would be [ '120', '12', '122', '-133']. In this case only 2 torsions phi
    and psi will be used to compute rmsd. The variable valid_torsion_count is
    used to count how many torsions we have to compare.  We can easily extend
    this to compare chi angles > 2.
    @param torsion_list1: The list of torsions for the first structure for a 
                          specific residue
    @type torsion_list1: list
    @param torsion_list2: The list of torsions for the second structure for a
                          specific residue
    @type torsion_list2: list
    @return : Returns the root mean square difference between all torsions that
              could be matched in the provided lists
    @rtype : float
    """
    torsions = zip(torsion_list1, torsion_list2)

    # Do the calculation with the mapping
    dist_squared = 0.0
    valid_torsion_count = 0
    for ref_index, test_index in torsions:
        try:
            # The max difference between two residues will be
            # 180 degrees 
            # If we had four dihedrals to compare i.e. phi/psi/chi1/chi2
            # the max would be sqrt(180*180*4/4) = 180
            # The min is 0
            d = ((ref_index - test_index + 180) % 360) - 180
        except TypeError:
            # This will be catching difference of two residues that do
            # not have the same number of dihedrals.  We just want to
            # pass and not increment the counter b/c we should at least
            # be able to compare phi/psi for GLY and phi/psi/chi1 for ALA
            pass
        else:
            dist_squared += d * d
            valid_torsion_count += 1

    if valid_torsion_count:
        res_rmsd = math.sqrt(dist_squared / float(valid_torsion_count))
    else:
        res_rmsd = 0.0

    return res_rmsd


def get_matched_residue_data(st1, st2, cutoff):
    """
    Gets all residues from st1 and st2 that have C-alphas within
    the supplied cutoff distance.

    """
    st1_indices = analyze.evaluate_asl(st1, '(atom.ptype " CA ")')
    st2_indices = analyze.evaluate_asl(st2, '(atom.ptype " CA ")')

    st1_residues = []
    st2_residues = []
    distance_data = []
    for i in st1_indices:
        st1_atom = st1.atom[i]
        # Set the minimum distance to the ignore cutoff
        distance = cutoff
        closest_atom = None

        # Find the closest atom from st2 to the st1_atom within dist criterion
        for j in st2_indices:
            st2_atom = st2.atom[j]
            new_dist = measure.measure_distance(st1_atom, st2_atom)
            if new_dist < distance:
                distance = new_dist
                closest_atom = st2_atom

        if closest_atom:
            st1_residues.append(st1_atom.getResidue())
            st2_residues.append(closest_atom.getResidue())
            distance_data.append(distance)

    return st1_residues, st2_residues, distance_data


def get_torsional_rmsd_data(dihedral_dict1, dihedral_dict2):
    """
    This will get the torsional rmsd for each residue calling the
    calculate_in_place_torsion_rmsd function for each residue in the residue
    lists.
    @param dihedral_dict1: Ordered dictionary of residue keys in the order the residues
                  were read from the structure
    @type L{collections.OrderedDict}
    @param dihedral_dict2: Ordered dictionary of residue keys in the order the residues
                  were read from the structure
    @type L{collections.OrderedDict}
    """

    rmsd_torsion_data = []
    for res1_info, res2_info in zip(dihedral_dict1, dihedral_dict2):
        #NOTE some rmsd's will be reported as 0.0049999999999954525 in the
        # example of 3MBX-mut.mae and 3MBX-orig.mae the following if statement
        # prints the dihedrals associated with one of these differences and
        # we see it is one decimal being rounded in the psi angle
        # [-137.32, 86.61, -168.98, -179.57] [-137.32, 86.6, -168.98, -179.57]

        torsion1 = dihedral_dict1[res1_info]
        torsion2 = dihedral_dict2[res2_info]
        if torsion1 is not None and torsion2 is not None:
            rmsd_torsions = calculate_in_place_torsion_rmsd(torsion1, torsion2)
            rmsd_torsion_data.append([res1_info, res2_info, rmsd_torsions])

    return rmsd_torsion_data


def get_residue_key(atomobj):
    """
    This will generate a key that is compatible with the protein.analysis
    dihedral dictionary.
    @param atomobj: An atomobject which we will use to extract the residue info
    @type atomobj: L{schrodinger.structure._StructureAtom}
    @return : Returns a key for the protein.analysis dictionary
    @rtype: string
    """
    resnum = atomobj.resnum
    reschain = atomobj.chain
    resinscode = atomobj.inscode
    pdbres = atomobj.pdbres
    key = "%s:%s%3i" % (reschain, pdbres, resnum)
    if resinscode != " ":
        key += ":%s" % resinscode

    return key


def get_ordered_dict(struct):
    """
    This will generate an ordered dictionary of the residues
    that match the CA atom type selection.
    @param struct: A structure to generate the residues dictionary
                   for using the CA atom selection
    @type struct: L{schrodinger.structure}
    @return : A ordered dictionary of residues
    @rtype : L{collections.orderedDict}
    """
    asl_CA = "atom.ptype CA"
    atom_list = analyze.evaluate_asl(struct, asl_CA)
    residues = collections.OrderedDict()
    for atom in atom_list:
        atomobj = struct.atom[atom]
        key = get_residue_key(atomobj)
        if key not in residues:
            resobj = atomobj.getResidue()
            residues[key] = get_res_info(resobj)

    return residues


def assign_dihedral_properties(struct):
    """
    This will build a dictionary of dihedral angles for each residue returning
    a list containing phi, psi, chi1, and chi2 which will be used to compute
    the per residue torsional rmsd between two structures. 
    @param struct: A structure object used to compute the torsionals of the 
                   residues.
    @type struct: L{schrodinger.structure}
    @return: An ordered dictionary containing lists of the dihedral angles
              as values and residue info strings as keys.
    @rtype: OrderedDict
    """
    residues = get_ordered_dict(struct)

    report_struct = analysis.Report(
        struct, sets_to_run=['SIDECHAIN DIHEDRALS', 'BACKBONE DIHEDRALS'])
    dihedral_dict = {}
    dihedrals = report_struct.get_set('BACKBONE DIHEDRALS')
    for point in dihedrals.points:
        try:
            key = point.descriptor
        except (ValueError, IndexError):
            print 'Skipping unknown residue'
        else:
            if key in residues:
                dihedral_dict[key] = point.values[:2]
            else:
                print "not in residues :, ", key

    dihedrals = report_struct.get_set('SIDECHAIN DIHEDRALS')
    for point in dihedrals.points:
        try:
            key = point.descriptor
        except (ValueError, IndexError):
            print 'Skipping unknown residue'
        else:
            if key in residues:
                dihedral_dict[key] = dihedral_dict[key] + point.values[:2]
            else:
                print "not in residues :, ", key

    final_dict = collections.OrderedDict()
    for res_key, res_info in residues.iteritems():
        final_dict[res_info] = dihedral_dict.get(res_key)

    return final_dict


def get_rmsd_data(st1, st2, opts_dist):
    """
    For each residue in st1, calculate the RMSD to the corresponding residue
    in st2.  If "opts_dist" is None, then residue equivalence is determined by
    the residue index in the CT, otherwise the distance method is used.
    """

    try:
        matching_atoms1 = set(analyze.evaluate_asl(st1, opts.asl))
        matching_atoms2 = set(analyze.evaluate_asl(st2, opts.asl))
    except:
        print 'Invalid ASL: %s' % opts.asl
        sys.exit(1)

    if opts_dist is not None:
        (residues1, residues2, distance_data) = get_matched_residue_data(
            st1, st2, opts_dist)
    else:
        residues1 = st1.residue
        residues2 = st2.residue
        # If there is no distance_data set all values to None for zip
        distance_data = [None] * len(residues1)

    rmsd_data = []
    for res1, res2, dist in zip(residues1, residues2, distance_data):
        at_list1 = [ai for ai in res1.getAtomIndices() if ai in matching_atoms1]
        at_list2 = [ai for ai in res2.getAtomIndices() if ai in matching_atoms2]

        if len(at_list1) == 0:
            # Skip reference structures not matching the ASL
            continue

        res1_info = get_res_info(res1)
        res2_info = get_res_info(res2)

        if len(at_list1) != len(at_list2):
            # If the atom list lengths don't match, then RMSD can't be calculated.
            # This will also happen if the res2 did not match the ASL.
            # We still need to put something into the rmsd_data list in order
            # for list lengths to be consistent between different structures.
            res_rmsd = None
        else:
            # Get the residue comparison
            res_rmsd = rmsd.calculate_in_place_rmsd(st1, at_list1, st2,
                                                    at_list2, True)

        if opts_dist is not None:
            rmsd_data.append([res1_info, res2_info, res_rmsd, dist])
        else:
            rmsd_data.append([res1_info, res2_info, res_rmsd])

    return rmsd_data


def calc_rmsd_for_st_pair(st1, st2, dist):
    """
    Calculate the per-residue RMSDs for the given structure pair.

    @type st1: structure.Structure
    @param st1: Reference structure.

    @type st2: structure.Structure
    @param st2: Comparison structure.

    @type dist: float
    @param dist: If specified, use the distance method instead of residue index
        method to determine similar residues between the structures. Residues
        farther than this distance are rejected.
    """

    rmsd_data = get_rmsd_data(st1, st2, opts.dist)

    if opts.torsional is not None:
        dihedral_dict1 = assign_dihedral_properties(st1)
        dihedral_dict2 = assign_dihedral_properties(st2)
        torsional_rmsd_data = get_torsional_rmsd_data(dihedral_dict1,
                                                      dihedral_dict2)
        for t_res1, t_res2, rmsd_tors in torsional_rmsd_data:
            for rmsd_entry in rmsd_data:
                res1, res2 = rmsd_entry[:2]
                if t_res1 == res1 and t_res2 == res2:
                    rmsd_entry.append(rmsd_tors)
                    break

    if opts.write_bfactor_struct and opts.torsional:
        for entry in torsional_rmsd_data:
            for res in st1.residue:
                res1_info = get_res_info(res)
                if res1_info == entry[0]:
                    # This is highly annoying there is a res.temperature_factor
                    # which gets the bfactor but not a set function!
                    for atom in res.getAtomIndices():
                        st1.atom[atom].temperature_factor = entry[2]
                    break
        st1.write(opts.write_bfactor_struct)
    return rmsd_data


if __name__ == "__main__":
    # Parse the command line, and check that the files are valid
    parser = get_parser()
    (opts, args) = parser.parse_args()

    for filename in args:
        if not os.path.isfile(filename):
            parser.error('Cannot find file: %s' % filename)

    ref_st = None
    ref_res_count = None
    rmsd_data_list = []

    # Iterate over all structures from the given files:
    for i, st in enumerate(structure.MultiFileStructureReader(args)):
        if i == 0:
            ref_st = st
            ref_res_count = len(st.residue)
            continue

        print 'Processing structure %i...' % i
        if opts.dist is None:
            if len(st.residue) != ref_res_count:
                parser.error(
                    'The length of residues do not match for the structure 1 '
                    'and structure %i' % i)

        rmsd_data = calc_rmsd_for_st_pair(ref_st, st, opts.dist)
        rmsd_data_list.append(rmsd_data)

    if not rmsd_data_list:
        parser.error('There must be at least 2 input structures')

    # Write the result to the CSV file:
    with open(opts.output, "wb") as ofile:
        writer = csv.writer(ofile)
        print 'Number of structure comparisons:', len(rmsd_data_list)
        print 'ASL filter: %s' % opts.asl
        if opts.dist is not None:
            print 'Cutoff for C-alphas: %f' % opts.dist

        if opts.write_matrix:
            # Generate a list of residue names in the reference structure:
            ref_residues = [data[0] for data in rmsd_data_list[0]]
            # Write this as the header row:
            writer.writerow([''] + ref_residues)

            # Each next row will be a list of RMSD for the comparison structures:
            for i, rmsd_data in enumerate(rmsd_data_list, start=1):
                row = ["Structure%i" % i]
                for res_data in rmsd_data:
                    rmsd = res_data[2]
                    if rmsd is None:
                        rmsd = ''
                    else:
                        rmsd = "%.4f" % rmsd
                    row.append(rmsd)
                writer.writerow(row)
        else:
            header = ['Residue 1', 'Residue 2', 'RMSD']
            if opts.dist is not None:
                header.append('Distance between C-alphas')
                if opts.torsional:
                    header.append('TORSIONAL RMSD')
            writer.writerow(header)

            for i, rmsd_data in enumerate(rmsd_data_list):
                if i != 0:
                    # Separator row:
                    writer.writerow([])

                for res_data in rmsd_data:
                    res1, res2, rmsd = res_data[:3]
                    if rmsd is None:
                        continue
                    row = [res1, res2, "%.4f" % rmsd]

                    if opts.dist is not None:
                        row.append("%.2f" % res_data[3])
                        if opts.torsional:
                            if len(res_data) == 5:
                                tors_rmsd = "%.4f" % res_data[4]
                            else:
                                tors_rmsd = ''
                            row.append(tors_rmsd)
                    writer.writerow(row)

        print "Written:", opts.output
