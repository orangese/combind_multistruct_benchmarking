__doc__ = """
Renumber residues based on a structural alignment to a template structure.

A script which renumbers residues based on a SKA structural alignment. In
"reference" mode, all residue get assigned the numbering of the first structure.
In "alignment" mode, all structures are renumbered based on a residue's position
in a multiple sequence alignment. If the -nosuper option is specified, the
coordinates are left unchanged, otherwise the structures returned are
superimposed using SKA.

NOTE:
This script changes non-standard amino-acid types to standard ones for SKA
compatibility e.g. HIE, HID, HIP to HIS, and CYX to CYS.

'alignment' mode is limited to 52 input structures and the renumbering is very
dependant on the sequences used (i.e. count, their length etc.).

Copyright Schrodinger, LLC. All rights reserved.
"""

###############################################################################
# Globals, constants and imports

import os
import sys
import argparse
import string

from schrodinger import structure
from schrodinger.structutils import analyze
from schrodinger.utils import cmdline, fileutils, log
from schrodinger.application import ska
from schrodinger.ui.sequencealignment import globals as seqdata

SKA_KEYS = {'MODE': 'align', 'ORDER': 'seed', 'USE_AUTOMATIC_SETTINGS': 'yes'}

INSCODES = [' '] + list(string.ascii_uppercase)
GAP_CHARS = ("-", ".")
REF_KEY = '1'

logger = log.get_output_logger("adjust_residue_numbering.py")
logger.setLevel(log.INFO)


###############################################################################
def get_initial_structures_from_file(infile, numstructs, template_id):
    """
    Get the structures that will be renumbered from a file, check to make sure
    they contain protein atoms and are standardized for ska. Return a template
    structure and a list of query structure objects.

    @param infile:  Input Maestro file.
    @type infile:  str

    @param numstructs:  Number of structures to use.
    @type numstructs:  int

    @param template_id:  Structure number to use as the template.
    @type template_id:  int

    @return: Template structure, list of query structures.
    @rtype: L{Structure<structure.Structure>}, list of
    L{Structure<structure.Structure>}.
    """

    template_st = None
    initial_structures = []

    for i, st in enumerate(structure.StructureReader(infile), start=1):
        if not analyze.evaluate_asl(st, '(protein)'):
            logger.error("Structure number %d (%s) does not contain any "
                         "protein atoms." % (i, st.title))
            sys.exit(1)

        st = standardize(st)

        if i == template_id:
            template_st = st
        else:
            initial_structures.append(st)

    if initial_structures:
        structure_count = len(initial_structures)
        if numstructs is None:
            numstructs = structure_count

        if numstructs < structure_count:
            initial_structures = initial_structures[:numstructs]
        elif numstructs > structure_count:
            logger.warning("Number of structures to be imported is larger "
                           "than the number of query structures in the input "
                           "file %s. All structures will be used." % infile)
    elif not initial_structures:
        logger.error("Insufficient number of structures in %s." % infile)
        sys.exit(1)

    # Make sure template_id is > 0 and < numstructs
    if template_st is None:
        logger.error("Template structure cannot be found.")
        sys.exit(1)

    return template_st, initial_structures


###############################################################################
def standardize(new_ct):
    """
    Rename protein residues so that the structure is compatible with ska.

    @param new_ct:  Original structure.
    @type new_ct:  L{Structure<structure.Structure>}

    @return:  Standardized structure for ska.
    @rtype:  L{Structure<structure.Structure>}
    """

    for res in new_ct.residue:

        # Derive a standard name from 1-letter code:
        rcode = res.getCode()

        if rcode != 'X':
            res.pdbres = seqdata.AMINO_ACIDS[rcode][0]

    return new_ct


###############################################################################
def renumber_non_std_residues(st, starting_resnum):
    """
    Renumber all non-standard residues in the given structure consecutively.

    @type st: L{structure.Structure}
    @param st: Structure to renumber

    @type starting_resnum: int
    @param starting_resnum: Residue number to start renumbering from.
    """
    new_resnum = starting_resnum
    for res in st.residue:
        if res.getCode() == 'X':
            res.resnum = new_resnum
            res.inscode = ' '
            new_resnum += 1


###############################################################################
def renumber_st_by_alignment(st, alignment):
    """
    Renumber the residues in the given structure based on the specified
    alignment sequence. So the first residue in the sequence will be
    numbered 1, second 2, etc.  Residue number is incremented with
    each gap in the sequence.

    @type st: L{structure.Structure}
    @param st: Structure to renumber

    @type alignment: str
    @param alignment: Sequence by which to renumber.
    """

    query_prop_st = st.copy()

    st = ska._standardize(st, reorder=True, rename=False)

    for key, value in query_prop_st.property.iteritems():
        st.property[key] = value

    std_residues = (res for res in structure.get_residues_unsorted(st)
                    if res.getCode() != 'X')

    prev_res_str = None
    new_resnum = 0
    for code in alignment:
        new_resnum += 1
        if code in GAP_CHARS:
            # There is no residue at this position
            continue

        res = std_residues.next()
        while str(res) == prev_res_str:
            # This residue may have been broken into 2, as in GLU163 in 1uy9
            # Use the same rensum/inscode
            res.resnum = new_resnum
            res.inscode = ' '
            res = std_residues.next()
        # The next residue found
        prev_res_str = str(res)
        assert res.getCode() == code

        res.resnum = new_resnum
        res.inscode = ' '

    # FIXME Why are all non-stadard residues numbered this way? This makes
    # sense for ligands, but not for residues within the protein.
    renumber_non_std_residues(st, new_resnum + 1)

    return st


###############################################################################
def renumber_by_alignment(ska_data, final_sts):
    """
    Renumber the residues in structures, based on the aligned sequences.
    Residues will get assigned resnums according to their position in the
    alignment.

    @param ska_data:  SKA data object
    @type ska_data:  ska._SkaData

    @param final_sts:  Dictionary of int key:  L{Structure<structure.Structure>}
    @type final_sts:  dict

    @return:  List of modified structures.
    @rtype:  list

    @return:  Renumbered structure.
    @rtype:  L{structure.Structure}
    """

    out_sts = []

    for key in sorted(ska_data.align.keys(), key=int):
        alignment = ska_data.align[key]
        st = final_sts[key]
        logger.debug("Structure title: %s" % st.title)
        logger.debug("Structure alignment: %s" % alignment)
        renum_st = renumber_st_by_alignment(st, alignment)
        logger.info('Renumbered structure "%s"' % st.title)
        out_sts.append(renum_st)

    logger.info("Renumbered %i structures" % len(out_sts))

    return out_sts


###############################################################################
def renumber_query_st(template_st, query_st, template_seq, query_seq):
    """
    Renumber the residues in <query_st> based on the positions of the residues
    in the template.

    @param template_st:  Template structure.
    @type template_st:  L{structure.Structure}

    @param query_st:  Query structure.
    @type query_st:  L{structure.Structure}

    @param template_seq:  Template sequence.
    @type template_seq:  str

    @param query_seq:  Query sequence.
    @type query_seq:  str

    @return:  Renumbered query structure.
    @rtype:  L{structure.Structure}
    """

    # Ensure that the residue order in the structures is consistent with the
    # alignment SKA produces.

    query_prop_st = query_st.copy()

    query_st = ska._standardize(query_st, reorder=True, rename=False)
    template_st = ska._standardize(template_st, reorder=True, rename=False)

    for key, value in query_prop_st.property.iteritems():
        query_st.property[key] = value

    # Determine what residue number to start the assignment with. It will
    # be set to the first residue number in the template minus the number
    # of leading gaps in the template.

    n_term_gaps = 0

    for seq in template_seq:
        if seq in GAP_CHARS:
            n_term_gaps += 1
        else:
            break

    template_residues = [
        res for res in structure.get_residues_unsorted(template_st)
        if res.getCode() != 'X'
    ]

    query_residues = (res for res in structure.get_residues_unsorted(query_st)
                      if res.getCode() != 'X')

    first_st_res_num = template_residues[0].resnum
    start_res_num = first_st_res_num - n_term_gaps

    prev_query_res_str = None

    new_resnum = start_res_num
    new_ins = ' '

    seq_i = -1
    template_i = 0

    for template_aa, query_aa in zip(template_seq, query_seq):
        seq_i += 1
        # Find the template residue corresponding to this seq position:
        if template_aa not in GAP_CHARS:
            template_res = template_residues[template_i]
            assert template_res.getCode() == template_aa
            template_i += 1
        else:
            template_res = None

        # Find the query residue corresponding to this seq position:
        if query_aa not in GAP_CHARS:
            query_res = query_residues.next()
            while str(query_res) == prev_query_res_str:
                # This residue may have been broken into 2, as in GLU163 in 1uy9
                # Use the same rensum/inscode
                query_res.resnum = new_resnum
                query_res.inscode = new_ins
                query_res = query_residues.next()
            # The next query residue found
            prev_query_res_str = str(query_res)
            if query_res.getCode() != query_aa:
                query_res.applyStyle()
                print 'ERROR: Code for position %i (residue %s %s) in query CT is "%s"; sequence has "%s"' % (
                    seq_i, prev_query_res_str, query_res.pdbres.strip(),
                    query_res.getCode(), query_aa)
                sys.exit(1)
        else:
            # There is no query residue at this position
            continue

        if seq_i < n_term_gaps:
            # If this query residue is at the beginning of the sequence, and
            # the template sequence has not yet started, base the residue
            # numbering on start_res_num.
            new_resnum = start_res_num + seq_i
            new_ins = ' '
        elif template_res is not None:
            # This query residue has an equivalent template residue
            new_resnum = template_res.resnum
            new_ins = template_res.inscode
        else:
            if all((code in GAP_CHARS for code in template_seq[seq_i:])):
                # If we reached the end of the template sequence, then keep
                # incrementing the residue number:
                new_resnum += 1
                new_ins = ' '
            else:
                # Template has a gap here; we need to generate modified residue
                # numbers (between the current residue and the next residue
                # number) by adding insertion codes.
                try:
                    new_ins = INSCODES[INSCODES.index(new_ins) + 1]
                except IndexError:
                    raise RuntimeError("Exhausted 26 insertion codes.")

        query_res.resnum = new_resnum
        query_res.inscode = new_ins

    # Now renumber the non-standard residues:
    # FIXME Why are all non-stadard residues numbered this way? This makes
    # sense for ligands, but not for residues within the protein.
    renumber_non_std_residues(query_st, new_resnum + 1)

    return query_st


###############################################################################
def renumber_by_reference(ska_data, final_sts):
    """
    Renumber the residue in structures, based on the residue numbering in the
    reference structure. "Equivalent" residues in each structure will get the
    same residue number.

    @param ska_data:  SKA data
    @type ska_data:  list

    @param final_sts:  Dictionary of int key:  L{Structure<structure.Structure>}
    @type final_sts:  dict

    @return:  List of modified structures.
    @rtype:  list
    """

    # Template CT is the same for all alignments; but sequences are different,
    # since pairwise alignments were done.
    template_st = final_sts[REF_KEY]
    out_sts = [template_st]

    logger.info("Template (not renumbering): %s" % template_st.title)

    for i, data in enumerate(ska_data):
        if data.align:
            assert len(data.align) == 2
            for key, _seq in data.align.iteritems():
                if key is REF_KEY:
                    template_seq = _seq
                else:
                    query_seq = _seq
                    query_st = final_sts[key]

            # Set the starting residue numbering to begin a the first residue
            # number minus the number of gaps preceding it in the alignment.

            logger.info("Processing query: %s" % query_st.title)
            logger.debug(" Template sequence: %s" % template_seq)
            logger.debug("    Query sequence: %s" % query_seq)
            try:
                query_st = renumber_query_st(template_st, query_st,
                                             template_seq, query_seq)
            except RuntimeError:
                logger.warning(
                    'WARNING: This query sequence has gaps greater '
                    'than 26 residues, when aligned to the template. Verify '
                    'that the structures are sufficiently similar.')
                logger.warning('         The structure will be skipped.')
                continue

            out_sts.append(query_st)

    logger.info("Renumbered %i of %i queries" % (len(out_sts) - 1,
                                                 len(ska_data)))

    return out_sts


###############################################################################
def parse_args():
    """
    Parse the command line options.

    @return:  All script arguments and options.
    @rtype:  class:`argparse.Namespace
    """

    parser = cmdline.create_argument_parser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("infile", help="Input Maestro file.")

    parser.add_argument("outfile", help="Output Maestro file.")

    parser.add_argument(
        "-number",
        type=int,
        help="Number of structures to be imported. If no value"
        " is given, all structures will be used.")

    parser.add_argument(
        "-template_id",
        type=int,
        default=1,
        help="The number of the structure to be used as the "
        "reference structure. If no value is given, the "
        "first structure in the file will be used.")

    parser.add_argument(
        "-renumber",
        choices=['reference', 'alignment'],
        default="reference",
        help="Renumbering scheme. In 'reference' mode, "
        "query structures are renumbered using the residue"
        " numbering in the reference structure. In "
        "'alignment' mode', all structures are renumbered "
        "based on a residue's position in a multiple "
        "sequence alignment. Default is 'reference'.")

    parser.add_argument(
        "-nosuper",
        action='store_true',
        default=False,
        help="Don't change the coordinates when renumbering.")

    parser.add_argument(
        "-debug", action='store_true', default=False, help="Enable debugging.")

    args = parser.parse_args()

    if not os.path.isfile(args.infile):
        print "\nInput file %s cannot be found.\n" % args.infile
        sys.exit(1)

    if not fileutils.is_maestro_file(args.infile):
        print "\nInput file %s does not appear to be in Maestro format.\n" % \
              args.infile
        sys.exit(1)

    if not fileutils.is_maestro_file(args.outfile):
        print "\nOutput file %s does not appear to be in Maestro format.\n" % \
              args.outfile
        sys.exit(1)

    return args


##############################################################################
def main():
    """
    Main body of the script.
    """

    mod_st_list = None

    structure_dict = {}

    cmd_args = parse_args()
    if cmd_args.debug:
        logger.setLevel(log.DEBUG)

    # Parse the input file and get the template structure and a list of query
    # structures.
    template, query_sts = get_initial_structures_from_file(
        cmd_args.infile, cmd_args.number, cmd_args.template_id)

    # Create tuples required for ska.multiple_align_ct or ska.pairwise_align_ct
    ska_template = (REF_KEY, template)
    ska_queries = [(str(j), st) for j, st in enumerate(query_sts, start=2)]

    # Create a list of structures with their original coordinates.
    safe_tuple = [(ska_tup[0], ska_tup[1].copy()) for ska_tup in ska_queries]

    # Determine the renumbering mode the alignment is running in.
    if cmd_args.renumber == 'reference':
        ska_data = ska.pairwise_align_ct(
            ska_template,
            ska_queries,
            keywords=SKA_KEYS,
            debug=cmd_args.debug,
            reorder=True)
    else:

        try:
            ska_data = ska.multiple_align_ct(
                ska_template,
                ska_queries,
                keywords=SKA_KEYS,
                debug=cmd_args.debug,
                reorder=True)
        except Exception as msg:
            print msg
            print "Multiple sequence alignment cannot be generated as one or " \
                  "more structures failed to align."
            sys.exit(1)

    structure_dict[ska_template[0]] = ska_template[1]

    # Create a list of structures either with the original coordinates or with
    # the newly aligned coordinates.
    if cmd_args.nosuper:
        for st_number, orig_st in safe_tuple:
            structure_dict[st_number] = orig_st
    else:
        for st_number, moved_st in ska_queries:
            structure_dict[st_number] = moved_st

    # Renumber the structures according to the user input.
    if ska_data is not None:
        if cmd_args.renumber == 'reference':
            logger.debug("Renumbering by reference")
            mod_st_list = renumber_by_reference(ska_data, structure_dict)
        if cmd_args.renumber == 'alignment':
            logger.debug("Renumbering by alignment")
            mod_st_list = renumber_by_alignment(ska_data, structure_dict)

    if mod_st_list:
        writer = structure.StructureWriter(cmd_args.outfile)

        for mod_st in mod_st_list:
            writer.append(mod_st)
        writer.close()
    else:
        logger.error("No structures returned.")
        sys.exit(1)


if __name__ == '__main__':
    cmdline.main_wrapper(main)
