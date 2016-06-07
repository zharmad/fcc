#!/usr/bin/env python
# -*- coding: UTF-8  -*-

"""
Calculates a matrix of fraction of common contacts between two or more structures.

Authors:
        RODRIGUES Joao
        TRELLET Mikael
        MELQUIOND Adrien
        CHEN Po-chia
"""

def parse_chain_symmetries(chsymmstr):
    """
    Chain symmetry parsing routine.
    Begin with a string describing the symmetry groups and methof of permutation.
    Example for a dimer, and cyclic trimer: "AB,cDEF"
    This function returns two list of lists:
    (1) the group(s) interpreted from the string, and
    (2) the list of all permutations for each sets.
    
    Note that the first permutation must be identity for compatibility.
    
    ~~TODO~~ Allow users to specify explicit permutations, e.g. "(ABCD,CDAB)"
    """
    import itertools

    # Parse string
    grps=[] ; ll=[]
    cyclic=[]; bCyclic=False
    #bExplicit=False
    for char in chsymmstr:
        if char == ',':
            #end group
            if len(ll)>0:
                grps.append(ll)
                cyclic.append(bCyclic)
            ll=[]
            bCyclic=False
        elif char == 'c':
            #specify cyclic for this group.
            bCyclic=True
        else:
            #add member to group
            asc=ord(char)-64
            ll.append(asc)
    if len(ll)>0:
        grps.append(ll)
        cyclic.append(bCyclic)

    #Generate permutations.
    #First generate permutations from the chain symmetries above.
    permutations=[]
    for grp, cyclic in zip(grps,cyclic):
        if cyclic:
            p = [ [grp[i-j] for i in range(len(grp))] for j in range(len(grp)) ]
        else:
            p = itertools.permutations(grp)
        ll = [i for i in iter(p)]
        permutations.append(ll)

    return grps, permutations

def map_segid(seg, in_map, out_map):
    """
    Mapping function that is used instead of, say, map()
    """
    for i in xrange(len(in_map)):
        if seg==in_map[i]:
            return out_map[i]
    return seg

def map_contacts(contacts, in_map, out_map):
    """
    Interpret the integer string containing information about the
    resid and segid of the contact pairs.
    Note: this must synchronise with the C code as well as others,
    Chosen C syntax is: %02d%02d%05d%05d
    """
    #~~CONTACT_SYNTAX~~
    #n_digits = int(math.log(contacts[i],10)). First 4 are segid digits.
    seg1_offset=int(1e12)
    seg2_offset=int(1e10)
    res1_offset=int(1e5)
    res2_offset=int(1e0)
    res_size=int(1e5)
    seg_size=int(1e2)

    l=[]
    for frame in contacts:
        ll=[]
        for pair in frame:
            # Here is the %02d magic with int.
            seg1=(pair%(seg1_offset*seg_size))/seg1_offset
            seg2=(pair%(seg2_offset*seg_size))/seg2_offset
            nseg1=map_segid(seg1, in_map, out_map)
            nseg2=map_segid(seg2, in_map, out_map)
            if nseg2 > nseg1:
                npair=pair + (nseg1-seg1)*seg1_offset + (nseg2-seg2)*seg2_offset
            else:
                #Keep the sorting such that ab==ba.
                npair= pair%(res2_offset*res_size)/res2_offset*res1_offset + \
                       pair%(res1_offset*res_size)/res1_offset*res2_offset + \
                       nseg2*seg1_offset + \
                       nseg1*seg2_offset

            ll.append(npair)
            #if pair != npair:
            #    print in_map, ", n:", len(in_map)
            #    print out_map, ", n:", len(out_map)
            #    print map_segid(seg1, in_map, out_map)
            #    print "%d %d %d %d" % (seg1, nseg1, seg2, nseg2)
            #    print "%d -> %d" % (pair, npair)
            #print "-- frame 1st & last:", ll[0], ll[-1]
        l.append(set(ll))
    return l

def generate_contact_symmetries(contacts, permutations):
    """
    Generate all allowed copies of contacts, based on the given list of permutations.
    """
    import itertools
    # Keep contacts as a 1D list to simplify matters.
    # Extract the null transformation.
    chsymm=[ permutations[i][0] for i in xrange(len(permutations)) ]
    in_flat = [i for sublist in chsymm for i in sublist]
    mapped_contacts = []
    map_list = []
    for out_tmp in itertools.product(*permutations):
        out_flat=[i for sublist in out_tmp for i in sublist]
        #print in_flat, "-->", out_flat
        # Keep contacts as a 1D list to match base code.
        mapped_contacts.extend(map_contacts(contacts, in_flat, out_flat))
        # Keep map_list as 2D, len(map_list) contains the number of images.
        map_list.append(out_flat)

    return mapped_contacts, map_list

# Contact Parsing routines
def parse_contact_file(f_list, ignore_chain):
    """Parses a list of contact files."""
    
    if ignore_chain:
        # ~~CONTACT_SYNTAX~~ Altered from old code
        #contacts = [ [ int(l[0:5]+l[6:-1]) for l in open(f)] for f in f_list if f.strip()]
        contacts = [ [ int(l[4:-1]) for l in open(f)] for f in f_list if f.strip()]
    else:
        contacts = [ set([ int(l) for l in open(f)]) for f in f_list if f.strip()]

    return contacts

# FCC Calculation Routine
def calculate_fcc(listA, listB):
    """
    Calculates the fraction of common elements between two lists
    taking into account chain IDs
    """
    
    cc = len(listA.intersection(listB))
    cc_v = len(listB.intersection(listA))
    
    return (cc, cc_v)
    
def calculate_fcc_nc(listA, listB):
    """
    Calculates the fraction of common elements between two lists
    not taking into account chain IDs. Much Slower.
    """
    
    largest,smallest = sorted([listA, listB], key=len)
    ncommon = len([ele for ele in largest if ele in smallest])
    return (ncommon, ncommon)

# Matrix Calculation

def calculate_pairwise_matrix(contacts, ignore_chain, n_images=1):
    """ Calculates a matrix of pairwise fraction of common contacts (FCC).
        Outputs numeric indexes.

        contacts: list_of_unique_pairs_of_residues [set/list]
        
        Returns pairwise matrix as an iterator, each entry in the form:
        FCC(cplx_1/cplx_2) FCC(cplx_2/cplx_1)
        
        When images are included, construct only the submatrix containing
        the fcc between the base member with all images.
        Note all fcc between images map down to this subset,
        so can be avoided. ('.' in matrix below)
        # \abgde
        #  \....
        #   \...
        #    \..
        #     \.        
    """
    #Second failsafe.
    if ignore_chain and n_images>1:
        sys.stderr.write(" ERROR: cannot both ignore all chains and utilise chain symmetries!\n")
        sys.exit(2)

    contact_lengths = []
    for c in contacts:
        try:
            ic = 1.0/len(c)
        except ZeroDivisionError:
            ic = 0
        contact_lengths.append(ic)
      
    if ignore_chain:
        calc_fcc = calculate_fcc_nc
    else:
        calc_fcc = calculate_fcc

    n_xval = len(contacts)/n_images
    n_yval = len(contacts)
    for i in xrange(n_xval):
        for k in xrange(i+1, n_yval):
            cc, cc_v = calc_fcc(contacts[i], contacts[k])
            fcc, fcc_v = cc*contact_lengths[i], cc*contact_lengths[k]
            yield (i+1, k+1, fcc, fcc_v)

def _output_fcc(output, values, f_buffer):

    buf = []
    for i in values:
        buf.append(i)
        if len(buf) == f_buffer:
            output( ''.join(["%s %s %1.3f %1.3f\n" %(i[0],i[1],i[2],i[3]) for i in buf]) )
            buf = []
    output( ''.join(["%s %s %1.3f %1.3f\n" %(i[0],i[1],i[2],i[3]) for i in buf]) )

def output_chain_permutations(filename, n_confs, maps):
    """
    Log the permutations used for future use in resetting PDB files.
    """
    fp=open(filename,'w')
    print >> fp, '# Conformers:', n_confs
    print >> fp, '# Permutations:', len(maps)
    print >> fp, '# list of chain permutes below'
    for m in maps:
        mtxt = [ chr(a+64) for a in m ]
        print >> fp, str(mtxt).strip('[]').replace(',','').replace('\'','')
    fp.close()    
    
if __name__ == '__main__':
    
    import optparse
    import sys
    from time import time, ctime
    import os
    
    USAGE = "%s <contacts file 1> <contacts file 2> ... [options]\n" %os.path.basename(sys.argv[0])
    
    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option('-o', '--output', dest="output_file", action='store', type='string',
                        default=sys.stdout,
                        help='Output File [default: STDOUT]')
    parser.add_option('-f', '--file', dest="input_file", action='store', type='string',
                        help='Input file (one contact file name per line)')
    parser.add_option('-b', '--buffer_size', dest="buffer_size", action='store', type='string',
                        default=50000, 
                        help='Buffer size for writing output. Number of lines to cache before writing to file [default: 50000]')
    parser.add_option('-i', '--ignore_chain', dest="ignore_chain_char", action='store_true',
                        help='Ignore chain character in residue code. Use for homomeric complexes.')
    parser.add_option('-s', '--symmetry', dest="chain_symm_list", action='store', type='string',
                        default='',
                        help='Account for symmetries in the input chains. Use -m in ./cluster_fcc.py afterwards.'
                             'Pairwise matrix becomes (N_contacts)x(N_contacts*N_permutations).'
                             'Example: "AB,CDE" for a complex between a homodimer and a homotrimer, or '
                             '"cABCD" for a cyclic permutation of a tetramer')

    (options, args) = parser.parse_args()
    
    if options.input_file:
        args = [name.strip() for name in open(options.input_file)]
    
    if len(args) < 2:
        sys.stderr.write("- Provide (at least) two structures to calculate a matrix. You provided %s.\n" %len(args))
        sys.stderr.write(USAGE)
        sys.exit(1)

    sys.stderr.write("+ BEGIN: %s\n" %ctime())
    if options.ignore_chain_char:
        sys.stderr.write("+ Ignoring chains. Expect a considerable slowdown!!\n")
        exclude_chains = True
    else:
        exclude_chains = False

    #Interpret chain syntax.
    if len(options.chain_symm_list)>0:
        bConsiderSymmetry = True
        sys.stderr.write("...parsing chain symmetries list.\n")
        chsymm, permutations = parse_chain_symmetries(options.chain_symm_list)
        sys.stderr.write("...interpreted %d symmetry groups(s),\n" % len(chsymm))
        sys.stderr.write("   with %s permutations.\n" % str([len(permutations[i]) for i in xrange(len(chsymm))]).strip('[]') )
    else:
        bConsiderSymmetry = False
        sys.stderr.write("...no chain symmetries given. Will not consider images.\n")

    if exclude_chains==True and bConsiderSymmetry==True:
        sys.stderr.write(" ERROR: cannot both ignore all chains and use chain symmetries!\n")
        sys.exit(2)

    t_init = time()
    sys.stderr.write("+ Parsing %i contact files\n" %len(args))

    c = parse_contact_file(args, exclude_chains)
    
    if bConsiderSymmetry==True:
        #Generate all permutations and write it for future retracing.
        c, all_perms = generate_contact_symmetries(c, permutations)
        output_chain_permutations(options.output_file+'.maps', len(c)/len(all_perms), all_perms)
        m = calculate_pairwise_matrix(c, exclude_chains, len(all_perms))
    else:
        m = calculate_pairwise_matrix(c, exclude_chains)

    if isinstance(options.output_file, str):
        f = open(options.output_file, 'w')
    else:
        f = options.output_file

    sys.stderr.write("+ Calculating Matrix\n") # Matrix is calculated when writing. Generator property.
    sys.stderr.write("+ Writing matrix to %s\n" %f.name)
    _output_fcc(f.write, m, options.buffer_size)
    
    if isinstance(options.output_file, str):
        f.close()
    t_elapsed = time()-t_init
    sys.stderr.write("+ END: %s [%6.2f seconds elapsed]\n" %(ctime(), t_elapsed))
