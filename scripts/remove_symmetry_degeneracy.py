#!/usr/bin/env python
# -*- coding: UTF-8  -*-

class stringList(object):
    """Defines a 'clusterable' Element"""

    __slots__ = ['name', 'lines']

    def __init__(self, name):
        self.name = name
        self.lines = []

    def add_line(self, line):
        self.lines.append(line)

    def num_lines(self):
        return len(self.lines)

    def get_lines(self):
        return self.lines


def read_permutations(filename):
    f = open(filename, 'r')
    perm=[]
    for line in f:
        if line[0]=='#' or line=='\n':
            continue
        perm.append( line.split() )
    f.close()
    return perm

def read_mapping_from_matrix(filename):
    f = open(filename, 'r')

    #Target conformer to match
    tconf=1

    #Detect symmetries and compress
    for line in f:
        pass
    last=line.split()
    n_conf=int(last[0])
    n_image=int(last[1])/n_conf
    sys.stderr.write("Found %d conformers with %d images each.\n" % (n_conf, n_image))
    f.seek(0)

    minimum=0.0
    minlist=[minimum]*n_conf
    #Default
    mapping=[0]*n_conf
    #Read in FCC values and pick the maximum.
    for line in f:
        i, j, dRM, dMR = line.split()
        i=int(i)
        #End here
        if i < tconf:
            continue
        if i > tconf:
            break
        j=int(j)
        avdist=0.5*(float(dRM)+float(dMR))
        #Obtain base conformer of image
        jres = ( (j-1) % n_conf)+1
        if jres == tconf:
            continue
        image= int( (j-1)/n_conf)
        #print "map: ", i, jres, image, avdist
        if minlist[jres-1] < avdist:
            minlist[jres-1]=avdist
            mapping[jres-1]=image
            #print "...updated."

    f.close()
    sys.stderr.write("= = Matching image for each conformer:\n")
    sys.stderr.write("%s\n" % str(mapping))    
    #print mapping
    return mapping, n_conf, n_image


def read_pdbs(filename):
    pdb_list = [name.strip() for name in open(filename) if name.strip()]
    return pdb_list


def map_pdb_chains(out_prefix, pdb_list, mapping, permutations, bReverse):
    """
    Note that we are trying not to depend on external PDB manipulation libraries
    like MDAnalysis or mdtraj. So the PDB writing is reimplemented here
    in a very skeleton fashion.
    """
    n_conf  = len(mapping)
    n_image = len(permutations)

    base_image=permutations[0]

    #Process each file in the PDB list and map the chains and segids.
    for i in xrange(n_conf):
        in_fp =open(pdb_list[i], 'r')
        #out_fp=open(out_prefix+str(i+1)+'.pdb', 'w')

        if bReverse:
            src_perm=permutations[mapping[i]]
            targ_perm=permutations[0]
        else:
            src_perm=permutations[0]
            targ_perm=permutations[mapping[i]]
        print "= = = Mapping from ", src_perm, " to ", targ_perm

        chainList = {}
        #bFirst=True
        for line in in_fp:
            if line[0:4]!='ATOM' and line[0:6]!='HETATM':
                #print >> out_fp, line[:-1]
                continue
            # Map chains here. Position 22 is the chain.
            #if bFirst:
            #    print line[21]
            #    bFirst=False
            try:
                ch_index=src_perm.index(line[21])
                new_chain=targ_perm[ch_index]
                new_line=line[0:21]+new_chain+line[22:]
            except:
                new_chain=line[21]
                new_line=line
                #print >> out_fp, line[:-1]
                #continue

            # Append input line to our collection
            # Create or Retrieve Elements
            if new_chain not in chainList:
                r = stringList(new_chain)
                chainList[new_chain] = r
            else:
                r = chainList[new_chain]
            #NB: Remove the new-line marker to prevent duplication.
            r.add_line(new_line[:-1])

            #new_chain=targ_perm[ch_index]
            #new_line=line[0:21]+new_chain+line[23:-1]
            #print >> out_fp, new_line

        in_fp.close()

        out_fp=open(out_prefix+str(i+1)+'.pdb', 'w')
        # Rerank chains to place them in new order.
        atom_id=1
        for key, obj in sorted(chainList.items()):
            for line in obj.get_lines():
                # Renumber atoms now.
                idstr='%5d' % atom_id
                line=line[0:6]+idstr+line[11:]
                print >> out_fp, line
                atom_id+=1
            print >> out_fp, "TER"
        print >> out_fp, "END"

        out_fp.close()

#def debug_mapping(mapping, permutations):



if __name__ == "__main__":

    import optparse
    import sys
    from time import time, ctime
    import os

    USAGE="%s -s <matrix file> -m <permutation map> -l <pdb list> [-o output_prefix]" %os.path.basename(sys.argv[0])

    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option('-o', '--output', dest="out_prefix", action='store', type='str',
                    default='out',
                    help='Output File [STDOUT]')
    parser.add_option('-s', '--matrix', dest="matrix_file", action='store', type='str',
                    default='',
                    help='Matrix file produced from calc_fcc_matrix.py')
    parser.add_option('-l', '--pdb_list', dest="pdb_list", action='store', type='str',
                    default='',
                    help='Starting list of PDB file.')
    parser.add_option('-m', '--perm_file', dest="perm_file", action='store', type='str',
                    default='',
                    help='Permutation map file produced from calc_fcc_matrix.py')

    (options, args) = parser.parse_args()

    if sys.version_info[0:2] < (2,6):
        cur_version = "%s.%s" %sys.version_info[0:2]
        sys.stderr.write("- Python version not supported (%s). Please use 2.5 or newer.\n" %cur_version )
        sys.exit(1)

    if options.perm_file=='' or options.pdb_list=='' or options.matrix_file=='':
        sys.stderr.write("USAGE: %s\n" %USAGE)
        sys.exit(1)


    # Read Matrix
    sys.stderr.write("+ BEGIN: %s\n" %ctime())
    t_init = time()

    try:
        permutations = read_permutations(options.perm_file)
    except IOError:
        sys.stderr.write("File not found: %s\n" % options.perm_file)
        sys.exit(1)

    try:
        mapping, n_conf, n_image = read_mapping_from_matrix(options.matrix_file)
    except IOError:
        sys.stderr.write("File not found: %s\n" % options.matrix_file)
        sys.exit(1)

    try:
        pdb_list = read_pdbs(options.pdb_list)
    except IOError:
        sys.stderr.write("File not found: %s\n" % options.pdb_list)
        sys.exit(1)

    #Take each PDB file and conduct a reverse chain mapping.
    if n_conf != len(pdb_list):
        sys.stderr.write(
                "ERROR: the number of conformers from the matrix (%d) and pdb list (%d) do not agree!\n"
                % (n_conf,len(pdb_list)) )
        sys.exit(2)
    if n_image != len(permutations):
        sys.stderr.write(
                "ERROR: the number of images from the matrix (%d) and permutations file (%d) do not agree!\n"
                % (n_image, len(permutations)) )
        sys.exit(2)

    #debug_mapping(mapping, permutations)

    map_pdb_chains(options.out_prefix, pdb_list, mapping, permutations, False)

    t_elapsed = time()-t_init
    sys.stderr.write( "+ END: %s [%3.2f seconds]\n" %(ctime(), t_elapsed))
