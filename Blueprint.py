#!/usr/bin/env  python

from Bio.PDB import PDBParser, PDBIO

import re
from itertools import groupby
from collections import defaultdict

# Classes for manipulation of Rosetta's Blueprint files.

class Segment:
    def __init__(self, id, sstype,bp_data = [ ] , residues = None):
        self.id = id
        self.sstype = sstype
        self.bp_data = bp_data
        self.residues = residues # and iterable of biopython's BIO.PDB.Residue.Residue

    def renumerate(self, new_start_nres):
        indexer = new_start_nres
        for res in self.residues:
            res.id = (res.id[0], indexer, res.id[2])

            indexer += 1

    def abego(self):
         abego_st = ''
         for i in range(len(self.bp_data)): 
            abego_st+=self.bp_data[i][2][-1]

         return abego_st

    def set_abego(self,abego_string):
         for i in range(len(self.bp_data)):
            self.bp_data[i][2] = self.bp_data[i][2][0] + abego_string[i]


class SSPair:
    def __init__(self, strand1_segment, strand2_segment, orientation , register_shift = 0, offset = None, hairpin=False):
        self.strands = (strand1_segment, strand2_segment)
        self.paired = { strand1_segment : strand2_segment,
                         strand2_segment : strand1_segment }
        self.orientation = orientation # 'P' for parallel 'A' for antiparallel
        self.register_shift = register_shift
        self.offset = offset
        self.hairpin = hairpin
        if self.hairpin and self.offset==None:
            self.offset=0

    def get_pair(self, strand):
        return self.paired[strand]

class FoldInfo:
    def __init__(self, sspairs = []):
        self._pairs = set(sspairs)
        self._paired  = None
        self._sspairs = None
        self._update_pair_info()

    def _update_pair_info(self):
        self._paired  = defaultdict(list)
        self._sspairs = defaultdict(list)
        for sspair in self.pairs:
            self._paired[sspair.strands[0]].append(sspair.strands[1])
            self._paired[sspair.strands[1]].append(sspair.strands[0])
            self._sspairs[sspair.strands[0]].append(sspair)
            self._sspairs[sspair.strands[1]].append(sspair)
    
    def pairs(self):
        return self._pairs

    def remove_pair(self, sspair):
        self._pairs.remove(sspair)
        self._update_pair_info()

    def add_pair(self, sspair):
        self._pairs.add(sspair)
        self._update_pair_info()


class HSSTriplet:
    def __init__(self, helix, strand1, strand2):
        self.helix = helix
        self.strand1 = strand1
        self.strand2 = strand2

class Blueprint:
    def __init__(self, blueprint_file = None, pdbfile = None, structure = None, segments = None, data=None):
        if pdbfile:
            self.structure = PDBParser().get_structure(pdbfile, pdbfile)
        else:
            self.structure = structure

        if segments:
            self.segments = segments
            self.bp_data = [ ]
            self.segment_dict = { }
            for seg in segments:
                self.bp_data += seg.bp_data
                self.segment_dict[seg.id] = seg

        if blueprint_file and not data:
            # read the blueprint file and initialize segments
            # if self.structure is available put the residues in the segments.

            #self.segments = [ ]
            foldinfo_register = ""
            hsstriplet_register = ""
            # be careful here. originally i used the first line, but recently i need the third line (i dont remember why)
	        # i meed the first line for reading a plain reference blueprint without abegos
            register = re.compile('^\s*(\d+)\s+(\w+)\s+(\w+)\s+(.+)')
            #register = re.compile('^\s*(\d+)\s+(\w+)\s+(\w-|\w)\s+(.+)')
            #register = re.compile('^\s*(\d+)\s+(\w+)\s+(\w\D)\s+(.+)')
            data = [ ]
            for line in open(blueprint_file):
                if line.startswith('FOLDINFO'):
                    foldinfo_register = line.strip()
                elif line.startswith('HSSTRIPLET'):
                    hsstriplet_register = line.strip()
                elif line.startswith('HSSTRIAD'):
                    hsstriplet_register = line.strip()
                elif line.startswith('SSPAIR'):
                    #r = re.compile("(\d+)-(\d+).(\w).")
                    r = re.compile("(\d+)-(\d+).(\w).([-]?\d+)")
                    self.sspairs = r.findall(line)
                elif line.startswith('HHPAIR') or line[0]=='#':
                    pass
                else:
                    r = register.split(line)
                    data.append( [int(r[1]), r[2], r[3], r[4]] )

        if blueprint_file or data:
            # group the tuples in lists by their secondary structure and initiliaze the segments
            # grab the residues from the structure if this is available
            # self.bp_data contains all blueprint residue data
            # self segment_dict is a dict of segments where the keys are the ID for the ss segment. For example
            # H3 means Helix 3.
            self.segments = [ ]
            self.bp_data = [ ]
            self.segment_dict = { }
            res_index = 0
            segment_count = { 'L' : 1, 'H': 1 , 'E': 1 }
            residues = list(self.structure.get_residues()) if self.structure else None
            for sstype, bp_data  in groupby(data, key = lambda x: x[2][0]):
                resdata = list(bp_data)
                self.bp_data += resdata
                id = sstype + str(segment_count[sstype])
                segment_count[sstype] += 1
                seg = None
                if self.structure:
                    segment_residues = [ ]
                    for data in resdata:
                        segment_residues.append(residues[res_index])
                        res_index += 1
                    seg = Segment(id, sstype, resdata, segment_residues)
                else:
                    seg = Segment(id, sstype,  resdata)
                # append the segment to the segment list
                self.segments.append(seg)
                # insert the segment to the segment dict
                self.segment_dict[id] = seg
                #use the segment_dict to fill foldinfo and hsstriplet
                ##  I AM GOING TO FINISH THIS LATER BECAUSE IT IS GOING TO BE TRICKY TO SET UP THE FOLDS WITH THE SWAPP
                ##  MEANWHILE I AM GOING TO MODIFY dump_blueprint to take the foldinfo and hss tripplet as arguments
                #get_fold_tokens  = re.compile('(\d+-\d+\.[AP]\.-?\d)')
                #fold_tokens = get_fold_tokens.findall(foldinfo_register)
                #for ft in fold_tokens:
                #    pass

    def topology(self):
        return  reduce(lambda x,y: x + '-' + y , [s.id for s in self.segments])

    def topology_lengths(self):
        topol1 = reduce(lambda x,y: x + '-' + y , [s.id for s in self.segments])
        elements = re.compile("[HEL]\d+")
        ss_lst = elements.findall(topol1)
        topol2='' ; topol3=''
        for ss in ss_lst:
            seg = self.segment_dict[ss]
            n = len(seg.bp_data)
            topol2+='%s%s-' %(ss[0],n)
            topol3+='%s[%s-%s]' %(ss[0],n,n)

        return topol2, topol3

    def ss_tag(self):
        H = 0
        E = 0
        for s in self.segments:
            if s.sstype == 'H':
                H += 1
            elif s.sstype == 'E':
                E += 1
            else:
                pass
        return "%dH%dE" % (H, E)

    def freeze_all(self):
        for res in self.bp_data:
            res[3] = '.'

    def remodel_all(self):
        for res in self.bp_data:
            res[3] = 'R'

    def remodel_segment(self, index = None, id = None, index_to_zero = False,loop_edge = False,edges = True, abego=False):
        res_for_remodel = []
        if index:
            for res in self.segments[index].bp_data:
                res_for_remodel.append(res)
        elif id :
            for res in self.segment_dict[id].bp_data:
                res_for_remodel.append(res)

#	if edges:
#		prev_res_num = self.segment_dict[id].bp_data[0][0]-1	
#		next_res_num = self.segment_dict[id].bp_data[-1][0]+1	
#		res_for_remodel.append(self.bp_data[prev_res_num-1])
#		res_for_remodel.append(self.bp_data[next_res_num-1])

        for res in res_for_remodel:
            if index_to_zero:
                res[0] = 0
            res[3] = 'R'
        if abego:
            if len(res[2]) == 2:
                ss = res[2][0]
                res[2]="%sX" %ss
        
        #if loop_edge:
        if edges:
           for i in range(1, len(self.segments) - 1):
                prev_seg = self.segments[i-1]
                seg = self.segments[i]
                next_seg = self.segments[i+1]
                #if seg.sstype == 'L' or edges:
                if seg.id == id:
                    if seg.bp_data[0][3] == 'R':
                        prev_seg.bp_data[-1][3] = 'R'
                    if seg.bp_data[-1][3] == 'R':
                        next_seg.bp_data[0][3] = 'R'



    def residue_segment(self, pos):
        its_segment=''
        for segment in self.segment_dict.keys():
            seg = self.segment_dict[segment]
            for res in seg.bp_data:
                if res[0] == pos:
                    its_segment = segment
                    break
            else:
                continue
            break

        return its_segment

    def segment_lengths(self):
        return reduce(lambda i,j: i + '-' + j , [s.sstype + str(len(s.bp_data)) for s in self.segments])


    def reindex_blueprint(self, start = 1, rebuild_index_to_zero = False):
        indexer = start
        for bp_data in self.bp_data:
            if rebuild_index_to_zero and bp_data[3] == 'R':
                bp_data[0] = 0
            else:
                bp_data[0] = indexer
                indexer += 1

    def set_segment_abego(self, segment=None, abegos=None):
        seg = self.segment_dict[segment]
        if len(abegos) == len(seg.bp_data):
            for abego,res in zip(abegos,seg.bp_data):
                res[2]='%s%s' %(res[2][0],abego)

    def segment_list(self):
        r = re.compile('([HEL]\d+)-?')
        seg_list = r.findall(self.topology())
        return seg_list

    def sequence(self):
        st=''
        for res in self.bp_data:
            st+=res[1]
        return st

    def secstruct(self):
        st=''
        for res in self.bp_data:
                st+=res[2][0]
        return st



    def dump_blueprint(self, filename, header_lines = []):
        '''header lines are for setting foldinfo, hsstriplet or any other register on the top of the blueprint.'''
        out = open(filename, 'w')
        for line in header_lines:
            line.strip() # avoid doble carriage return
            out.write(line + '\n')
        for r in self.bp_data:
            out.write("%d    %s    %s    %s\n" % tuple(r))
        out.close()

    def dump_pdb(self, filename):
        io = PDBIO()
        io.set_structure(self.structure)
        io.save(filename)

    def swapp_segments(self, index1, index2):
        '''This function swaps the segments, reindexes the blueprint and PDB file
        and set for remodelling the segments directly conected to the swapped segments.
        The rest of the structure is set frozen.'''
        #freeze the structure and delete the residues conected to the swapped segments for remodel
        #add to the blueprint the corresponding insertions for the deleted residues.
        self.freeze_all()
        self.remodel_segment(index1 - 1, index_to_zero = True)
        self.remodel_segment(index1 + 1, index_to_zero = True)
        self.remodel_segment(index2 - 1, index_to_zero = True)
        self.remodel_segment(index2 + 1, index_to_zero = True)

        #wapp the self.segments
        self.segments[index1], self.segments[index2] = self.segments[index2], self.segments[index1]
        #renumerate the blueprint and the residues
        indexer = 1
        residues_to_detach = set()
        for segment in self.segments:
            for i in range(0,len(segment.bp_data)):
                if segment.bp_data[i][0]  == 0:
                    residues_to_detach.add(segment.residues[i])
                    continue
                segment.bp_data[i][0] = indexer
                id = segment.residues[i].id
                segment.residues[i].id = (id[0], indexer, id[2])
                indexer += 1

        # detach the residues of the residues directly connected to the swapp
        # this is done to avoid clashes during the remodelling
        for res in residues_to_detach:
            p = res.get_parent()
            p.detach_child(res.id)


        # sort the residues in the structure accoriding to the new indexing
        for chain in self.structure.get_chains():
            chain.child_list = sorted(chain.child_list, key = lambda r: r.id[1]) 

        #now that the elements have been reindexed self.bp_data and self.residues must be updated
        self.bp_data = reduce(lambda x,y :  x + y, [s.bp_data for s in self.segments])
        self.residues = reduce(lambda x,y :  x + y, [s.residues for s in self.segments])
        


if __name__ == '__main__':
    blue = Blueprint('f066/nmr_01.blueprint','f066/f066.NMR.2l69_01.pdb')

