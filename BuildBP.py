import itertools
import re
import sys
import os
import copy
from Blueprint import Blueprint 
from argparse import ArgumentParser
import numpy as np


def WriteResfile(**kwargs):
	bluefile = kwargs.get('blueprint')
	resfile_name = kwargs.get('resfile_name')

        resfile = open(resfile_name,'w')
        resfile.write('AUTO\nstart\n')

	blue = Blueprint(bluefile)
	blue.reindex_blueprint(start=1)

	# hairpin
        hairpinloops = ['L2','L8','L4']
        for loop in hairpinloops:
                seg = blue.segment_dict[loop]
                pos = int(seg.bp_data[0][0])
                abego = ''
                for i in range(len(seg.bp_data)): abego+=seg.bp_data[i][2][-1]
                if abego=='EA':
                        resfile.write('%i A PIKAA G\n' %(pos))
                        resfile.write('%i A PIKAA N\n' %(pos+1))
                elif abego=='GG':
                        resfile.write('%i A PIKAA ND\n' %(pos))
                        resfile.write('%i A PIKAA G\n' %(pos+1))


        # arch
        archloops = ['L3','L5']
        for loop in archloops:
                seg = blue.segment_dict[loop]
                pos = int(seg.bp_data[0][0])
                abego = ''
                for i in range(len(seg.bp_data)): abego+=seg.bp_data[i][2][-1]
                if abego=='BB':
			resfile.write('%i A PIKAA GILV\n' %(pos-1))
                        resfile.write('%i A PIKAA EGLPST\n' %(pos))
                        resfile.write('%i A PIKAA P\n' %(pos+1))
			resfile.write('%i A PIKAA GV\n' %(pos+2))
                elif abego=='BBB':
			resfile.write('%i A PIKAA KEV\n' %(pos-1))
                        resfile.write('%i A PIKAA PSL\n' %(pos))
                        resfile.write('%i A PIKAA PEK\n' %(pos+1))
			resfile.write('%i A PIKAA P\n' %(pos+2))
			resfile.write('%i A PIKAA TGIV\n' %(pos+3))
                elif abego=='BAA':
			resfile.write('%i A PIKAA DGPS\n' %(pos))
                        resfile.write('%i A PIKAA PE\n' %(pos+1))
			resfile.write('%i A PIKAA NGET\n' %(pos+2))
                elif abego=='BABB':
			resfile.write('%i A PIKAA DKT\n' %(pos-1))
                        resfile.write('%i A PIKAA SAGLTHP\n' %(pos))
                        resfile.write('%i A PIKAA DEPKNTS\n' %(pos+1))
                        resfile.write('%i A PIKAA EKRSNQ\n' %(pos+2))
                        resfile.write('%i A PIKAA DPS\n' %(pos+3))
                        #resfile.write('%i A PIKAA G\n' %(pos+4))			
                elif abego=='ABB':
			resfile.write('%i A PIKAA ALVSTI\n' %(pos-1))
			resfile.write('%i A PIKAA DEKN\n' %(pos))
                        #resfile.write('%i A PIKAA EKT\n' %(pos+1))
                        resfile.write('%i A PIKAA P\n' %(pos+2))
                        resfile.write('%i A PIKAA GILV\n' %(pos+3))
                elif abego=='BBG':
                        #resfile.write('%i A PIKAA PA\n' %(pos+1))
                        #resfile.write('%i A PIKAA GN\n' %(pos+2))
                        #resfile.write('%i A PIKAA ST\n' %(pos+3))

			# expand profile including info from BGB loop
			resfile.write('%i A PIKAA PN\n' %(pos))
                        resfile.write('%i A PIKAA PRAKE\n' %(pos+1))
                        resfile.write('%i A PIKAA GN\n' %(pos+2))
                        resfile.write('%i A PIKAA STED\n' %(pos+3))

                elif abego=='BGB':
			resfile.write('%i A PIKAA EKRST\n' %(pos))
                        resfile.write('%i A PIKAA G\n' %(pos+1))
                        resfile.write('%i A PIKAA ED\n' %(pos+2))
                elif abego=='BBGB':
                        resfile.write('%i A PIKAA DEKPQRSTGA\n' %(pos))
                        resfile.write('%i A PIKAA AKQP\n' %(pos+1))
			resfile.write('%i A PIKAA G\n' %(pos+2))
			resfile.write('%i A PIKAA ADEGQST\n' %(pos+3))
                elif abego=='BBGBB':
			resfile.write('%i A PIKAA VILMFYW\n' %(pos-1))
                        resfile.write('%i A PIKAA P\n' %(pos))
                        #resfile.write('%i A PIKAA AKQPRST\n' %(pos+1))
                        resfile.write('%i A PIKAA G\n' %(pos+2))
                        #resfile.write('%i A PIKAA ADEGQST\n' %(pos+3))
			#resfile.write('%i A PIKAA VILMFYWH\n' %(pos+4))
			resfile.write('%i A PIKAA H\n' %(pos+5))
                elif abego=='ABABB':
                        resfile.write('%i A PIKAA N\n' %(pos-1))
                        #resfile.write('%i A PIKAA AL\n' %(pos))
                        resfile.write('%i A PIKAA ST\n' %(pos+1))
                        resfile.write('%i A PIKAA DSGH\n' %(pos+2))
                        #resfile.write('%i A PIKAA EKQST\n' %(pos+3))
                        resfile.write('%i A PIKAA P\n' %(pos+4))
			#resfile.write('%i A PIKAA G\n' %(pos+5))

			
	resfile.close()

	aalist=""
	for line in open(resfile_name):
		if 'PIKAA' in line:
			pos = line.split()[0]
			aalist+='%s,' %pos

	return aalist[:-1]


def Bulged(strand):
        flag=False ; bulgepos=None
        for i in range(1,len(strand.bp_data)-1):
                prev_res = strand.bp_data[i-1]
                res = strand.bp_data[i]
                next_res = strand.bp_data[i+1]
                if 'EA' == res[2] and 'EB' == prev_res[2] and 'EB' == next_res[2]:
                    flag=True
                    bulgepos = res[0]
        return bulgepos

def EBulged(strand):
        flag=False ; bulgepos=None
        for i in range(1,len(strand.bp_data)-1):
                prev_res = strand.bp_data[i-1]
                res = strand.bp_data[i]
                next_res = strand.bp_data[i+1]
                if 'EE' == res[2] and 'EB' == prev_res[2] and 'EB' == next_res[2]:
                    flag=True
                    bulgepos = res[0]
        return bulgepos

def MakePlainBlueprint(ss,comb,bluefile):
        c = comb     
        out_file = open(bluefile,'w')
        struct = {'H':'HA', 'E':'EB', 'L':'LA'}
        total_length=sum(c)
        k=0
        curr_ss = ss[k]
        for i in range(1,total_length+1):
                if i>sum(c[:k+1]):
                        k+=1
                        curr_ss = ss[k]
                out_file.write('0  V  %s  R\n' %(curr_ss))
        out_file.close()

#------------
def MakeRefBlueprint(ss,comb,bulges,**kwargs):	
	refbluefile = kwargs.get('refblue')
	e_bulges = kwargs.get('e_bulges',None)
        c = comb     
        out_file = open(refbluefile,'w')
        struct = {'H':'HA', 'E':'EB', 'L':'LA'}
        total_length=sum(c)
        k=0
        curr_ss = ss[k]
        for i in range(1,total_length+1):
                if i>sum(c[:k+1]):
                        k+=1
                        curr_ss = ss[k]
                out_file.write('0  V  %s  R\n' %(curr_ss))
        out_file.close()

        # Put bulges
        blue = Blueprint(refbluefile)
       	bluelist = []
       	for seg in blue.segments:
       	        if seg.id in bulges.keys():
       	                for j, res in enumerate(seg.bp_data):
       	                        if j == bulges[seg.id]-1:
       	                                #res[2] = 'EE' ; res[1] = 'G'
					res[2] = 'EA'
				
		if e_bulges:
	                if seg.id in e_bulges.keys():
        	                for j, res in enumerate(seg.bp_data):
                	                if j == e_bulges[seg.id]-1:
                        	                res[2] = 'EE' ; res[1] = 'G'

       	blue.dump_blueprint(refbluefile)
	#os.chdir('../')

#---------------
def Shift(**kwargs):
        refbluefile = kwargs.get('refblue')
        seg1 = kwargs.get('seg1')
        seg2 = kwargs.get('seg2')
        refblue = Blueprint(refbluefile)
        nseg1 = len( refblue.segment_dict[seg1].bp_data )
        nseg2 = len( refblue.segment_dict[seg2].bp_data )
        shift = nseg1 - nseg2
        return shift
#---------------
def MakeFirstBlueprint(**kwargs):
	tail = [[0, 'V', 'L', 'R']]
	refbluefile = kwargs.get('refblue')
	segments = kwargs.get('segments')
	newbluefile = kwargs.get('newblue')
	ss_pairing = kwargs.get('ss_pairing')
	hs_pairing = kwargs.get('hs_pairing')
	hh_pairing = kwargs.get('hh_pairing')
	seg_adapt = kwargs.get('adapt',None)
	seg_abego = kwargs.get('seg_abego')

	refblue = Blueprint(refbluefile)
	shift=0
	if seg_adapt:
	        seg1 = seg_adapt.keys()[0] # strand to adapt
        	seg2 = seg_adapt[seg1] # referemce strand
	        nseg1 = len( refblue.segment_dict[seg1].bp_data )
	        nseg2 = len( refblue.segment_dict[seg2].bp_data )
        	shift = nseg1 - nseg2

        bp_data_new = []
        for seg in segments:
                if shift == 0:
                        for res in refblue.segment_dict[seg].bp_data:
                                bp_data_new.append(res)
                elif shift > 0:
                        if seg==seg1:
                                for k in range(0,nseg1-shift):
                                        bp_data_new.append(refblue.segment_dict[seg1].bp_data[k])
                        else:
                                for res in refblue.segment_dict[seg].bp_data:
                                        bp_data_new.append(res)
                elif shift < 0:
                        if seg==seg2:
                                for k in range(shift,nseg2):
                                        bp_data_new.append(refblue.segment_dict[seg2].bp_data[k])
                                        bp_data_new.append(refblue.segment_dict[seg1].bp_data[k])
                        else:
                                for res in refblue.segment_dict[seg].bp_data:
                                        bp_data_new.append(res)

	refblue.bp_data = tail + bp_data_new + tail

	# write strand pairing	
	ss_line = "SSPAIR "
	# Get first strand to set "start"
	segments.sort()
	npairs=0
	for i,key in enumerate(ss_pairing.keys()):
            for st in ss_pairing[key]:
		s1=key
		s2 = st[:st.find('.')]
		orient = st[-1]
		# compute effective length of strands for calculating register shift
		#------
		if Bulged(refblue.segment_dict[s1]):
			n1 = len( refblue.segment_dict[s1].bp_data ) - 1
		else:
			n1 = len( refblue.segment_dict[s1].bp_data )
                if Bulged(refblue.segment_dict[s2]):
                        n2 = len( refblue.segment_dict[s2].bp_data ) - 1
                else:
                        n2 = len( refblue.segment_dict[s2].bp_data )
		#------
		if orient == 'A':
		#	shift = n1-n2
			shift=99
		
		if int(s1[1:]) < int(s2[1:]):
			seg1_tag=int(s1[1:])
			seg2_tag=int(s2[1:])
		else:
                        seg1_tag=int(s2[1:])
                        seg2_tag=int(s1[1:])

		if npairs==0:
			ss_line += '%s-%s.%s.%s' %(seg1_tag,seg2_tag,orient,shift)	
		else:
			ss_line += ';%s-%s.%s.%s' %(seg1_tag,seg2_tag,orient,shift)	
		npairs+=1

	header = [ss_line]
	if hs_pairing != None:
		header.append( hs_pairing )

	if hh_pairing !=None:
		header.append( hh_pairing )

	refblue.dump_blueprint(newbluefile,header_lines=header)

        blue0_top = refblue.topology() # for abego conversion later
        blue0_cp = copy.deepcopy(refblue)

        # abego loop
        # Adapt global abego motif to current stage
        #---------------------------------
        if seg_abego != None:
                new_abego={};conversor={}
                r=re.compile('[HEL]')
                top= blue0_cp.topology()
                #curr_ss = r.findall(top)
                counter=0
		curr_ss=[]
                for i,seg in enumerate(segments):
                        ss = seg[0]
			curr_ss.append(ss)
			if ss=='L':
				newindex = curr_ss.count(ss)+1
			else:
				newindex = curr_ss.count(ss)

			conversor[seg]='%s%s' %(ss,newindex)

                for seg in seg_abego:
                        new_abego[conversor[seg]] = seg_abego[seg]


                blue = Blueprint('%s' %(newbluefile))
                if new_abego != None:
                        for seg in new_abego.keys():
                                abego=new_abego[seg]
                                for i,res in enumerate(blue.segment_dict[seg].bp_data):
                                        res[2]+=abego[i]

                blue.dump_blueprint(newbluefile,header_lines=header)
	
	# Abego 'B' for building strand
	os.system("sed  's/ E / EB/' %s > %s.b" %(newbluefile,newbluefile))
				
			
#-------------------
def AddSegmentToBlueprint(**kwargs):
	tail = [[0, 'V', 'L', 'R']]
	refbluefile = kwargs.get('refblue')
	segments = kwargs.get('segments')
	blue0file= kwargs.get('blue0')
	newbluefile = kwargs.get('newblue')
	append = kwargs.get('append')
	insert_between_first = kwargs.get('insert_between_first')
	insert_between_last = kwargs.get('insert_between_last')
	ss_pairing = kwargs.get('ss_pairing')
	ss_pairing_shift = kwargs.get('ss_pairing_shift')
	hs_pairing = kwargs.get('hs_pairing')
	hh_pairing = kwargs.get('hh_pairing')
	seg_abego = kwargs.get('seg_abego')
	specific_abego = kwargs.get('specific_abego')
	insert = kwargs.get('insert')
	only_remodel = kwargs.get('only_remodel',False)

	blue0 = Blueprint(blue0file)
	blue0.reindex_blueprint(start=1)
	blue0.freeze_all()

	refblue = Blueprint(refbluefile)
	bp_data_new = []
	for seg in segments:
		for res in refblue.segment_dict[seg].bp_data:
			bp_data_new.append(res)		

        if append==True:
                blue0.bp_data[-2][3]='R'
                blue0.bp_data[-1][3]='R' ; blue0.bp_data[-1][2]='L'
                bp_data_new = blue0.bp_data + bp_data_new[1:] + tail

        elif append==False:
                blue0.bp_data[0][3]='R' #; blue0.bp_data[0][2]='L'
                blue0.bp_data[1][3]='R'
		if segments[-1][0] != 'L':
			blue0.bp_data[0][2]= bp_data_new[-1][2]
                bp_data_new =  tail + bp_data_new[:-1] + blue0.bp_data

	elif insert_between_first:
		# replace original residues by insert ones. so we can use original bp as blue0
		index1 = insert_between_first - 1
		index2 = insert_between_last - 1
		if index1>1:  # we need at least one point fixed (not R)
			blue0.bp_data[index1][3]='R' #; blue0.bp_data[index1][2]='LX'
		
		bp_data_new = blue0.bp_data[:index1+1] + bp_data_new # + blue0.bp_data[index2:] 
		shift = index2-index1-1
		for k in range(index2,len(blue0.bp_data)):
			blue0.bp_data[k][0]-=shift
		blue0.bp_data[index2][3]='R'
		bp_data_new = bp_data_new + blue0.bp_data[index2:]

        else:
                bp_data_new = blue0.bp_data		

	blue0_top = blue0.topology() # for abego conversion later
	blue0_cp = copy.deepcopy(blue0)

        # if we dont add residues then we dont need to replace bp_data. it is already update with R
        if only_remodel==False or isinstance(only_remodel,list):
                blue0.bp_data = bp_data_new

        # write strand pairing
	if ss_pairing != None:
         ss_line = "SSPAIR "
  	 keys = ss_pairing.keys()
	 keys.sort()
	 npairs=0
         for i,key in enumerate(keys):
	    for st in ss_pairing[key]:
		s1=key
	        s2 = st[:st.find('.')]
                orient = st[-1]

                # compute effective length of strands for calculating register shift
                #------
                if Bulged(refblue.segment_dict[s1]):
                        n1 = len( refblue.segment_dict[s1].bp_data ) - 1
                else:
                        n1 = len( refblue.segment_dict[s1].bp_data )
                if Bulged(refblue.segment_dict[s2]):
                        n2 = len( refblue.segment_dict[s2].bp_data ) - 1
                else:
                        n2 = len( refblue.segment_dict[s2].bp_data )
                #------

                #if orient == 'A':
                #shift = n1-n2
		shift=99	
		if npairs==0:
                        ss_line += '%s-%s.%s.%s' %(int(s1[1:]),int(s2[1:]),orient,shift)
                else:
                        ss_line += ';%s-%s.%s.%s' %(int(s1[1:]),int(s2[1:]),orient,shift)
		npairs+=1

	header=[]
	if ss_pairing_shift != None:
	        header.append( ss_pairing_shift )
	else:
		header.append( ss_line )
        if hs_pairing != None:
                header.append( hs_pairing )

        if hh_pairing !=None:
                header.append( hh_pairing )

	blue0.dump_blueprint(newbluefile, header_lines=header)


        # insertion
        if insert:
                seg = insert.keys()[0]
                shift = insert[seg]
                blue = Blueprint('%s' %(newbluefile))
                blue_aux = copy.deepcopy(blue)
        #       print blue.segment_dict[seg].bp_data
                blue_aux.reindex_blueprint()
                insert_index = blue_aux.segment_dict[seg].bp_data[-1][0] # we want to insert just before the 1st residue of the next seg, so that we insert at the very end

                for k in range(shift):
                        blue.bp_data.insert(insert_index, [0, 'V', 'E', 'R'])
        #       print blue.segment_dict[seg].bp_data

		# FOR NEW CONSTRAINTS WE NEED 7 RESIDUES TO USE BEND CONSTRAINTS (5 + 2 INSERTED)
		#for i in range(1,5+1):
		#	blue.segment_dict[seg].bp_data[-i][3]='R'

                blue.dump_blueprint(newbluefile,header_lines=header)



	# abego loop
	# Adapt global abego motif to current stage
	#---------------------------------
	if seg_abego != None:
		new_abego={};conversor={}
		r=re.compile('[HEL]')
		top= blue0_cp.topology()
		curr_ss = r.findall(top)
		counter=0
		for seg in segments:
			ss = seg[0]
			if append:
				if ss=='L':
					newindex = curr_ss.count(ss)
				else:
					newindex = curr_ss.count(ss) + 1
					curr_ss.append(ss)
			elif append==False:
				curr_ss.insert(counter,ss)
				newindex = curr_ss[:counter+1].count(ss)
				if ss=='L': newindex+=1 # because Nter is loop L1 (and is not added through segments)
				counter+=1
			elif insert_between_first:	
				r=re.compile('[HEL]\d') ; seg_list = r.findall(top)
				flag=False
				for k,sg in enumerate(seg_list): # find segment whhere insert_first is located
					for res in blue0_cp.segment_dict[sg].bp_data:
						if res[0]==insert_between_first+1:
							flag=True
							break
					if flag:
						break
	
				first_seg=k
				newindex = curr_ss[:first_seg+1+counter].count(ss)
				counter+=1

			conversor[seg]='%s%s' %(ss,newindex)

		for seg in seg_abego:
			new_abego[conversor[seg]] = seg_abego[seg]

					
        	blue = Blueprint('%s' %(newbluefile))
        	if new_abego != None:
        	        for seg in new_abego.keys():
        	        	abego=new_abego[seg]
        	        	for i,res in enumerate(blue.segment_dict[seg].bp_data):
        	                	res[2]+=abego[i]

        	blue.dump_blueprint(newbluefile,header_lines=header)
	
	if specific_abego: # only abegos for specific positions
		blue = Blueprint('%s' %(newbluefile))
		for seg in specific_abego:
			pos,letter =  specific_abego[seg]
			blue.segment_dict[seg].bp_data[pos][2]+=letter

		blue.dump_blueprint(newbluefile,header_lines=header)

	# only remodel.
	blue = Blueprint('%s' %(newbluefile))
        if isinstance(only_remodel,list):
                for seg in only_remodel:
                        blue.remodel_segment(id=seg)

		blue.dump_blueprint(newbluefile,header_lines=header)
		
	# Abego 'B' for building strand
	os.system("sed  's/ E / EB/' %s > %s.b" %(newbluefile,newbluefile))

#-------------------
def write_dummy_pdb(filename):
    dummy_pdb = 'ATOM      1  N   GLY A   1       0.346   1.210   0.744  1.00  0.00 \n\
ATOM      2  CA  GLY A   1       1.687   1.135   0.174  1.00  0.00 \n\
ATOM      3  C   GLY A   1       2.383   2.488   0.222  1.00  0.00 \n\
ATOM      4  O   GLY A   1       2.996   2.918  -0.752  1.00  0.00 \n\
ATOM      5 1H   GLY A   1      -0.448   0.981   0.185  1.00  0.00 \n\
ATOM      6 2H   GLY A   1       0.109   0.647   1.536  1.00  0.00 \n\
ATOM      7 3H   GLY A   1       0.005   2.079   1.101  1.00  0.00 \n\
ATOM      8 1HA  GLY A   1       2.277   0.413   0.741  1.00  0.00 \n\
ATOM      9 2HA  GLY A   1       1.615   0.810  -0.864  1.00  0.00 \n\
ATOM     10  N   GLY A  2       2.284   3.158   1.368  1.00  0.00 \n\
ATOM     12  CA  GLY A  2       3.918   4.450   2.676  1.00  0.00 \n\
ATOM     13  C   GLY A  2       3.551   4.379   3.850  1.00  0.00 \n\
ATOM     14  O   GLY A  2       1.859   5.564   1.752  1.00  0.00 \n\
ATOM     15 1H   GLY A  2       1.061   5.930   0.512  1.00  0.00 \n\
ATOM     16 2H   GLY A  2      -0.012   6.930   0.747  1.00  0.00 \n\
ATOM     17 3H   GLY A  2      -0.863   7.182  -0.405  1.00  0.00 \n\
ATOM     18 1HA  GLY A  2      -0.564   8.037  -1.404  1.00  0.00 \n\
ATOM     19 2HA  GLY A  2       0.540   8.749  -1.379  1.00  0.00 \n'
    out = open(filename, 'w')
    out.write(dummy_pdb)


def XMLReplaceXXXYYY(**kwargs):
        xml_lines = kwargs.get('xml_lines')
        identifier = kwargs.get('identifier')
        xxx = kwargs.get('xxx')
        yyy = kwargs.get('yyy')

        for i,line in enumerate(xml_lines):
                if identifier in line:
                        if xxx:
                                line = line.replace('xxx','%s' %(xxx))
				xml_lines[i] = line

                        if yyy:
                                line = line.replace('yyy','%s' %(yyy))
				xml_lines[i] = line
	#return xml_lines

def XMLReplaceTagsValues(**kwargs):
        xml_lines = kwargs.get('xml_lines')
        identifier = kwargs.get('identifier')
        tags = kwargs.get('tags')
        values = kwargs.get('values')

        for i,line in enumerate(xml_lines):
                if identifier in line:
			for tag,value in zip(tags,values):
				line = line.replace('%s' %(tag),'%s' %(value))
				xml_lines[i] = line


def XMLReplaceString(**kwargs):
        xml_in = kwargs.get('xml_in')
        xml_out = kwargs.get('xml_out')
        identifier = kwargs.get('identifier')
        string = kwargs.get('string')
        filein = open(xml_in,'r')
        fileout = open(xml_out,'w')
        for line in filein:
                if identifier in line:
                        line = line.replace('xxx','%s' %(string))
                fileout.write(line)
        filein.close()
        fileout.close()

def XMLAtomsDistance(**kwargs):
        xml_in = kwargs.get('xml_in')
        xml_out = kwargs.get('xml_out')
        identifier = kwargs.get('identifier')
        res_i = kwargs.get('res_i')
        res_j = kwargs.get('res_j')

        filein = open(xml_in,'r')
        fileout = open(xml_out,'w')
        for line in filein:
                if identifier in line:
                        line = line.replace('xxx','%s' %(res_i))
                        line = line.replace('yyy','%s' %(res_j))
                fileout.write(line)
        filein.close()
        fileout.close()

def XMLSegmentAtomDistance(**kwargs):
        xml_in = kwargs.get('xml_in')
        xml_out = kwargs.get('xml_out')
        bluefile = kwargs.get('blue')
        identifier = kwargs.get('identifier')
        compound_name = kwargs.get('compound_name')
        segment = kwargs.get('segment')
        resid = kwargs.get('resid')
        blue = Blueprint(bluefile)
        blue.reindex_blueprint(start=1)
        seg = blue.segment_dict[segment]
        fileout = open(xml_out,'w')
        for line in open(xml_in):
                if identifier in line:
                        fileout.write(line)
                        for j,res in enumerate(seg.bp_data):
                                line2 = '\t\t<AtomicDistance name="cn%i" residue1=%i residue2=%i atomname1=O atomname2=N distance=4.0 confidence=1 />\n' %(j,resid,res[0]) # 1:bulge ; 2:loop

                                fileout.write(line2)
                        # write compount statement
                        fileout.write( '\n\t\t<CompoundStatement name=%s >\n' %(compound_name) )
                        for k in range(j+1):
                                fileout.write( '\t\t\t<OR filter_name=cn%i />\n' %(k) )
                        fileout.write( '\t\t</CompoundStatement>' )
                else:
                        fileout.write(line)

def AmbiguousConstraints(list1,list2):
	st='AmbiguousConstraint\n'
	for res1 in list1:
		for res2 in list2:
			st += "AtomPair N %i O %i BOUNDED 3.5 4.5 0.5\n" %(res1,res2)
	st+='END_AMBIGUOUS\n'	
	return st

def ReplaceLine(line,word1,word2):
	if word1 in line:
		line.replace(word1,word2)

def MotifTopology(topol,motif_filename):
	elements = re.compile("[HEL]")
	ss = elements.findall(topol)
	relengths = re.compile("(\d+)-(\d+)")
	relengths2 = re.compile("(\d+),(\d+)") # for specific lengths, not ranges
	lengths = relengths.findall(topol)
	lengths2 = relengths2.findall(topol)

	## for interpreting specific lengths
	index=-1
	specific=[]
	for st in topol:
		if st=='[':
			index+=1
		if st==',':
			specific.append(index)

	comb=[] ; j=0 ; k=0
	for i in range(len(ss)):
		if i in specific:
			frag = lengths2[j]
			comb.append( [int(l) for l in frag] )		
			j+=1
		else:
			fragment = lengths[k]
			comb.append(range(int(fragment[0]),int(fragment[1])+1))
			k+=1

	# index ss elements for reading motif
	nl=0; nh=0; ne=0
	s_index=[]
	for s in ss:
		if s=='L':
			nl+=1
			s_index.append(s+'%s' %(nl))
		if s=='H':
			nh+=1
			s_index.append(s+'%s' %(nh))
		if s=='E':
			ne+=1
			s_index.append(s+'%s' %(ne))

	#----------------
	# read motifile
	#----------------
	motifile = open(motif_filename)
	header = motifile.readline()
	motif = motifile.readline()
	dic_motif={}
	for a,b in zip(header.split(),motif.split()):
		dic_motif[a] = b

	dic_abego={}
	for k,s in enumerate(s_index):
		if s in dic_motif.keys():
			motif_value = dic_motif[s]
			if motif_value.isdigit():
				comb[k]=[int(motif_value)]
			else:
				comb[k]=[len(motif_value)]
				dic_abego[s] = motif_value


	combinations = list(itertools.product(*comb))
	print 'Number of combinations: %s' %(len(combinations))

	return ss,combinations, dic_abego

def GetCombinations(topol):
	elements = re.compile("[HEL]")
	ss = elements.findall(topol)
	relengths = re.compile("(\d+)-(\d+)")
	relengths2 = re.compile("(\d+),(\d+)") # for specific lengths, not ranges
	lengths = relengths.findall(topol)
	lengths2 = relengths2.findall(topol)

	## for interpreting specific lengths
	index=-1
	specific=[]
	for st in topol:
		if st=='[':
			index+=1
		if st==',':
			specific.append(index)

	comb=[] ; j=0 ; k=0
	for i in range(len(ss)):
		if i in specific:
			frag = lengths2[j]
			comb.append( [int(l) for l in frag] )		
			j+=1
		else:
			fragment = lengths[k]
			comb.append(range(int(fragment[0]),int(fragment[1])+1))
			k+=1

	combinations = list(itertools.product(*comb))
	print 'Number of combinations: %s' %(len(combinations))

	return ss,combinations

def HBondConstraints(donor,acceptor): # return string for cst file
    hb_ang = np.deg2rad(160.)
    hb_ang_tol=np.deg2rad(20.0)
    st = "AtomPair N %i O %i HARMONIC 3.0 0.5\n" %(donor,acceptor) 
    st+= "Angle N %i H %i O %i HARMONIC %3.1f %3.1f\n" %(donor,donor,acceptor,hb_ang,hb_ang_tol)
    st+= "Angle H %i O %i C %i HARMONIC %3.1f %3.1f\n" %(donor,acceptor,acceptor,hb_ang,hb_ang_tol)
    return st	

def CircularHBondConstraints(donor,acceptor): # return string for cst file
    hb_ang = np.deg2rad(160.)
    hb_ang_tol=np.deg2rad(20.0)
    st = "AtomPair N %i O %i HARMONIC 3.0 0.5\n" %(donor,acceptor)
    #st+= "Angle N %i H %i O %i CIRCULARHARMONIC %3.1f %3.1f\n" %(donor,donor,acceptor,hb_ang,hb_ang_tol)
    st+= "Angle H %i O %i C %i CIRCULARHARMONIC %3.1f %3.1f\n" %(donor,acceptor,acceptor,hb_ang,hb_ang_tol)
    return st

def BulgeCircularHBondConstraints(donor,acceptor): # return string for cst file
    hb_ang = np.deg2rad(160.)
    hb_ang_tol=np.deg2rad(20.0)
    # bulge position
    st = "AtomPair N %i O %i HARMONIC 2.9 0.25\n" %(donor,acceptor)
    st+= "Angle N %i H %i O %i CIRCULARHARMONIC %3.1f %3.1f\n" %(donor,donor,acceptor,hb_ang,hb_ang_tol)
    st+= "Angle H %i O %i C %i CIRCULARHARMONIC %3.1f %3.1f\n" %(donor,acceptor,acceptor,hb_ang,hb_ang_tol)
    # bulge position +1
    st+= "AtomPair N %i O %i HARMONIC 3.4 0.3\n" %(donor+1,acceptor)
    st+= "Angle N %i H %i O %i CIRCULARHARMONIC %3.1f %3.1f\n" %(donor,donor,acceptor,hb_ang,hb_ang_tol)
    st+= "Angle H %i O %i C %i CIRCULARHARMONIC %3.1f %3.1f\n" %(donor,acceptor,acceptor,hb_ang,hb_ang_tol)
    return st

def PairConstraints(a,b,value,tol): # return string for cst file
    tag = 'pair-%s.%s' %(a,b)
    st = "AtomPair CA %i CA %i BOUNDED %3.1f %3.1f %3.1f 0.5 %s\n" %(a,b,value-tol,value+tol,tol/2,tag) 
    return st

def HarmonicPairConstraints(a,b,value,sd): # return string for cst file
    st = "AtomPair CA %i CA %i HARMONIC %3.1f %3.1f\n" %(a,b,value,sd)
    return st


def AngleConstraints(a,b,c,value,tol,tag): # return string for cst file
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)
    st = "Angle CA %i CA %i CA %i BOUNDED %3.1f %3.1f 0.5 %s\n" %(a,b,c,ang-ang_tol,ang+ang_tol,tag)
    return st

def HarmonicAngleConstraints(a,b,c,value,tol): # return string for cst file
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)
    st = "Angle CA %i CA %i CA %i HARMONIC %3.1f %3.1f\n" %(a,b,c,ang,ang_tol)
    return st


def CircularHarmonicAngleConstraints(a,b,c,value,tol): # return string for cst file
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)
    st = "Angle CA %i CA %i CA %i CIRCULARHARMONIC %3.1f %3.1f\n" %(a,b,c,ang,ang_tol)
    return st


def CstTypeAngleConstraints(a,b,c,value,tol,cst_type): # return string for cst file
    st=''	
    if cst_type=='harmonic':
	    ang = np.deg2rad(value)
	    ang_tol = np.deg2rad(tol)
	    st = "Angle CA %i CA %i CA %i HARMONIC %3.1f %3.1f\n" %(a,b,c,ang,ang_tol)
    elif cst_type=='bounded':
            ang = np.deg2rad(value)
            ang_tol = np.deg2rad(tol)
	    sd = ang_tol/2
	    tag="ang_%i.%i.%i" %(a,b,c)
	    st = "Angle CA %i CA %i CA %i BOUNDED %3.2f %3.2f %3.2f 0.5 %s\n" %(a,b,c,ang-ang_tol,ang+ang_tol,sd,tag)
    return st

def DihedralConstraints(a,b,value,tol,tag): # return string for cst file
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)    
    sd = ang_tol/2
    st = "Dihedral CB %i CA %i CA %i CB %i BOUNDED %3.2f %3.2f %3.2f 0.5 %s\n" %(a,a,b,b,ang-ang_tol,ang+ang_tol,sd,tag)
    return st

def CaDihedralConstraints(a,b,c,d,value,tol): # return string for cst file
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)
    sd = ang_tol/2
    tag="ang_%i.%i.%i.%i" %(a,b,c,d)
    st = "Dihedral CA %i CA %i CA %i CA %i BOUNDED %3.2f %3.2f %3.2f 0.5 %s\n" %(a,b,c,d,ang-ang_tol,ang+ang_tol,sd,tag)
    return st

def CircularHarmonicCaDihedralConstraints(a,b,c,d,value,tol): # return string for cst file
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)
    st = "Dihedral CA %i CA %i CA %i CA %i CIRCULARHARMONIC %3.1f %3.1f\n" %(a,b,c,d,ang,ang_tol)
    return st

def HarmonicDihedralConstraints(a,b,value,tol): # return string for cst file
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)
    st = "Dihedral CB %i CA %i CA %i CB %i HARMONIC %3.1f %3.1f\n" %(a,a,b,b,ang,ang_tol)
    return st

def CircularHarmonicDihedralConstraints(a,b,value,tol): # return string for cst file
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)
    st = "Dihedral CB %i CA %i CA %i CB %i CIRCULARHARMONIC %3.1f %3.1f\n" %(a,a,b,b,ang,ang_tol)
    return st

def HelixConstraints(bluefile):
	blue = Blueprint(bluefile)
	blue.reindex_blueprint(start=1)

	h1 = blue.segment_dict['H1']
	h1npos = int(h1.bp_data[0][0])
	h1cpos = int(h1.bp_data[-1][0])

	h3 = blue.segment_dict['H3']
	h3npos = int(h3.bp_data[0][0])
	h3cpos = int(h3.bp_data[-1][0])

	e3 = blue.segment_dict['E3']
	e3npos = int(e3.bp_data[0][0])
	e3cpos = int(e3.bp_data[-1][0])

	e4 = blue.segment_dict['E4']
	e4npos = int(e4.bp_data[0][0])
	e4cpos = int(e4.bp_data[-1][0])

        lines=''                

	# Nter E3
	atomlist1=[]
	atomlist1.append( blue.bp_data[e3npos-1][0] )
	atomlist1.append( blue.bp_data[e3npos-2][0] )

	# cter H1
	nh1 = len(h1.bp_data)
	atomlist2=[]
	for i in range(0,4):
        	index = (h1cpos-1) - i
        	atomlist2.append( blue.bp_data[index][0] )
	
	value = 8.0 ; tol=2.0 ## PARAMETERS
	for a in atomlist1:
		for b in atomlist2[:2]:
			st = PairConstraints(a,b,value,tol) ; lines+=st # C-ter H1 attached to E3 (Nter)

	# nter H1
	nh1 = len(h1.bp_data)
	atomlist1=[]
	for i in range(0,1):
        	index = (h1npos-1) +  i
        	atomlist1.append( blue.bp_data[index][0] )


	ne3 = len(e3.bp_data)
	atomlist2=[]
	for i in range(0,1):
        	index = (e3cpos-1) - i
        	atomlist2.append( blue.bp_data[index][0] )
        	index = (e4npos+1) + i
        	atomlist2.append( blue.bp_data[index][0] )	

	value = 8.0 ; tol=2.0 # PARAMETERS
	for a in atomlist1:
		for b in atomlist2:
			st = PairConstraints(a,b,value,tol) ; lines+=st # Nter of H1 attached to Cter E3 and Nter E4

	
	return lines	


def AmbiguousCst(cst_lst):
        header = 'AmbiguousConstraint\n'
        for cst in cst_lst:
                header += cst
        header += 'END_AMBIGUOUS\n'
        return header


def MultiCst(cst_lst):
        header = 'MultiConstraint\n'
        for cst in cst_lst:
                header += cst
        header += 'END_MULTI\n'
        return header

def ConstraintsStrandCurvature(**kwargs):
	segment = kwargs.get("strand")
	positions = kwargs.get("positions")
	bend = float( kwargs.get("bend") )
	bend_tol = float( kwargs.get("bend_tol") )
	bend_bulge = float( kwargs.get("bend_bulge") )
        twist = float( kwargs.get("twist") )
        twist_tol = float( kwargs.get("twist_tol") )
	bluefile = kwargs.get("bluefile")

	blue = Blueprint(bluefile)
	blue.reindex_blueprint(start=1)
	seg = blue.segment_dict[segment]
	blue = Blueprint(bluefile)
	cst_st=''
	# bending
	if bend:
	   if positions == None:
		positions = range( 2,len(seg.bp_data)-2)
           for i in positions:
		pos = seg.bp_data[i][0]
		if seg.bp_data[i][2] == 'EA': # bulge
			st = AngleConstraints(pos-2,pos,pos+2,180-bend_bulge,bend_tol,"bend%s.%i" %(segment,pos))
		else: # non-bulged
			st = AngleConstraints(pos-2,pos,pos+2,180-bend,bend_tol,"bend%s.%i" %(segment,pos))
		cst_st += st
	  
	   pos1 = seg.bp_data[0][0]
	   pos2 = seg.bp_data[-1][0]
	   if len(seg.bp_data) % 2 ==0:
		cen1 = pos1 + len(seg.bp_data)/2
		cen2 = pos1 + len(seg.bp_data)/2 + 1
		st = AngleConstraints(pos1,cen1,pos2,180-bend_bulge,bend_tol,"edgebend%s.%i" %(segment,cen1))
		st = AngleConstraints(pos1,cen2,pos2,180-bend_bulge,bend_tol,"edgebend%s.%i" %(segment,cen2))

	   else:
		cen = pos1 + len(seg.bp_data)/2
		st = AngleConstraints(pos1,cen,pos2,180-bend_bulge,bend_tol,"edgebend%s.%i" %(segment,cen))
	   cst_st += st

	   st = AngleConstraints(pos-2,pos,pos+2,180-bend_bulge,bend_tol,"bend%s.%i" %(segment,pos))
	# twisting
	if twist:
	   if positions == None:
		positions = range( len(seg.bp_data)-2)
	   for i in positions:
		pos1 = seg.bp_data[i][0]
		pos2 = pos1+2
                st = DihedralConstraints(pos1,pos2,twist,twist_tol,'dih%s.%i' %(segment,pos1))
		cst_st += st
	return cst_st


def RegularStrandCurvature(**kwargs):
	segment = kwargs.get("strand")

	level = kwargs.get("level") # 1 or 2

	bend = kwargs.get("global_bend",None)
	bend_tol = kwargs.get("global_bend_tol",10.0) 
        twist = kwargs.get("global_twist",None )
        twist_tol = kwargs.get("global_twist_tol",5.0 )

        bend_area = kwargs.get("bend_area_value",None )
        bend_area_value = kwargs.get("bend_area_value",None )
        bend_area_value_tol = kwargs.get("bend_area_value_tol",5.0)  # n-ter, middle, c-ter

        twist_area = kwargs.get("twist_area")  # n-ter, middle, c-ter
        twist_area_value = kwargs.get("twist_area_value",None )
        twist_area_value_tol = kwargs.get("twist_area_value_tol",5.0 )

	bend_positions = kwargs.get("bend_positions",None )
	bend_positions_value = kwargs.get("bend_positions_value",None )
	bend_positions_value_tol = kwargs.get("bend_positions_value_tol",None )

	constraint_type = kwargs.get("constraint_type","harmonic" )

	blue = kwargs.get("blueprint")
	blue.reindex_blueprint(start=1)
	seg = blue.segment_dict[segment]
	cst_st=''
	#################
	n = len(seg.bp_data)        

	############
	# BENDING
	############

	#----------------------
	# Set all triads for bend calculation
	#----------------------
	bend_triads = []
	step=level*2
	for i in range(0,n-step):
		 if i+step*2 < n:
			pos1 = i
			pos2 = i + step*1
			pos3 = i + step*2
		        bend_triads.append([pos1,pos2,pos3])

        #----------------------
        # Define bend areas
        #----------------------
	# By positions. This is used especially for positions paired to bulges, where we expect higher bend than the rest of triads
	bend_positions_triads=[]
	if bend_positions:
		for relpos in bend_positions:
			if relpos < 0: # when giving rel position from C-terminal (for E3)
				pos = n+relpos
			else:
				pos=relpos
			for triad in bend_triads:
				if pos==triad[1]: # central position of triad
					bend_positions_triads.append( triad )

	# By areas					
	bend_area_triads=[]
        if bend_area =='n-term':
                bend_area_triads = [ bend_triads[0] ]
        elif bend_area =='c-term':
                bend_area_triads = [ bend_triads[-1] ]
        elif bend_area =='center':
                index =  len(bend_triads)/2 - 1
                if len(bend_triads) % 2 == 0:
                        bend_area_triads = [ pair for k in bend_triads[index:index+2] ]
                else:
                        bend_area_triads = [ bend_triads[index] ]

        ############
        # TWIST
        ############

      	#----------------------
       	# Set all pairs for twist calculation
	#----------------------
	twist_pairs = []
	step=level*2
	for i in range(0,n-step):
		if i+step < n:
			pos1 = i
			pos2 = i + step*1
			twist_pairs.append([pos1,pos2])       

	#----------------------
	# Define twist areas
	#----------------------
	twist_area_pairs=[]
	if twist_area =='n-term':
		twist_area_pairs = [ twist_pairs[0] ]
	elif twist_area =='c-term':
		twist_area_pairs = [ twist_pairs[-1] ]
        elif twist_area =='center':
		index =  len(twist_pairs)/2 - 1
		if len(twist_pairs) % 2 == 0:
			twist_area_pairs = [ pair for pair in twist_pairs[index:index+2] ]
		else:
			twist_area_pairs = [ twist_pairs[index] ]
				
	#########################
	# INTRODUCE CONSTRAINTS
	#########################
    	# After indentifying combination of positions for bending and twist... Put constraints
	# Calculate Bends
	for triad in bend_triads:
		a,b,c = triad	
		pos1 = seg.bp_data[b][0]
		pos2 = seg.bp_data[a][0]
		pos3 = seg.bp_data[c][0]         
		if bend_positions and triad in bend_positions_triads:
			st = CstTypeAngleConstraints(pos2,pos1,pos3,180-bend_positions_value,bend_positions_value_tol,constraint_type) ; cst_st += st
		elif bend_area_value and triad in bend_area_value_triads and triad not in bend_positions_triads:
			st = CstTypeAngleConstraints(pos2,pos1,pos3,180-bend_area_value,bend_area_value_tol,constraint_type) ; cst_st += st
		elif bend:
			st = CstTypeAngleConstraints(pos2,pos1,pos3,180-bend,bend_tol,constraint_type) ; cst_st += st

	# Calculate Twists
	e_kink = EBulged(seg)
	for pair in twist_pairs:
		a,b = pair
		pos1 = seg.bp_data[a][0]
		pos2 = seg.bp_data[b][0]
		# discard dyads with gly kink
		if e_kink == pos1 or e_kink == pos2:
			continue
		if twist_area_value and pair in twist_area_pairs:
			st = DihedralConstraints(pos1,pos2,twist_area_value,twist_area_value_tol,'dih%i.%i' %(pos1,pos2)) ; cst_st += st
		elif twist:
			st = DihedralConstraints(pos1,pos2,twist,twist_tol,'dih%i.%i' %(pos1,pos2)) ; cst_st += st 

	return cst_st

def BulgedStrandCurvature(**kwargs):
	# Remember to use the blueprint bp.b with abegos, because this is used for checking whether the strand is bulged.
	level = kwargs.get("level") # 1 or 2
	segment = kwargs.get("strand")
	bend = kwargs.get("bend",None)
	bend_tol = kwargs.get("bend_tol",5.0) 
        twist = kwargs.get("twist",None)
        twist_tol = kwargs.get("twist_tol",5.0)
	blue = kwargs.get("blueprint")
	constraint_type = kwargs.get("constraint_type","harmonic" )

	seg = blue.segment_dict[segment]
	nres = len(seg.bp_data)
	cst_st=''

	pos1 = Bulged(seg)
	pos2 = pos1+1
	bulgepos_fromN = pos1-seg.bp_data[0][0]+1
	bulgepos_fromC = pos1-seg.bp_data[-1][0]-1
	
	#-------------
	# Bending
	#-------------
	if bend:
		if level==1:
			if nres > 6 and ( bulgepos_fromN >=3 and bulgepos_fromC <=-4 ):
				st = CstTypeAngleConstraints(pos1-2,pos1,pos1+3,180-bend,bend_tol,constraint_type) ; cst_st += st
				st = CstTypeAngleConstraints(pos2-3,pos2,pos2+2,180-bend,bend_tol,constraint_type) ; cst_st += st
				
		if level==2:
			if nres >=10 and ( bulgepos_fromN >=5 and bulgepos_fromC <=-6 ):
				st = CstTypeAngleConstraints(pos1-4,pos1,pos1+5,180-bend,bend_tol,constraint_type) ; cst_st += st
				st = CstTypeAngleConstraints(pos2-5,pos2,pos2+4,180-bend,bend_tol,constraint_type) ; cst_st += st
				#if EBulged(seg) and ( EBulged(seg)-pos1 ) == 2: # ABE bulge
				#	pos3 = EBulged(seg)
					#st = CstTypeAngleConstraints(pos3-6,pos3,pos3+3,180-bend,bend_tol,constraint_type) ; cst_st += st
				#	bend1 = 180-bend
				#	st = CstTypeAngleConstraints(pos2-5,pos2,pos2+4,bend1,bend_tol,constraint_type) ; cst_st += st
				#	st = CstTypeAngleConstraints(pos2,pos2-5,pos2+4,(180-bend1)/2.0,bend_tol,constraint_type) ; cst_st += st
				#	st = CstTypeAngleConstraints(pos2-5,pos2+4,pos2-5,(180-bend1)/2.0,bend_tol,constraint_type) ; cst_st += st
				#else:
                                #	st = CstTypeAngleConstraints(pos1-4,pos1,pos1+5,180-bend,bend_tol,constraint_type) ; cst_st += st
	                        #        st = CstTypeAngleConstraints(pos2-5,pos2,pos2+4,180-bend,bend_tol,constraint_type) ; cst_st += st	

        #----------------------
        # Twist
        #----------------------
	if twist:
		# Set all pairs
	        twist_pairs = []
		ini=seg.bp_data[0][0]
		end=seg.bp_data[-1][0]
	        step=level*2
		p1=pos1-1 - 2*(level-1) # left
		p2=pos1+2 + 2*(level-1) # right
		twist_pairs.append([p1,p2])

		# ignore twist with gy kinks as they dont have CB
		e_kink = EBulged(seg)
		for pair in twist_pairs:
			a,b=pair
			if e_kink and e_kink in pair:
                        	continue
			else:
				st = DihedralConstraints(a,b,twist,twist_tol,'dih%i.%i' %(a,b)) ; cst_st += st

	return cst_st


def ProlineAbegoB_PhiPsiConstraints(pos):
        phi_min = -75.0 - 5.0 
        phi_max = -75.0 + 5.0
        psi_min = 150.0 - 5.0
        psi_max = 150.0 + 5.0

	st=''      
	st+="Dihedral C %i N %i CA %i C %i BOUNDED %3.1f %3.1f 0.1 phi_%s\n" %(pos-1,pos,pos,pos,np.deg2rad(phi_min),np.deg2rad(phi_max),pos)
	st+="Dihedral N %i CA %i C %i N %i BOUNDED %3.1f %3.1f 0.1 psi_%s\n" %(pos,pos,pos,pos+1,np.deg2rad(psi_min),np.deg2rad(psi_max),pos)

	# Complement with omega constraint
        #st+="Dihedral CA %i C %i N %i CA %i BOUNDED %3.2f %3.2f 0.1 omega_%s\n" %(pos-1,pos-1,pos,pos,np.deg2rad(170.0),np.deg2rad(180.0),pos)
        #st+="Dihedral CA %i C %i N %i CA %i BOUNDED %3.2f %3.2f 0.1 omega_%s\n" %(pos,pos,pos+1,pos+1,np.deg2rad(170.0),np.deg2rad(180.0),pos+1)

	return st


def GlycineExtended_PhiPsiConstraints(pos):
        phi_min = -180.0 - 20.0
        phi_max = -180.0 + 20.0
        psi_min = 180.0 - 20.0
        psi_max = 180.0 + 20.0

        st=''
        st+="Dihedral C %i N %i CA %i C %i BOUNDED %3.1f %3.1f 0.3 phi_%s\n" %(pos-1,pos,pos,pos,np.deg2rad(phi_min),np.deg2rad(phi_max),pos)
        st+="Dihedral N %i CA %i C %i N %i BOUNDED %3.1f %3.1f 0.3 psi_%s\n" %(pos,pos,pos,pos+1,np.deg2rad(psi_min),np.deg2rad(psi_max),pos)

	return st

def AbegoConstraintSegment(segments,blue):
	blue.reindex_blueprint(start=1)
	cst=''
	for segment in segments:
		seg = blue.segment_dict[segment]
		for res in seg.bp_data:
			pos = res[0]
			if len(blue.bp_data[pos-1][2])==2: # abego specified.
				st = AbegoPhiPsiConstraints(pos,blue)
				cst+=st
	return cst

def HairpinPairingResidues(blue,segment1,segment2,shift):
    # segment1 is the first strand in the hairpin according to sequence
    s1 = blue.segment_dict[segment1]
    s2 = blue.segment_dict[segment2]
    map1={} ; map2={}
    pairs=[]

    if Bulged(s1) and Bulged(s2):
        b1pos = Bulged(s1)
        shift1 = int(s1.bp_data[-1][0]) - b1pos
        b2pos = Bulged(s2)
        shift2 = b2pos - int(s2.bp_data[0][0])
	# shifts from the hairpin
	count1=-1
	count2=-1
        while count2 < len(s2.bp_data)-1:
	    if count2 < shift2:
		count2+=1
	    elif count2==shift2: # bulge s2 found
		count2+=2
	    elif count2 > shift2: # ignore position bpos+1 for pairing
		count2+=1
	
	    if count1 < shift1-2:
		count1+=1
            elif count1==shift1-2: # bulge s1 found
		count1+=2
            elif count1 > shift1-2: # bulge s1 found
		count1+=1	    

            pos2 = int(s2.bp_data[count2][0])
            pos1 = int(s1.bp_data[-1-count1][0])	    

            pairs.append([pos1,pos2])

	    if count1 == len(s1.bp_data)-1:
		break

    elif Bulged(s1):
	b1pos = Bulged(s1)
	shift1 = int(s1.bp_data[-1][0]) - b1pos
	count1=-1
	count2=-1
	while count2 < len(s2.bp_data)-1:
            count2+=1
            if count1 < shift1-2:
                count1+=1		
            elif count1==shift1-2: # bulge s1 found
                count1+=2
            elif count1 > shift1-2: # bulge s1 found
                count1+=1
            pos2 = int(s2.bp_data[count2][0])
            pos1 = int(s1.bp_data[-1-count1][0])

            pairs.append([pos1,pos2])	    

            if count1 == len(s1.bp_data)-1:
                break
	
    elif Bulged(s2):
        b2pos = Bulged(s2)
        shift2 = b2pos - int(s2.bp_data[0][0])
        count1=-1
        count2=-1
        while count2 < len(s2.bp_data)-1:
	    count1+=1
            if count2 < shift2:
                count2+=1
            elif count2==shift2: # bulge s2 found
                count2+=2
            elif count2 > shift2: # ignore position bpos+1 for pairing
                count2+=1
            pos2 = int(s2.bp_data[count2][0])
            pos1 = int(s1.bp_data[-1-count1][0])

            pairs.append([pos1,pos2])

            if count1 == len(s1.bp_data)-1:
                break	

#    elif Bulged(s1): 
#        b1pos = Bulged(s1)
#        shift = int(s1.bp_data[-1][0]) - b1pos
#        for i in range(len(s2.bp_data)):
#            pos2 = int(s2.bp_data[i][0]) # Antiparallel
#            if i < shift:
#                pos1 = int(s1.bp_data[-1-i][0])
#            elif i+1 < len(s1.bp_data):
#                pos1 = int(s1.bp_data[-1-i-1][0])
#            else:
#                break
#            pairs.append([pos1,pos2])
    # The bulgepos+1 is the one paired to the second strand

#    elif Bulged(s2):
#        b2pos = Bulged(s2)
#        shift = b2pos - int(s2.bp_data[0][0])
#        for i in range(len(s1.bp_data)):
#            pos1 = int(s1.bp_data[-1-i][0])
#            if i < shift:
#                pos2 = int(s2.bp_data[i][0])
#            elif i+1 < len(s2.bp_data):
#                pos2 = int(s2.bp_data[i+1][0])
#            else:
#                break
#            pairs.append([pos1,pos2])




    elif Bulged(s1)==None and Bulged(s2)==None:
	 if shift<0:
         	for i in range(len(s2.bp_data)):
	            pos2 = int(s2.bp_data[i][0]) # Antiparallel     
        	    if i < (len(s1.bp_data)+shift):
                	pos1 = int(s1.bp_data[-1+shift-i][0])
	            else:
        	        break
		    #print pos1,pos2
	            pairs.append([pos1,pos2])
         elif shift>=0:
		i=0
                for j in range(shift,len(s2.bp_data)):
                    pos2 = int(s2.bp_data[j][0]) # Antiparallel
                    if i < len(s1.bp_data):
                        pos1 = int(s1.bp_data[-1-i][0])
                    else:
                        break
                    pairs.append([pos1,pos2])
		    i+=1



    # correct for shift
    # this simple solution for bulges only works for abs(shift) = 1 # for larger shifts it will be needed a list of up/down.


    if Bulged(s1)==None and Bulged(s2)==None:
	pass
    else:
     for pair in pairs:
        if shift < 0:
                count=1
                while count<=abs(shift):
                        if Bulged(s1) and ( pair[0] == Bulged(s1)+2 ):
                                pair[0]-=2
                        else:
                                pair[0]-=1
                        count+=1
        elif shift > 0:
		count=1
		while count<=abs(shift):
			if Bulged(s2) and ( pair[1] == Bulged(s2) ):
				pair[1]+=2
			else:
				pair[1]+=1
			count+=1
			



    return pairs

def BulgedStrandOrientations(**kwargs):
	# returns list of indices within the segment.
        strand = kwargs.get('strand')
        blue = kwargs.get('blueprint')
        seg = blue.segment_dict[strand]
        downside=[] ; upside=[] # indices
        if isinstance(Bulged(seg),list):
                b1pos = Bulged(seg)[0]
                b2pos = Bulged(seg)[1]
                shift1 = int(seg.bp_data[-1][0]) - b1pos
                shift2 = int(seg.bp_data[-1][0]) - b2pos

                start = seg.bp_data[0][0]
                bulge_index1 = b1pos - start
                # left side
                for i in range(bulge_index1,-1,-2):
                    downside.append(i)
                for i in range(bulge_index1-1,-1,-2):
                    upside.append(i)

                k=bulge_index1+1
                c=0
                while seg.bp_data[k][0] <= seg.bp_data[-1][0]:
                        if seg.bp_data[k][0] != b2pos:
                                if c % 2==0:
                                        downside.append(k)
                                else:
                                        upside.append(k)
                                c+=1
                        else:
                                downside.append(k)
                                c=0
                        k+=1
                        if k==len(seg.bp_data): break

        else:
		bulgepos=Bulged(seg)
                start = seg.bp_data[0][0]
                bulge_index = bulgepos - start
                # left side
                for i in range(bulge_index,-1,-2):
                    downside.append(i)
                for i in range(bulge_index-1,-1,-2):
                    upside.append(i)
                # right side
                for i in range(bulge_index+1,len(seg.bp_data),2):
                    downside.append(i)
                for i in range(bulge_index+2,len(seg.bp_data),2):
                    upside.append(i)

        upside.sort()
        downside.sort()

        return upside, downside


def HbondsBulgedStrand(**kwargs):
        strand1 = kwargs.get('strand1')
        strand2 = kwargs.get('strand2')
        blue = kwargs.get('blueprint')
	offset = kwargs.get('offset_s1',0)
	shift = kwargs.get('shift',0)

        pair1 = HairpinPairingResidues(blue,strand1,strand2,shift)
	#pair1 = pair1[offset:] # it is more convenient just to change the pairings.
	# in regular hairpin it is easier to change the inipos, which is the reference for hbond. Whereas here the reference is the bulge position

	#print pair1

        seg1 = blue.segment_dict[strand1]
        seg2 = blue.segment_dict[strand2]

	if Bulged(seg1) and Bulged(seg2):
		b1pos = Bulged(seg1)
		b2pos = Bulged(seg2)
		seg=seg1
        else:
		if Bulged(seg1):
                	seg = seg1
	        elif Bulged(seg2):
        	        seg = seg2
	        else:
        	        print 'Warning: None of the strands is bulged'		
		b1pos = Bulged(seg)
		b2pos = None

        hblist=[]
        # to the left of the bulge
        i=0
        pos = b1pos-2
        while pos >= seg.bp_data[0][0]:
            hblist.append(pos)
	    pos -= 2	

	hblist.append(b1pos)

        pos = b1pos+3
        # to the right of the bulge
        while pos <= seg.bp_data[-1][0]:
            hblist.append(pos)
	    pos+=2

	#print pair1
	#print hblist
        hbpairs=[]
        for pair in pair1:
            hbpair = list( set(hblist).intersection(set(pair)) )
            if len(hbpair)==1:
                hbpairs.append(pair)
	#print hbpairs

        cst_st = ''
        for hbpair in hbpairs:
                pos1,pos2 = hbpair
                # Add hbond for bulge position (which is not C-alpha paired)
                if b1pos in hbpair: # position b1pos+1 is in the pairing list (not the bulgepos, so we add the hbond here)
                        # identify paired position
                        paired_pos = list( set(hbpair).difference(set([b1pos])) )[0]
			st = BulgeCircularHBondConstraints(b1pos,paired_pos) ; cst_st += st
			st = CircularHBondConstraints(paired_pos,b1pos+1) ; cst_st += st
                elif b2pos != None and b2pos in hbpair: # position b1pos+1 is in the pairing list (not the bulgepos, so we add the hbond here)
                        # identify paired position
                        paired_pos = list( set(hbpair).difference(set([b2pos])) )[0]
			st = BulgeCircularHBondConstraints(b2pos,paired_pos) ; cst_st += st
			st = CircularHBondConstraints(paired_pos,b2pos+1) ; cst_st += st
		else:
	                st = CircularHBondConstraints(pos1,pos2) ; cst_st += st
	                st = CircularHBondConstraints(pos2,pos1) ; cst_st += st


	return cst_st




def HbondingResiduesInBulgedStrand(**kwargs):
        seg = kwargs.get('strand_segment')
        blue = kwargs.get('blueprint')

        b1pos = Bulged(seg)

        hblist=[]
        # to the left of the bulge
        i=0
        pos = b1pos
        while pos > seg.bp_data[0][0]:
            pos = (b1pos-2)-2*i
            i+=1
            hblist.append(pos)

        pos = b1pos
        # to the right of the bulge
        i=0
        while pos < seg.bp_data[-1][0]:
            pos = (b1pos+1)+2*i
            i+=1
            hblist.append(pos)


        return hblist


def AllSheetSegmentPairs(blue):
        cst_st=''
        pair1 = HairpinPairingResidues(blue,'E1','E2')
        pair2 = HairpinPairingResidues(blue,'E2','E3')
        pair3 = HairpinPairingResidues(blue,'E3','E4')

        pairs=[]
        pairs.extend(pair1)
        pairs.extend(pair2)
        pairs.extend(pair3)		

	dic_pairs={}
	for pair in pairs:
		a,b = pair
		seg_a = blue.residue_segment(a)
		seg_b = blue.residue_segment(b)
		dic_pairs.setdefault(a,{})
		dic_pairs.setdefault(b,{})
		dic_pairs[a][seg_b]=b
		dic_pairs[b][seg_a]=a
	
	return dic_pairs		

	
def HbondsRegularHairpin(**kwargs):
        strand1 = kwargs.get('strand1')
        strand2 = kwargs.get('strand2')
        blue = kwargs.get('blueprint')
        offset = kwargs.get('offset_s1',0)
        n_hbpairs = kwargs.get('n_hbpairs','all')
        shift = kwargs.get('shift',0)
	exclude_hbs = kwargs.get('exclude_hbs',[])

        pair1 = HairpinPairingResidues(blue,strand1,strand2,shift)
	#print pair1

        seg1 = blue.segment_dict[strand1]
        seg2 = blue.segment_dict[strand2]


        # to the left of the bulge
        inipos = seg1.bp_data[-1][0] - offset # start counting from hairpin loop
        #if shift <=0:
        #       inipos = seg1.bp_data[-1+shift][0] - offset

        pos=inipos
        hblist=[pos]
        # all residues of seg1 are hbonded to seg2

        while pos >= seg1.bp_data[0][0]+2:
            pos -= 2
            if pos != Bulged(seg1): # in case this hbond pattern is used with a bulged hairpin.
                    hblist.append(pos)

        hbpairs=[]
        for pair in pair1:
            hbpair = list( set(hblist).intersection(set(pair)) )
            if len(hbpair)==1:
                hbpairs.append(pair)

        #print hbpairs

        if n_hbpairs == 'all':
                nhbonds_max = len(hbpairs)
        else:
                nhbonds_max = n_hbpairs


        cst_st = ''
        exclude_hbs_positive = [len(hbpairs)+k for k in exclude_hbs]
        for i,hbpair in enumerate(hbpairs):
                if len(exclude_hbs) > 0:
                        if ( (i in exclude_hbs) or (i in exclude_hbs_positive) ):
                                continue

                pos1,pos2 = hbpair
                if blue.bp_data[pos1-1][1] != 'P':
                        st = CircularHBondConstraints(pos1,pos2) ; cst_st += st
                if blue.bp_data[pos2-1][1] != 'P':
                        st = CircularHBondConstraints(pos2,pos1) ; cst_st += st
                if i+1==nhbonds_max:
                        break

        return cst_st




	
def FlatSheetConstraints(blue):
        cst_st=''
        pair1 = HairpinPairingResidues(blue,'E1','E2')
        pair2 = HairpinPairingResidues(blue,'E2','E3')
        pair3 = HairpinPairingResidues(blue,'E3','E4')
        # Get Triads
        pairs=[]
        pairs.extend(pair1)
        pairs.extend(pair2)
        pairs.extend(pair3)
        triads=[]
        for i,p1 in enumerate(pairs):
            for j,p2 in enumerate(pairs):
                if j>i:
                    diff = set(p1).difference(set(p2))
                    if len(diff) == 1:
                        t = list( set(p1).union(set(p2)) )
                        triads.append(t)

        quarts=[]
        for i,p1 in enumerate(triads):
            for j,p2 in enumerate(triads):
                if j>i:
                    diff = set(p1).difference(set(p2))
                    if len(diff) == 1:
                        t = list( set(p1).union(set(p2)) )
                        quarts.append(t)


        # Set Constraints
        for quart in quarts:
            quart.sort()
            a,b,c,d = quart
            st = HarmonicAngleConstraints(a,b,c,170,5.0) ; cst_st += st
            st = HarmonicAngleConstraints(a,c,d,170,5.0) ; cst_st += st

        return cst_st


def Alternate(a):
    if 1 in a:
        direction=1
        if a.index(1) %2 == 0:
           even=True
        else:
           even=False
    else:
        direction=-1
        if a.index(-1) %2 == 0:
           even=True
        else:
           even=False


    for i in range(len(a)):
     if i%2==0:
        if even:
           a[i]=direction
        else:
           a[i]=-direction
     else:
        if even:
           a[i]=-direction
        else:
           a[i]=direction



def L6Loops(orient_l6):
        a,b = orient_l6
        abego=[]
        if a==1 and b==1:
		#l6=['ABABB','BBG','ABB']
		l6=['ABABB']
        elif a==1 and b==-1:
                l6=['BBGB','BABB']
        elif a==-1 and b==1:
                l6=['BABB','BBGB']
	elif a==-1 and b==-1:
		print ' down-down loop: ignore this topology'
		l6=[]

	return l6
	




def OrientSegmentsFromLoops(segment1,segment2,loop,sidech):
	loop_orient = {'ABB':'11','BBG':'11','BABB':'01', 'BBGB':'10','BBABB':'11','BBGBB':'11','ABABB':'11'}

	if loop_orient[loop]=='11':
		sidech[segment1][-1]=1
		sidech[segment2][0]=1

        elif loop_orient[loop]=='00':
                sidech[segment1][-1]=-1
                sidech[segment2][0]=-1

        elif loop_orient[loop]=='10':
                sidech[segment1][-1]=1
                sidech[segment2][0]=-1

        elif loop_orient[loop]=='01':
                sidech[segment1][-1]=-1
                sidech[segment2][0]=1

