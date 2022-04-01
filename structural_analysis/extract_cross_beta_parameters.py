import numpy as np
from math import ceil,atan2
import matplotlib
import matplotlib.pyplot as plt
import glob,os,copy
import scipy
from scipy.spatial.transform import Rotation as R
pyrosetta.init("-mute all")
from xyzMath3 import *

def ScaleXYZVector(xyzvector,factor):
    xyzvector.x*=factor
    xyzvector.y*=factor
    xyzvector.z*=factor

def MidPosition(pose,p1,p2):
    mid_xyz= pose.residue(p1).atom("CA").xyz() + pose.residue(p2).atom("CA").xyz()
    ScaleXYZVector(mid_xyz, 1/2)
    return mid_xyz


def MidPositionList(pose,poslist):
    cm=numeric.xyzVector_double_t()
    for pos in poslist:
        cm+=pose.residue(pos).atom("CA").xyz()
    ScaleXYZVector(cm, 1/len(poslist))
    return cm

def MakeBlueprintFromPose(pose):
        DSSP=protocols.moves.DsspMover()
        DSSP.apply(pose)    # populates the pose's Pose.secstruct
        ss = pose.secstruct()
        seq = pose.sequence()
        # Get abegos
        abego = core.sequence.get_abego(pose) # list from 1 to nres
        abego = list(abego)
        # build list
        res_list=[]
        for k,ss_res in enumerate(ss):
                res_list.append( [k+1,seq[k],ss_res+abego[k],'.'] )

        blue = Blueprint(data=res_list)
        return blue


def Angle(v1,v2):
    scalar = v1.dot(v2) / ( np.linalg.norm(v1) * np.linalg.norm(v2) )
    return np.rad2deg( np.arccos(scalar) )

def StrandDirectionTriad(pose,triad):
        a,b,c = triad
        cen1 = np.array( pose.residue( a ).xyz("CA") + pose.residue( b ).xyz("CA") ) / 2.0
        cen2 = np.array( pose.residue( b ).xyz("CA") + pose.residue( c ).xyz("CA") ) / 2.0
        #cen1 = pose.residue( a ).xyz("CA") + pose.residue( b ).xyz("CA")
        #cen2 = pose.residue( b ).xyz("CA") + pose.residue( c ).xyz("CA")
        #ScaleXYZVector(cen1,0.5)
        #ScaleXYZVector(cen2,0.5)
        v = cen2-cen1
        return v

def CheckStrandBreak(pose,blue,seg):
    v_map={}
    for k in range(0,len(seg.bp_data)-3,3):
        triad=[seg.bp_data[k][0],seg.bp_data[k+1][0],seg.bp_data[k+2][0]]
        v = StrandDirectionTriad(pose,triad)
        v_map[triad[0]]=v

    # check directions
    keys=list(v_map.keys())
    keys.sort()
    start_v = v_map[keys[0]]
    start_v2 = v_map[keys[0]]

    break_position=0
    for key in keys[1:]:
        dot_prod = np.dot(start_v,v_map[key])
        if dot_prod<0:
            break_position=key
            break

    return break_position


def StrandPairing(pose,blue,pairing="A"):
    dssp=core.scoring.dssp.Dssp(pose)
    pairset=core.scoring.dssp.StrandPairingSet(pose)
    strand_dic_pairs={} ; strand_dic_a_pairs={} ; strand_dic_p_pairs={} ; sspair_lst=[] ; sspair_info={}
    for i in range( 1, pairset.size()+1 ):
        sp=pairset.strand_pairing(i)
        if pairing=="A" and not sp.antiparallel():
            continue
        elif pairing=="P" and sp.antiparallel():
            continue
        start=sp.begin1()
        end=sp.end1()
        segment = blue.residue_segment(start)
        partner = sp.get_pair(start)
        partner2 = sp.get_pair(end)
        segment2 = blue.residue_segment(partner)

        # save residue ranges for strand pairing.
        range1=[start,end]
        range2=[partner,partner2]
        st_pair="%s-%s" %(segment,segment2)
        sspair_info[st_pair] = [range1,range2]


        if "E" in segment and "E" in segment2:
            sspair_lst.append([segment,segment2])

            strand_dic_pairs.setdefault(segment,[])
            strand_dic_pairs[segment].append(segment2)
            strand_dic_pairs.setdefault(segment2,[])
            strand_dic_pairs[segment2].append(segment)

    # organise strands by beta sheets
    sspairs=copy.deepcopy(sspair_lst) # to avoid modifying the original list below.
    sheets={}
    if len(sspair_lst) > 0:
        n_sheets=1
        sheets[n_sheets]=sspairs[0]
        for sspair in sspairs[1:]:
            strand1,strand2 = sspair
            flag=True
            for sheet_number in sheets.keys():
                sheet_strands = sheets[sheet_number]
                if strand1 in sheet_strands and strand2 not in sheet_strands:
                    sheets[sheet_number].append(strand2)
                    flag=False
                elif strand2 in sheet_strands and strand1 not in sheet_strands:
                    sheets[sheet_number].append(strand1)
                    flag=False
            if flag: # make new sheet
                n_sheets+=1
                sheets[n_sheets]=sspair

    return sspair_lst,strand_dic_pairs,sspair_info


def CleanBlueprintByStrandBreaks(pose,blue):
    #################################
    # Now redefine blueprint strands according to breaks (loops recognized  as strands)
    # this will help to better identify beta sheets, and avoid recognize sandwiches as a single sheet due to a long curved strand.
    #################################
    ini_strand_list=[k for k in blue.segment_list() if 'E' in k]
    strand_breaks=[]
    for strand in ini_strand_list:
        # check strand breaks to correctly assign beta sheets
        seg=blue.segment_dict[strand]
        if seg.length() >= 10:
            break_pos = CheckStrandBreak(pose,blue,seg)
            if break_pos !=0:
                strand_breaks.append(break_pos)
        if seg.length()<=3:
            for pos in seg.bp_data:
                abego=pos[2][1]
                pos[2]="L%s" %abego

    for break_position in strand_breaks:
        blue.bp_data[break_position-1][2]="LX"
        blue.bp_data[break_position-2][2]="LX"
        blue.bp_data[break_position][2]="LX"

    blue.refresh()

def dihedral_xyz(p1,p2,p3,p4):
    a = p2 - p1
    b = p3 - p2
    c = p4 - p3

    x = -np.dot(a/np.linalg.norm(a),c/np.linalg.norm(c)) + np.dot(a/np.linalg.norm(a),b/np.linalg.norm(b)) * np.dot(b/np.linalg.norm(b),c/np.linalg.norm(c))
    cross_bc = np.cross(b/np.linalg.norm(b),c/np.linalg.norm(c))
    y = np.dot(a/np.linalg.norm(a),cross_bc)

    return np.rad2deg( atan2(y,x) )


def SimpleCrossBetaParamsOnPairing(pose,blue,cb_segments,sspairing):
    seg1 = blue.segment_dict[cb_segments[0]] # b-arch1
    seg2 = blue.segment_dict[cb_segments[1]] # b-arch1
    seg3 = blue.segment_dict[cb_segments[2]] # b-arch2
    seg4 = blue.segment_dict[cb_segments[3]] # b-arch2

    # cross beta pairing 1
    s1_range,s3_range = sspairing["%s-%s" %(cb_segments[0],cb_segments[2])]
    # cross beta pairing 2
    s2_range,s4_range = sspairing["%s-%s" %(cb_segments[1],cb_segments[3])]

    # sspairing info ensure correct ranges for calcularing twist.
    e1_start_ss,e1_end_ss = s1_range
    e3_end_ss,e3_start_ss = s3_range # anti
    e2_start_ss,e2_end_ss = s2_range
    e4_end_ss,e4_start_ss = s4_range # anti

    #print(s1_range)
    #print(s3_range)
    #print(s2_range)
    #print(s4_range)
    #print(sspairing)

    mid_e1_e3_n = MidPosition(pose,e1_start_ss,e3_end_ss)
    mid_e1_e3_c = MidPosition(pose,e1_end_ss,e3_start_ss)
    mid_e2_e4_n = MidPosition(pose,e2_start_ss,e4_end_ss)
    mid_e2_e4_c = MidPosition(pose,e2_end_ss,e4_start_ss)

    arch1_dist = (mid_e1_e3_c-mid_e2_e4_n).length()
    arch2_dist = (mid_e1_e3_n-mid_e2_e4_c).length()

    ################################
    # twist two cross beta pairs
    ################################

    p1=pose.residue(e1_end_ss-2).atom("CA").xyz()
    p2=pose.residue(e1_end_ss).atom("CA").xyz()
    p3=pose.residue(e2_start_ss).atom("CA").xyz()
    p4=pose.residue(e2_start_ss+2).atom("CA").xyz()
    arch1_twist = dihedral_xyz(p1,p2,p3,p4)

    p1=pose.residue(e3_end_ss-2).atom("CA").xyz()
    p2=pose.residue(e3_end_ss).atom("CA").xyz()
    p3=pose.residue(e4_start_ss).atom("CA").xyz()
    p4=pose.residue(e4_start_ss+2).atom("CA").xyz()
    arch2_twist = dihedral_xyz(p1,p2,p3,p4)


    p1=pose.residue(seg1.bp_data[-3][0]).atom("CA").xyz()
    p2=pose.residue(seg1.bp_data[-1][0]).atom("CA").xyz()
    p3=pose.residue(seg2.bp_data[0][0]).atom("CA").xyz()
    p4=pose.residue(seg2.bp_data[2][0]).atom("CA").xyz()
    arch1_twist2 = dihedral_xyz(p1,p2,p3,p4)
    arch1_dist2 = p2.distance(p3)
    
    

    p1=pose.residue(seg3.bp_data[-3][0]).atom("CA").xyz()
    p2=pose.residue(seg3.bp_data[-1][0]).atom("CA").xyz()
    p3=pose.residue(seg4.bp_data[0][0]).atom("CA").xyz()
    p4=pose.residue(seg4.bp_data[2][0]).atom("CA").xyz()
    arch2_twist2 = dihedral_xyz(p1,p2,p3,p4)
    arch2_dist2 = p2.distance(p3)


    # lateral beta-arch distance
    # arch1    
    v1 = Vec(pose.residue(seg2.bp_data[0][0]).atom("CA").xyz() - pose.residue(seg1.bp_data[-1][0]).atom("CA").xyz())

    #triad=[seg.bp_data[-3][0],seg.bp_data[-2][0],seg.bp_data[-1][0]]
    #s1 = Vec(StrandDirectionTriad(pose,triad))
    #cross_s1v1 = s1.cross(v1)
    #cross_s1v1.normalize()
    #arch1_lat_dist = v1.dot(cross_s1v1)
    

    s1mid = Vec(MidPositionList(pose,list(range(e1_start_ss,e1_end_ss+1))))
    s3mid = Vec(MidPositionList(pose,list(range(e3_start_ss,e3_end_ss+1))))  
    v31 = s3mid-s1mid
    v31.normalize()
    arch1_lat_dist = v1.dot(v31)
    
    # arch2    
    v1 = Vec(pose.residue(seg4.bp_data[0][0]).atom("CA").xyz() - pose.residue(seg3.bp_data[-1][0]).atom("CA").xyz())
    #s2mid = Vec(MidPositionList(pose,list(range(e2_start_ss,e2_end_ss+1))))
    #s4mid = Vec(MidPositionList(pose,list(range(e4_start_ss,e4_end_ss+1))))  
    #v42 = s4mid-s2mid
    #v42.normalize()
    arch2_lat_dist = v1.dot(v31)
        
    

    ################################
    # twist of beta arches
    ################################

    sandwich_twist = dihedral_xyz(mid_e1_e3_n,mid_e1_e3_c,mid_e2_e4_n,mid_e2_e4_c)

    ################################
    # distance between two cross beta pairs
    ################################

    cm13=numeric.xyzVector_double_t()
    for pos in s1_range+s3_range:
        cm13+=pose.residue(pos).atom("CA").xyz()
    ScaleXYZVector(cm13, 1/len(s1_range+s3_range))

    cm24=numeric.xyzVector_double_t()
    for pos in s2_range+s4_range:
        cm24+=pose.residue(pos).atom("CA").xyz()
    ScaleXYZVector(cm24, 1/len(s2_range+s4_range))

    dist_cen = (cm13-cm24).length()

    ################################
    
    # 
    
    seg1 = blue.segment_dict[cb_segments[0]] # b-arch1
    seg2 = blue.segment_dict[cb_segments[1]] # b-arch1
    seg3 = blue.segment_dict[cb_segments[2]] # b-arch2
    seg4 = blue.segment_dict[cb_segments[3]] # b-arch2    
    

    #return sandwich_twist,arch1_twist,arch2_twist,arch1_twist2,arch2_twist2,arch1_dist2,arch2_dist2,dist_cen
    return arch1_twist2,arch2_twist2,arch1_dist2,arch2_dist2,arch1_lat_dist,arch2_lat_dist



def SimpleCrossBetaRotationalParameters(pose,blue,cb_segments,sspairing):
    seg1 = blue.segment_dict[cb_segments[0]] # b-arch1
    seg2 = blue.segment_dict[cb_segments[1]] # b-arch1
    seg3 = blue.segment_dict[cb_segments[2]] # b-arch2
    seg4 = blue.segment_dict[cb_segments[3]] # b-arch2

    # cross beta pairing 1
    s1_range,s3_range = sspairing["%s-%s" %(cb_segments[0],cb_segments[2])]
    # cross beta pairing 2
    s2_range,s4_range = sspairing["%s-%s" %(cb_segments[1],cb_segments[3])]

    # sspairing info ensure correct ranges for calcularing twist.
    e1_start_ss,e1_end_ss = s1_range
    e3_end_ss,e3_start_ss = s3_range # anti
    e2_start_ss,e2_end_ss = s2_range
    e4_end_ss,e4_start_ss = s4_range # anti

    


    ################################
    # twist two cross beta pairs
    ################################
    p1 = MidPosition(pose,e1_start_ss,e3_end_ss)
    p2 = MidPositionList(pose,[e1_start_ss,e1_end_ss,e3_start_ss,e3_end_ss])
    p3 = MidPositionList(pose,[e2_start_ss,e2_end_ss,e4_start_ss,e4_end_ss])
    p4 = MidPosition(pose,e2_end_ss,e4_start_ss)

    twist = dihedral_xyz(p1,p2,p3,p4)

    ################################


    ################################
    # tilt two cross beta pairs
    ################################

    p1 = MidPosition(pose,e2_end_ss,e4_start_ss)
    p2 = MidPositionList(pose,[e2_start_ss,e4_start_ss,e2_end_ss,e4_end_ss])
    p3 = MidPosition(pose,e2_start_ss,e2_end_ss)
    p4 = MidPositionList(pose,[e1_start_ss,e1_end_ss,e3_start_ss,e3_end_ss])

    st=''
    for i,p in enumerate([p1,p2,p3,p4]):
        st+='pseudoatom p%i, pos=[%.1f, %.1f, %.1f]\n' %(i,p.x,p.y,p.z)

    print('tilt atoms:')
    print(st)
    tilt = dihedral_xyz(p1,p2,p3,p4)

    ################################



    ################################
    # roll two cross beta pairs
    ################################
    p1 = MidPosition(pose,e2_end_ss,e4_start_ss)
    p2 = MidPositionList(pose,[e1_start_ss,e3_end_ss,e2_end_ss,e4_start_ss])
    p3 = MidPositionList(pose,[e1_end_ss,e3_start_ss,e2_start_ss,e4_end_ss])
    p4 = MidPosition(pose,e1_end_ss,e3_start_ss)


    roll = dihedral_xyz(p1,p2,p3,p4)

    ################################


    return twist,tilt,roll


def CrossBetaPositions(pose,blue,cb_segments,sspairing):
    seg1 = blue.segment_dict[cb_segments[0]] # b-arch1
    seg2 = blue.segment_dict[cb_segments[1]] # b-arch1
    seg3 = blue.segment_dict[cb_segments[2]] # b-arch2
    seg4 = blue.segment_dict[cb_segments[3]] # b-arch2

    # cross beta pairing 1
    s1_range,s3_range = sspairing["%s-%s" %(cb_segments[0],cb_segments[2])]
    # cross beta pairing 2
    s2_range,s4_range = sspairing["%s-%s" %(cb_segments[1],cb_segments[3])]

    # sspairing info ensure correct ranges for calcularing twist.
    e1_start_ss,e1_end_ss = s1_range
    e3_end_ss,e3_start_ss = s3_range # anti
    e2_start_ss,e2_end_ss = s2_range
    e4_end_ss,e4_start_ss = s4_range # anti

    return s1_range,s3_range,s2_range,s4_range
    
    
def cross_beta_euler_rotational_parameters(p,res1,res2):
    # start transforming stub1 to the reference axis
    stub1 = stub(res1[1],res1[2],res1[0])
    stub2 = stub( Vec(0,0,0),Vec(1,0,0),Vec(0,1,0) )
    xform = stub2 * ~stub1 # stub2 upstream (reference), stub1 downstream
    axis, ang, cen = xform.rotation_axis_center()
    
    for ir in range(1, p.size() + 1):
        for ia in range(1, p.residue(ir).natoms() + 1):
            aid = pyrosetta.rosetta.core.id.AtomID(ia, ir)
            oldxyz = Vec(p.xyz(aid))
            newxyz = xform * oldxyz
            p.set_xyz(aid, newxyz.to_rosetta())    
            
    # reorient stubs after the first transformation
    res1_new = [] ; res2_new=[]
    for pos1,pos2 in zip(res1,res2):
        res1_new.append(xform * pos1)
        res2_new.append(xform * pos2)              

    # find transformation from stub1 to stub2
    stub1 = stub(res2_new[1],res2_new[2],res2_new[0])
    stub2 = stub(res1_new[1],res1_new[2],res1_new[0])

    xform = stub2 * ~stub1 # stub2 upstream (reference), stub1 downstream
    axis, ang, cen = xform.rotation_axis_center()
    rotmat = rotation_matrix(axis,ang)
        
    for ir in range(1, p.size() + 1):
        for ia in range(1, p.residue(ir).natoms() + 1):
            aid = pyrosetta.rosetta.core.id.AtomID(ia, ir)
            oldxyz = Vec(p.xyz(aid))
            newxyz = xform * oldxyz
            p.set_xyz(aid, newxyz.to_rosetta())            

    distance = xform.t.length()

    # convert rotation matrix from xyzMath format to scipy
    rotmat_spy = np.zeros((3,3))
    rotmat_spy[0,0] = rotmat.colx().x
    rotmat_spy[0,1] = rotmat.coly().x
    rotmat_spy[0,2] = rotmat.colz().x
    rotmat_spy[1,0] = rotmat.colx().y
    rotmat_spy[1,1] = rotmat.coly().y
    rotmat_spy[1,2] = rotmat.colz().y
    rotmat_spy[2,0] = rotmat.colx().z
    rotmat_spy[2,1] = rotmat.coly().z
    rotmat_spy[2,2] = rotmat.colz().z
            
            
    r=R.from_matrix(rotmat_spy)
    return distance,r.as_euler('xyz',degrees=True)


def cross_beta_rotational_parameters_intermediate_frame(p,res1,res2):
    # start transforming stub1 to the reference axis
    stub1 = stub(res1[1],res1[2],res1[0])
    stub2 = stub( Vec(0,0,0),Vec(1,0,0),Vec(0,1,0) )
    xform = stub2 * ~stub1 # stub2 upstream (reference), stub1 downstream
    axis, ang, cen = xform.rotation_axis_center()
    
    p.dump_pdb("test.pdb")
    
    for ir in range(1, p.size() + 1):
        for ia in range(1, p.residue(ir).natoms() + 1):
            aid = pyrosetta.rosetta.core.id.AtomID(ia, ir)
            oldxyz = Vec(p.xyz(aid))
            newxyz = xform * oldxyz
            p.set_xyz(aid, newxyz.to_rosetta())    
            
    p.dump_pdb("test_start.pdb")        
    
    # reorient stubs after the first transformation
    res1_new = [] ; res2_new=[]
    for pos1,pos2 in zip(res1,res2):
        res1_new.append(xform * pos1)
        res2_new.append(xform * pos2)              
        
    res1_new=res1 ; res2_new=res2
    
    # find transformation from stub1 to stub2
    stub1 = stub(res2_new[1],res2_new[2],res2_new[0])
    stub2 = stub(res1_new[1],res1_new[2],res1_new[0])

    xform = stub2 * ~stub1 # stub2 upstream (reference), stub1 downstream
    axis, ang, cen = xform.rotation_axis_center()
    distance = xform.t.length()
    rotmat = rotation_matrix(axis,ang)
    ang_deg_mid = np.rad2deg(ang)/2.0
    ang_mid = np.deg2rad(ang_deg_mid)
    tr_mid = xform.t*0.5
    
    xform_mid = Xform(R=rotation_matrix_degrees(axis,ang_deg_mid),t=tr_mid)
    
          
    for ir in range(1, p.size() + 1):
        for ia in range(1, p.residue(ir).natoms() + 1):
            aid = pyrosetta.rosetta.core.id.AtomID(ia, ir)
            oldxyz = Vec(p.xyz(aid))
            #newxyz = xform_mid * oldxyz
            newxyz = xform * oldxyz
            p.set_xyz(aid, newxyz.to_rosetta())            

    #distance = xform_mid.t.length()


    # reorient stubs after the second transformation
    res1_mid = [] ; res2_mid=[]
    for pos1,pos2 in zip(res1_new,res2_new):
        res1_mid.append(xform_mid * pos1)
        res2_mid.append(xform_mid * pos2)  

#    ref_stub = stub(res2_mid[1],res2_mid[2],res2_mid[0])
#    ref_stub = stub(res1_mid[1],res1_mid[2],res1_mid[0])
#    ref_stub = stub(res1_mid[1],res1_mid[2],res1_mid[0])
    ref_stub = stub(res1[1],res1[2],res1[0])    
    r1 = ref_stub.R
    e1 = Vec(r1.colx().x,r1.coly().x,r1.colz().x)
    e2 = Vec(r1.colx().y,r1.coly().y,r1.colz().y)
    e3 = Vec(r1.colx().z,r1.coly().z,r1.colz().z)
    
    param1 = (ang*e1).dot(axis)
    param2 = (ang*e2).dot(axis)
    param3 = (ang*e3).dot(axis)        
    
    
    
    # convert rotation matrix from xyzMath format to scipy
    
    return distance,np.array([param1,param2,param3])

def strand_pair_central_residues_from_range(a,b,c,d,n):
    length = b-a+1
    center=math.ceil(length/2)
    all_residues=list(range(a,b+1)) + list(range(c,d+1))
    residues=[]

    if length > n:
        if n%2 == 0:        
            for k in range(a+center-math.ceil(n/2),a+center+math.ceil(n/2)):
                residues.append(k)
            for k in range(d-center+math.ceil(n/2),d-center-math.ceil(n/2),-1):
                residues.append(k)                
                

        else:
            for k in range(a+center-math.ceil((n-1)/2),a+center+math.ceil((n+1)/2)):
                residues.append(k)
            for k in range(d-center+math.ceil((n-1)/2),d-center-math.ceil((n+1)/2),-1):
                residues.append(k)                
    else:
        residues=all_residues
        
    return residues    
     
     
     

###########################################################    
   

align_templates=True
data2=[]
for pdbname in glob.glob('*.pdb'): # iterate over the pdbs of the present directory
#  try:
    pose=pyrosetta.pose_from_pdb(pdbname)
    blue=MakeBlueprintFromPose(pose)
    CleanBlueprintByStrandBreaks(pose,blue)
    sspair_lst,strand_dic_pairs,sspair_info = StrandPairing(pose,blue,"A")

    topol = blue.topology()
    strand_segments = [ s for s in topol.split('-') if "E" in s]

    start=0
    hairpins=0 ; arches=0
    pairs_arch=[]
    for s1,s2 in zip(strand_segments[start:-1],strand_segments[start+1:]):
        pair=[s1,s2]
        if pair in sspair_lst:
                hairpins+=1
        else:
            arches+=1
            pairs_arch.append(pair)

    #print pairs_arch
    all_cross_beta_segments=[]
    all_cross_beta_positions=[]
    for a1,arch1 in enumerate(pairs_arch):
        for a2,arch2 in enumerate(pairs_arch):
            if a2>a1:
                cb_pair1=[arch1[0],arch2[0]]
                cb_pair2=[arch1[1],arch2[1]]
                if cb_pair1 in sspair_lst and cb_pair2 in sspair_lst:
                    #print(cb_pair1)
                    #print(cb_pair2)
                    #cb_segments = cb_pair1+cb_pair2
                    #cb_segments.sort()
                    #cb_pairs=[cb_pair1,cb_pair2]
                    
                    cb_segments=arch1+arch2
                    
                    ### version 1 ###
                    #sandwich_twist,arch1_twist,arch2_twist,arch1_twist2,arch2_twist2,arch1_dist,arch2_dist,dist_cen = SimpleCrossBetaParamsOnPairing(pose,blue,cb_segments,sspair_info)
                    #data2.append([sandwich_twist,arch1_twist,arch2_twist,arch1_twist2,arch2_twist2,arch1_dist,arch2_dist,dist_cen])
                    
                    arch1_twist2,arch2_twist2,arch1_dist2,arch2_dist2,arch1_lat_dist,arch2_lat_dist =  SimpleCrossBetaParamsOnPairing(pose,blue,cb_segments,sspair_info)                   
                    
                    
                    
                    
                    #print ("%s %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f" %(pdbname,sandwich_twist,arch1_twist,arch2_twist,arch1_twist2,arch2_twist2,arch1_dist,arch2_dist,dist_cen))
                    
                    
                    ### version 2 ###
                    #twist,tilt,roll = SimpleCrossBetaRotationalParameters(pose,blue,cb_segments,sspair_info)
                    #data2.append([twist,tilt,roll])
                    #print ("%s %.2f %.2f %.2f" %(pdbname,twist,tilt,roll))                    
                    s1_range,s3_range,s2_range,s4_range = CrossBetaPositions(pose,blue,cb_segments,sspair_info)
                    #print(s1_range,s3_range,s2_range,s4_range)
                    
                    
                    # define positions
                    e1_start_ss,e1_end_ss = s1_range
                    e3_end_ss,e3_start_ss = s3_range # anti
                    e2_start_ss,e2_end_ss = s2_range
                    e4_end_ss,e4_start_ss = s4_range # anti
                    #print(s1_range,s3_range,s2_range,s4_range)

                    # need at least 5-residue strands in the cross beta motif for making the alignment with 5-res template used for transforms.
                    if ( e1_end_ss - e1_start_ss + 1 ) < 5 or ( e2_end_ss - e2_start_ss + 1 ) < 5 or ( e3_end_ss - e3_start_ss + 1 ) < 5 or ( e4_end_ss - e4_start_ss + 1 ) < 5: continue

                    # build cb-pose
                    cb_sel = pyrosetta.rosetta.utility.vector1_unsigned_long()
                    #for i in list(range(e1_start_ss,e1_end_ss+1))+list(range(e3_start_ss,e3_end_ss+1))+list(range(e2_start_ss,e2_end_ss+1))+list(range(e4_start_ss,e4_end_ss+1)):                 
                       #cb_sel.append(i)
                                              
                       
                    # alternative selection for later superposition to reference 5-residue strand dimer.   
                    e24_central = strand_pair_central_residues_from_range(e2_start_ss,e2_end_ss,e4_start_ss,e4_end_ss,5)
                    # find 5 central residues
                    e13_central = strand_pair_central_residues_from_range(e1_start_ss,e1_end_ss,e3_start_ss,e3_end_ss,5)
                    e13_central.sort()
                    e24_central.sort()                                                     


                    # build sheet 1 pose
                    cb_pose1 = pyrosetta.rosetta.core.pose.Pose()
                    cb_sel1 = pyrosetta.rosetta.utility.vector1_unsigned_long()
                    for i in e13_central:
                        cb_sel1.append(i)
                        
                    n=len(cb_sel1)        
                    ft = pyrosetta.rosetta.core.kinematics.FoldTree(n)
                    pyrosetta.rosetta.core.pose.create_subpose(pose, cb_sel1, ft, cb_pose1)
                    
                    # build sheet 2 pose
                    cb_pose2 = pyrosetta.rosetta.core.pose.Pose()
                    cb_sel2 = pyrosetta.rosetta.utility.vector1_unsigned_long()
                    for i in e24_central:
                        cb_sel2.append(i)        

                    n=len(cb_sel2)
                    ft = pyrosetta.rosetta.core.kinematics.FoldTree(n)
                    pyrosetta.rosetta.core.pose.create_subpose(pose, cb_sel2, ft, cb_pose2)
                    
                    # read reference strand dimers with two different pleatings
                    refpose1a = pyrosetta.pose_from_pdb('ref_strand_dimer_6_res-5A.pdb')
                    refpose1b = pyrosetta.pose_from_pdb('ref_strand_dimer_6_res-5B.pdb')

                    refpose2a = refpose1a.clone()
                    refpose2b = refpose1b.clone()

                    sel = pyrosetta.rosetta.utility.vector1_unsigned_long()
                    for i in range(1,11):
                        sel.append(i)


                    # now identify the best superposition between reference dimers and target cb residues.
                    # sheet 1
                    superimp1 = pyrosetta.rosetta.protocols.simple_moves.SuperimposeMover(cb_pose1,1,10,1,10,True)
                    superimp1.apply(refpose1a)
                    superimp1.apply(refpose1b)

                    rms1a = protocols.toolbox.pose_manipulation.superimpose_pose_on_subset_CA(refpose1a,cb_pose1,sel)
                    rms1b = protocols.toolbox.pose_manipulation.superimpose_pose_on_subset_CA(refpose1b,cb_pose1,sel)

                    if rms1a < rms1b:
                       merged_pose = refpose1a.clone()

                    else:
                       merged_pose = refpose1b.clone()

                    # sheet 2
                    superimp2 = pyrosetta.rosetta.protocols.simple_moves.SuperimposeMover(cb_pose2,1,10,1,10,True)
                    superimp2.apply(refpose2a)                  
                    superimp2.apply(refpose2b)

                    rms2a = protocols.toolbox.pose_manipulation.superimpose_pose_on_subset_CA(refpose2a, cb_pose2,sel)
                    rms2b = protocols.toolbox.pose_manipulation.superimpose_pose_on_subset_CA(refpose2b, cb_pose2,sel)
                    
                    if rms2a < rms2b:
                       pyrosetta.rosetta.core.pose.append_pose_to_pose(merged_pose, refpose2a, False)      
                    else:
                       pyrosetta.rosetta.core.pose.append_pose_to_pose(merged_pose, refpose2b, False)


                    if not align_templates:
                        merged_pose = cb_pose1.clone()
                        pyrosetta.rosetta.core.pose.append_pose_to_pose(merged_pose, cb_pose2, False)
    
                    # the two best fits are merged in a single pose.
                    
                    # positions for building the stubs
                    e1_start_ss,e1_end_ss = 1,5
                    e3_start_ss,e3_end_ss = 6,10
                    e2_start_ss,e2_end_ss = 11,15
                    e4_start_ss,e4_end_ss = 16,20

                    p1 = Vec(MidPosition(merged_pose,e1_start_ss,e3_end_ss))
                    #p11 = Vec(MidPosition(merged_pose,e1_end_ss,e3_start_ss))
                    #p2 = Vec(MidPositionList(merged_pose,[e1_start_ss,e1_end_ss,e3_start_ss,e3_end_ss]))
                    #p3 = Vec(MidPositionList(merged_pose,[e2_start_ss,e2_end_ss,e4_start_ss,e4_end_ss]))
                    p4 = Vec(MidPosition(merged_pose,e2_end_ss,e4_start_ss))
                    #p41 = Vec(MidPosition(merged_pose,e2_start_ss,e4_end_ss))
                    
                    ### alternate definition
                    p2 = Vec(MidPositionList(merged_pose,list(range(e1_start_ss,e1_end_ss+1))+list(range(e3_start_ss,e3_end_ss+1))))
                    p3 = Vec(MidPositionList(merged_pose,list(range(e2_start_ss,e2_end_ss+1))+list(range(e4_start_ss,e4_end_ss+1))))
                    
                    # sheet 1 vector normal to plane
                    #s1mid = Vec(MidPositionList(merged_pose,[e1_start_ss,e1_end_ss]))
                    #s3mid = Vec(MidPositionList(merged_pose,[e3_start_ss,e3_end_ss]))
                    s1mid = Vec(MidPositionList(merged_pose,list(range(e1_start_ss,e1_end_ss+1))))
                    s3mid = Vec(MidPositionList(merged_pose,list(range(e3_start_ss,e3_end_ss+1))))  
                    s1 = p1-p2
                    #s11 = p2-p11
                    s13_plane = s1.cross(s3mid-s1mid)
                    #s13_plane2 = s11.cross(s3mid-s1mid)
                    #s13_plane = p3-p2
                    p3b = p2 + s13_plane.normalized()
                    #p3b = p2 + (s13_plane+s13_plane2).normalized()
                    
                    
                    # sheet 2 vector normal to plane                    
                    #s2mid = Vec(MidPositionList(merged_pose,[e2_start_ss,e2_end_ss]))
                    #s4mid = Vec(MidPositionList(merged_pose,[e4_start_ss,e4_end_ss]))
                    s2mid = Vec(MidPositionList(merged_pose,list(range(e2_start_ss,e2_end_ss+1))))
                    s4mid = Vec(MidPositionList(merged_pose,list(range(e4_start_ss,e4_end_ss+1))))
                    s2 = p4-p3
                    #s21 = p3-p41
                    s24_plane = s2.cross(s4mid-s2mid)
                    #s24_plane2 = s21.cross(s4mid-s2mid)
                    p3c = p3 - s24_plane.normalized() # negative to be oriented to the equivalent side as the first sheet.                    
                    #p3c = p3 - (s24_plane+s24_plane2).normalized() # negative to be oriented to the equivalent side as the first sheet.                    
                    # alternative
                    #s24_plane = s13_plane
                    #p3c = p3 + s24_plane.normalized()
                    
                    
                    # stubs for the two sheets
                    stub1_res = [p1,p2,p3b]
                    stub2_res = [p4,p3,p3c]
                    
                    # parameters
                    p=merged_pose.clone()
                    dist,euler_angles = cross_beta_euler_rotational_parameters(p,stub1_res,stub2_res)
                    twist,roll,tilt = euler_angles
                    
                    p2=merged_pose.clone()
                    dist2,mid_angles = cross_beta_rotational_parameters_intermediate_frame(p2,stub1_res,stub2_res)
                    #mid_ang1, mid_ang2, mid_ang3 = mid_angles
                    mid_ang1, mid_ang2, mid_ang3 = np.rad2deg(mid_angles)
                    #print(pdbname,'%.2f %.2f %.2f %.2f' %(dist,twist,roll,tilt))
                    
                    cbeta_cum_lat_dist = arch1_lat_dist-arch2_lat_dist
                    
                    print(pdbname,'%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f' %(dist,twist,roll,tilt,arch1_twist2,arch2_twist2,arch1_lat_dist,arch2_lat_dist,cbeta_cum_lat_dist))
#  except:
#     pass
