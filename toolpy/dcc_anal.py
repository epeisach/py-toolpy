#!/usr/bin/env /apps/python-2.6.1/bin/python
#=================================================================
#   This module is to analyse the data items calculated from dcc
#  (an independent one)
#=================================================================
import os,sys, math
import cifparse as cif, util


def process(args):
    
   # print(args)
    flist,file, vlist = '','',[]
    for i,x in enumerate(args):
        if '-file' in x :
            file=args[i+1]
        if '-matt' in x :
           get_matt_coeff(args[i+1],1)
           sys.exit()
        if '-list' in x :
            flist=args[i+1]
        elif '-anal' in x :
            print('doing analysis')
        
    if file:
        vlist= get_list(file)
        return
    elif flist:
        all_value = get_all_item(flist)
        print ('The output file is %s' %all_value)

##########################################################
def get_matt_coeff(file, id):
    '''calculate Matthew_coeff: file is in pdb format;
       id=0, do not consider MTRIX records;  id==1 consider
    '''

    fp=open(file, 'r')
    
    nop, nmat=0,1
    cell=[1,1,1,1,1,1]
    sym='X'
    res,atom=[],[]
    amass, hmass = 0,0

    for x in fp:
        if 'REMARK 290 ' in x[:12] and '555' in x and ',' in x[23:32]:
            nop=nop+1
        elif 'SEQRES' in x[:6] :
            t=x[17:79].split()
            res.extend(t)
        elif 'MTRIX3' in x[:6]  and '1' not in x[55:].strip() :
            nmat  = nmat +1
        elif 'CRYST1' in x[:6]:
            c=x[7:54].split()
            cell=[float(y) for y in c]
            sym=x[54:65].strip()
#            break

        elif ('ATOM' in x[:4] or 'HETA' in x[:4]) and 'HOH' not in x[15:20] :
            atom.append(x)
            amass = amass + atom_mass(x[76:78].strip())*float(x[54:60])
            #amass = amass + atom_mass(x[76:78].strip())

     #   elif 'ENDMDL' in x[:6]:
     #       break
    fp.close()

    cell_vol=cell_volume(cell)
    nsym=sg_nsym(sym)
    
    if nsym ==-1:
        print('Error: space group (%s) is not in the list (%s)' %(sym, file))
        nsym=nop
#----------

    rmass=0
    for x in res:
        resm=residue_mass(x)
        if resm<1:
            m1=non_standard_res(x,atom)
            rmass = rmass + m1
        else:
           # print(x, resm)
            rmass = rmass + resm 

    res_mass=rmass 
    atm_mass=amass 
   
#    res_mass, atm_mass=rmass, amass
#    print(cell_vol, atm_mass,res_mass, nsym, nmat)
    amatt, asolv = calc_matt(cell_vol, atm_mass, nsym, nmat, id)
    rmatt, rsolv = calc_matt(cell_vol, res_mass, nsym, nmat, id)
    print('%s :nsym, nmat, matt, solv = %2d %2d %5.2f  %6.2f | \
%5.2f  %6.2f ' %(file,nsym, nmat, amatt,asolv, rmatt,rsolv))        

##########################################################
def non_standard_res(res,atom):

    amass, n =0, 0
    for i, x in enumerate(atom):
        if res == x[17:20].strip():
            if i>0 and x[22:26] != atom[i-1][22:26] and n>0 : break
            n=n+1
            amass = amass + atom_mass(x[76:78].strip())
           # print("%2d" %i,res,"%6.1f" %amass, x[22:26], atom[i-1][22:26], x)
            
     
    return amass
##########################################################
def calc_matt(cell_vol, mass, nsym, nmat, id) :

    mass = nsym*mass
    if nmat>0 and id>0 : mass=mass*nmat

    matt, solv = 0, 1
    if(cell_vol >10 and mass>10):
        matt = cell_vol/mass;
        solv =(1.0-1.239/matt)*100.0;

    return matt, solv

##########################################################


def get_all_item(flist):

    head='''data_dcc_summary
#
loop_
_pdbx_density.pdbid
_pdbx_density.ls_d_res_high
_pdbx_density.diff_R_work_report_calc
_pdbx_density.R_value_R_work_report
_pdbx_density.R_value_R_free_report
_pdbx_density.R_value_R_work
_pdbx_density.R_value_R_free
_pdbx_density.correlation_coeff_Fo_to_Fc
_pdbx_density.real_space_R
_pdbx_density.correlation
_pdbx_density.fom
_pdbx_density.Biso_mean
_pdbx_density.B_wilson
_pdbx_density.Matthew_coeff
_pdbx_density.solvent_content
_pdbx_density.K_solvent
_pdbx_density.B_solvent
_pdbx_density.tls_group_number
_pdbx_density.ncs_group_number
_pdbx_density.mtrix_number
_pdbx_density.partial_B_value_correction_success
_pdbx_density.iso_B_value_type
_pdbx_density.translational_pseudo_symmetry
_pdbx_density.reflns_twin
_pdbx_density.twin_by_xtriage
_pdbx_density.twin_by_refmac
_pdbx_density.Padilla-Yeates_L2_mean
_pdbx_density.Z_score_L_test
_pdbx_density.anisotropy
_pdbx_density.I_over_sigI_resh
_pdbx_density.I_over_sigI_diff
_pdbx_density.program
_pdbx_density.space_group_name_H-M
'''



    fout=flist + '_all.cif'
    fp=open(flist, 'r')
    fw=open(fout, 'w')
    fw.write(head)
    for x in fp:
        file='sum-all/pdb' + x.lower().strip() + '.ent_rcc_sum.cif'
        #file='test/pdb' + x.lower().strip() + '.ent_rcc_sum.cif'
        vlist= get_list(file)
        for y in vlist:
            fw.write('%s ' %y)
        fw.write('\n')

    fp.close()
    fw.close()
        
    return fout    
##########################################################
def get_list(file):

    vlist=[]

    if not util.check_file(200, file):
        print('Error: file (%s) do not exist' %file)
        return vlist

    flist=open(file,'r').readlines()
    pdbid='XXXX'
    for x in flist:
        if 'data_' in x:
            pdbid=x.split('_')[1].strip()
            break
            
    items,values = cif.cifparse(flist, '_pdbx_density.')

    sym=cif.parse_values(items,values,'_pdbx_density.space_group_name_H-M')
    res=cif.parse_values(items,values,'_pdbx_density.ls_d_res_high');
    rw=cif.parse_values(items,values,'_pdbx_density.R_value_R_work');
    rf=cif.parse_values(items,values,'_pdbx_density.R_value_R_free');
    biso=cif.parse_values(items,values,'_pdbx_density.Biso_mean');
    bwil=cif.parse_values(items,values,'_pdbx_density.B_wilson');
    l2=cif.parse_values(items,values,'_pdbx_density.Padilla-Yeates_L2_mean');
    z=cif.parse_values(items,values,'_pdbx_density.Z_score_L_test');
    fom=cif.parse_values(items,values,'_pdbx_density.fom');
    isig=cif.parse_values(items,values,'_pdbx_density.I_over_sigI_resh');
    isigd=cif.parse_values(items,values,'_pdbx_density.I_over_sigI_diff');
    pst=cif.parse_values(items,values,'_pdbx_density.translational_pseudo_symmetry');
    bsol=cif.parse_values(items,values,'_pdbx_density.B_solvent');
    ksol=cif.parse_values(items,values,'_pdbx_density.K_solvent');
    tlst=cif.parse_values(items,values,'_pdbx_density.partial_B_value_correction_success');
    ntls=cif.parse_values(items,values,'_pdbx_density.tls_group_number');
    nncs=cif.parse_values(items,values,'_pdbx_density.ncs_group_number');
    nmtx=cif.parse_values(items,values,'_pdbx_density.mtrix_number');
    matt=cif.parse_values(items,values,'_pdbx_density.Matthew_coeff');
    solv=cif.parse_values(items,values,'_pdbx_density.solvent_content');
    dpix=cif.parse_values(items,values,'_pdbx_density.Cruickshank_dpi_xyz');
    rtwin=cif.parse_values(items,values,'_pdbx_density.reflns_twin');
    xtwin=cif.parse_values(items,values,'_pdbx_density.twin_by_xtriage');
    tmp=cif.parse_values(items,values,'_pdbx_density.iso_B_value_type');
    if tmp : btype=tmp[0][0]
    ctwin_t=cif.parse_values(items,values,'_pdbx_density.twin_operator');
    ctwin='N'
    if '2:' in ctwin_t[0] : ctwin='Y'
    anis=cif.parse_values(items,values,'_pdbx_density.anisotropy');


# looped 
    items,values = cif.cifparse(flist, '_pdbx_density_corr.')
    prog=cif.parse_values(items,values,'_pdbx_density_corr.program');
    resh=cif.parse_values(items,values,'_pdbx_density_corr.ls_d_res_high');
    rwork=cif.parse_values(items,values,'_pdbx_density_corr.ls_R_factor_R_work');
    rfree=cif.parse_values(items,values,'_pdbx_density_corr.ls_R_factor_R_free');
    fcc=cif.parse_values(items,values,'_pdbx_density_corr.correlation_coeff_Fo_to_Fc');
    rsr=cif.parse_values(items,values,'_pdbx_density_corr.real_space_R');
    dcc=cif.parse_values(items,values,'_pdbx_density_corr.correlation');
    detail=cif.parse_values(items,values,'_pdbx_density_corr.details');     
    

    nr, nc=0, 0
    for i, x in enumerate(detail):
        if 'Best' in x :
            nc=i
            break
        
    rprog, cprog = prog[nr].replace(' ', ''), prog[nc]
    crw, crf, fcc, rsr, dcc=rwork[nc], rfree[nc], fcc[nc],rsr[nc],dcc[nc]
    rw_crw='?'
    if util.is_number(rw[0]) and util.is_number(crw):
        t=int (1000*(float(rw[0]) -float(crw)))
        rw_crw='%d' %(t)
    all=[pdbid, res,rw_crw, rw, rf, crw, crf, fcc, rsr, dcc, fom, biso,bwil,
         matt,solv, ksol, bsol,  ntls,nncs, nmtx,
         tlst,btype, pst, rtwin,xtwin, ctwin,  l2, z,anis, isig,isigd, rprog, sym]

    
    all_new=[]
    for x in all:
        t=x
        if not x :
           t='?'
        else:
            if (type(x)==list ):t=x[0]

        y=t.replace(' ', '_')
        if util.is_number(y) and '.' in y:
            y='%.2f' %float(y)
        
        all_new.append(y)
    
    return all_new


##########################################################
def residue_mass(resname):
    '''get mass of residue; 
       AAD: the residue mass of modified aa; DD is modified NA.
    '''
    resd={
"GLY" : 57.05 ,
"ALA" : 71.08 ,
"VAL" : 99.13 ,
"LEU" : 113.16 ,
"ILE" : 113.16 ,
"SER" : 87.08 ,
"THR" : 101.10 ,
"CYS" : 103.14 ,
"PRO" : 97.12 ,
"PHE" : 147.18 ,
"TYR" : 163.18 ,
"TRP" : 186.21 ,
"HIS" : 138.15 ,
"ASP" : 110.05 ,
"ASN" : 114.10 ,
"GLU" : 122.06 ,
"GLN" : 128.13 ,
"MET" : 131.19 ,
"LYS" : 129.18 ,
"ARG" : 157.19 ,
"DA" : 328.20 ,
"DT" : 303.19 ,
"DG" : 344.20 ,
"DC" : 304.18 ,
"A" : 328.20 ,
"T" : 303.19 ,
"G" : 344.20 ,
"C" : 304.18 ,
"U" : 305.16  
} 
    mass=0
    if resname in resd:
        mass=resd[resname]
#    elif len(resname)<3:
#        mass=310

    return mass    
##########################################################
def sg_nsym(sg):
    '''contains space group names and the number of operators.
    '''
    
    sym={
      "A 1" :          2,
      "A 1 2 1" :     4,
      "A 2" :         4,
      "B 1 1 2" :     4,
      "B 2" :         4,
      "B 2 21 2" :    8,
      "C 2" :         4,
      "C 1 2 1" :     4,
      "C 21" :        4,
      "C 1 21 1" :    4,
      "C 2(A 112)" :  4,
      "C 2 2 2" :     8,
      "C 2 2 21" :    8,
      "C 4 21 2" :   16,
      "F 2 2 2" :    16,
      "F 2 3" :      48,
      "F 4 2 2" :    32,
      "F 4 3 2" :    96,
      "F 41 3 2" :   96,
      "I 1 2 1" :     4,
      "I 1 21 1" :    4,
      "I 2" :         4,
      "I 2 2 2" :     8,
      "I 2 3" :      24,
      "I 21" :        4,
      "I 21 3" :     24,
      "I 4" :         8,
      "I 21 21 21" :  8,
      "I 4 2 2" :    16,
      "I 4 3 2" :    48,
      "I 41" :        8,
      "I 41 2 2" :   16,
      "I 41 3 2" :   48,
      "P 1" :         1,
      "P -1" :        2,
      "P 2" :         2,
      "P 1 2 1" :     2,
      "P 1 1 2" :     2,
      "P 2 2 2" :     4,
      "P 2 3" :      12,
      "P 2 2 21" :    4,
      "P 2 21 21" :   4,
      "P 2 21 2" :    4,
      "P 21" :        2,
      "P 1 21 1" :    2,
      "P 1 1 21" :    2,
      "P 21(C)" :     2,
      "P 21 2 21" :   4,
      "P 21 3" :     12,
      "P 21 21 2" :   4,
      "P 21 21 2 A" : 4,
      "P 21 21 21" :  4,
      "P 3" :         3,
      "P 3 1 2" :     6,
      "P 3 2 1" :     6,
      "P 31" :        3,
      "P 31 1 2" :    6,
      "P 31 2 1" :    6,
      "P 32" :        3,
      "P 32 1 2" :    6,
      "P 32 2 1" :    6,
      "P 4" :         4,
      "P 4 2 2" :     8,
      "P 4 3 2" :    24,
      "P 4 21 2" :    8,
      "P 41" :        4,
      "P 41 2 2" :    8,
      "P 41 3 2" :   24,
      "P 41 21 2" :   8,
      "P 42" :        4,
      "P 42 2 2" :    8,
      "P 42 3 2" :   24,
      "P 42 21 2" :   8,
      "P 43" :        4,
      "P 43 2 2" :    8,
      "P 43 3 2" :   24,
      "P 43 21 2" :   8,
      "P 6" :         6,
      "P 6 2 2" :    12,
      "P 61" :        6,
      "P 61 2 2" :   12,
      "P 62" :        6,
      "P 62 2 2" :   12,
      "P 63" :        6,
      "P 63 2 2" :   12,
      "P 64" :        6,
      "P 64 2 2" :   12,
      "P 65" :        6,
      "P 65 2 2" :   12,
      "H 3" :         9,
      "R 3" :         3,
      "H 3 2" :      18,
      "R 3 2" :       6 ,
      "I 41/A" :      16 ,
      "P 1 21/C 1" :   4 ,
      "I -4 C 2" :  16 
}

    if sg in sym:
        nsym=sym[sg]
    else:
        nsym=-1
    return nsym

##########################################################
def cell_volume(cell):
       
    alpha = 3.14159 * cell[3]/180;
    beta  = 3.14159 * cell[4]/180;
    gamma = 3.14159 * cell[5]/180;

    cell_volume = cell[0] * cell[1] * cell[2] * \
math.sqrt(1.- math.cos(alpha)*math.cos(alpha) -  \
math.cos(beta)*math.cos(beta) - math.cos(gamma)*math.cos(gamma) \
     + 2.0 * math.cos(alpha) * math.cos(beta) * math.cos(gamma));
    return cell_volume;
##########################################################
def atom_mass(atom):
    '''selected mass of atom (only for estimating matthew coeff.)
    '''
    at={
"H" : 1.0079 ,
"D" : 1.0079 ,
"C" : 12.011 ,
"N" : 14.0067 ,
"O" : 15.9994 ,
"F" : 18.9984 ,
"MG" : 24.305 ,
"P" : 30.97376 ,
"S" : 32.06 , 
"ZN" : 65.38 ,
"SE" : 78.96 ,
"BR" : 79.904 ,
"X" : 40
}
    mass=40
    if atom in at: mass=at[atom]

    return mass  

##########################################################

if __name__ == '__main__':
    process(sys.argv)
    





    
