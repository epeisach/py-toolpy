# =================================================================
#  Calculate the Matthew coefficient and solvent content.
#
# =================================================================
import math
import util
from cifparse import cif2pdb


##################################################################
def matt_coeff(infile, outfile):
    """calculate Matthew_coeff and solven content: file is in pdb format;"""

    if not util.check_file(100, infile):
        print("Error: file (%s) not exist" % infile)
        return

    file = infile
    if util.is_cif(infile):
        file = cif2pdb(infile)

    fp = open(file, "r")

    cell = [1, 1, 1, 1, 1, 1]
    spt, nop, nmat, sym, res, atom, res = 1, 0, 1, "X", [], [], []
    rmass, amass, armass = 0, 0, 0
    hetres, aname, rest = [], [], ""
    for x in fp:
        if "REMARK 290 " in x[:12] and "555" in x and "," in x[23:32]:
            nop = nop + 1
        elif "SEQRES" in x[:6]:
            t = x[17:79].split()
            res.extend(t)
        #        elif 'SPLIT' in x[:6] :
        #            t=x[6:].split()
        #            spt=len(t)

        elif "MTRIX3" in x[:6] and "1" not in x[55:].strip():
            nmat = nmat + 1
        elif "CRYST1" in x[:6]:
            c = x[7:54].split()
            cell = [float(y) for y in c]
            sym = x[54:65].strip()

        elif "ATOM" in x[:4] or "HETA" in x[:4]:
            if "HOH" in x[17:20] or "DOD" in x[17:20]:
                continue

            atom.append(x)
            occ = float(x[54:60])
            amass = amass + atom_mass(x[76:78].strip()) * occ
            t = x[17:27]  # comp_ch_res_int
            aname.append(x[76:78].strip())

            if t != rest:
                comp = t[:3].strip()
                restmp = residue_mass(comp)
                if restmp < 1:
                    hetres.append(comp)
                else:
                    armass = armass + residue_mass(comp)

                rest = t
                aname = []

        elif "ENDMDL" in x[:6]:
            break

    fp.close()

    cell_vol = cell_volume(cell)
    nsym = sg_nsym(sym)
    if nsym == -1:
        print("Error: space group (%s) is not in the list (%s)" % (sym, file))
        nsym = nop
    # ----------

    for x in hetres:
        armass = armass + non_standard_res(x, atom)

    for x in res:
        resm = residue_mass(x)
        if resm < 1:
            m1 = non_standard_res(x, atom)
            rmass = rmass + m1
        else:
            rmass = rmass + resm

    amatt, asolv = calc_matt(cell_vol, amass, nsym, nmat, spt)  # by atom, occ
    rmatt, rsolv = calc_matt(cell_vol, rmass, nsym, nmat, spt)  # by SEQRES
    armatt, arsolv = calc_matt(cell_vol, armass, nsym, nmat, spt)  # residue

    matt, solv = -1, -1

    if 2.0 < rmatt < 5:
        matt, solv = rmatt, rsolv
    elif 2.0 < armatt < 5:
        matt, solv = armatt, arsolv
    elif 2.0 < amatt < 5:
        matt, solv = amatt, asolv
    else:
        matt, solv = armatt, arsolv
        print("Warning: packing problem (%s), Matthew_coeff=%.2f; Solvent=%.2f" % (file, matt, solv))

    if util.is_cif(infile):
        util.delete_file(file)

    print("%s : split nsym, nmat= %2d %2d %2d" % (file, spt, nsym, nmat))
    print("By ATOM:    matt= %6.2f , solvent= %6.2f " % (amatt, asolv))
    print("By SEQRES:  matt= %6.2f , solvent= %6.2f " % (rmatt, rsolv))
    print("By residue: matt= %6.2f , solvent= %6.2f " % (armatt, arsolv))
    print("Possible:   matt= %6.2f , solvent= %6.2f " % (matt, solv))

    print("\nmass_total_atom=%.1f ;  cell_vol=%.1f" % (amass, cell_vol))
    error = "?"
    if matt > 8.7 or matt < 1.5:
        error = "Warning: Matthew_coefficient(%.2f) is abnormal. Possible incomplete content of ASU (or a split entry)." % matt
    if matt == 0.0 and solv == 1.0:
        error = "?"  # space group problem

    if outfile:
        fw = open(outfile, "w")
        ss = """data_matt
#
_packing.Matthew_coefficient  %6.2f
_packing.solvent_content     %6.2f
_packing.error  "%s"\n
        """ % (
            matt,
            solv,
            error,
        )
        fw.write(ss)
        fw.close()
        print("The output file = %s\n" % outfile)
    return matt, solv


##########################################################
def non_standard_res(res, atom):
    """mass of none standard residues are calculated by each atoms
    in the atom list.
    """

    amass, n = 0, 0
    for i, x in enumerate(atom):
        if res == x[17:20].strip():
            if i > 0 and x[22:26] != atom[i - 1][22:26] and n > 0:
                break
            n = n + 1
            amass = amass + atom_mass(x[76:78].strip())

    return amass


##########################################################
def calc_matt(cell_vol, mass, nsym, nmat, spt):
    """cell_vol: cell volume; mass: mass in ASU;
    nsym: symmetry operators;  nmat: number of matrix; spt: number of split
    """

    if mass < 1 or cell_vol < 2:
        return 0, 0

    mass = spt * nsym * mass
    if nmat > 0:
        mass = mass * nmat

    matt, solv = 0, 1
    if cell_vol > 10 and mass > 10:
        matt = cell_vol / mass
        solv = (1.0 - 1.239 / matt) * 100.0

    return matt, solv


##########################################################
def residue_mass(resname):
    """get mass of residue;
    AAD: the residue mass of modified aa; DD is modified NA.
    """
    resd = {
        "GLY": 57.05,
        "ALA": 71.08,
        "VAL": 99.13,
        "LEU": 113.16,
        "ILE": 113.16,
        "SER": 87.08,
        "THR": 101.10,
        "CYS": 103.14,
        "PRO": 97.12,
        "PHE": 147.18,
        "TYR": 163.18,
        "TRP": 186.21,
        "HIS": 138.15,
        "ASP": 110.05,
        "ASN": 114.10,
        "GLU": 122.06,
        "GLN": 128.13,
        "MET": 131.19,
        "LYS": 129.18,
        "ARG": 157.19,
        "DA": 328.20,
        "DT": 303.19,
        "DG": 344.20,
        "DC": 304.18,
        "A": 328.20,
        "T": 303.19,
        "G": 344.20,
        "C": 304.18,
        "U": 305.16,
        "MG": 24.305,
        "ZN": 65.38,
    }
    mass = 0
    if resname in resd:
        mass = resd[resname]
    #    elif len(resname)<3:
    #        mass=310

    return mass


##########################################################
def sg_nsym(sg):
    """contains space group names and the number of operators."""

    sym = {
        "A 1": 2,
        "A 1 2 1": 4,
        "A 2": 4,
        "B 1 1 2": 4,
        "B 2": 4,
        "B 2 21 2": 8,
        "C 2": 4,
        "C 1 2 1": 4,
        "C 21": 4,
        "C 1 21 1": 4,
        "C 2(A 112)": 4,
        "C 2 2 2": 8,
        "C 2 2 21": 8,
        "C 4 21 2": 16,
        "F 2 2 2": 16,
        "F 2 3": 48,
        "F 4 2 2": 32,
        "F 4 3 2": 96,
        "F 41 3 2": 96,
        "I 1 2 1": 4,
        "I 1 21 1": 4,
        "I 2": 4,
        "I 2 2 2": 8,
        "I 2 3": 24,
        "I 21": 4,
        "I 21 3": 24,
        "I 4": 8,
        "I 21 21 21": 8,
        "I 4 2 2": 16,
        "I 4 3 2": 48,
        "I 41": 8,
        "I 41 2 2": 16,
        "I 41 3 2": 48,
        "P 1": 1,
        "P -1": 2,
        "P 2": 2,
        "P 1 2 1": 2,
        "P 1 1 2": 2,
        "P 2 2 2": 4,
        "P 2 3": 12,
        "P 2 2 21": 4,
        "P 2 21 2": 4,
        "P 2 21 21": 4,
        "P 21": 2,
        "P 1 21 1": 2,
        "P 1 1 21": 2,
        "P 21(C)": 2,
        "P 21 2 21": 4,
        "P 21 3": 12,
        "P 21 21 2": 4,
        "P 21 21 2 A": 4,
        "P 21 21 21": 4,
        "P 3": 3,
        "P 3 1 2": 6,
        "P 3 2 1": 6,
        "P 31": 3,
        "P 31 1 2": 6,
        "P 31 2 1": 6,
        "P 32": 3,
        "P 32 1 2": 6,
        "P 32 2 1": 6,
        "P 4": 4,
        "P 4 2 2": 8,
        "P 4 3 2": 24,
        "P 4 21 2": 8,
        "P 41": 4,
        "P 41 2 2": 8,
        "P 41 3 2": 24,
        "P 41 21 2": 8,
        "P 42": 4,
        "P 42 2 2": 8,
        "P 42 3 2": 24,
        "P 42 21 2": 8,
        "P 43": 4,
        "P 43 2 2": 8,
        "P 43 3 2": 24,
        "P 43 21 2": 8,
        "P 6": 6,
        "P 6 2 2": 12,
        "P 61": 6,
        "P 61 2 2": 12,
        "P 62": 6,
        "P 62 2 2": 12,
        "P 63": 6,
        "P 63 2 2": 12,
        "P 64": 6,
        "P 64 2 2": 12,
        "P 65": 6,
        "P 65 2 2": 12,
        "H 3": 9,
        "R 3": 3,
        "H 3 2": 18,
        "R 3 2": 6,
        "I 41/A": 16,
        "P 1 21/C 1": 4,
        "P 1 21/N 1": 4,
        "I -4 C 2": 16,
    }

    if sg in sym:
        nsym = sym[sg]
    else:
        nsym = -1
    return nsym


##########################################################
def cell_volume(cell):
    """volume of cell."""

    alpha = 3.14159 * cell[3] / 180
    beta = 3.14159 * cell[4] / 180
    gamma = 3.14159 * cell[5] / 180

    cell_volume = (
        cell[0]
        * cell[1]
        * cell[2]
        * math.sqrt(
            1.0 - math.cos(alpha) * math.cos(alpha) - math.cos(beta) * math.cos(beta) - math.cos(gamma) * math.cos(gamma) + 2.0 * math.cos(alpha) * math.cos(beta) * math.cos(gamma)
        )
    )

    return cell_volume


##########################################################
def atom_mass(atom):
    """selected mass of atom (only for estimating matthew coeff.)"""
    at = {"H": 1.0079, "D": 1.0079, "C": 12.011, "N": 14.0067, "O": 15.9994, "F": 18.9984, "MG": 24.305, "P": 30.97376, "S": 32.06, "ZN": 65.38, "SE": 78.96, "BR": 79.904, "X": 40}
    mass = 40
    if atom in at:
        mass = at[atom]

    return mass
