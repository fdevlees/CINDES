#!/bin/env python
''' this module contains all functions related to the NWChem program '''
import re
import time
import string
import random
import os
multiplicity = {1: 'singlet', 2: 'doublet', 3: 'triplet', 4: 'quartet'}

from job import BaseJob


class NWChemJob(BaseJob):
    extension = ''
    script = 'ID_NWChem'
    cmd = 'nwchem'

    def write(self):
        pass
    # -----

    def termination(self, logpath):
        with open(logpath, 'r') as fid:
            text = fid.readlines()[-4:]
            if re.search('Zhang', ''.join(text)):
                fid.close()
                return 1
            elif re.search('IGNORE', ''.join(text)):
                fid.close()
                return 2
            else:
                return 3
    # -----

    def errortermination(self, debug=False):
        from CINDES.cclib.parser.nwchemparser import NWChem
        mymol = NWChem(self.logpath).parse()
        from CINDES.utils import utils
        t = utils.PeriodicTable()

        if hasattr(mymol, 'atomcoords'):
            coords = map(list, mymol.atomcoords[-1])
            #print coords
            for sym, xyz in zip(mymol.atomnos, coords):
                xyz.insert(0, t.element[sym])
            print "atomcoords and added elements:"
            for item in mymol.atomcoords[-1]:
                print ' '.join(map(str, item))
            if debug == True:
                import submitter
                if hasattr(mymol, 'optdone'):
                    if mymol.optdone == False:
                        print "Optimizations not converged!"
                    elif mymol.optdone == True:
                        print "Optimization is converged!"
                fid = open(self.filepath, 'r')  # change .log in .com extension and read input file
                multcharge = re.compile(' zmatrix')  # a regex for the mult charge line
                newfile = []
                once = 0  # only find that line once
                for line in fid:  # copy file exept for the zmat found in the inputfile
                    if multcharge.match(line) and once == 0:  # when found
                        print "match!"
                        once += 1
                        while True:
                            line = next(fid)  # take al new lines
                            if line == ' end\n':  # end of zmat
                                # now instead of this zmat that is now completely skipped place in newfile
                                # the last coordinates of the crashed run
                                coords = filter(lambda x: not x[0] is None, coords)
                                xyz_f = [" " + item[0] +
                                         ' '.join(map("{:12.6f}".format, item[1:])) + '\n' for item in coords]
                                print xyz_f
                                newfile.extend(xyz_f)
                                break
                    else:
                        newfile.append(line)  # copy that line because it is not the zmat found in the inputfile
                print "=" * 20
                newfilepath = "{}/{}zzz".format(self.path, self.name)
                open(newfilepath, 'w').writelines(newfile)
                print "newfile written in: ", newfilepath
                errorjob = NWChemJob(newfilepath, self.calc)
                errorjob.submit()
                self.errorpath = "{}.log".format(newfilepath)
                return newfilepath
        else:
            print "has no coords in file"
            return False


def write_geom(mol, fid, geom=None):
    # set attributes:
    if geom is None:
        Azmat = 'zmat'
        Axyz = 'xyz'
    else:
        Azmat = 'zmat{}'.format(str(geom))
        Axyz = 'xyz{}'.format(str(geom))

    if hasattr(mol, Azmat):
        zmat = getattr(mol, Azmat)
        fid.write("geometry\n zmatrix\n")
        for i in range(len(zmat)):
            fid.write("  ")
            for item in zmat[i]:
                fid.writelines("%s " % item)
            fid.write("\n")
        fid.write(" end\nend\n")
    elif hasattr(mol, Axyz):
        xyz = getattr(mol, Axyz)
        fid.write("geometry\n")
        fid.write(xyz)
        fid.write("end\n")
        #fid.write('\n')
    else:
        print "attribute not found:", Axyz
        raise AttributeError
    return


def write_subjob(fid, job):
    hotline = job['hotline']
    func, basis = [w for w in hotline.split() if '/' in w][0].split('/')
    if func in ['rhf', 'uhf']:
        theory = 'scf'
    elif func in ['mp2', 'scf', 'ccsd', 'ccsd(t)']:
        theory = func
    else:
        theory = 'dft'
        if func[0] == 'u':  # so functional is defined as ub3lyp or ub3pw91
            odft = True
            func = func[1:]  # i.e. remove the u from the functional
        else:
            odft = False
    fid.write("charge {charge}\n".format(charge=job['charge']))
    # basisset/functional
    fid.write("basis\n * library {}\nend\n".format(basis))
    if theory == 'scf':
        fid.write('scf\n {}\n maxiter 50\n {}\nend\n'.format(
            multiplicity[int(job['mult'])], func))
    if theory == 'dft':
        fid.write('dft\n')
        fid.write(' iterations 100\n')
        if odft:
            fid.write(' odft\n')  # makes it an open-shell unrestricted calculation
        # for quadratic convergence:
        if 'scf=xqc' in hotline:
            fid.write(' cgmin\n')
        if 'fukui' in hotline:
            fid.write(' fukui\n')
        fid.write(' mult {}\n'.format(job['mult']))
        #if func=='b3lyp':xc='vwn_1_rpa 0.19 lyp 0.81 HFexch 0.20 slater 0.80 becke88 nonlocal 0.72'
        if func == 'b3lyp':
            xc = 'b3lyp'
        elif func == 'b3pw91':
            xc = 'acm'
        elif func == 'bp86':
            xc = 'becke88 perdew86'
        elif func == 'b3p86':
            xc = 'vwn_1_rpa 1.00 perdew86 0.81 HFexch 0.20 slater 0.80 becke88 nonlocal 0.72'
        elif func == 'cam-b3lyp':
            xc = 'xcamb88 1.00 lyp 0.81 vwn_5 0.19 hfexch 1.00 cam 0.33 cam_alpha 0.19 cam_beta 0.46'
        elif func == 'blyp':
            xc = 'becke88 lyp'
        else:
            raise NameError('no valid functional')
        fid.write(' xc {}\n'.format(xc))
        fid.write('end\n')
    if 'cosmo' in hotline:
        dielec = filter(lambda x: 'cosmo' in x, hotline.split())[0].split('=')[1]
        fid.write("cosmo\n {}\nend\n".format(dielec))
    if 'smd' in hotline:
        solvent = filter(lambda x: 'cosmo' in x, hotline.split())[0].split('=')[1]
        fid.write("cosmo\n do_cosmo_smd true\n solvent {}\nend\n".format(solvent))
    if 'opt' in hotline:
        fid.write('driver\n maxiter 100\nend\n')
    if 'pop=npa' in hotline:
        fid.write('property\n nbofile\nend\n')

    fid.write("task {}".format(theory))
    # functional mult
    if 'opt' in hotline:
        fid.write(" optimize")
    if 'freq' in hotline:
        fid.write(" frequencies")
    fid.write("\n\n")
    return


#---- START FILEWRITER2 THIS ONLY FOR MAKERS TRY TO MAKE THIS ONE UNIVERSAL ----#
def filewriter(mol, calc, pos=None):
    '''    This function creates a file with the geometry contained in zmat
    The name of the file contains the index in the name
    '''
    # ------------
    # this function uses globals: identify, path
    # ------------
    index = mol.index
    paras = calc
    jobs = paras['jobs']

    if pos:
        name = "{0}{1}_{2}".format(paras['identify'], str(index), str(pos))
        filename = name
        filepath = '{0}/{1}/{2}'.format(paras['path'], str(index), filename)
        geom = pos
        Job = NWChemJob(filepath, calc)
        Job.pos = pos
    else:
        name = paras['identify'] + str(index)
        filename = name
        filepath = paras['path'] + '/' + filename
        try:
            geom = calc['geom']
        except KeyError:
            geom = None
        Job = NWChemJob(filepath, calc)
    mol.addjob(Job)

    fid = open(filepath, 'w')

    # JOB 1
    # extract info
    job1 = jobs[0]

    # write info
    fid.write("echo\nstart {filename}\n".format(filename=filename))
    #fid.write("memory 1500 mb\n")
    fid.write("memory total 800 stack 200 heap 200 global 400 mb\n")

    # set scratchdir:
    #hash = ''.join(random.choice(string.ascii_lowercase) for _ in range(4))
    #spath = "/scratch/leuven/100/vsc10010/REDOX/CALC/{}".format(filename.rsplit('_',1)[1])
    #try: 
    #    os.makedirs(spath)
    #except OSError:
    #    if not os.path.isdir(spath):
    #        raise

    #fid.write('scratch_dir ' + spath + '\n')
    fid.write("title \"{filename}\"\n".format(filename=filename))
    # here the zmat
    write_geom(mol, fid, geom)
    write_subjob(fid, job1)

    # THE OTHER JOBS
    for i, job in enumerate(paras['jobs'][1:]):
        # write title
        fid.write("title \"{filename}\"\n".format(filename=filename))
        write_subjob(fid, job)

    fid.close()
    return


def get_paths(mols):
    ''' get all paths that need to be examined later '''
    files = []
    for molecule in mols:
        files.extend(get_molpaths(molecule))
    return files


def get_molpaths(mol, fileparameters):
    paths = []
    for job in mol.jobs:
        log = get_logpath(job)
        paths.append(log)
    return paths
