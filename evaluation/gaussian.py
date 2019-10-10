#!/bin/env python
''' this module contains all functions related to the Gaussian09 program '''
import re
from job import BaseJob


class GaussianJob(BaseJob):
    extension = '.log'
    script = 'ID_gauss'
    cmd = 'g09'
    assign_geom = True

    def __init__(self, filepath, calc=dict(), worker=False):
        super(GaussianJob, self).__init__(filepath, calc=calc, worker=worker)

        # for gaussian geometries are automatically asigned to attributes unless this variable is set to False

        self.nNormalTermination = 0
        for job in calc['jobs']:
            if 'freq' in job['hotline']:
                self.nNormalTermination += 1
            self.nNormalTermination += 1
        return

    def write(self):
        ''' overwrites the standard BaseJob write method '''
        pass

    def getlog(self):
        return self.name + '.log'

    def parse(self, path=None):
        if path is None:
            path = self.logpath
        from CINDES.evaluation.cclib.parser.gaussianparser import Gaussian
        try:
            mymol = Gaussian(self.logpath).parse()
        except Exception as e:
            print "cclib read error with:", path
            errorfile = 'error_{}.log'.format(self.name)
            print "see {} for more details".format(errorfile)
            with open(errorfile, 'w') as f:
                f.write(e)
            return False
        return mymol

    # ----- might as well be a static method
    def termination(self, logpath, raise_errors=True):
        try:
            with open(logpath, 'r') as fid:
                text = ''.join(fid.readlines()[-3:])
                if (re.search('Normal termination', text)
                  and not (re.search('Initial command', text)
                  and not re.search('Link1', text))):
                    ret = 1
                elif re.search('IGNORE', text):
                    ret = 2
                elif re.search('open-new-file', text):
                    ret = 3
                else:
                    ret = 0
        except IOError:
            ret = 0
        return ret
    # -----

    def errortermination(self, debug=False):
        import time
        mymol = self.parse(self.logpath)
        if not mymol:
            print "did not manage to read logfile without Normal Termination"
            return False
        from CINDES.utils import utils
        t = utils.PeriodicTable()
        if hasattr(mymol, 'atomcoords'):
            coords = map(list, mymol.atomcoords[-1])
            #print coords
            for sym, xyz in zip(mymol.atomnos, coords):
                xyz.insert(0, t.element[sym])
            print "atomcoords and added elements:"
            # for item in mymol.atomcoords[-1]:
            #    print ' '.join(map(str,item))
            if debug == True:
                if hasattr(mymol, 'optdone'):
                    if mymol.optdone == False:
                        print "Optimizations not converged!"
                    elif mymol.optdone == True:
                        print "Optimization is converged!"
                fid = open(self.filepath, 'r')  # change .log in .com extension and read input file
                multcharge = re.compile('^\-?[01]\s[12]')  # a regex for the mult charge line
                newfile = []
                once = 0  # only find that line once
                for line in fid:  # copy file exept for the zmat found in the inputfile
                    if multcharge.match(line) and once == 0:  # when found
                        once += 1
                        newfile.append(line)  # the line with the match itself has to be included in the newfile
                        while True:
                            line = next(fid)  # take al new lines
                            if line == '\n':  # end of zmat
                                # now instead of this zmat that is now completely skipped place in newfile
                                # the last coordinates of the crashed run
                                xyz_f = [item[0] + ' '.join(map("{:12.6f}".format, item[1:])) + '\n' for item in coords]
                                newfile.extend(xyz_f)
                                newfile.extend(['\n'])
                                break
                    elif self.worker and 'nprocshared' in line:
                        # change %nprocshared=* to %nprocsared=28
                        newfile.append('%nprocshared=28\n')
                        #raise NotImplementedError('and what if there is not a nprocshared in file??!')
                    else:
                        newfile.append(line)  # copy that line because it is not the zmat found in the inputfile
                newfilepath = "{}/{}zzz.com".format(self.path, self.name)
                self.errorfile = newfilepath
                open(newfilepath, 'w').writelines(newfile)
                print "newfile written in: ", newfilepath
                errorjob = GaussianJob(newfilepath, self.calc)

                errorjob.submit()
                self.errorpath = "{}/{}zzz.log".format(self.path, self.name)
                return True
        else:
            return False

    def rewrite(self, nprocs=None):
        with open(self.filepath, 'r') as f:
             oldfile = f.readlines()
        newfile = []
        for line in oldfile:
            if not nprocs is None and 'nprocshared' in line:
                newfile.append('%nprocshared=28\n')
            else:
                newfile.append(line)
        with open(self.filepath, 'w') as f:
            f.writelines(newfile)
        print "jobfile rewritten"
        return


def writegeom(mol, fid, geom=None):
    print "geom", geom,
    # set attributes
    if geom is None:
        Azmat = 'zmat'
        Axyz = 'xyz'
    elif geom.startswith('ZMAT'):
        Azmat = geom
    else:
        Azmat = 'zmat{}'.format(str(geom))
        Axyz = 'xyz{}'.format(str(geom))


    # write geom. zmat has priority over xyz!
    if hasattr(mol, Azmat):
        zmat = getattr(mol, Azmat)
        # here the zmat
        for i in range(len(zmat)):
            for item in zmat[i]:
                fid.writelines("%s " % item)
            fid.write("\n")
        fid.write("\n")
    elif hasattr(mol, Axyz):
        xyz = getattr(mol, Axyz)
        fid.write(xyz)
        fid.write("\n")
    else:
        print "attribute not found:", Axyz
        raise AttributeError
    return


#---- START FILEWRITER2 THIS ONLY FOR MAKERS TRY TO MAKE THIS ONE UNIVERSAL ----#
def filewriter(mol, calc, pos=None):
    '''    This function creates a file with the geometry contained in zmat
    The name of the file contains the index in the name
    '''
    # ------------
    # this function uses globals: identify, path

    # and paras:
    # -path
    # -identify
    # -nprocs
    # -jobs
    # ------------
    index = mol.index
    paras = calc
    jobs = paras['jobs']

    if not pos is None: # NB it can be 0!
        name = "{0}{1}_{2}".format(paras['identify'], str(index), str(pos))
        filename = "{}.com".format(name)
        filepath = '{0}/{1}/{2}'.format(paras['path'], str(index), filename)
        geom = '{}_{:d}'.format(calc['geom'], pos)
        print "geom in filewriter:", geom
        Job = GaussianJob(filepath, calc)
        Job.pos = pos
    else:
        name = paras['identify'] + str(index)
        filename = "{}.com".format(name)
        filepath = paras['path'] + '/' + filename
        try:
            geom = calc['geom']
        except KeyError:
            geom = None
        Job = GaussianJob(filepath, calc)


    mol.addjob(Job)


    fid = open(filepath, 'w')

    mem = calc.get('mem', 1800)
    memoryline = "%mem={:d}MB\n".format(mem)

    # JOB 1
    job1 = jobs[0]
    fid.write("%chk=" + name + ".chk\n")
    fid.write(memoryline)
    if not paras['nprocs'] == 1:
        fid.write("%nprocshared=" + str(paras['nprocs']) + "\n")
    fid.write(job1['hotline'])  # first gaussianline
    fid.write("\n\n")
    fid.write(str(mol) + "\n\n")
    fid.write("{} {}\n".format(job1['charge'], job1['mult']))

    # write geom:
    writegeom(mol, fid, geom)

    # write extra inputs:
    write_extra_lines(fid, job1, name)

    fid.write('\n')
    # THE OTHER JOBS
    # for i, (charge, mult, line) in enumerate(paras['gaussianlines'][1:]):
    for i, job in enumerate(jobs[1:]):
        fid.write("--link1--\n")
        fid.write("%chk=" + name + ".chk\n")
        fid.write(memoryline)
        if not paras['nprocs'] == 1:
            fid.write("%nprocshared=" + str(paras['nprocs']) + "\n")
        fid.write(job['hotline'])
        fid.write("\n\n")
        fid.write(str(index) + " {}th calc\n\n".format(i + 2))
        if not 'allcheck' in job['hotline']:
            fid.write("{} {}\n\n".format(job['charge'], job['mult']))
        write_extra_lines(fid, job, name)
        fid.write("\n")
    fid.close()
    return Job

def write_extra_lines(fid, job, name):
    p_solvent = re.compile('scrf.*read')
    if 'rdfreq' in job['hotline']:
        fid.write(' {}\n'.format(calc['rdfreq']))
    if any(item in job['hotline'] for item in ['out=wfx', 'output=wfx']):
        fid.write('{}.wfx\n'.format(name))
    if any(item in job['hotline'] for item in ['out=wfn', 'output=wfn']):
        fid.write('{}.fwn\n'.format(name))
    if p_solvent.match(job['hotline']):
        fid.write('{}\n'.format(calc['solvent_info']))
    return


def get_paths(mols):
    ''' get all paths that need to be examined later
    this could be a population method '''
    files = []
    for molecule in mols:
        files.extend(get_molpaths(molecule))
    return files


def get_molpaths(mol):
    paths = []
    for job in mol.jobs:
        log = get_logpath(job)
        paths.append(log)
    return paths


def get_logpath(job):
    log = job.name + '.log'
    return log
