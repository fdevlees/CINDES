#!/bin/env python
'''this is a module that facilitates fancy printing options

   It's main part consist of building a decorator that logs the start and the end
   of a function to the screen. This is done by preceding the function definition by the line:
   @log_io()

   possible argument: signator='<character>'
   character has to have a length of 1

'''

import numpy as np
import logging
#logging.basicConfig(format='%(levelname)s:%(message)s')
#logging.basicConfig(format='%(message)s')

#logger = logging.getLogger()
#handler = logging.StreamHandler()
#logger.addHandler(handler)
#handler.setFormatter(logging.Formatter('%(message)s'))
print = logging.info


def logpopulation(mols_todo, mols_nodo):
    logpop = []
    p=logpop.append
    p("|      NEW POPULATION CONSTRUCTED:")
    p("|   mols_todo:")
    if mols_todo:
        for mol in mols_todo:
            p("|      {}".format(mol))
    else:
        p("|      -")
    p("|   mols_nodo:")
    if mols_nodo:
        for mol in mols_nodo:
            p("|      {}".format(mol))
    else:
        p("|      -")
    logging.info('\n'.join(logpop))
    return


def pre(f,**kwargs):
    printname("BEFORE: " + f.__name__,**kwargs)
    return

def post(f, time=None, **kwargs):
    if time:
        name = ' AFTER: {} | consumed {} s'.format(f.__name__, time)
    else:
        name = ' AFTER: ' + f.__name__
    printname(name,**kwargs)
    return

def printname(name,signator='*',space='\n'*1):
    message = []
    maxl=100
    l = len(name)
    rest = (l + 8)
    sides = (maxl-rest)/2
    side = sides*' '
    message.append(space + side + rest*signator)
    message.append(side + 2*signator + '  ' + str(name) + '  ' + 2*signator)
    message.append(side + rest*signator + space)
    print("\n".join(message))


from functools import wraps
def log_io(pre=pre,post=post,signator='*', print_time=False):
    assert len(signator)==1
    def decorate(f):
        @wraps(f)
        def wrapped(*args,**kwargs):
            if print_time:
                import time
                start = time.clock()
            pre(f,signator=signator)
            r = f(*args,**kwargs)
            if print_time:
                consumed = time.clock() - start
                post(f,signator=signator, time=consumed)
            else:
                post(f,signator=signator)
            return r
        return wrapped
    return decorate

def print_title(message, outline='c',signator='=',newlines=False):
    lines = message.split('\n')
    l = max([ len(line) for line in lines ])
    maxl=100
    out = []
    if newlines:out.append('')
    if outline=='l':
        out.append(signator*(l+4))
        for line in lines:
            out.append('  '+line+'  ')
        out.append(signator*(l+4))
    elif outline=='c':
        rest = ( l + 8 )
        sides = (maxl - rest)/2
        side = sides*' '
        if rest/2 == float(rest)/2:
            out.append(side + rest*signator)
        else:
            out.append(side + (rest+1)*signator)
        for line in lines:
            m = len(line)
            nspace = 2 + (l - m)/2
            if m/2 == float(m)/2:
               out.append(side + 2*signator + nspace*' ' + str(line) + nspace*' ' + 2*signator)
            else:
               out.append(side + 2*signator + nspace*' ' + str(line) + (nspace+1)*' ' + 2*signator)
        if rest/2 == float(rest)/2:
            out.append(side + rest*signator)
        else:
            out.append(side + (rest+1)*signator)
    if newlines:out.append('')
    print("\n".join(out))
    return

def sprint(n,*args,**kwargs):
    '''tries to prints the first n items of iterable objects'''
    def printitem(item):
        if isinstance(item,float):
          try:
            if float(int(item))==item and item<10000.0: #for whole numbers as floats that are not too large
                print('{:4.0f}.'.format(item),end='')
            else:
                print('{:10.4e}'.format(item),end='')
          except ValueError:
            print('{:4.0f}.'.format(0.0),end='')
        else:
            print(item, end='')
    #if kwargs:
    #    dictlist = kwargs.items()
    #    args = tuple(dictlist) + args
    minlen = min([ len(arg) for arg in args ])
    for i in range(min((n,minlen))):
        for item in args:
            try:
                a = item[i]
            except IndexError:
                a = None
                break
            finally:
                if isinstance(a,list) or isinstance(a,np.ndarray):
                    for item in a:
                        printitem(item)
                else:
                    printitem(a)
        print()
    return


def complexprint(data, func=lambda arg: arg, strfunc= lambda *s: s):
    ''' this function prints all first level items of a collection(list/tuple/dict) on one line
    with the option of applying a function on each first level item.
    a strfunc argument is available to be applied on each string encountered.

    example 1.
    - I have a large nested list. with a lot of strings / tuples / nested list and i want to print(each item on one line
    - I want also all strings that contain underscores to be splitted up:
    complexprint(data, strfunc= lambda x:x.split('_')

    example 2.
    - I have a large nested list
    - I want of each first dimension item only the minimum configuration:
    complexprint(data, func = lambda y: min(y, key= lambda x:x[1]) )

    example 3.
    - to do

    '''
    def printje(item, level=1):
        if level>10:
            print("type was:", type(item))
            print("max recursion reached")
            raise SystemExit('stop')
        if type(item) in [list,tuple,dict]:
            for it in item:
                printje(it,level=level+1)
        else:
            if type(item)==str:
                for it in strfunc(item): #string function. 
                    print(it, end='')
            else:
                print(item, end='')
        return

    for totalcycle in data:
        # totalcycle = 0func(totalcycle) # a function on zeroth level
        for siterun in totalcycle:
            item = func(siterun) # a function on first level
            printje(item)
            print()
    else:
        print("&"*20)
    return

import sys

def dump(obj, nested_level=0, ret=[], n=10):
    spacing = '   '
    if type(obj) == dict:
        ret.append( '%s{' % ((nested_level) * spacing) )
        for k, v in list(obj.items()):
            if hasattr(v, '__iter__'):
                ret.append('%s%s:' % ((nested_level + 1) * spacing, k) )
                dump(v, nested_level + 1, ret)
            else:
                ret.append( '%s%s: %s' % ((nested_level + 1) * spacing, k, v))
        ret.append( '%s}' % (nested_level * spacing))
    elif type(obj) == list or type(obj) == tuple or 'numpy' in str(type(obj)):
        ret.append( '%s[' % ((nested_level) * spacing))
        nprint = min( len(obj), n )
        for v in obj[: nprint ]:
            if hasattr(v, '__iter__'):
                dump(v, nested_level + 1, ret)
            else:
                ret.append( '%s%s' % ((nested_level + 1) * spacing, v))
        if nprint==n:
            ret.append( '%s%s' % ((nested_level + 1) * spacing, '... ({})'.format(len(obj)) ))
        ret.append('%s]' % ((nested_level) * spacing))
    else:
        ret.append( '%s%s' % (nested_level * spacing, obj))
    return '\n'.join(ret)



if __name__=='__main__':
    @log_io(signator='=')
    def A(x,y):
       print('in function A')
       return x+y
    print(A(3,4))

    print_title("C I N D E S\nAn Inverse Molecular Design Program\nwritten by Jos L. Teunissen",newlines=True)
    print_title("this would be a subtitle\nbut it is just a test",outline='l',signator='^')
