#!/usr/bin/env python

## simple script to "translate" the histograms into massfitter compatible format
## not sure if this is the smartest way to go, but...
## it should also replace the QCD templates

import sys
import ROOT
import os
import re
import optparse

sys.path.append('../scripts/')
sys.path.append('scripts/')
import myRootFunctions as mRF


proc_dict = {'VV':'vv',
             'W+jets':                   'wjets',
             'QCD':                      'qcd',
             'Z#rightarrow ll':          'zjets',
             'Single top':               'singletop',
             't#bar{t}':                 'ttbar',
             't#bar{t}161.5':            'ttbar',
             't#bar{t}163.5':            'ttbar',
             't#bar{t}166.5':            'ttbar',
             't#bar{t}169.5':            'ttbar',
             't#bar{t} 172.5':           'ttbar',
             't#bar{t}172.5':            'ttbar',
             't#bar{t}175.5':            'ttbar',
             't#bar{t}178.5':            'ttbar',
             't#bar{t}181.5':            'ttbar',
             't#bar{t}184.5':            'ttbar',
             't#bar{t}systtunep11':      'ttbarsystp11',
             't#bar{t}systtunep11tev':   'ttbarsystp11tev',
             't#bar{t}systtunep11nocr':  'ttbarsystp11nocr',
             't#bar{t}systtunep11mpihi': 'ttbarsystp11mpihi',
             't#bar{t}systq2down':       'ttbarsystq2down',
             't#bar{t}systq2up':         'ttbarsystq2up',
             't#bar{t}systmepsdown':     'ttbarsystmepsdown',
             't#bar{t}systmepsup':       'ttbarsystmepsup',
             't#bar{t}systpowhegpy':     'ttbarsystpowhegpy',
             't#bar{t}systpowheghw':     'ttbarsystpowheghw',
             'other t#bar{t}':           'ttbarV',
             'data':                     'data'}





def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inputfile' , dest='inputfile'  , help='Name of the input file containing the histograms.' , default=None)
    parser.add_option('-o', '--outputfile', dest='outputfile' , help='Name of the output file containing the histograms.', default='lxyhist.root')
    parser.add_option('-b', '--batch'     , dest='batch'      , help='Use this flag to run root in batch mode.'          , default=False, action='store_true')

    (opt, args) = parser.parse_args()

    if opt.inputfile is None:
        parser.error('No input file defined!')
    inputfile = opt.inputfile
    outputfile = opt.outputfile
    batch = opt.batch

    if batch:
        sys.argv.append( '-b' )
        ROOT.gROOT.SetBatch()

    i_file = mRF.openTFile(inputfile)

    channels = ['e','mu']
    procs =  [k for k in mRF.GetKeyNames(i_file,'')]

    hists = {}
    for proc in procs:
        keys = [k for k in mRF.GetKeyNames(i_file,proc)]

        for key in keys:
            ## Select jetlxy or secvtxmass histograms for e and mu channels:
            try: chan, var = key.split('_', 1)
            except ValueError: continue
            if not (chan == 'e'      or chan == 'mu'):         continue
            if not (var  == 'jetlxy' or var  == 'secvtxmass'): continue

            hist = mRF.getHist(i_file, proc + '/' + key)

            if var == 'jetlxy':
                var = var.replace('jet','max')
            if var == 'secvtxmass':
                var = var.replace('secvtx','jet')

            ## Extract generated top mass
            if not 'ttbar' in proc_dict[proc]:
                mass = '0'
            else:
                try:
                    mass = proc.split('}')[-1].replace('.','')
                except IndexError:
                    mass = '1725'
            if mass == '' or 'syst' in mass:
                mass = '1725'

            name = chan+'_'+var+'_'+proc_dict[proc]+'_'+mass
            # print name
            hists[name] = hist.Clone() ## FIXME (name.replace(' ',''))

    o_file = ROOT.TFile.Open(outputfile,'recreate')
    for chan in channels:
        o_file.mkdir(chan)

    for h in hists:
        channel = h.split('_')[0]
        if channel == 'e':    o_file.cd(channel)
        elif channel == 'mu': o_file.cd(channel)
        else: continue
        hists[h].Write()
    o_file.Close()

if __name__ == '__main__':
    main()
    print 'done'

