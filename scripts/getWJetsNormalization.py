#!/usr/bin/env python




import os
import sys
import optparse
import subprocess
from pprint import pprint
from array import array
import re

import ROOT
from ROOT import *

sys.path.append('../scripts/')
sys.path.append('scripts/')
import myRootFunctions as mRF


def getFiles(path):
    l = os.listdir(path)
    files = [f for f in l if f.endswith('.root')]
    return files

def getSamples(path):
    files = getFiles(path)
    ## FIXME use regular expressions
    # samples = list(set([re.sub(r'_', '', f) for f in files]))
    samples = []
    for f in files:
        li = f.split('_')
        li.pop(-1)
        name = '_'.join(li)
        samples.append(name)
    # samples = list(set(['_'.join([f.split('_').pop()]) for f in files]))
    return list(set(samples))


def mergeFiles(inputdir, outputdir, samples, force):
    print 'merging the following samples:'
    for sample in samples:
        ## hadd the files
        if force:
            cmd = 'hadd -f '+outputdir+'/'+sample+'.root '+' '.join([inputdir+'/'+f for f in samples[sample]])
        else:
            cmd = 'hadd '+outputdir+'/'+sample+'.root '+' '.join([inputdir+'/'+f for f in samples[sample]])
        print cmd
        os.system(cmd)


def makeDirectory(path):
    if not os.path.exists(path):
        print 'creating output directory '+path
        os.makedirs(path)
    else:
        print 'directory '+path+' already exists'


def getNormalization(f,chan,tag,verbose=0):
    h_wjets  = mRF.getHist(f,'W+jets/'+tag)
    h_dyll   = mRF.getHist(f,'Z#rightarrow ll/'+tag)
    h_stop   = mRF.getHist(f,'Single top/'+tag)
    h_ttbar  = mRF.getHist(f, 't#bar{t}172.5/'+tag)
    h_qcd    = mRF.getHist(f, 'QCD/'+tag)
    h_vv     = mRF.getHist(f, 'VV/'+tag)
    h_ttv    = mRF.getHist(f, 'other t#bar{t}/'+tag)
    h_others = h_qcd.Clone("others")
    h_others.Add(h_vv)
    h_others.Add(h_ttv)
    h_data   = mRF.getHist(f,'data/'+tag)
    mcsum = h_wjets.Integral()+h_stop.Integral()+h_dyll.Integral()+h_qcd.Integral()+h_vv.Integral()+h_ttv.Integral()+h_ttbar.Integral()
    if verbose > 0:
        print '========================================================================='
        print ' %-25s ' % tag
        print '-----------------------------'
        print ' %16s | %8.0f    (total yields)' % ('Single Top', h_stop.Integral())
        print ' %16s | %8.0f ' % ('Drell-Yan', h_dyll.Integral())
        print ' %16s | %8.0f ' % ('QCD',   h_qcd.Integral())
        print ' %16s | %8.0f ' % ('VV',    h_vv.Integral())
        print ' %16s | %8.0f ' % ('ttW/Z', h_ttv.Integral())
        print ' %16s | %8.0f ' % ('ttbar', h_ttbar.Integral())
        print ' %16s | %8.0f ' % ('W+jets', h_wjets.Integral())
        print '-----------------------------'
        print ' %16s | %8.0f ' % ('Sum MC', mcsum)
        print ' %16s | %8.0f ' % ('Data', h_data.Integral())

    ## sum single top to it
    # st_samples = ['T_s-channel','T_t-channel','T_tW-channel','Tbar_s-channel','Tbar_t-channel','Tbar_tW-channel']
    # st_weight = 1
    # for st_sample in st_samples:
    #     f_st = mRF.openTFile(path+'/'+st_sample+'.root')
    #     h_st = mRF.getHist(f_st,'hist_charge'+chan)
    #     for key in topinfo.samples.keys():
    #         if st_sample in key:
    #             nevents = topinfo.samples[key][0]
    #             xsec = topinfo.samples[key][1]
    #             st_weight = topinfo.lumi[chan.replace('_','')] * xsec / nevents
    #     h_st.Scale(st_weight)
    #     h_wjets.Add(h_st)


    ## Yields
    n_wjets_plus   = h_wjets.GetBinContent(3)
    n_wjets_minus  = h_wjets.GetBinContent(1)
    n_data_plus    = h_data.GetBinContent(3)
    n_data_minus   = h_data.GetBinContent(1)

    n_dyll_plus    = h_dyll.GetBinContent(3)
    n_dyll_minus   = h_dyll.GetBinContent(1)
    n_stop_plus    = h_stop.GetBinContent(3)
    n_stop_minus   = h_stop.GetBinContent(1)
    n_others_plus  = h_others.GetBinContent(3)
    n_others_minus = h_others.GetBinContent(1)
    n_ttbar_plus   = h_ttbar.GetBinContent(3)
    n_ttbar_minus  = h_ttbar.GetBinContent(1)

    ## Asymmetries
    a_wjets  = float(n_wjets_plus-n_wjets_minus)/float(n_wjets_plus+n_wjets_minus)
    a_data   = float(n_data_plus-n_data_minus)/float(n_data_plus+n_data_minus)
    a_dyll   = float(n_dyll_plus-n_dyll_minus)/float(n_dyll_plus+n_dyll_minus)
    a_stop   = float(n_stop_plus-n_stop_minus)/float(n_stop_plus+n_stop_minus)
    a_others = float(n_others_plus-n_others_minus)/float(n_others_plus+n_others_minus)
    a_ttbar  = float(n_ttbar_plus-n_ttbar_minus)/float(n_ttbar_plus+n_ttbar_minus)

    ## Scale factors
    n_corr = float(n_data_plus-n_data_minus) / a_wjets

    n_uncorr = n_wjets_plus+n_wjets_minus
    # n_uncorr = h_wjets.Integral(0,4) ## Juerg?
    weight = float(h_wjets.Integral())/n_uncorr ## is always == 1.0?

    ## We're scaling actually the DIFFERENCE, not the asymmetry
    w_sf   = (n_data_plus-n_data_minus) / (n_wjets_plus-n_wjets_minus) ## == n_corr/n_uncorr
    w_frac = n_corr / (n_data_plus+n_data_minus)  ## Equal to a_data/a_wjets

    ## Errors
    n_wjets_plus_err   = h_wjets.GetBinError(3)
    n_wjets_minus_err  = h_wjets.GetBinError(1)
    n_data_plus_err    = h_data.GetBinError(3)
    n_data_minus_err   = h_data.GetBinError(1)
    n_dyll_plus_err    = h_dyll.GetBinError(3)
    n_dyll_minus_err   = h_dyll.GetBinError(1)
    n_stop_plus_err    = h_stop.GetBinError(3)
    n_stop_minus_err   = h_stop.GetBinError(1)
    n_others_plus_err  = h_others.GetBinError(3)
    n_others_minus_err = h_others.GetBinError(1)
    n_ttbar_plus_err   = h_ttbar.GetBinError(3)
    n_ttbar_minus_err  = h_ttbar.GetBinError(1)

    sum_wjets_err   = sqrt(pow(n_wjets_plus_err,2)+pow(n_wjets_minus_err,2))
    sum_data_err    = sqrt(pow(n_data_plus_err,2)+pow(n_data_minus_err,2))
    sum_dyll_err    = sqrt(pow(n_dyll_plus_err,2)+pow(n_dyll_minus_err,2))
    sum_stop_err    = sqrt(pow(n_stop_plus_err,2)+pow(n_stop_minus_err,2))
    sum_others_err  = sqrt(pow(n_others_plus_err,2)+pow(n_others_minus_err,2))
    sum_ttbar_err   = sqrt(pow(n_ttbar_plus_err,2)+pow(n_ttbar_minus_err,2))

    a_wjets_err   = a_wjets   * sqrt(pow((sum_wjets_err  /float(n_wjets_plus-n_wjets_minus)),2)    +pow((sum_wjets_err  /float(n_wjets_plus+n_wjets_minus)),2))
    a_data_err    = a_data    * sqrt(pow((sum_data_err/float(n_data_plus-n_data_minus)),2)+pow((sum_data_err/float(n_data_plus+n_data_minus)),2))
    a_dyll_err    = a_dyll    * sqrt(pow((sum_dyll_err /float(n_dyll_plus-n_dyll_minus)),2)  +pow((sum_dyll_err /float(n_dyll_plus+n_dyll_minus)),2))
    a_stop_err    = a_stop    * sqrt(pow((sum_stop_err /float(n_stop_plus-n_stop_minus)),2)  +pow((sum_stop_err /float(n_stop_plus+n_stop_minus)),2))
    a_others_err  = a_others  * sqrt(pow((sum_others_err /float(n_others_plus-n_others_minus)),2)  +pow((sum_others_err /float(n_others_plus+n_others_minus)),2))
    a_ttbar_err   = a_ttbar   * sqrt(pow((sum_ttbar_err /float(n_ttbar_plus-n_ttbar_minus)),2)  +pow((sum_ttbar_err /float(n_ttbar_plus+n_ttbar_minus)),2))

    n_corr_err = n_corr * sqrt(pow(sum_data_err/float(n_data_plus-n_data_minus),2)+pow(a_wjets_err/a_wjets,2))
    n_uncorr_err = sqrt(h_wjets.Integral())*weight
    w_sf_err = w_sf * sqrt(pow(n_uncorr_err/n_uncorr,2)+pow(n_corr_err/n_corr,2))

    if verbose == 0:
        line = '%-16s : %5.2f+-%4.2f | %9.0f | %9.0f | %7.4f+-%6.4f | %7.4f+-%6.4f | %5.2f%% '
        print line % (tag[:-7], w_sf, w_sf_err, n_data_plus-n_data_minus, n_wjets_plus-n_wjets_minus, a_data, a_data_err, a_wjets, a_wjets_err, (h_wjets.Integral()/mcsum)*100)
    if verbose > 0:
        print '-------------------------------------------------------------------'
        print ' sample  |  N+       |  N-       |  N+ - N-  |     Asymmetry      |'
        print '-------------------------------------------------------------------'
        print ' Wjets   | %9.0f | %9.0f | %9.0f | %7.4f +- %7.4f |' % ( n_wjets_plus,  n_wjets_minus,  (n_wjets_plus -n_wjets_minus),  a_wjets, a_wjets_err)
        print ' ttbar   | %9.0f | %9.0f | %9.0f | %7.4f +- %7.4f |' % ( n_ttbar_plus,  n_ttbar_minus,  (n_ttbar_plus-n_ttbar_minus),   a_ttbar,  a_ttbar_err)
        print ' s-top   | %9.0f | %9.0f | %9.0f | %7.4f +- %7.4f |' % ( n_stop_plus,   n_stop_minus,   (n_stop_plus-n_stop_minus),     a_stop,  a_stop_err)
        print ' DY+jets | %9.0f | %9.0f | %9.0f | %7.4f +- %7.4f |' % ( n_dyll_plus,   n_dyll_minus,   (n_dyll_plus-n_dyll_minus),     a_dyll,  a_dyll_err)
        print ' others  | %9.0f | %9.0f | %9.0f | %7.4f +- %7.4f |' % ( n_others_plus, n_others_minus, (n_others_plus-n_others_minus), a_others,  a_others_err)
        print ' Data    | %9.0f | %9.0f | %9.0f | %7.4f +- %7.4f |' % ( n_data_plus,   n_data_minus,   (n_data_plus-n_data_minus),     a_data, a_data_err)
        print '-------------------------------------------------------------------'
        print ' W_sf    |    %6.3f +- %6.3f ' % (w_sf, w_sf_err)
        # print '----------------------------------------------------------------'
        # print ' Ncorr   | (%5.2f +- %5.2f) x 10^3' % (n_corr/1e3, n_corr_err/1e3)
        print '----------------------------------------------------------------'
        print ' This SF corresponds to an estimated fraction of Wjets of %5.2f%%' % (w_frac*100)
        print '    (uncorrected fraction on MC was:                      %5.2f%%' % ((h_wjets.Integral()/mcsum)*100)
        print '----------------------------------------------------------------'

    n_wjets = [n_wjets_plus, n_wjets_minus, n_wjets_plus_err, n_wjets_minus_err]
    n_data = [n_data_plus, n_data_minus, n_data_plus_err, n_data_minus_err]
    a_wjets = [a_wjets, a_wjets_err]
    a_data = [a_data, a_data_err]
    n_corr = [n_corr, n_corr_err]
    n_uncorr = [n_uncorr, n_uncorr]

    return [w_sf, w_sf_err]

    # writeToFile(tag,n_data, n_wjets, a_wjets, a_data, n_corr, n_uncorr, sf)

def writeToFile(tag,n_data, n_mc, a_mc, a_data, n_corr, n_uncorr, sf):
    outfile = open(tag+'.txt','w')
    outfile.write(' & inclusive & muon channel & electron channel \\\\ \n')
    outfile.write('\\hline \n')
    outfile.write('$N^{+}_{data}$ & '+str("%.0f" % n_data[0])+' & '+str("%.0f" % mu_n_data[0])+' & '+str("%.0f" % e_n_data[0])+'\\\\ \n')
    outfile.write('$N^{-}_{data}$ & '+str("%.0f" % n_data[1])+' & '+str("%.0f" % mu_n_data[1])+' & '+str("%.0f" % e_n_data[1])+'\\\\ \n')
    outfile.write('$N^{+}_{MC}$  & '+str("%.0f" % n_mc[0])+' & '+str("%.0f" % mu_n_mc[0])+'     & '+str("%.0f" % e_n_mc[0])+'\\\\ \n')
    outfile.write('$N^{-}_{MC}$  & '+str("%.0f" % n_mc[1])+' & '+str("%.0f" % mu_n_mc[1])+'      & '+str("%.0f" % e_n_mc[1])+'\\\\ \n')
    outfile.write('\\hline \n')
    outfile.write('$A_{data}$ & $'+str("%.3f" % a_data[0])+' \pm '+str("%.3f" % a_data[1])+'$ & $'+str("%.3f" % mu_a_data[0])+'\pm'+str("%.3f" % mu_a_data[1])+'$&$'+str("%.3f" % e_a_data[0])+'\pm'+str("%.3f" % e_a_data[1])+'$\\\\ \n')
    outfile.write('$A_{MC}$ & $'+str("%.3f" % a_mc[0])+' \pm '+str("%.3f" % a_mc[1])+'$ & $'+str("%.3f" % mu_a_mc[0])+'\pm'+str("%.3f" % mu_a_mc[1])+'$&$'+str("%.3f" % e_a_mc[0])+'\pm'+str("%.3f" % e_a_mc[1])+'$\\\\ \n')
    outfile.write('\\hline \n')
    outfile.write('$N_{MC}$ & '+str("%.0f" % n_uncorr[0])+' & '+str("%.0f" % mu_n_uncorr[0])+' & '+str("%.0f" % e_n_uncorr[0])+'\\\\ \n')
    outfile.write('\\hline \n')
    outfile.write('\\hline \n')
    outfile.write('$N_{MC}^{corr}$ & '+str("%.0f" % n_corr[0])+' & '+str("%.0f" % mu_n_corr[0])+' & '+str("%.0f" % e_n_corr[0])+' \\\\ \n')
    outfile.write('W+jets SF & $'+str("%.2f" % sf[0])+'\pm '+str("%.2f" % sf[1])+'$ & $'+str("%.2f" % mu_sf[0])+'\pm'+str("%.2f" % mu_sf[1])+'$ & $'+str("%.2f" % e_sf[0])+'\pm'+str("%.2f" % e_sf[1])+'$\\\\ \n')
    outfile.write('\\hline \n')
    outfile.close()




def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inputfile', dest='inputfile' , help='Name of the input fileectory containing the histograms.' , default=None)
    parser.add_option('-b', '--batch'    , dest='batch'     , help='Use this flag to run root in batch mode.'                , action='store_true', default=False)
    parser.add_option('-v', '--verbose'  , dest='verbose'   , help='Verbose level. [default = %default]'                     , action='store', type='int', default=0)

    (opt, args) = parser.parse_args()

    if opt.inputfile is None:
        parser.error('No input file defined!')

    inputfile = opt.inputfile
    batch = opt.batch

    if batch:
        ## run ROOT in batch mode
        sys.argv.append( '-b' )
        ROOT.gROOT.SetBatch()

    sfs = {}
    i_file = mRF.openTFile(inputfile)

    channels = ['e','mu']
    print 'channels: ',channels
    levels = set([k.split('_',1)[-1] for k in mRF.GetKeyNames(i_file,'data') if 'charge' in k and not 'chargeCheck' in k and not len(k.split('_')) == 2])
    print 'levels: ',levels
    if opt.verbose==0:
        print '----------------------------------------------------------------------------------------------------'
        print '                 :     W_sf    | diff data | diff wjets|     A data      |     A wjets     | % Wjets'
        print '----------------------------------------------------------------------------------------------------'
    for chan in channels:
        for level in levels:
            sf = getNormalization(i_file,chan,chan+'_'+level, verbose=opt.verbose)
            sfs[chan+'_'+level.replace('charge','')] = sf
        if opt.verbose==0: print '----------------------------------------------------------------------------------------------------'

    # pprint(sfs)

    ## write weights to a root file
    sf_hist = TH1F('wjetssf','wjetssf',len(sfs),0,len(sfs))
    i=0
    for l in sfs:
        i+=1
        sf_hist.GetXaxis().SetBinLabel(i,l)
        sf_hist.SetBinContent(i,sfs[l][0])
        sf_hist.SetBinError(i,sfs[l][1])
    sf_out = mRF.openTFile('data/weights/top_wjetssf.root','RECREATE')
    sf_hist.Write()
    sf_out.Close()


    # print '* Inclusive normalization:'
    # n_data, n_mc, a_mc, a_data, n_corr, n_uncorr, sf = getNormalization(inpath,'')
    # print '* e+jets normalization:'
    # e_n_data, e_n_mc, e_a_mc, e_a_data, e_n_corr, e_n_uncorr, e_sf = getNormalization(inputfile,'e')
    # print '* mu+jets normalization:'
    # mu_n_data, mu_n_mc, mu_a_mc, mu_a_data, mu_n_corr, mu_n_uncorr, mu_sf = getNormalization(inputfile,'mu')



    # outfile = open('wjetsNormalization.txt','w')
    # outfile.write(' & inclusive & muon channel & electron channel \\\\ \n')
    # outfile.write('\\hline \n')
    # outfile.write('$N^{+}_{data}$ & '+str("%.0f" % n_data[0])+' & '+str("%.0f" % mu_n_data[0])+' & '+str("%.0f" % e_n_data[0])+'\\\\ \n')
    # outfile.write('$N^{-}_{data}$ & '+str("%.0f" % n_data[1])+' & '+str("%.0f" % mu_n_data[1])+' & '+str("%.0f" % e_n_data[1])+'\\\\ \n')
    # outfile.write('$N^{+}_{MC}$  & '+str("%.0f" % n_mc[0])+' & '+str("%.0f" % mu_n_mc[0])+'     & '+str("%.0f" % e_n_mc[0])+'\\\\ \n')
    # outfile.write('$N^{-}_{MC}$  & '+str("%.0f" % n_mc[1])+' & '+str("%.0f" % mu_n_mc[1])+'      & '+str("%.0f" % e_n_mc[1])+'\\\\ \n')
    # outfile.write('\\hline \n')
    # outfile.write('$A_{data}$ & $'+str("%.3f" % a_data[0])+' \pm '+str("%.3f" % a_data[1])+'$ & $'+str("%.3f" % mu_a_data[0])+'\pm'+str("%.3f" % mu_a_data[1])+'$&$'+str("%.3f" % e_a_data[0])+'\pm'+str("%.3f" % e_a_data[1])+'$\\\\ \n')
    # outfile.write('$A_{MC}$ & $'+str("%.3f" % a_mc[0])+' \pm '+str("%.3f" % a_mc[1])+'$ & $'+str("%.3f" % mu_a_mc[0])+'\pm'+str("%.3f" % mu_a_mc[1])+'$&$'+str("%.3f" % e_a_mc[0])+'\pm'+str("%.3f" % e_a_mc[1])+'$\\\\ \n')
    # outfile.write('\\hline \n')
    # outfile.write('$N_{MC}$ & '+str("%.0f" % n_uncorr[0])+' & '+str("%.0f" % mu_n_uncorr[0])+' & '+str("%.0f" % e_n_uncorr[0])+'\\\\ \n')
    # outfile.write('\\hline \n')
    # outfile.write('\\hline \n')
    # outfile.write('$N_{MC}^{corr}$ & '+str("%.0f" % n_corr[0])+' & '+str("%.0f" % mu_n_corr[0])+' & '+str("%.0f" % e_n_corr[0])+' \\\\ \n')
    # outfile.write('W+jets SF & $'+str("%.2f" % sf[0])+'\pm '+str("%.2f" % sf[1])+'$ & $'+str("%.2f" % mu_sf[0])+'\pm'+str("%.2f" % mu_sf[1])+'$ & $'+str("%.2f" % e_sf[0])+'\pm'+str("%.2f" % e_sf[1])+'$\\\\ \n')
    # outfile.write('\\hline \n')
    # outfile.close()




if __name__ == '__main__':
    main()
