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


def getNormalization(f,chan,tag):

    h_wjets = mRF.getHist(f,'W+jets/'+tag)
    h_data  = mRF.getHist(f,'data/'+tag)

    print tag,h_wjets.Integral(),h_data.Integral()



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


    n_mc_plus    = h_wjets.GetBinContent(3)
    n_mc_minus   = h_wjets.GetBinContent(1)
    n_data_plus  = h_data.GetBinContent(3)
    n_data_minus = h_data.GetBinContent(1)
    a_mc   = float(n_mc_plus-n_mc_minus)/float(n_mc_plus+n_mc_minus)
    n_corr = float(n_data_plus-n_data_minus) / a_mc

    a_data = float(n_data_plus-n_data_minus)/float(n_data_plus+n_data_minus)

    ## error
    n_mc_plus_err    = h_wjets.GetBinError(3)
    n_mc_minus_err   = h_wjets.GetBinError(1)
    n_data_plus_err  = h_data.GetBinError(3)
    n_data_minus_err = h_data.GetBinError(1)

    sum_mc_err   = sqrt(pow(n_mc_plus_err,2)+pow(n_mc_minus_err,2))
    sum_data_err = sqrt(pow(n_data_plus_err,2)+pow(n_data_minus_err,2))

    a_mc_err   = a_mc * sqrt(pow((sum_mc_err/float(n_mc_plus-n_mc_minus)),2)+pow((sum_mc_err/float(n_mc_plus+n_mc_minus)),2))
    n_corr_err = n_corr * sqrt(pow(sum_data_err/float(n_data_plus-n_data_minus),2)+pow(a_mc_err/a_mc,2))

    a_data_err   = a_data * sqrt(pow((sum_data_err/float(n_data_plus-n_data_minus)),2)+pow((sum_data_err/float(n_data_plus+n_data_minus)),2))



    ##
    n_uncorr = h_wjets.Integral(0,4)
    weight = float(h_wjets.Integral())/n_uncorr
    n_uncorr_err = sqrt(h_wjets.Integral())*weight

    w_sf = n_corr/n_uncorr
    w_sf_err = w_sf * sqrt(pow(n_uncorr_err/n_uncorr,2)+pow(n_corr_err/n_corr,2))

    w_frac = n_corr / (n_data_plus+n_data_minus)

    print '-----------------------------------------------------'
    print 'N+(MC)  : %10.2f +- %5.2f' % (n_mc_plus, n_mc_plus_err)
    print 'N-(MC)  : %10.2f +- %5.2f' % (n_mc_minus, n_mc_minus_err)
    print 'N+(data): %10.0f +- %5.2f' % (n_data_plus, n_data_plus_err)
    print 'N-(data): %10.0f +- %5.2f' % (n_data_minus, n_data_minus_err)
    print 'A (MC)  :  %9.4f +- %6.4f' % (a_mc, a_mc_err)
    print 'A (data):  %9.4f +- %6.4f' % (a_data, a_data_err)
    print '----------------------------------------------------'
    print 'Ncorr   : %10.2f +- %5.2f' % (n_corr, n_corr_err)
    print 'Nuncorr : %10.2f +- %5.2f' % (n_uncorr, n_uncorr_err)
    print '----------------------------------------------------'
    print 'W_sf    :  %9.4f +- %9.4f ' % (w_sf, w_sf_err)
    print '===================================================='
    print 'Fraction: ', w_frac
    print '----------------------------------------------------'

    n_mc = [n_mc_plus, n_mc_minus, n_mc_plus_err, n_mc_minus_err]
    n_data = [n_data_plus, n_data_minus, n_data_plus_err, n_data_minus_err]
    a_mc = [a_mc, a_mc_err]
    a_data = [a_data, a_data_err]
    n_corr = [n_corr, n_corr_err]
    n_uncorr = [n_uncorr, n_uncorr]
    sf = [w_sf, w_sf_err]

    return sf

    # writeToFile(tag,n_data, n_mc, a_mc, a_data, n_corr, n_uncorr, sf)

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
    for chan in channels:
        for level in levels:
            sf = getNormalization(i_file,chan,chan+'_'+level)
            sfs[chan+'_'+level.replace('charge','')] = sf

    pprint(sfs)

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
