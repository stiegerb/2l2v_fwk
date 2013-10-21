###############################################
## Instructions to run the l+jets analysis:

## Create the raw histograms:
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/topljets_samples.json -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/TopAnalysis/top_5311/raw/  -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False  -s 8nh

## Also for the additional mass points and the systematic samples:
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/topljets_syst_samples.json  -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/TopAnalysis/top_5311/raw/  -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False"  -s 8nh

## Prepare the histograms and make some control plots:
## First merge the output files:
runPlotter --iLumi 19683 --inDir ~/work/TopAnalysis/top_5311/raw/ --json data/topljets_samples.json      --outFile ~/work/TopAnalysis/top_5311/plotter_wjets_raw.root --forceMerge --noPlots
runPlotter --iLumi 19683 --inDir ~/work/TopAnalysis/top_5311/raw/ --json data/topljets_syst_samples.json --outFile ~/work/TopAnalysis/top_5311/plotter_syst_raw.root  --forceMerge --noPlots
## Then make the plots (e.g. only the met ones)
runPlotter --iLumi 19683 --inDir ~/work/TopAnalysis/top_5311/raw/ --json data/topljets_samples.json  --outFile ~/work/TopAnalysis/top_5311/plotter_wjets_raw.root --useMerged --noLog --plotExt .pdf --outDir plots/ --only met
runPlotter --iLumi 19683 --inDir ~/work/TopAnalysis/top_5311/raw/ --json data/topljets_samples.json  --outFile ~/work/TopAnalysis/top_5311/plotter_wjets_raw.root --useMerged --noLog --plotExt .png --outDir plots/ --only met

## Run the W+jets background estimate:
./scripts/getWJetsNormalization.py -b -i ~/work/TopAnalysis/top_5311/plotter_wjets_raw.root

## Re-run the W+jets MC samples and apply the background scale factors:
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/topljets_samples.json  -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/TopAnalysis/top_5311/raw/  -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False @weightsFile='data/weights/top_wjetssf.root'" -s 8nh -t WJets
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/topljets_samples.json  -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/TopAnalysis/top_5311/raw/  -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False @weightsFile='data/weights/top_wjetssf.root'" -s 8nh -t W1Jets
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/topljets_samples.json  -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/TopAnalysis/top_5311/raw/  -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False @weightsFile='data/weights/top_wjetssf.root'" -s 8nh -t W2Jets
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/topljets_samples.json  -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/TopAnalysis/top_5311/raw/  -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False @weightsFile='data/weights/top_wjetssf.root'" -s 8nh -t W3Jets
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/topljets_samples.json  -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/TopAnalysis/top_5311/raw/  -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False @weightsFile='data/weights/top_wjetssf.root'" -s 8nh -t W4Jets

## Run the QCD background estimates:
./scripts/getQCDNormalization.py -b -i ~/work/TopAnalysis/top_5311/plotter_wjets_scaled.root

## Re-run the QCD MC samples:
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/topljets_samples.json  -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/TopAnalysis/top_5311/raw/  -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False @weightsFile='data/weights/top_qcdsf.root'" -s 8nh -t QCD

## Prepare the histograms and make control plots:
runPlotter --iLumi 19683 --inDir ~/work/TopAnalysis/top_5311/raw/ --json data/topljets_samples.json  --outFile ~/work/TopAnalysis/top_5311/plotter_wjets.root --forceMerge --noPlots

runPlotter --iLumi 19683 --inDir ~/work/TopAnalysis/top_5311/raw/ --json data/topljets_samples.json  --outFile ~/work/TopAnalysis/top_5311/plotter_wjets.root --useMerged --noLog --plotExt .pdf --outDir plots/ --only met
runPlotter --iLumi 19683 --inDir ~/work/TopAnalysis/top_5311/raw/ --json data/topljets_samples.json  --outFile ~/work/TopAnalysis/top_5311/plotter_wjets.root --useMerged --noLog --plotExt .png --outDir plots/ --only met

######################################################
## Prepare the histograms for the mass extraction:

./scripts/prepareHistograms.py -i ~/work/TopAnalysis/top_5311/plotter_wjets_raw.root -o lxyhist_raw.root ## why 'raw' here?
./scripts/prepareHistograms.py -i ~/work/TopAnalysis/top_5311/plotter_syst_raw.root -o lxyhist_syst.root

hadd -f lxyhist.root lxyhist_raw.root lxyhist_syst.root

FitSecVtxDistributions templ="/afs/cern.ch/user/j/jueugste/cmssw/topmass/CMSSW_5_3_11/src/UserCode/llvv_fwk/lxyhist.root"

## Finally, do the mass fits:
## for each mass point, submit to batch with wrapper script from juerg
./scripts/massfitter.py toys -n 100 -s 10000 -i FitSecVtxDistributionsWS.root -o test -m 1635 -c e -b
## with merged rootfile
./scripts/massfitter.py calib -i FitSecVtxDistributionsWS.root -r test/bias_hists_e.root -c e
./scripts/massfitter.py fit -i FitSecVtxDistributionsWS.root -r bias_summary.root -c e
