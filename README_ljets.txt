instructions to run the l+jets analysis:

**create the raw histograms:

runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/topljets_samples.json -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/top_5311/raw/  -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False  -s 8nh



** also for the additional mass points and teh systematic samples:

runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/topljets_syst_samples.json  -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/top_5311/raw/  -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False"  -s 8nh




**prepare the histograms and make some control plots:

runPlotter --iLumi 19683 --inDir ~/work/top_5311/raw/  --json data/topljets_samples.json   --outFile ~/work/top_5311/plotter_wjets_raw.root   --forceMerge --noLog



**run the W+jets background estimate:

./scripts/getWJetsNormalization.py -b -i ~/work/top_5311/plotter_wjets_raw.root



**re-run the W+jets MC samples and apply the background scale factors:

runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/topljets_samples.json  -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/top_5311/raw/  -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False @weightsFile='data/weights/top_wjetssf.root'" -s 8nh -t WJets

runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/topljets_samples.json  -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/top_5311/raw/  -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False @weightsFile='data/weights/top_wjetssf.root'" -s 8nh -t W1Jets

runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/topljets_samples.json  -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/top_5311/raw/  -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False @weightsFile='data/weights/top_wjetssf.root'" -s 8nh -t W2Jets

runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/topljets_samples.json  -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/top_5311/raw/  -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False @weightsFile='data/weights/top_wjetssf.root'" -s 8nh -t W3Jets

runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/topljets_samples.json  -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/top_5311/raw/  -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False @weightsFile='data/weights/top_wjetssf.root'" -s 8nh -t W4Jets



**run the QCD background estimates:



**re-run the QCD MC samples:


**prepare the histograms and make some control plots:

runPlotter --iLumi 19683 --inDir ~/work/top_5311/raw/  --json data/topljets_samples.json   --outFile ~/work/top_5311/plotter_wjets_raw.root   --forceMerge --noLog

runPlotter --iLumi 19683 --inDir ~/work/top_5311/raw/  --json data/topljets_syst_samples.json   --outFile ~/work/top_5311/plotter_syst_raw.root   --forceMerge --noLog


** prepare the histograms for the mass extraction:

./scripts/prepareHistograms.py -i ~/work/top_5311/plotter_wjets_raw.root -o lxyhist_raw.root

./scripts/prepareHistograms.py -i ~/work/top_5311/plotter_syst_raw.root -o lxyhist_syst.root

hadd -f lxyhist.root lxyhist_raw.root lxyhist_syst.root



FitSecVtxDistributions templ="/afs/cern.ch/user/j/jueugste/cmssw/topmass/CMSSW_5_3_11/src/UserCode/llvv_fwk/lxyhist.root"
