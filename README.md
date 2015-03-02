Everything starts here
====

Log into your github account and click "Fork" at the top-right corner of this
page to get your own copy of this repository.

After that, open a Terminal and log into gridui:

    ssh -Y gridui.ifca.es -o ServerAliveInterval=240

Set the CMS environment:

    source /cvmfs/cms.cern.ch/cmsset_default.sh

And the architecture:

    export SCRAM_ARCH=slc6_amd64_gcc481

Now we can choose our favorite CMSSW release.

    cmsrel CMSSW_7_2_0
    cd CMSSW_7_2_0/src
    cmsenv


Get the material
====

    git clone https://github.com/your_github_username/MuonTreeAnalyzer.git MuonTreeAnalyzer


File description
====

*muonAnalyzer contains the main analysis functions.
 
*RunPROOF_muonAnalyzer executes muonAnalyzer using PROOF. Edit it if you need to
modify the location of the source rootfiles or add new samples.

*runAll_muonAnalyzer.sh is just a quick executable.

The rest of the files shouldn't be modified.

Let's run the code
====

Edit the file runAll_muonAnalyzer.sh:

     vim runAll_muonAnalyzer.sh

Comment/uncomment the necessary lines depending on the samples you want to run.

Now load root and PAF:

    source /gpfs/csic_projects/cms/sw/ROOT/current/root/bin/thisroot.sh
    export PAFPATH=/gpfs/csic_projects/cms/PROOF/paf/
    export PATH=$PAFPATH/bin:$PATH

And run the code:

    ./runAll_muonAnalyzer.sh


It is commit time
====

First get the latest changes in the repository, if any:

    git pull https://github.com/juanracasti/MuonTreeAnalyzer.git

And then commit your changes:

    git status
    git add <filepattern>
    git commit -m 'Modified'
    git push

Now open a new tab in your browser with your copy of the repository:

    https://github.com/your_github_username/MuonTreeAnalyzer

Open a new pull request with your commits and a brief description




