Lads, lads, lads. Hereby a quick introduction on how to work with github.
--------------------------------------------------------------------------------------------------------------------------------------------------------------
Getting started:
To get started, you first have to clone the repository. The link to clone can be found in the github website. Use git clone UrlHere . Keep in mind that it will be 
cloned into the folder where you opened git bash.
The master is the branch where all completed code ends up, so when working on code. ALWAYS use your own branch. A branch can be made by using
git checkout -b BranchName. This command is a combination of 2 commands, where checkout is used to change between branches and -b is the same a using
git branch BranchName, which makes a branch but leaves you in the original branch instead of directly checking out to the new branch.

Work, work, work, work, work:
In the cloned folder, make a folder to store your work. And just make some code and save it in this folder. Using git bash you can than see the changes made using 
git status, check everything changed to see if its not something you do not want. Then use git add . to put all changes in the queue for a commit. Afterwards
git commit -m "SomeMessage" will make a commit. This will all still be on your computer itself. If you want it to be on the net to for instance make a pull request,
use git push. Making commits is a very handy way to keep an overview of the progress of your work and using github it can be easily seen where some mistakes were 
made and a commit can be reversed if needed. 

Finished branch:
When your work is finished, you want to merge it with the master. This can be done by going to the web client and issueing a pull request. DO NOT TRY TO MERGE UNFINISHED
WORK INTO THE MASTER, THIS WILL RESULT IN A TELLY. When the pull request is accepted it will be merged to the master and will appear in everyones work. 
After this, use git pull origin master, to set your computer up to date with the remote master. This will load some merges made by other people since the last time 
you pulled or cloned. Then start a new branch using git checkout -b BranchName again. Make sure you make a new branch FROM THE MASTER.

Important lines:
git status
git add .
git commit -m "CommitMessage"
git push
git pull origin master
git checkout
git checkout -b "BranchName"

XXX Tom
