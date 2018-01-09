ConPair-master tracked files
Bruce Wollison 20180109
ConPair is a concordance/contamination detection tool (https://www.ncbi.nlm.nih.gov/pubmed/27354699) available in repo on git: https://github.com/nygenome/conpair
Current version using 0.1
While this tool is useful for concordance/contamination testing, the original code requires 3 seperate command calls resulting in 2 output files. Since we plan to use it in the same workflow everytime, and also want to hook it up to our DB to start automatically, we created a small workflow management system that starts with querying the Riker DB to get the samples that need to be run, finishing with organizing the results into a report. This required refactoring the 3 ConPair scripts, since they weren't previous callable as python modules. The 3 tracked files are the refactored ConPair scripts, while the rest of the ConPair repo can substitute all other files those scripts call and should be placed in the same folder.
 
