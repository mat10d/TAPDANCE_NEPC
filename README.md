# TAPDANCEV2

Porting of TAPDANCE modified version from Tim Starr. General instructions on command line operability included below.

### Set up:

#### Download mysql (8.0.19) through homebrew and set up password (password in config.pl is "tapdance2password")

	brew install mysql

	mysqladmin -u root password ‘tapdance2password’

	mysql -u root -p
	
	mysql > SET GLOBAL sql_mode=(SELECT REPLACE(@@sql_mode,'ONLY_FULL_GROUP_BY',''));
	mysql > SET GLOBAL local_infile = 'ON';
	mysql > CREATE DATABASE tapdance;
	mysql > QUIT;
	
There is no need to redownload mysql and setup the mysql password each time. For any subsequent work, proceed with below instructions.
	
#### Set up the Perl5 Database Interface driver for the MySQL database (DBD-mysql-4.050)

	cpanm DBD::mysql
	cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)
	cpan DBD::mysql

#### Specify path to bowtie (using bowtie 1.2.3--newest bowtie 1 update)

Need to add the mm10 genome folder in the indexes subfolder of downloaded bowtie-1.2.3 folder (http://bowtie-bio.sourceforge.net/index.shtml).

	PATH=$PATH:/Users/matteo/Desktop/CIS_results/bowtie-1.2.3
	
### Running TAPDANCE:

Get the seqs.tab data [here]([url](https://www.dropbox.com/scl/fi/fqlqlj6et1ezymaqqsm9g/seqs.tab?rlkey=pdizvin409j0783d57nbkf3j1&dl=0)).

Need to switch into working directory (that contains config.pl, /data, /lib folders).
	
	cd Desktop/CIS_results/TAPDANCE/final_tapdance_run
	
To map all the fragments to reads, run the first perl script.
	
	perl lib/TAPDANCE.pl
	
To run CIS identification for various metadata files, there is no need to remap fragments, only the following command is required. Please edit the $library_percent variable in config.pl if you want to change the insertion threshold.

	perl lib/TAP2.pl

Primary results will be contained in results/[name of the metadata descriptor in the metadata.tab file]. Ensure that you save and export and rename results folders, as risk of results being overwritten by re-runs.
