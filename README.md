# TAPDANCEV2

Porting of TAPDANCE modified version from Tim Starr. General instructions on command line operability included below.

### Set up:

##### Download mysql (8.0.19) through homebrew and set up password (password in config.pl is "tapdance2password")

	brew install mysql

	mysqladmin -u root password ‘tapdance2password’

	mysql -u root -p
	
	mysql > SET GLOBAL sql_mode=(SELECT REPLACE(@@sql_mode,'ONLY_FULL_GROUP_BY',''));
	mysql > SET GLOBAL local_infile = 'ON';
	mysql > CREATE DATABASE tapdance;
	mysql > QUIT;
	
##### Set up the Perl5 Database Interface driver for the MySQL database (DBD-mysql-4.050)

	cpanm DBD::mysql
	cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)
	cpan DBD::mysql

##### Specify path to bowtie (using bowtie 1.1.2)

Need to add the mm9 genome in the indexes subfolder of downloaded bowtie-1.1.2 folder.

	PATH=$PATH:/Users/matteodibernardo/Desktop/bowtie-1.1.2/
	
### Running TAPDANCE:

Need to switch into working directory (that contains config.pl, /data, /lib folders).
	
	perl lib/TAPDANCE.pl
	perl lib/TAP2.pl

Primary results will be contained in results/all/cis_all-nr-colon-0.0001
