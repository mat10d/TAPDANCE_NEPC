# TAPDANCEV2
Porting of TAPDANCE modified version from Tim Starr

## Download mysql (8.0.19) through homebrew and set up password

brew install mysql

mysqladmin -u root password ‘tapdance2password’

mysql -u root -p

	SET GLOBAL sql_mode=(SELECT REPLACE(@@sql_mode,'ONLY_FULL_GROUP_BY',''));

	SET GLOBAL local_infile = 'ON';

	CREATE DATABASE tapdance;

	QUIT;
	
## Set up the Perl5 Database Interface driver for the MySQL database (DBD-mysql-4.050)

cpanm DBD::mysql

cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)

cpan DBD::mysql

## Specify path to bowtie (using 1.1.2)

PATH=$PATH:/Users/matteodibernardo/Desktop/bowtie-1.1.2/
