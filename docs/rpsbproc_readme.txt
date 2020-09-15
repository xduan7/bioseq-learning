Post-RPSBLAST Processing Utility v0.1      README file, revised 11 September 2016
================================================================================
https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/rpsbproc/README


CONTENTS:

1. SUMMARY
  1.1 Standalone RPS-BLAST
  1.2 rpsbproc command line utility
2. DOWNLOAD THE SOFTWARE PACKAGE
3. DOWNLOAD DATA FILES
4. RUN THE PROGRAMS
5. PARSE RPSBPROC OUTPUT
6. SAMPLE OUTPUT FROM RPSBPROC

-----------------------------
1. SUMMARY
-----------------------------
The rpsbproc command line utility is an addition to the standalone version of
Reverse Position-Specific BLAST (RPS-BLAST), also known as CD-Search (Conserved
Domain Search).

RPS-BLAST is used to identify conserved domains, or functional units, within a
query sequence. RPS-BLAST searches a protein sequence (or a protein translation
of a nucleotide sequence) against a database of profiles that represent
conserved domains. This is the opposite of PSI-BLAST, which searches a profile
against a database of protein sequences, hence the term 'Reverse'.

The next two sections provide a brief description of the output produced by
Standalone RPS-BLAST, and the output produced by the rpsbproc command line
utility, which post-processes the RPS-BLAST output.

-----------------------------
1.1 Standalone RPS-BLAST
-----------------------------

The standalone RPS-BLAST is packaged with the BLAST executables (available
on the NCBI FTP site at https://ftp.ncbi.nih.gov/blast/executables/LATEST/),
and is also available as part of the NCBI toolkit distribution
(see https://ftp.ncbi.nih.gov/toolbox). Additional details are provided in
https://ftp.ncbi.nih.gov/blast/documents/rpsblast.html.

For each query sequence, standalone RPS-BLAST lists the conserved domain
models that scored above a certain threshold (default set to an evalue of 10),
sorted by scores. The information provided for each hit includes the
conserved domain's PSSMID, a set of scores (e-value, bitscore, etc)
and the actual sequence alignment between the conserved domain and
the query sequence.

The actual program that runs RPS-BLAST locally, "rpsblast," has existed since
early Fall of 2000. In October 2008, with the release of BLAST 2.2.18,
the RPS-BLAST executable was split into two separate programs: "rpsblast" for
protein queries and "rpstblastn" for nucleotide queries. Both programs search
a query sequence against a database of profiles, and using a BLAST-like
algorithm, look for alignments between the query sequences and the profiles of
conserved domains. They produce results such as those described in the
preceding paragraph.

-----------------------------
1.2 rpsbproc command line utility
-----------------------------

The rpsbproc command line utility is available, as of December 2014,
from the Conserved Domain Database (CDD) FTP site:
https://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/

It post-processes the results of local RPS-BLAST searches in order to provide a
non-redundant view of the search results, and to provide additional annotation
on query sequences, such as domain superfamilies and functional sites, similar to
the annotation provided by the corresponding web services (e.g., the NCBI Batch
CD-Search web service at
http://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi).

Specifically, the rpsbproc utility reads the output of rpsblast/rpstblastn,
fills in domain superfamily and functional site information for each region of
the sequence, re-sorts the hits by a different standard, and calculates a
set of non-redundent representative hits. In this way, it turns the
raw alignments into domain/site annotations on the query sequence at different
redundancy level, basically produce the same data as web-based Batch CD-Search
service does. The annotation data is presented in tab-delimited tables to be processed
either programatically or manually with a spreadsheet (see details below).

See the CDD and CD-Search help document for additional details about
superfamilies, conserved sites, and more:

http://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml

-------------------------------------------------
2. DOWNLOAD THE SOFTWARE PACKAGES
-------------------------------------------------
The rpsbproc package is available at:
	https://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/

Rpsbproc binaries are build for Win-x64 and linux-x64, they are directly executable,
no installation required. Rpsblast and rpstblastn also have pre-built executables for
multiple systems that are available at:
	https://ftp.ncbi.nih.gov/blast/executables/LATEST/


Download and untar the tarball:

	$ wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/RpsbProc-x64-linux.tar.gz
	$ tar -xzf RpsbProc-x64-linux.tar.gz
	$ ls -l
		drwxr-xr-x 3 shlu structure    4096 Sep 26 21:12 RpsbProc-x64-linux
		-rw-r--r-- 1 shlu structure 5368503 Sep 26 21:12 RpsbProc-x64-linux.tar.gz


The RpsbProc-x64-linux directory should have the following layout:

	.
	|-- README
	|-- rpsbproc
	|-- rpsbproc.ini
	`-- utils
		|-- getblbins.sh
		`-- getcdddata.sh



For those who need (or desire) to build these utilities locally can download the source code
tarballs. The rpsbproc program, as well as rpsblast and rpstblastn programs, are NCBI C++ toolkit
applications and therefore require NCBI C++ toolkit to build. The blast source archive is
actually a subset of NCBI C++ toolkit, with all blast programs added in. It is recommended
to download it and insert the rpsbproc source code in as an additional application and build all programs
together. Please download the source archive https://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/RpsbProc-src.tar.gz
and follow instructions in the README file there.


Rpsblast and rpstblastn also have pre-built executables for multiple systems that are available at:
	https://ftp.ncbi.nih.gov/blast/executables/LATEST/

Download the one for your platform and untar it. At the time or writing, the latest package for linux is ncbi-blast-2.9.0+-x64-linux.tar.gz.

	$ cd RpsbProc-x64-linux
	$ wget https://ftp.ncbi.nih.gov/blast/executables/LATEST/ncbi-blast-2.9.0+-x64-linux.tar.gz
	$ tar -xzf ncbi-blast-2.9.0+-x64-linux.tar.gz

Copy the rpsblast and rpstblastn programs to current directory:

	$ cp ncbi-blast-2.9.0+/bin/rpsblast ncbi-blast-2.9.0+/bin/rpstblastn .
	$ ls -l

		-rwxr-xr-x 1 shlu structure     32331 Sep 26 13:21 README
		drwxr-xr-x 4 shlu structure      4096 Mar 11  2019 ncbi-blast-2.9.0+
		-rw-r--r-- 1 shlu structure 243817734 Mar 11  2019 ncbi-blast-2.9.0+-x64-linux.tar.gz
		-rwxr-xr-x 1 shlu structure  43489856 Sep 26 22:26 rpsblast
		-rwxr-xr-x 1 shlu structure  15739416 Sep 26 17:04 rpsbproc
		-rwxr-xr-x 1 shlu structure       240 Sep 13  2018 rpsbproc.ini
		-rwxr-xr-x 1 shlu structure  43485504 Sep 26 22:26 rpstblastn
		drwxr-xr-x 2 shlu structure      4096 Sep 16 16:50 utils

On linux platform, these steps can be done automatically by running utils/getblbins.sh (Get Blast Binaries). The script
automatically retrieves the latest linux binary and create symbolic links from the current directory to rpsblast and rpstblastn:

	$ utils/getblbins.sh
	$ ls -l
		-rwxr-xr-x 1 shlu structure    32331 Sep 26 13:21 README
		drwxr-xr-x 3 shlu structure     4096 Sep 26 23:29 ncbi-blast-2.9.0+-x64-linux
		lrwxrwxrwx 1 shlu structure       60 Sep 26 23:29 rpsblast -> ./ncbi-blast-2.9.0+-x64-linux/ncbi-blast-2.9.0+/bin/rpsblast
		-rwxr-xr-x 1 shlu structure 15739416 Sep 26 17:04 rpsbproc
		-rwxr-xr-x 1 shlu structure      240 Sep 13  2018 rpsbproc.ini
		lrwxrwxrwx 1 shlu structure       62 Sep 26 23:29 rpstblastn -> ./ncbi-blast-2.9.0+-x64-linux/ncbi-blast-2.9.0+/bin/rpstblastn
		drwxr-xr-x 2 shlu structure     4096 Sep 26 23:10 utils

--------------------------------------
3. DOWNLOAD DATA FILES
--------------------------------------

The rpsblast and rpstblastn executables need domain databases in order to work.
Download pre-formatted search databases from:

	https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian

Create db directory and download database files
	$ mkdir ./db
	$ wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz && tar -xzf Cdd_LE.tar.gz -C ./db && rm -f Cdd_LE.tar.gz
	$ wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_NCBI_LE && tar -xzf Cdd_NCBI_LE -C ./db && rm -f Cdd_NCBI_LE
	$ wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cog_LE.tar.gz && tar -xzf Cog_LE.tar.gz -C ./db && rm -f Cog_LE.tar.gz
	$ wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Kog_LE.tar.gz && tar -xzf Kog_LE.tar.gz -C ./db && rm -f Kog_LE.tar.gz
	$ wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Prk_LE.tar.gz && tar -xzf Prk_LE.tar.gz -C ./db && rm -f Prk_LE.tar.gz
	$ wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Smart_LE.tar.gz && tar -xzf Smart_LE.tar.gz -C ./db && rm -f Smart_LE.tar.gz
	$ wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Tigr_LE.tar.gz && tar -xzf Tigr_LE.tar.gz -C ./db && rm -f Tigr_LE.tar.gz


The rpsbproc program will need a set of domain-annotation data files to function.
Download these files from:

	https://ftp.ncbi.nih.gov/pub/mmdb/cdd/

Create data directory and download domain-annotation files

	$ mkdir ./data
	$ wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz -O ./data/cddid.tbl.gz && gzip -d ./data/cddid.tbl.gz
	$ wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdtrack.txt -O ./data/cdtrack.txt
	$ wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/family_superfamily_links -O ./data/family_superfamily_links
	$ wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddannot.dat.gz -O ./data/cddannot.dat.gz && gzip -d ./data/cddannot.dat.gz
	$ wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddannot_generic.dat.gz -O ./data/cddannot_generic.dat.gz && gzip -d ./data/cddannot_generic.dat.gz
	$ wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/bitscore_specific.txt -O ./data/bitscore_specific.txt

On linux platform, these steps can be done automatically by running utils/getcdddata.sh (Get CDD Data). The script
will obtain all databases and datafiles, decompress them if necessary, and put in ./data and ./db by default.

	$ ls -l
		-rwxr-xr-x 1 shlu structure    32291 Sep 18 16:41 README
		drwxr-xr-x 2 shlu structure     4096 Sep 18 17:09 data
		drwxr-xr-x 2 shlu structure     4096 Sep 18 17:11 db
		drwxr-xr-x 2 shlu structure     4096 Sep 18 16:53 example
		drwxr-xr-x 3 shlu structure     4096 Sep 18 17:06 linux
		drwxr-xr-x 2 shlu structure     4096 Sep 18 17:10 logs
		-rwxr-xr-x 1 shlu structure       27 Sep 13 09:18 precedences_sample.txt
		lrwxrwxrwx 1 shlu structure       38 Sep 18 17:07 rpsblast -> ./linux/ncbi-blast-2.9.0+/bin/rpsblast
		lrwxrwxrwx 1 shlu structure       40 Sep 18 17:07 rpstblastn -> ./linux/ncbi-blast-2.9.0+/bin/rpstblastn
		-rwxr-xr-x 1 shlu structure 15670008 Sep 18 13:40 sparclbl
		-rwxr-xr-x 1 shlu structure      446 Sep 18 16:25 sparclbl.ini
		drwxr-xr-x 2 shlu structure     4096 Sep 16 16:50 utils

	$ ls -l data/
		-rw-r--r-- 1 shlu structure  1014866 Sep 18 17:09 bitscore_specific.txt
		-rw-r--r-- 1 shlu structure  3875945 Sep 18 17:09 cddannot.dat
		-rw-r--r-- 1 shlu structure  2455749 Sep 18 17:09 cddannot_generic.dat
		-rw-r--r-- 1 shlu structure 26712116 Sep 18 17:09 cddid.tbl
		-rw-r--r-- 1 shlu structure  1324235 Sep 18 17:09 cdtrack.txt
		-rw-r--r-- 1 shlu structure  1674528 Sep 18 17:09 family_superfamily_links
		-rw-r--r-- 1 shlu structure 29836008 Sep 18 17:09 specific_arch.txt
		-rw-r--r-- 1 shlu structure 17694472 Sep 18 17:09 superfamily_arch.txt


The README file on the CDD FTP site provides information about the
little-endian databases and the scope of data contained in each one.
It also provides information about customizing search databases:

     https://ftp.ncbi.nih.gov/pub/mmdb/cdd/README

Please note that we can no longer maintain pre-computed search databases for
big_endian architectures. Search databases distributed via the big_endian
FTP directory are outdated. (As mentioned in the descriptive file
(https://ftp.ncbi.nih.gov/blast/documents/rpsblast.html) for the RPS-BLAST
executable, RPS-BLAST uses a BLAST database, but also has some other files that
contain a precomputed lookup table for the profiles to allow the search to
proceed faster.  Unfortunately it was not possible to make this lookup table
architecture independent (like the BLAST databases themselves) and one cannot
take an RPS-BLAST database prepared on a big-endian system (e.g., Solaris Sparc)
and run it on a small-endian system (e.g., NT). The RPS-BLAST database must be
prepared again for the small-endian system.)


Now the local RPSBLAST system is ready to run.

--------------------------------------
4. RUN THE PROGRAMS
--------------------------------------
If the above steps are followed to prepare the package, the current directory
should look like:
	.
	|-- README
	|-- data
	|   |-- bitscore_specific.txt
	|   |-- cddannot.dat
	|   |-- cddannot_generic.dat
	|   |-- cddid.tbl
	|   |-- cdtrack.txt
	|   `-- family_superfamily_links
	|-- db
	|   |-- Cdd.aux
	|   |-- Cdd.freq
	|  ......
	|-- rpsblast
	|-- rpsbproc
	|-- rpstblastn
	`-- utils
		|-- getblbins.sh
		`-- getcdddata.sh
**IMPORTANT**: The filenames under ./data are important! make sure they are correct
after downloading from our ftp site and rename them if not match.

Both rpsblast and rpstblastn support query files with multiple FASTA sequences
and can deliver the results in many formats (run them with the "-help" for
detailed usage information), but rpsbproc only processes the asn1 formatted
Blast Archive object (by specifying "-outfmt 11" when running rpsblast and rpstblastn)
because it contains all the information necessary for rpsbproc to function. For backward
compatibility, this version still supports xml formatted BlastOutput object (by specifying "-outfmt 5"),
but it is deprecated and maybe removed in future versions. Data generated
by rpsblast/rpstblastn can be saved in a file then processed by rpsbproc:

	[~/localrpsb]$ rpsblast -query multi_protein.FASTA -db ./db/Cdd -evalue 0.01 -outfmt 11 -out multi_protein.asn
	[~/localrpsb]$ rpsbproc -i multi_protein.asn -o multi_protein.out -e 0.01 -m rep

Or be piped directly to rpsbproc:

	[~/localrpsb]$ rpsblast -query multi_protein.FASTA -db ./db/Cdd -evalue 0.01 -outfmt 11 | rpsbproc -i multi_protein.asn -o multi_protein.out -e 0.01 -m rep

Command line options supported by rpsbproc:

-i<in-filename>
	Input file that contains the XML data generated by rpsblast /rpstblastn.
	By default the program reads data from standard input.

-o<out-filename>
	File to which the processed results should be written. By default the
	output is directed to standard output.

-e<evalue-cutoff>
	Only processes hits with evalues equal or better (smaller in value)
	than the designated. Default value is 0.01. This option is convenient
	in that the rpsblast needs only be run once with a higher evalue
	(such as 1.0) and the results can be processed by rpsbproc with
	difference evalue-cutoffs.

-m<data-mode>
	Select redundancy level of domain hit data. Valid options are "rep"
	(concise), "std"(standard), and "full" (all hits). Please refer to http://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#RPSBResultsLevelsOfDetail
	for detailed information regarding the three redundancy levels.
	The default value is "rep".

-t<target-data>
	Select desired (target) data. Valid options are "doms" for
	domain/superfamily hits, "feats" for functional sites and
	structural motif annotations, or "both". Default is "both".

-d<datapath>
	Location of data files. By default, the program looks for a directory named 'data' the following order:
		current directory
		program directory (the directory where the rpsbproc program resides)

-q
	Quiet mode -- do not display program information and processing progress
	to stderr.

-h
	Display a brief usage screen, ignore other switches.

--------------------------------------
5. PARSE RPSBPROC OUTPUT
--------------------------------------

The rpsbproc executable presents the data in a tab-delimited table-like format.
It is designed such that it tries to accommondate both human read and program
processing. A sample output file is included at the bottom of this document.

There are two sections in an output file. The first section displays the program
information, parameters used for data processing, and a "template" to explain
the format and content of each column of the data tables. This section is
intended for human read. All lines in this section start with a '#' so programs
can treat them as "comment" lines that can be safely ignored.

The second section contains the real data intended to be programmatically
processed. All columns are delimited with a 'tab' character ('\t'). The data
section always starts with a "DATA" token and ends with an "ENDDATA" token.
In between, there can be several sessions, each of which start with a SESSION
token and end with an ENDSESSION token. A session is one XML object produced by
rpsblast or rpstblastn. Despite the number of FASTA sequences in a query file,
each run of rpsblast/rpstblastn produces only one XML object. Besides the
SESSION token, the line also contains several columns to display additional
information about the session:

	SESSION	<session-ordinal>	<program>	<database>	<score-matrix>	<evalue-threshold>

Each session is given an ordinal number to distinguish it from other sessions,
followed by the program/version, database, score matrix, and e-value threshold
used for the rpsblast/rpstblastn session. The session end line contains the
ENDSESSION token and the session ordinal number to indicate the end of the
particular session:

	ENDSESSION	<session-ordinal>

In each session, there can be multiple queries, starting with a QUERY token
line and ending with an ENDQUERY token line:

	QUERY	<query-id>	<seq-type>	<seq-length>	<definition-line>
	ENDQUERY

     The <query-id> is an ID assigned automatically by rpsblast/rpstblastn during
     a session, so is only unique within the session.
     The <seq-type> is either "Peptide" or "Nucleotide".
     The <seq-length> is the length of the query sequence (number of amino acid
     or nucleotide residues).
     The <definition-line> is the exact copy of the definition line of the FASTA
     sequence.

Within each QUERY section, there are three optional sub-sections:

	1. Domain/superfamily hits. If requested (with command line option
	"-t doms" or "-t both"), this section will appear, starting with a
	DOMAINS token and ending with an ENDDOMAINS token:

		DOMAINS
		ENDDOMAINS

		Within the domain hit section, each hit is repesented by one line
		with the following columns:

			<session-ordinal>	<query-id[readingframe]>	<hit-type>	<PSSM-ID>	<from>	<to>	<E-Value>	<bitscore>	<accession>	<short-name>	<incomplete>	<superfamily PSSM-ID>

		The session ordinal number and the query id are included in each
		line to make it unique across the whole rpsbproc session.
		The "[readingframe]" part following the query id only appears for
		nucleotide queries. Other columns are:

			<hit-type>: "Specific", "Non-Specific", "Superfamily" or
			"Multidom". Refer to
			http://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#RPSB_hit_types
			for more information.

			<PSSM-ID>: the PSSMID of the domain or superfamily

			<from>: Domain start position in query sequence coordinates

			<to>: Domain end (inclusive) position in query sequence coordinates

			<E-Value>: The E-value of this hit

			<bitscore>: The bitscore of this hit

			<accession>: The accession of the domain or superfamily

			<short-name>: The short-name of the domain or superfamily

			<incomplete>: Indicate if there are more than 20% missing
			from N- or C- terminal compare to the original domain.
			Values are:
			"-" (no more than 20% shorter on either terminals),
			"N" (N-terminal has 20% or more missing),
			"C"(C-terminal has 20% or more missing) or
			"NC" (both terminal has 20% or more missing)

			<superfamily PSSM-ID>: If applicable, this field contains
			the PSSMID of the superfamily to which the hit domain
			belongs. Otherwise, it will contain a dash '-'.

	2. Site annotations. If requested (with command line option "-t doms"
	or "-t both"), this section will appear, starting with a SITES token
	and ending with an ENDSITES token:

			SITES
			ENDSITES

		Each functional site is represented by one line, which contains
		the following columns:
			<session-ordinal>	<query-id[readingframe]>	<annot-type>	<title>	<residue(coordinates)>	<complete-size>	<mapped-size>	<source-domain>

		The first two columns are the same as in domain hit lines.
		Other columns are:

			<annot-type>: "Specific" (site-annotations found from
			specific domain hits) or "Generic" (site-annotations
			found from non-specific/superfamily hits).

			<title>: Description title of the functional site

			<residue(coordinates)>: one-letter residue symbols and
			the coordinate on the query sequence, each residue is
			separated by a comma (,). For protein queries,
			coordinate is just one number; for nucleotide queries,
			the residue symbol is still the translated amino acid,
			the coordinates, however, is a three-position range for
			the codon on query sequence, such as I44-46,F47-49.
			For minus strand, the range is like I6926-6924 with
			the bigger number comes first to indicate the reversed
			direction.

			<complete-size>: The number of residues in the original
			defined functional site.

			<mapped-size>: The number of residues of the functional
			site mapped to the query sequence.

			<source-domain>: The PSSMID of the domain on which the
			site is defined. For specific sites, it's the
			specific hit domain; for generic sites, it's the
			root domain of the domain hierarchy where the
			non-specific hit domain belongs.

	3. Structural motifs. This section appears with site annotations,
	   starting with a MOTIFS token and ending with an ENDMOTIFS token:

			MOTIFS
			ENDMOTIFS

		Each structural motif is represented by one line, which contains
		the following columns:
			<session-ordinal>	<query-id[readingframe]>	<title>	<from>	<to>	<source-domain>

		The first two columns are the same as in domain hit lines.
		Other columns are:

		<title>: The description title of the motif

		<from>,<to>: The range of the motif in query sequence
		coordinates. Structural motifs are defined as ranges instead of
		multiple residues

		<source-domain>: Some structural motifs are defined on conserved
		domains, whose PSSMID will appear here; other structural motifs
		are not defined on CDs thus do not have a source domain.

--------------------------------------
6. SAMPLE OUTPUT FROM RPSBPROC
--------------------------------------

#Post-RPSBLAST Processing Utility v0.1
#Config file:	rpsbproc.ini
#Input data file:	test/sample.xml
#Output data file:	test/sample.out
#E-Value cutoff:	0.01
#Redundancy:	Concise
#Data requested:	Domain hits and site annotations
#Output format -- tab-delimited table
#DATA
#SESSION	<session-ordinal>	<program>	<database>	<score-matrix>	<evalue-threshold>
#QUERY	<query-id>	<seq-type>	<seq-length>	<definition-line>
#DOMAINS
#<session-ordinal>	<query-id[readingframe]>	<hit-type>	<PSSM-ID>	<from>	<to>	<E-Value>	<bitscore>	<accession>	<short-name>	<incomplete>	<superfamily PSSM-ID>
#more such lines......
#ENDDOMAINS
#SITES
#<session-ordinal>	<query-id[readingframe]>	<annot-type>	<title>	<residue(coordinates)>	<complete-size>	<mapped-size>	<source-domain>
#more such lines......
#ENDSITES
#MOTIFS
#<session-ordinal>	<query-id[readingframe]>	<title>	<from>	<to>	<source-domain>
#more such lines......
#ENDMOTIFS
#ENDQUERY	<query-id>
#more query sections..
#ENDSESSION	<session-ordinal>
#more session sections..
#ENDDATA
#=====================================================================
DATA
SESSION	1	RPSBLAST 2.7.1+	/blast/db/blast/cdd	BLOSUM62	0.01
QUERY	Query_1	Peptide	430	gi|27085270|gb|AAN85445.1| SMAD8 protein [Mus musculus]
DOMAINS
1	Query_1	Specific	199814	13	136	4.19971e-86	260.898	cd10490	MH1_SMAD_1_5_9	-	260161
1	Query_1	Superfamily	260162	230	430	1.84736e-150	428.529	cl00056	MH2	-	-
ENDDOMAINS
SITES
1	Query_1	Specific	DNA binding site	K36,S40,S74,L75,R78,L79,Q80,V81,S82,K85,H104	11	11	199814
1	Query_1	Specific	Zn binding site	C68,C113,C125,H130	4	4	199814
1	Query_1	Generic	trimer interface	Y242,E250,T251,Q253,F264,T265,D266,P267,S268,R273,G277,L278,L279,S280,V282,E289,R293,H294,G296,Y369,T372,T395,E402,H404,H406	25	25	199819
ENDSITES
ENDQUERY
ENDSESSION	1
SESSION	2	RPSTBLASTN 2.7.1+	/blast/db/blast/cdd	BLOSUM62	0.01
QUERY	Query_1	Nucleotide	30180	gi|342851972|Borrelia afzelii PKo plasmid lp32-10, complete sequence
DOMAINS
2	Query_1[1]	Specific	251280	16624	17493	6.0196e-138	438.343	pfam02414	Borrelia_orfA	-	45
2	Query_1[1]	Specific	251598	17812	18234	8.70594e-72	241.506	pfam02890	DUF226	-	45
2	Query_1[1]	Specific	251280	11503	12075	3.14064e-70	243.047	pfam02414	Borrelia_orfA	N	45
2	Query_1[1]	Specific	250787	19276	19530	2.61136e-43	157.532	pfam01672	Plasmid_parti	-	45
2	Query_1[1]	Specific	250787	1846	2118	1.14753e-30	120.939	pfam01672	Plasmid_parti	-	45
2	Query_1[1]	Non-Specific	114930	13324	13443	6.37628e-11	64.3142	pfam06238	Borrelia_lipo_2	N	45
2	Query_1[1]	Non-Specific	251598	28285	28458	5.49995e-10	62.003	pfam02890	DUF226	N	45
2	Query_1[1]	Superfamily	265491	5188	5487	7.8357e-10	60.8474	cl14958	HSDR_N	-	-
2	Query_1[1]	Non-Specific	251280	2527	3081	7.9483e-10	63.1586	pfam02414	Borrelia_orfA	N	45
2	Query_1[1]	Non-Specific	251598	664	753	1.69066e-09	60.4622	pfam02890	DUF226	N	45
2	Query_1[1]	Non-Specific	115671	12772	12897	1.47407e-07	53.5286	pfam07032	DUF1322	C	45
2	Query_1[1]	Superfamily	272166	23854	24534	2.07183e-06	51.6026	cl21746	EpsG	-	-
2	Query_1[1]	Non-Specific	251280	19072	19767	8.60803e-05	47.7506	pfam02414	Borrelia_orfA	-	45
2	Query_1[1]	Non-Specific	226114	6079	6273	0.000112918	45.0542	COG3586	COG3586	N	45
2	Query_1[1]	Superfamily	272166	10507	10767	0.00105216	43.5134	cl21746	EpsG	N	-
2	Query_1[1]	Non-Specific	250680	4720	4908	0.0065359	39.2762	pfam01519	DUF16	N	45
2	Query_1[1]	Non-Specific	251280	14236	14874	0.00734381	41.5874	pfam02414	Borrelia_orfA	-	45
2	Query_1[1]	Multidom	224396	19069	20253	5.10176e-10	64.6994	COG1479	COG1479	-	45
2	Query_1[1]	Multidom	224396	2551	3045	0.00197166	43.8986	COG1479	COG1479	N	45
2	Query_1[2]	Specific	238997	1088	1300	5.6938e-09	62.7734	cd02042	ParA	N	271875
2	Query_1[2]	Specific	115964	12329	12748	7.71117e-61	234.958	pfam07341	DUF1473	-	45
2	Query_1[2]	Specific	251280	11207	11920	7.7439e-50	198.364	pfam02414	Borrelia_orfA	C	45
2	Query_1[2]	Specific	251666	13499	13789	5.87401e-44	179.104	pfam02999	Borrelia_orfD	-	45
2	Query_1[2]	Superfamily	271929	23390	23665	1.29102e-36	154.451	cl21509	ApoLp-III_like	N	-
2	Query_1[2]	Non-Specific	115940	12071	12298	5.44539e-28	125.946	pfam07316	DUF1463	N	45
2	Query_1[2]	Superfamily	271875	2	193	8.34729e-17	88.9669	cl21455	ABC_ATPase	NC	-
2	Query_1[2]	Specific	259034	26684	27478	3.42711e-10	67.0106	pfam14897	EpsG	-	272166
2	Query_1[2]	Non-Specific	251280	9587	10507	9.60884e-07	55.4546	pfam02414	Borrelia_orfA	-	45
2	Query_1[2]	Superfamily	271898	28811	29050	1.72577e-06	54.6842	cl21478	Mt_ATP-synt_B	N	-
2	Query_1[2]	Non-Specific	255394	8819	9874	0.000107423	48.521	pfam09491	RE_AlwI	-	45
2	Query_1[2]	Non-Specific	165315	7280	8293	0.000263094	47.3654	PHA03016	PHA03016	C	45
2	Query_1[2]	Specific	259034	6230	6739	0.000434934	46.595	pfam14897	EpsG	N	272166
2	Query_1[2]	Non-Specific	254343	22463	22690	0.00318818	43.8986	pfam07667	DUF1600	-	45
2	Query_1[2]	Multidom	224396	1268	2461	0.00256734	44.2838	COG1479	COG1479	-	45
2	Query_1[3]	Specific	238997	18237	18350	3.05374e-09	63.929	cd02042	ParA	C	271875
2	Query_1[3]	Specific	238997	18585	18761	3.53342e-06	53.5286	cd02042	ParA	N	271875
2	Query_1[3]	Specific	238997	756	872	3.56046e-05	50.447	cd02042	ParA	C	271875
2	Query_1[3]	Superfamily	271929	22950	23405	7.61085e-74	278.1	cl21509	ApoLp-III_like	C	-
2	Query_1[3]	Specific	251859	14796	15101	2.51807e-34	147.132	pfam03304	Mlp	-	45
2	Query_1[3]	Specific	114930	13113	13331	3.6503e-32	139.813	pfam06238	Borrelia_lipo_2	C	45
2	Query_1[3]	Specific	251598	28044	28352	4.54523e-32	139.428	pfam02890	DUF226	C	45
2	Query_1[3]	Non-Specific	253346	25707	26294	4.30364e-19	96.2857	pfam05714	Borrelia_lipo_1	-	45
2	Query_1[3]	Superfamily	271875	186	398	3.20379e-16	87.0409	cl21455	ABC_ATPase	N	-
2	Query_1[3]	Specific	259034	15258	16181	1.16935e-09	65.0846	pfam14897	EpsG	-	272166
2	Query_1[3]	Non-Specific	115671	12891	13001	3.26518e-09	63.5438	pfam07032	DUF1322	N	45
2	Query_1[3]	Non-Specific	251280	2688	3440	1.27865e-05	51.6026	pfam02414	Borrelia_orfA	-	45
2	Query_1[3]	Non-Specific	165315	3090	4169	0.000204623	47.7506	PHA03016	PHA03016	-	45
2	Query_1[3]	Specific	214445	7947	8393	0.000774987	45.8246	MTH00166	ND6	-	271904
2	Query_1[3]	Superfamily	272166	22341	22871	0.0021348	44.2838	cl21746	EpsG	N	-
2	Query_1[3]	Superfamily	272166	21237	21503	0.0033904	43.8986	cl21746	EpsG	N	-
2	Query_1[3]	Superfamily	271981	27096	27542	0.00351032	43.5134	cl21561	7tm_4	C	-
2	Query_1[3]	Superfamily	272166	2118	2525	0.00431283	43.5134	cl21746	EpsG	NC	-
2	Query_1[3]	Multidom	224113	18228	18956	1.9077e-52	207.223	COG1192	Soj	-	45
2	Query_1[3]	Multidom	224396	3153	4286	1.05473e-05	51.9878	COG1479	COG1479	-	45
2	Query_1[3]	Multidom	233080	186	425	9.75408e-19	95.1301	TIGR00665	DnaB	N	45
2	Query_1[3]	Multidom	224113	747	1037	6.54658e-18	92.4337	COG1192	Soj	C	45
2	Query_1[3]	Multidom	234229	2682	3614	9.29802e-05	48.9062	TIGR03490	Mycoplas_LppA	C	45
2	Query_1[-1]	Specific	259034	16834	17766	5.48379e-16	86.2705	pfam14897	EpsG	-	272166
2	Query_1[-1]	Specific	259034	1708	2721	6.3267e-12	72.7886	pfam14897	EpsG	-	272166
2	Query_1[-1]	Specific	259034	11143	12060	1.43775e-11	71.633	pfam14897	EpsG	-	272166
2	Query_1[-1]	Specific	259034	1237	2169	4.78903e-11	69.707	pfam14897	EpsG	-	272166
2	Query_1[-1]	Specific	259034	4948	5904	6.88935e-11	69.3218	pfam14897	EpsG	-	272166
2	Query_1[-1]	Specific	259034	19141	20061	2.92585e-09	63.929	pfam14897	EpsG	-	272166
2	Query_1[-1]	Superfamily	268042	16477	17055	1.54365e-08	61.6178	cl19689	Oxidored_q1	C	-
2	Query_1[-1]	Specific	259034	14056	15099	6.28352e-08	59.3066	pfam14897	EpsG	-	272166
2	Query_1[-1]	Specific	259034	7075	7911	1.24159e-05	51.6026	pfam14897	EpsG	-	272166
2	Query_1[-1]	Specific	259034	4357	5307	4.14665e-05	50.0618	pfam14897	EpsG	-	272166
2	Query_1[-1]	Specific	259034	418	753	0.000131278	48.521	pfam14897	EpsG	NC	272166
2	Query_1[-1]	Specific	259034	11827	12375	0.000214713	47.7506	pfam14897	EpsG	N	272166
2	Query_1[-1]	Superfamily	268042	2515	3084	0.00139548	45.0542	cl19689	Oxidored_q1	N	-
2	Query_1[-1]	Specific	224196	9919	10461	0.0016561	44.669	COG1277	NosY	N	271894
2	Query_1[-1]	Superfamily	272166	13687	14139	0.00528465	43.1282	cl21746	EpsG	N	-
2	Query_1[-1]	Specific	224196	28582	29196	0.00701699	42.743	COG1277	NosY	C	271894
2	Query_1[-1]	Multidom	177158	16618	17994	1.32692e-11	71.633	MTH00095	ND5	-	45
2	Query_1[-2]	Specific	252574	15600	16097	6.58907e-31	135.576	pfam04404	ERF	-	45
2	Query_1[-2]	Specific	259034	2433	3449	3.40237e-20	100.138	pfam14897	EpsG	-	272166
2	Query_1[-2]	Superfamily	271437	6582	7034	5.74264e-18	92.8189	cl00213	DNA_BRE_C	-	-
2	Query_1[-2]	Specific	259034	3150	4181	1.70967e-17	91.2781	pfam14897	EpsG	-	272166
2	Query_1[-2]	Specific	259034	14103	15083	1.19188e-10	68.5514	pfam14897	EpsG	-	272166
2	Query_1[-2]	Non-Specific	177131	16935	17807	3.49055e-09	63.5438	MTH00059	ND2	-	45
2	Query_1[-2]	Superfamily	272166	23262	23762	4.64315e-07	56.6102	cl21746	EpsG	-	-
2	Query_1[-2]	Non-Specific	258475	19716	20594	3.00969e-06	53.9138	pfam14296	O-ag_pol_Wzy	C	45
2	Query_1[-2]	Specific	259034	11562	12074	4.96021e-05	49.6766	pfam14897	EpsG	N	272166
2	Query_1[-2]	Non-Specific	258475	1365	2786	0.000230787	47.7506	pfam14296	O-ag_pol_Wzy	-	45
2	Query_1[-2]	Specific	224196	13686	14462	0.000365584	46.9802	COG1277	NosY	-	271894
2	Query_1[-2]	Non-Specific	251280	15249	15647	0.000802318	45.8246	pfam02414	Borrelia_orfA	NC	45
2	Query_1[-2]	Superfamily	268329	6942	7256	0.00128791	45.0542	cl19976	7tm_7	N	-
2	Query_1[-2]	Non-Specific	255908	5496	6002	0.00217536	44.2838	pfam10325	7TM_GPCR_Srz	C	45
2	Query_1[-2]	Superfamily	272166	24501	25472	0.00323109	43.8986	cl21746	EpsG	-	-
2	Query_1[-2]	Specific	257675	18849	19202	0.00348228	43.5134	pfam13346	ABC2_membrane_5	N	271894
2	Query_1[-2]	Non-Specific	254938	4668	4838	0.00467304	43.1282	pfam08634	Pet127	C	45
2	Query_1[-2]	Superfamily	271981	3834	4589	0.00705389	42.743	cl21561	7tm_4	-	-
2	Query_1[-2]	Superfamily	272166	20268	20975	0.0085514	42.3578	cl21746	EpsG	N	-
2	Query_1[-2]	Multidom	177157	2418	3455	1.25459e-11	71.633	MTH00094	ND4	-	45
2	Query_1[-2]	Multidom	233789	6630	7079	3.05713e-12	73.559	TIGR02225	recomb_XerD	N	45
2	Query_1[-2]	Multidom	234043	19629	20063	0.000122134	48.521	TIGR02871	spore_ylbJ	N	45
2	Query_1[-3]	Specific	259034	8756	9649	2.00834e-11	70.8626	pfam14897	EpsG	-	272166
2	Query_1[-3]	Specific	259034	9863	10546	1.90816e-10	67.781	pfam14897	EpsG	-	272166
2	Query_1[-3]	Specific	259034	2414	3292	1.94348e-09	64.3142	pfam14897	EpsG	C	272166
2	Query_1[-3]	Specific	259034	16976	17938	2.43849e-08	60.8474	pfam14897	EpsG	-	272166
2	Query_1[-3]	Specific	259034	572	1597	4.1738e-08	60.077	pfam14897	EpsG	-	272166
2	Query_1[-3]	Specific	259034	22880	23857	8.63418e-07	55.4546	pfam14897	EpsG	-	272166
2	Query_1[-3]	Specific	259034	1937	2830	1.38229e-06	55.0694	pfam14897	EpsG	-	272166
2	Query_1[-3]	Superfamily	272166	7235	7831	4.64133e-06	53.1434	cl21746	EpsG	-	-
2	Query_1[-3]	Specific	259034	18086	18964	1.02965e-05	51.9878	pfam14897	EpsG	-	272166
2	Query_1[-3]	Specific	259034	269	772	3.32135e-05	50.447	pfam14897	EpsG	N	272166
2	Query_1[-3]	Non-Specific	255908	11594	11911	5.9978e-05	49.6766	pfam10325	7TM_GPCR_Srz	C	45
2	Query_1[-3]	Non-Specific	258475	2822	4171	0.000109453	48.521	pfam14296	O-ag_pol_Wzy	-	45
2	Query_1[-3]	Superfamily	271677	7460	8404	0.00019345	47.7506	cl09928	Molybdopterin-Binding	N	-
2	Query_1[-3]	Specific	259034	28613	29173	0.000235147	47.3654	pfam14897	EpsG	C	272166
2	Query_1[-3]	Specific	259034	11357	12313	0.000404641	46.595	pfam14897	EpsG	-	272166
2	Query_1[-3]	Superfamily	272166	19073	19435	0.00068706	45.8246	cl21746	EpsG	N	-
2	Query_1[-3]	Non-Specific	258475	20189	20956	0.000787441	45.8246	pfam14296	O-ag_pol_Wzy	C	45
2	Query_1[-3]	Superfamily	272166	4046	5119	0.00207293	44.2838	cl21746	EpsG	-	-
2	Query_1[-3]	Non-Specific	258475	5624	6067	0.00354807	43.5134	pfam14296	O-ag_pol_Wzy	C	45
ENDDOMAINS
SITES
2	Query_1[2]	Specific	Magnesium ion binding site	D1121-1123	2	1	238997
2	Query_1[2]	Generic	ATP binding site	D53-55,Y56-58	9	2	238540
2	Query_1[2]	Generic	Walker B motif	I41-43,I44-46,F47-49,I50-52,D53-55	5	5	238540
2	Query_1[3]	Specific	P-loop	K18258-18260,G18261-18263,G18264-18266,V18267-18269,G18270-18272,K18273-18275,S18276-18278	7	7	238997
2	Query_1[3]	Specific	Magnesium ion binding site	S18276-18278	2	1	238997
2	Query_1[3]	Specific	Magnesium ion binding site	D18603-18605	2	1	238997
2	Query_1[3]	Specific	P-loop	K777-779,G780-782,G783-785,V786-788,S789-791,K792-794,S795-797	7	7	238997
2	Query_1[3]	Specific	Magnesium ion binding site	S795-797	2	1	238997
2	Query_1[-2]	Generic	active site	K6929-6927,H6656-6654,R6647-6645,H6596-6594	5	4	271175
2	Query_1[-2]	Generic	DNA binding site	R6971-6969,K6929-6927,I6926-6924,T6719-6717,Q6713-6711,H6656-6654	7	6	271175
2	Query_1[-3]	Generic	molybdopterin cofactor binding site 	K8215-8213,S8209-8207,I8206-8204,N8203-8201,I8143-8141,N8140-8138,N8134-8132,N8125-8123,E8122-8120,I8038-8036,I8035-8033,L7999-7997,I7627-7625,*7624-7622,K7621-7619,F7582-7580,I7579-7577,*7576-7574,F7567-7565,A7501-7499,K7498-7496,S7483-7481	25	22	238218
ENDSITES
ENDQUERY
ENDSESSION	2
ENDDATA

#Total Blastout object processed	2


================================================================================

  Aron Marchler-Bauer, Shennan Lu, Renata Geer, revised 11 September 2016