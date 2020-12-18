#!/usr/bin/env python

"""
The main class object for executing KMC kmer funcs
KMC GitHub: https://github.com/refresh-bio/KMC
KMC Paper: https://academic.oup.com/bioinformatics/article/33/17/2759/3796399
KMC Suppl: https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/33/17/10.1093_bioinformatics_btx304/3/btx304_supplementary.pdf?Expires=1611337056&Signature=4iJG6giNcZdiDuHwljf-SbELTll74FtIj3YIFvfESeZC~m39EZPJdSXfqAJStvCr5SmH9lHGRCdJGHBLseX~ZunAgFZBFFHikmODBI14Kq84ctkQMihTvBzU1rme~S6MpXcC1Erxavl~ckAEnE7jfwIRJbtm4bSkTk-sEcZKIHqR3H0SwhdN0zMmhqMFkwn~jvNo5Rd~yPFwq8aXtE2CBrMORgVUsu~ACFKnl7sWB2FLtsZ2zp~ENuVz28mYwZWkyFXnUlRq2sKHenjWdw4BChI~QDf5EULM2oXgx4dSDlLTaIaJjsZBGl8tKQGw5Ohz48YEbqO82pn14AyeJkaRsQ__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA
KMC Docs: https://github.com/refresh-bio/KMC/blob/master/API.pdf
KMC-tools Docs: https://github.com/refresh-bio/KMC/blob/master/kmc_tools.pdf
"""


import subprocess



class Kcount:
    """
    Calls the kmc counting functions to create a database.
    """
    def __init__(self, files, outdir, kmersize):
        self.files = files
        self.outdir = outdir
        self.kmersize = kmersize


    def call_kmc_on_files(self):
        """
        Calls kmc count on all files
        """
        for filename in self.files:
            outname = filename.split(".fastq")[0]
            self.call_kmc_count(filename, outname)


    def call_kmc_count(self, filename, outname):
        """
        takes a fastq file and calls KMC count on it.
        """
        # create command: 'kmc -k17 inputfastq outname outdir'
        cmd = [
            "kmc", "-k{}".format(self.kmersize), 
            filename,
            outname,
            self.outdir,
        ]

        # call subprocess on the command
        output = subprocess.check_output(cmd)

        # check whether it was successful or error
        # ...output



fileslist = [
    "/home/deren/Documents/ipyrad/isolation/reftest_fastqs/1A_0_R1_.fastq.gz",
    "/home/deren/Documents/ipyrad/isolation/reftest_fastqs/1B_0_R1_.fastq.gz",
]

counter = Kcount(fileslist, "./", 17)
counter.call_kmc_on_files()
