#!/usr/bin/env python

"""
Initialize a project with input data (fastq/a or SRA/URLs).
This will create the JSON project file and check that the 
paths to the data files are valid. Optionally, if you have not 
already done so, you can implement read-trimming/filtering using
the default settings with 'fastp', where the new trimmed files will
be written to the working dir, and a reference to these files saved
in the project JSON file.

Usage:
------

kmerkit init --name test --workdir /tmp/kmers --trim ...

"""

import os
import sklearn
from loguru import logger
import kmerkit
from kmerkit.utils import KmerkitError
from kmerkit.kschema import Project, Kinit



def init_project(name, workdir, fastq_dict, force=False):
    """
    Initialize a Project JSON file from input data.
    """
    return ProjectInit(name, workdir, fastq_dict, force)


class ProjectInit:
    """
    Initialize a Project JSON file from parsed input data.
    """
    def __init__(self, name, workdir, fastq_dict, force=False):

        self.fastq_dict = fastq_dict

        # expand workdir path
        workdir = os.path.realpath(os.path.expanduser(workdir))

        # get JSON path
        json_path = os.path.join(workdir, name + ".json")

        # bail out if json exists and not force
        if os.path.exists(json_path) and not force:
            msg = f"Project file exists ({json_path}), use force to overwrite."
            logger.error(msg)
            raise KmerkitError(msg)

        # check sample names for bad characters
        self.check_fastq_dict()

        # create workdir if it does not exist
        os.makedirs(workdir, exist_ok=True)

        # inherit parent class to type check against Project JSON schema
        proj = Project(
            name=name,
            workdir=workdir,
            versions={
                'kmerkit': kmerkit.__version__,
                'kmc': '3.1.1',
                'gemma': '0.9.83',
                'sklearn': sklearn.__version__,
            },
            kinit=Kinit(data=fastq_dict, commands={})
        )

        # save project file
        with open(json_path, 'w') as out:
            out.write(proj.json(indent=4))
        logger.info(f"Initialized new kmerkit project: {json_path}")



    def check_fastq_dict(self):
        """
        Expand file paths and check that they exist. Also,
        """
        okeys = list(self.fastq_dict.keys())
        for okey in okeys:
            val = self.fastq_dict[okey]

            # get new key without any strange characters
            newkey = (okey
                .replace("-", "_")
                .replace("@", "_")
                .replace(" ", "_")
            )

            # check type of val
            assert isinstance(val, list), (
                "Filepaths in fastq_dict must be stored as list objects.\n"
                "Example: fdict = {'a': ['a_R1.fastq', 'a_R2.fastq'], 'b'...}"
            )

            # expand each path in val and check exists
            file_list = []
            for path in val:
                fullpath = os.path.realpath(os.path.expanduser(path))
                assert os.path.exists(fullpath), (
                    f"file {fullpath} in fastq_dict cannot be found.")
                file_list.append(fullpath)

            # remove old key (in case of strange characters) and store new.
            del self.fastq_dict[okey]
            self.fastq_dict[newkey] = sorted(file_list)



    # def expand_filenames(self):
    #     """
    #     Allows for selecting multiple input files using wildcard
    #     operators like "./fastqs/*.fastq.gz" to select all fastq.gz
    #     files in the folder fastqs.
    #     """
    #     if isinstance(self.fastq_path, (str, bytes)):
    #         self.files = glob.glob(self.fastq_path)

    #         # raise an exception if no files were found
    #         if not any(self.files):
    #             msg = f"no fastq files found at: {self.files}"
    #             logger.error(msg)
    #             raise KmerkitError(msg)

    #         # sort the input files
    #         self.files = sorted(self.files)

    #     # report on found files
    #     logger.debug(f"found {len(self.files)} input files")



    # def check_samples(self):
    #     """
    #     Gets sample names from the input files and checks that all 
    #     have the same style of suffix (e.g., .fastq.gz).
    #     """
    #     # split file names to keep what comes before 'name_split'
    #     sample_names = [
    #         os.path.basename(i.split(self.name_split)[0]) for i in self.files
    #     ]

    #     # check that all sample_names are unique
    #     if len(set(sample_names)) != len(sample_names):
            
    #         # if not, then check each occurs 2X (PE reads)
    #         if not all([sample_names.count(i) == 2 for i in sample_names]):
    #             raise KmerkitError(
    #                 "Sample names are not unique, or in sets of 2 (PE). "
    #                 "You may need to try a different name_split setting."
    #             )                
    #         logger.debug("detected PE data")

    #     # store dict mapping names and files (or file pairs for PE)
    #     for sname, file in zip(sample_names, self.files):

    #         # names to input fastqs
    #         if sname in self.names_to_infiles:
    #             self.names_to_infiles[sname].append(file)
    #         else:
    #             self.names_to_infiles[sname] = [file]


if __name__ == "__main__":

    import kmerkit

    FILES = "~/Documents/kmerkit/data/amaranths/hybridus_*.fastq.gz"
    FASTQ_DICT = kmerkit.utils.get_fastq_dict_from_path(FILES)

    # init and save project  
    PROJ = init_project(
        name='test', workdir='/tmp', fastq_dict=FASTQ_DICT, force=True)

    # # load from schema with schema type-checking
    # PROJ = Project.parse_file("/tmp/test.json")
    # print(PROJ.json(indent=4))
