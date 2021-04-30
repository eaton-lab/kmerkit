## TODO


from loguru import logger
from kmerkit.kschema import Project, KtreeBase, KtreeParams, KtreeData


class Ktree:
    def __init__(self, json_file, min_z_score):

        # store Project as a dict        
        self.json_file = json_file
        self.project = Project.parse_file(self.json_file).dict()

        # store user input as pydantic object and as a dict
        self._params = KtreeParams(min_z_score=min_z_score)
        self.params = self._params.dict()


    def run(self):

        # run analysis to get results.
        # ...

        # create a pydantic data object to store outputs
        data = KtreeData(
            data_in="test.file", 
            data_out=['test1.txt', 'text2.txt']
        )

        # create final ktree pydantic object
        # add ktree object to project dict
        self.project["ktree"] = KtreeBase(params=self._params, data=data)

        # opena json file for writing and 
        with open(self.json_file, 'w') as out:
            json_string = Project(**self.project).json(indent=4)
            out.write(json_string)
        logger.info(json_string)



if __name__ == "__main__":


    import kmerkit
    kmerkit.set_loglevel("INFO")
    JSON = "/tmp/test.json"
    ktr = Ktree(JSON, min_z_score=2)
    ktr.run()

