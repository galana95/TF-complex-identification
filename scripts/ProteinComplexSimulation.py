import pandas as pd

class ProteinComplexes():
    def __init__(self,simulation_data):
        self.simulation_data = simulation_data
        self.COMPLEX_THRESHOLDS={1:100,2:100,3:66,4:75,5:60,6:50,7:42,8:37,9:33,10:30,11:27,12:25}
        self.proteincomplex_data=simulation_data['List of proteins']
