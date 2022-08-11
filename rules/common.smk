import os 
import pandas as pd

class RunInfo:
    def __init__(self, config):
        self.metadata = pd.read_csv(config['metadata'], sep='\t')
        self.metadata = self.metadata[self.metadata['isabl_aliquot_id'] == config['aliquot_id']]
        aliquot_ids = set(self.metadata.query('sample_category == "TUMOR"')['isabl_aliquot_id'].unique())
        self.metadata = self.metadata[self.metadata['isabl_application'].isin(['WGS-SOMATICCALLING', 'WGS-REMIXT-POSTPROCESS'])]
        for app in self.metadata['isabl_application'].unique():
            app_samples = set(self.metadata.query('isabl_application == @app')['isabl_aliquot_id'].unique())
            if len(aliquot_ids.difference(app_samples)) > 0:
                raise Exception(f'missing samples {aliquot_ids.difference(app_samples)}')
        self.paths = self.metadata.set_index(['isabl_aliquot_id', 'isabl_application', 'result_type'])['result_filepath'].to_dict()
        self.aliquot_ids = list(sorted(self.metadata['isabl_aliquot_id'].unique()))

runinfo = RunInfo(config)
