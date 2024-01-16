import gumbo_client
from gumbo_client.utils import NameMappingUtils

from depmap_omics_upload.gp_delivery_qc.pydantic_models import GumboClientConfig


class GumboClientWrapper:
    def __init__(self, gumbo_client_config: GumboClientConfig):
        self.pr_table_name = gumbo_client_config.pr_table_name
        self.pr_table_index = gumbo_client_config.pr_table_index
        self.seq_table_name = gumbo_client_config.seq_table_name
        self.seq_table_index = gumbo_client_config.seq_table_index

        if gumbo_client_config.env == "staging":
            # hard-coded for now, waiting for independent staging envs to be enabled
            self.client = gumbo_client.Client(
                config_dir="~/.config/gumbo-staging",
                username=gumbo_client_config.username,
            )
        else:
            self.client = gumbo_client.Client(username=gumbo_client_config.username)

        self.mapping_utils = NameMappingUtils()

    def close_gumbo_client(self):
        self.client.close()

    def read_pr_table(self):
        pr_table = self.client.get(self.pr_table_name)
        pr_table_camel = self.mapping_utils.rename_columns(self.pr_table_name, pr_table)
        pr_table_camel = pr_table_camel.set_index(self.pr_table_index)
        return pr_table_camel

    def read_seq_table(self):
        seq_table = self.client.get(self.seq_table_name)
        seq_table_camel = self.mapping_utils.rename_columns(
            self.seq_table_name, seq_table
        )
        seq_table_camel = seq_table_camel.set_index(self.seq_table_index)
        return seq_table_camel
