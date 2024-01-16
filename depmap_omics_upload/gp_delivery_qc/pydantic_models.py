from typing import Literal

from pydantic import BaseModel


class Workspace(BaseModel):
    source: str
    name: str
    file_type: str


class GumboClientConfig(BaseModel):
    env: str = Literal["staging", "production"]
    username: str
    pr_table_name: str
    pr_table_index: str
    seq_table_name: str
    seq_table_index: str


class GPDeliveryQCConfig(BaseModel):
    workspace: Workspace
    id_col_name: str
    profile_id_col_name: str
    data_type: str = Literal["wgs", "rna"]
    bam_col_name: str
    col_name_map: dict
    min_file_size: int
    gumbo_client_config: GumboClientConfig
