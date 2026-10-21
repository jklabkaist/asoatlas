__version__ = "v0.0.2"

import os, yaml
from jklib.config import Config


git_dir = os.path.dirname(os.path.abspath(__file__))

with open("%s/path.yaml" %(git_dir), "r") as ymlfile:
    CONFIG = yaml.safe_load(ymlfile)

_cfg = Config(CONFIG)
DB = _cfg.DB
REFERENCE = _cfg.Reference
DATA = _cfg.Data
