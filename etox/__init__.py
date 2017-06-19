import os
import sys
import json

from eTOX_ALLIES import __rootpath__

json_settings = os.path.join(__rootpath__, 'data/settings.json')
settings = {}
with open(json_settings) as st:
    settings = json.load(st)