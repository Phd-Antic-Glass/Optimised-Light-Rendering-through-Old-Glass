import os

import mitsuba
mitsuba.set_variant("scalar_rgb")

import Caustic_collection_scene as s

from mitsuba.core import Thread
from mitsuba.core.xml import load_file
from mitsuba.core.xml import load_dict

import Utility as u


## render parameters
resolution = 1024

# path to the scene XML file
filename = './Caustic_collection.xml'

# Add the scene directory to the FileResolver's search path
Thread.thread().file_resolver().append(os.path.dirname(filename))

## rendering
scene_dict = s.build_gpu_scene()
output_folder = "CausticCollection/"

# Fermat's principle (ours)
name = output_folder+"Caustics_collection_1Mspp_path"
u.render(load_dict(scene_dict), name)
