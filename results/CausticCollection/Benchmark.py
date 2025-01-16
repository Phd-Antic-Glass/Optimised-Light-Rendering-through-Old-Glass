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
timeout = 10
spp = 10000

# path to the scene XML file
filename = './Caustic_collection.xml'

# Add the scene directory to the FileResolver's search path
Thread.thread().file_resolver().append(os.path.dirname(filename))

## rendering
scene_dict = s.build_heightfield_scene()
output_folder = "CausticCollection/"

# Fermat's principle (ours)
s.FNEE_integrator["timeout"] = float(timeout)
scene_dict["integrator"] = load_dict(s.FNEE_integrator)
name = output_folder+"Caustics_collection_" +str(int(timeout))+"_path_fnee"
u.render(load_dict(scene_dict), name)

# SMS
s.SMS_integrator["timeout"] = timeout
scene_dict["integrator"] = load_dict(s.SMS_integrator)
s.SMS_integrator["timeout"] = float(timeout)
name = output_folder+"Caustics_collection_" +str(timeout)+"_path_sms_ms"
u.render(load_dict(scene_dict), name)

# SMS with same initialisation strategy than our Fermat based solver (for fair comparison of the root finder)
s.FNEE_integrator_SMS["timeout"] = float(timeout)
scene_dict["integrator"] = load_dict(s.FNEE_integrator_SMS)
name = output_folder+"Caustics_collection_" +str(int(timeout))+"_FNEE_SMS"
u.render(load_dict(scene_dict), name)
