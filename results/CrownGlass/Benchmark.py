import os
import sys
import Caustic_collection_scene as s

##################
# MITSUBA IMPORT #
##################
import mitsuba
mitsuba.set_variant("scalar_rgb")

# print("SET VARIANT : gpu_rgb")
# mitsuba.set_variant("gpu_rgb")

from mitsuba.core import Thread
from mitsuba.core.xml import load_file
from mitsuba.core.xml import load_dict

##########
# RENDER #
##########
def render(scene, name, jpg=False):
    # Get the scene's sensor (if many, can pick one by specifying the index)
    sensor = scene.sensors()[0]

    # Call the scene's integrator to render the loaded scene with the desired sensor
    scene.integrator().render(scene, sensor)

    # The rendered data is stored in the film
    film = sensor.film()

    # Write out data as high dynamic range OpenEXR file
    film.set_destination_file(name)
    film.develop()

    if(jpg):
        # Write out a tone-mapped JPG of the same rendering
        from mitsuba.core import Bitmap, Struct
        img = film.bitmap(raw=True).convert(
            Bitmap.PixelFormat.RGB, Struct.Type.UInt8, srgb_gamma=True)
        img.write(name + ".jpg")


#################
# SCENE LOADING #
#################
resolution = 1024
timeout = 1000
spp = 10000

# Absolute or relative path to the scene XML file
filename = './Caustic_collection.xml'

# Add the scene directory to the FileResolver's search path
Thread.thread().file_resolver().append(os.path.dirname(filename))
        
############
# # RENDER #
############
#CPU
scene_dict = s.build_heightfield_scene()

# base precision
s.FNEE_integrator["timeout"] = timeout
s.FNEE_integrator["uniqueness_threshold"] = 1e-5
s.FNEE_integrator["max_trials"] = 10000
scene_dict["integrator"] = load_dict(s.FNEE_integrator)
name = "result/Crown_zoom_" +str(timeout)+"_path_fnee_base_precision"
render(load_dict(scene_dict), name)

# high precision + more pdf trials (much slower)
scene_dict = s.build_heightfield_scene()

s.FNEE_integrator["timeout"] = timeout
s.FNEE_integrator["uniqueness_threshold"] = 5e-7
s.FNEE_integrator["max_trials"] = 1000000
scene_dict["integrator"] = load_dict(s.FNEE_integrator)
name = "result/Crown_zoom_" +str(timeout)+"_path_fnee_high_precision"
render(load_dict(scene_dict), name)

