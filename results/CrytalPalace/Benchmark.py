import os
import sys
import palace_scene as s

##################
# MITSUBA IMPORT #
##################
import mitsuba
mitsuba.set_variant("scalar_rgb")

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
resolution = 512
timeout = 100
spp = 20000

# Absolute or relative path to the scene XML file
filename = './Caustic_collection.xml'

# Add the scene directory to the FileResolver's search path
Thread.thread().file_resolver().append(os.path.dirname(filename))
        
############
# # RENDER #
############
def render_scene(crop=""):
    scene_dict = s.build_heightfield_scene()
    scene_dict["camera_1"]["film"]["width"] = resolution
    scene_dict["camera_1"]["film"]["height"] = resolution
    scene_dict["camera_1"]["sampler"]["sample_count"] = spp

    if(crop == "_BLUE"):
        scene_dict["camera_1"]["film"]["width"] = 1024
        scene_dict["camera_1"]["film"]["height"] = 1024
        scene_dict["camera_1"]["film"]["crop_offset_x"] = 554
        scene_dict["camera_1"]["film"]["crop_offset_y"] = 252
        scene_dict["camera_1"]["film"]["crop_height"] = 128
        scene_dict["camera_1"]["film"]["crop_width"] = 128
    if(crop == "_RED"):
        scene_dict["camera_1"]["film"]["width"] = 1024
        scene_dict["camera_1"]["film"]["height"] = 1024
        scene_dict["camera_1"]["film"]["crop_offset_x"] = 413
        scene_dict["camera_1"]["film"]["crop_offset_y"] = 839
        scene_dict["camera_1"]["film"]["crop_height"] = 128
        scene_dict["camera_1"]["film"]["crop_width"] = 128
            #BLUE
#            "crop_offset_x": 554,
#            "crop_offset_y": 252,
#            "crop_height": 128,
#            "crop_width": 128,
            #RED
#            "crop_offset_x": 413,
#            "crop_offset_y": 839,
#            "crop_height": 128,
#            "crop_width": 128,
    #  FNEE
    s.FNEE_integrator["timeout"] = timeout
    s.FNEE_integrator["use_SMS"] = False
    scene_dict["integrator"] = load_dict(s.FNEE_integrator)
    name = "result/palace_" +str(timeout)+"_path_fnee"+crop
    scene = load_dict(scene_dict)
    render(scene, name)

    # FNEE_SMS
    s.FNEE_integrator["timeout"] = timeout
    s.FNEE_integrator["use_SMS"] = True
    scene_dict["integrator"] = load_dict(s.FNEE_integrator)
    name = "result/palace_" +str(timeout)+"_path_fnee_SMS"+crop
    scene = load_dict(scene_dict)
    render(scene, name)

    # SMS 
    s.SMS_integrator["timeout"] = timeout
    scene_dict["integrator"] = load_dict(s.SMS_integrator)
    name = "result/palace_" +str(timeout)+"_path_sms_ms"+crop
    scene = load_dict(scene_dict)
    render(scene, name)

render_scene()
render_scene("_BLUE")
render_scene("_RED")

