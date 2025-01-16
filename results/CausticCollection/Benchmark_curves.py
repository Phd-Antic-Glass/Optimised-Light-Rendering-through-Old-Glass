import os
import sys
import Caustic_collection_scene as s

##################
# MITSUBA IMPORT #
##################
import mitsuba
mitsuba.set_variant("scalar_rgb")

from mitsuba.core import Thread
from mitsuba.core.xml import load_file
from mitsuba.core.xml import load_dict

import numpy as np
from matplotlib import pyplot as plt
import Utility as u


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
timeout = np.arange(10, 110, 10, dtype=np.float32)
# timeout = np.arange(1, 11, 1, dtype=np.float32)
# timeout = [10]
print(timeout)
spp = 10000

# Absolute or relative path to the scene XML file
filename = './Caustic_collection.xml'

# Add the scene directory to the FileResolver's search path
Thread.thread().file_resolver().append(os.path.dirname(filename))
############
# # RENDER #
############

scene_dict = s.build_heightfield_scene()
output_folder = "CausticCollection/"
mode = input("recompute images ? (y/n)")
if(mode == "y"):
    for i in range(0,len(timeout)):
        # figure Caustic collection
        # FNEE
        s.FNEE_integrator["timeout"] = float(timeout[i])
        scene_dict["integrator"] = load_dict(s.FNEE_integrator)
        name = output_folder+"Caustics_collection_" +str(int(timeout[i]))+"_path_fnee"
        u.render(load_dict(scene_dict), name)

        ## SMS
        #  s.SMS_integrator["timeout"] = timeout
        #  scene_dict["integrator"] = load_dict(s.SMS_integrator)
        #  name = "Caustics_collection_" +str(int(timeout))+"_path_sms_ms"
        #  u.render(load_dict(scene_dict), name)

        # FNEE_SMS
        s.FNEE_integrator_SMS["timeout"] = float(timeout[i])
        scene_dict["integrator"] = load_dict(s.FNEE_integrator_SMS)
        name = output_folder+"Caustics_collection_" +str(int(timeout[i]))+"_FNEE_SMS"
        u.render(load_dict(scene_dict), name)
        
X = timeout
Y_fnee = []
Y_sms = []
folder=output_folder
img_ref = np.asarray(u.exr_to_srgb(folder + "ref/Caustics_collection_1Mspp_path.exr"))
# compute mse for each images
for i in range(0,len(X)):
    img_sms = np.asarray(u.exr_to_srgb(output_folder+"Caustics_collection_" +str(int(timeout[i]))+"_FNEE_SMS.exr"))
    img_fnee = np.asarray(u.exr_to_srgb(output_folder+"Caustics_collection_" +str(int(timeout[i]))+"_path_fnee.exr"))

    img_sms = np.clip(img_sms, 0 , np.max(img_ref)*1.0)
    img_fnee = np.clip(img_fnee, 0 , np.max(img_ref)*1.0)

    # img_fnee = np.asarray(u.exr_to_srgb(folder + "/Caustic_panel_" +str(int(timeout[i]))+"_path_sms_ms.exr"))
    # img_sms = np.asarray(u.exr_to_srgb(folder + "/Caustic_panel_" +str(int(timeout[i]))+"_path_fnee.exr"))
    # print(img_fnee)

    rmse_fnee = np.mean(np.square(img_fnee - img_ref))
    Y_fnee.append(rmse_fnee)
    rmse_sms = np.mean(np.square(img_sms - img_ref))
    Y_sms.append(rmse_sms)

print(Y_fnee)
fig, ax = plt.subplots()
ax.plot(X,Y_fnee,"b", label="Fermat")
ax.plot(X,Y_sms,"g", label="SMS")
ax.legend(loc="upper right")
# ax.set_xlim(0,0)
ax.set_ylim(0,)

ax.set_xlabel("time (s)")
ax.set_ylabel("MSE")
# ax.set_
plt.savefig(output_folder+"convergence_plot.pdf")
#  ax.set_title(output_folder)
plt.show()


