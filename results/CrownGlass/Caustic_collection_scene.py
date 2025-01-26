import mitsuba
mitsuba.set_variant("scalar_rgb")
from mitsuba.core import ScalarTransform4f


##############
# INTEGRATOR #
##############
FNEE_integrator = {
    "type": "path_fnee",
    "max_depth": 8,
    "samples_per_pass": 1,
    "timeout": 2,
    "biased": False,
    "twostage": False,
    "halfvector_constraints": True,
    "mnee_init": False,
    "max_iterations": 20,
    "solver_threshold": 0.00001,
    "uniqueness_threshold": 0.000001,
    "max_trials": 1000,
    "bounces": 2,
    "use_SMS": False,
    "crop_caustic": False
}

SMS_integrator = {
    "type": "path_sms_ms",
    "max_depth": 8,
    "samples_per_pass": 1,
    "timeout": 2,
    "biased": False,
    "twostage": False,
    "halfvector_constraints": True,
    "mnee_init": False,
    "max_iterations": 20,
    "solver_threshold": 0.00001,
    "uniqueness_threshold": 0.000001,
    "max_trials": 1000,
    "bounces": 2,
}

SMS_AG_integrator = {
    "type": "path_sms_ms_AG",
    "max_depth": 8,
    "samples_per_pass": 1,
    "timeout": 2,
    "biased": False,
    "twostage": False,
    "halfvector_constraints": True,
    "mnee_init": False,
    "max_iterations": 20,
    "solver_threshold": 0.00001,
    "uniqueness_threshold": 0.000001,
    "max_trials": 1000,
    "bounces": 2,
}

path_integrator = {
    "type": "path",
    "max_depth": 8,
    "samples_per_pass": 1,
    "timeout": 2,
}


###############
# BSDF_WINDOW #
###############
bsdf_window = {
        "type": "dielectric",
        "int_ior": 1.54
        }

################
# BSDF_DIFFUSE #
################
bsdf_diffuse = {
        "type": "diffuse",
        "reflectance": {
            "type": "rgb",
            "value": [0.4, 0.4, 0.4],
            }
        }

###################
# HEIGHFIELD_DICT #
###################
Texture_heightfield = {
        "type": "bitmap",
        "id": "id_texture_simplex_1",
        "filename": "textures/Crown.exr",
        "raw": True
        }

##############
# PANEL MESH #
##############

panel_mesh = {
        "type": "obj",
        "filename": "meshes/caustic_panel_thin.obj",
        "bsdf": bsdf_window,
        "caustic_caster_single": True,
        "caustic_caster_multi" : True
        }

panel_manchon_mesh = {
        "type": "obj",
        "filename": "meshes/Panel_manchon.obj",
        "bsdf": bsdf_window,
        "caustic_caster_multi" : True,
        #        "caustic_caster_single": True,
        }

##############
# SCENE_DICT #
##############
scene_dict = {
    "type": "scene",

    "Floor": {
        "type": "obj",
        "filename": "meshes/floor.obj",
        "bsdf": bsdf_diffuse,
        "caustic_receiver": True,
    },

    "Frame": {
        "type": "obj",
        "filename": "meshes/frame.obj",
        "bsdf": bsdf_diffuse,
        "caustic_receiver": True,
    },

    "Lead": {
        "type": "obj",
        "filename": "meshes/lead.obj",
        "bsdf": bsdf_diffuse,
        "caustic_receiver": True,
    },

    "Cube": {
        "type": "obj",
        "filename": "meshes/cube.obj",
        "bsdf": bsdf_diffuse,
        "caustic_receiver": True,
    },

    "light_1": {
        "type": "sphere",
        "radius": 0.05,
        "emitter": {
            "type": "area",
            "radiance": 50000.0,
            "caustic_emitter_multi": True,
            "caustic_emitter_single": True
        },
        "to_world": ScalarTransform4f.translate([0, 2, 5.0])
    },

    "camera_1": {
        "type": "perspective",
        "fov": 40,
        "fov_axis": "smaller",
        "to_world": ScalarTransform4f.look_at(origin=[5.4398980140686035, 0.26405173540115356, 2.311628818511963],
                                              target=[
            4.942975044250488 , -0.5388131737709045, 1.9822726249694824],
            up=[-0.17331193387508392, -0.2800687253475189, 0.9442057609558105]),
        "film": {
            "type": "hdrfilm",
            "rfilter": {"type": "box"},
            "width": 1024,
            "height": 1024,
            "crop_offset_x": 256,
            "crop_offset_y": 750,
            "crop_height": 64,
            "crop_width": 64,
        },
        "sampler": {
            "type": "independent",
            "sample_count": 100000,
        }
    },

    # "Texture_cive_1": Texture_cive_1,

    "integrator": path_integrator
}


def add_window_to_scene(scene_dict, window_name, id, to_world_1, to_world_2, H1, H2, texture_id_1, texture_id_2, L=1, W=1):
    Heighfield_dict = {
        "type": "heightfield",
        "to_world": ScalarTransform4f.translate([2.9, -3.26, 1.0]) * ScalarTransform4f.rotate([0, 0, 1], -90.0),
        "L": L,
        "W": W,
        "H": 0.00,
        "flip_normals": False,
        "caustic_bouncer": True,
        "caustic_caster_single": False,
        "caustic_caster_multi": True,
        "double_refraction_caster": True,
        "self_id": 1,
        "pair_id": 2,
        "partition_number": 1.0,
        "bsdf": bsdf_window,
        "heightmap": Texture_heightfield
    }

    Heighfield_dict["H"] = H1
    Heighfield_dict["to_world"] = to_world_1
    Heighfield_dict["flip_normals"] = False
    Heighfield_dict["self_id"] = id * 2
    Heighfield_dict["pair_id"] = id * 2 + 1
    Heighfield_dict["heightmap"] = texture_id_1
    Heighfield_dict["caustic_bouncer"] = False
    Heighfield_dict["caustic_caster_multi"] = True
    Heighfield_dict["caustic_caster_single"] = True
    scene_dict[window_name + "_1"] = Heighfield_dict.copy()

    Heighfield_dict["H"] = H2
    Heighfield_dict["self_id"] = id * 2 + 1
    Heighfield_dict["pair_id"] = id * 2
    Heighfield_dict["to_world"] = to_world_2
    Heighfield_dict["flip_normals"] = True
    Heighfield_dict["heightmap"] = texture_id_2
    Heighfield_dict["caustic_bouncer"] = True
    Heighfield_dict["caustic_caster_single"] = False
    Heighfield_dict["caustic_caster_multi"] = False

    scene_dict[window_name + "_2"] = Heighfield_dict.copy()

def build_heightfield_scene():
    Texture_heightfield = {
        "type": "bitmap",
        "id": "id_texture_simplex_1",
        "filename": "textures/Crown.exr",
        "raw": True
        }
    scene_heightfield_dict = scene_dict.copy()
    e = 0.04
    T1 = ScalarTransform4f.translate(
                [1.9 , -3.3, 0.001 ]) * ScalarTransform4f.rotate([0, 0, 1], -90.0)
    T2 = ScalarTransform4f.translate(
                [1.9 , -3.3 + e, 0.001 ]) * ScalarTransform4f.rotate([0, 0, 1], -90.0)
    add_window_to_scene(scene_heightfield_dict, "H1", 1,
            T1, T2, 0.00, 0.04, Texture_heightfield, Texture_heightfield, L=2, W=2)
    return scene_heightfield_dict

