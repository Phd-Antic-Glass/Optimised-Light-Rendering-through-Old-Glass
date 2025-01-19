import mitsuba
mitsuba.set_variant("scalar_rgb")
from mitsuba.core.xml import load_dict
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
    "solver_threshold": 0.0001,
    "uniqueness_threshold": 0.000001,
    "max_trials": 1000,
    "bounces": 2,
    "use_SMS": False
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
    "solver_threshold": 0.0001,
    "uniqueness_threshold": 0.000001,
    "max_trials": 1000,
    "bounces": 2,
}

FNEE_integrator_SMS = {
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
    "use_SMS": True
}

path_integrator = {
    "type": "path",
    "max_depth": 8,
    "samples_per_pass": 1,
    "timeout": 10,
}


###############
# BSDF_WINDOW #
###############
bsdf_window = load_dict({
    "type": "dielectric",
    "int_ior": 1.54
})

################
# BSDF_DIFFUSE #
################
bsdf_diffuse = load_dict({
    "type": "diffuse",
    "reflectance": {
        "type": "rgb",
        "value": [0.4, 0.4, 0.4],
    }
})

###################
# HEIGHFIELD_DICT #
###################
Texture_cive_1 = load_dict({
    "type": "bitmap",
    "id": "id_texture_cive_1",
    "filename": "textures/Heightmaps/Cive_1.exr",
    "raw": True
})

Texture_manchon_1 = load_dict({
    "type": "bitmap",
    "id": "id_texture_manchon_1",
    "filename": "textures/Heightmaps/manchon_1_1024.exr",
    "raw": True
})

Texture_cive_2 = load_dict({
    "type": "bitmap",
    "id": "id_texture_cive_2",
    "filename": "textures/Heightmaps/Cive_2.exr",
    "raw": True
})

Texture_flat = load_dict({
    "type": "bitmap",
    "id": "id_texture_flat",
    "filename": "textures/Heightmaps/flat.exr",
    "raw": True
})

Texture_simplex = load_dict({
    "type": "bitmap",
    "id": "id_texture_simplex_1",
    "filename": "textures/Heightmaps/Octave_1_LPF.exr",
    "raw": True
})

##############
# PANEL MESH #
##############
panel_mesh = load_dict({
    "type": "obj",
    "filename": "meshes/caustic_panel_thin.obj",
    "bsdf": bsdf_window,
        "caustic_caster_single": True,
    "caustic_caster_multi" : True
})

panel_manchon_mesh = load_dict({
    "type": "obj",
    "filename": "meshes/Panel_manchon.obj",
    "bsdf": bsdf_window,
    "caustic_caster_multi" : True,
#        "caustic_caster_single": True,
})

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

    "Plomb": {
        "type": "obj",
        "filename": "meshes/plomb.obj",
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
            4.942975044250488 + 0.1, -0.5388131737709045, 1.9822726249694824],
            up=[-0.17331193387508392, -0.2800687253475189, 0.9442057609558105]),
        "film": {
            "type": "hdrfilm",
            "rfilter": {"type": "box"},
            "width": 1024,
            "height": 1024,
        },
        "sampler": {
            "type": "independent",
            "sample_count": 100000,
        }
    },

    # "Texture_cive_1": Texture_cive_1,

    "integrator": load_dict(path_integrator)
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
        "caustic_caster_multi": False,
        "double_refraction_caster": True,
        "self_id": 1,
        "pair_id": 2,
        "partition_number": 1.0,
        "bsdf": bsdf_window,
        "heightmap": Texture_cive_1
    }

    Heighfield_dict["H"] = H1
    Heighfield_dict["to_world"] = to_world_1
    Heighfield_dict["flip_normals"] = False
    Heighfield_dict["self_id"] = id * 2
    Heighfield_dict["pair_id"] = id * 2 + 1
    Heighfield_dict["heightmap"] = texture_id_1
    Heighfield_dict["caustic_bouncer"] = True
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
    Heighfield_dict["caustic_caster_single"] = True
    Heighfield_dict["caustic_caster_multi"] = True

    scene_dict[window_name + "_2"] = Heighfield_dict.copy()

def build_heightfield_scene():
    scene_heightfield_dict = scene_dict.copy()
    Textures = [[Texture_flat, Texture_cive_1], [Texture_flat, Texture_manchon_1], [
        Texture_flat, Texture_cive_2], [Texture_flat, Texture_simplex]]

#      H = [[0, 0.0], [0, 0.0], [0, 0.0], [0, 0.0]] # flat window
    H = [[0, 0.01], [0, 0.001], [0, 0.001], [0, 0.001]] # displaced windows

    e = 0.02 # glass panel thickness

    ## 4x4 panels
    for i in range(2):
        for j in range(2):
            index = 2*i + j
            T1 = ScalarTransform4f.translate(
                [1.9 + i, -3.3, 0.0 + j]) * ScalarTransform4f.rotate([0, 0, 1], -90.0)
            T2 = ScalarTransform4f.translate(
                [1.9 + i, -3.3 + e, 0.0 + j]) * ScalarTransform4f.rotate([0, 0, 1], -90.0)
            add_window_to_scene(scene_heightfield_dict, "H" + str(i) + str(j), index,
                                T1, T2, H[index][0], H[index][1], Textures[index][0], Textures[index][1], L=1, W=1)

    return scene_heightfield_dict

def build_gpu_scene():
    scene_gpu_dict = scene_dict.copy()
    scene_gpu_dict["integrator"] = path_integrator
    scene_gpu_dict["panel"] = panel_mesh
    return scene_gpu_dict

def build_mesh_scene():
    scene_gpu_dict = scene_dict.copy()
    scene_gpu_dict["integrator"] = path_integrator
    scene_gpu_dict["panel"] = panel_mesh
    return scene_gpu_dict

def build_manchon_mesh_scene():
    scene_gpu_dict = scene_dict.copy()
    scene_gpu_dict["integrator"] = path_integrator
    scene_gpu_dict["panel"] = panel_manchon_mesh
    return scene_gpu_dict

#########
# SCENE #
#########
