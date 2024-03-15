"""
FILE: fw_model_landmark.py
LAST MODIFIED: 15-03-2024
DESCRIPTION: Functions for creating evaluator functions for anatomic landmarks
on fieldwork models

You can get a list of all landmarks by

fw_model_landmarks.landmarkNames()

To evaluate a landmark (e.g. femoral head centre), we first generate a function
to evaluate that landmark:

femur_HC_evaluator = fw_model_landmarks.makeLandmarkEvaluator('femur-HC', gf)

where gf is the geometric_field of a femur. Then to evaluate the the landmark,
we call the generated function with the gf parameters:

femur_HC_coords = femur_HC_evaluator(gf.field_parameters)

It is implemented this way so that landmark coordinates can be evaluated quickly
during fitting optimisations by simpling providing the gf parameters. The
femur_HC_evaluator function can now be called in your script with new gf
parameters to get the new landmark coordinates if the parameters change (e.g.
if the gf has been fitted or transformed).

===============================================================================
This file is part of GIAS3. (https://github.com/musculoskeletal/gias3)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import numpy as np

from gias3.common import geoprimitives
from gias3.fieldwork.field import geometric_field


def _make_evaluator(node_index):
    def evaluator(mesh_params):
        return mesh_params[:, node_index].squeeze()
    return evaluator


# ===========================================================================#
# femur landmark variables

FEMUR_LEFT_MEC_NODE = 630   # was 633
FEMUR_RIGHT_MEC_NODE = 628
FEMUR_LEFT_LEC_NODE = 550   # was 546
FEMUR_RIGHT_LEC_NODE = 550
FEMUR_LEFT_GT_NODE = 172
FEMUR_RIGHT_GT_NODE = 174
FEMUR_LEFT_LT_NODE = 276
FEMUR_RIGHT_LT_NODE = 276

FEMUR_SUB_TROCHANTER_NODES = [259, 260, 261, 262, 263, 291, 292, 293, 294, 307, 308, 309]
FEMUR_MID_SHAFT_NODES = [323, 334, 335, 336, 337, 346, 347, 348, 319, 320, 321, 322, 323]
FEMUR_CONDYLE_ALIGNMENT_NODES = [546, 633]

FEMUR_HEAD_ELEMENT = 0
FEMUR_HEAD_ELEMENTS = [0, 1, 2, 3, 4, 5, 6, 7]
FEMUR_NECK_ELEMENTS = [13, 14, 15, 16]
FEMUR_SHAFT_ELEMENTS = [23, 24, 25, 40, 41, 42]
FEMUR_NECK_LONG_ELEMENTS = [1, 2, 4, 6, 13, 14, 15, 16]
FEMUR_LATERAL_CONDYLE_ELEMENTS = [43, 44, 45, 46, 47]
FEMUR_MEDIAL_CONDYLE_ELEMENTS = [48, 49, 50, 51, 52, 53]
FEMUR_GREATER_TROCHANTER_ELEMENTS = [8, 9, 10, 11, 12]
FEMUR_PROXIMAL_ELEMENTS = [*FEMUR_HEAD_ELEMENTS, *FEMUR_GREATER_TROCHANTER_ELEMENTS, 17, 18, 19, 20, 21, 22]
FEMUR_DISTAL_ELEMENTS = [*FEMUR_LATERAL_CONDYLE_ELEMENTS, *FEMUR_MEDIAL_CONDYLE_ELEMENTS, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39]
FEMUR_LATERAL_EPICONDYLE_ELEMENTS = [45, ]
FEMUR_MEDIAL_EPICONDYLE_ELEMENTS = [53, ]


def make_evaluator_femur_head_centre(gf, flattened=False, side='left'):
    if flattened:
        headNodes = []
        for e in FEMUR_HEAD_ELEMENTS:
            headNodes += list(gf.ensemble_field_function.mapper._element_to_ensemble_map[e].keys())
            headNodes = list(set(headNodes))
    else:
        headNodes = list(gf.ensemble_field_function.mapper._element_to_ensemble_map[FEMUR_HEAD_ELEMENT].keys())

    def eval_femur_head_centre(mesh_params):
        return geoprimitives.fitSphereAnalytic(mesh_params[:, headNodes].squeeze().T)[0]

    return eval_femur_head_centre


def make_evaluator_femur_medial_epicondyle(gf, side='left'):
    if side == 'left':
        return _make_evaluator(FEMUR_LEFT_MEC_NODE)
    elif side == 'right':
        return _make_evaluator(FEMUR_RIGHT_MEC_NODE)


def make_evaluator_femur_lateral_epicondyle(gf, side='left'):
    if side == 'left':
        return _make_evaluator(FEMUR_LEFT_LEC_NODE)
    elif side == 'right':
        return _make_evaluator(FEMUR_RIGHT_LEC_NODE)


def make_evaluator_femur_greater_trochanter(gf, side='left'):
    if side == 'left':
        return _make_evaluator(FEMUR_LEFT_GT_NODE)
    elif side == 'right':
        return _make_evaluator(FEMUR_RIGHT_GT_NODE)


def make_evaluator_femur_lesser_trochanter(gf, side='left'):
    if side == 'left':
        return _make_evaluator(FEMUR_LEFT_LT_NODE)
    elif side == 'right':
        return _make_evaluator(FEMUR_RIGHT_LT_NODE)


def make_evaluator_femur_knee_centre(gf, side='left'):
    evalMC = make_evaluator_femur_medial_epicondyle(gf, side=side)
    evalLC = make_evaluator_femur_lateral_epicondyle(gf, side=side)

    def eval_femur_knee_centre(mesh_params):
        mc = evalMC(mesh_params)
        lc = evalLC(mesh_params)
        return (mc + lc) * 0.5

    return eval_femur_knee_centre


# ===========================================================================#
# hemi pelvis landmark variables

HEMI_PELVIS_LEFT_ASIS_NODE = 464
HEMI_PELVIS_RIGHT_ASIS_NODE = 466
HEMI_PELVIS_LEFT_PSIS_NODE = 384
HEMI_PELVIS_RIGHT_PSIS_NODE = 384
HEMI_PELVIS_LEFT_PS_NODE = 90
HEMI_PELVIS_RIGHT_PS_NODE = 92
HEMI_PELVIS_LEFT_PT_NODE = 102
HEMI_PELVIS_RIGHT_PT_NODE = 103
HEMI_PELVIS_LEFT_IS_NODE = 233
HEMI_PELVIS_RIGHT_IS_NODE = 233
HEMI_PELVIS_LEFT_IT_NODE = 26
HEMI_PELVIS_RIGHT_IT_NODE = 25
HEMI_PELVIS_LEFT_AN_NODE = 287
HEMI_PELVIS_RIGHT_AN_NODE = 289

HEMI_PELVIS_ACETABULUM_ELEMENTS = [36, 35, 38, 39, 40, 41, 42]


def make_evaluator_hemi_pelvis_acetabular_centre(gf):
    m = gf.ensemble_field_function.mapper._element_to_ensemble_map
    acNodes = []
    for e in HEMI_PELVIS_ACETABULUM_ELEMENTS:
        acNodes.extend(list(m[e].keys()))

    def eval_acetabular_centre(mesh_params):
        return geoprimitives.fitSphereAnalytic(mesh_params[:, acNodes].squeeze().T)[0]

    return eval_acetabular_centre


def make_evaluator_hemi_pelvis_asis(gf, side='left'):
    if side == 'left':
        return _make_evaluator(HEMI_PELVIS_LEFT_ASIS_NODE)
    elif side == 'right':
        return _make_evaluator(HEMI_PELVIS_RIGHT_ASIS_NODE)


def make_evaluator_hemi_pelvis_psis(gf, side='left'):
    if side == 'left':
        return _make_evaluator(HEMI_PELVIS_LEFT_PSIS_NODE)
    elif side == 'right':
        return _make_evaluator(HEMI_PELVIS_RIGHT_PSIS_NODE)


def make_evaluator_hemi_pelvis_pubis_symphysis(gf, side='left'):
    if side == 'left':
        return _make_evaluator(HEMI_PELVIS_LEFT_PS_NODE)
    elif side == 'right':
        return _make_evaluator(HEMI_PELVIS_RIGHT_PS_NODE)


def make_evaluator_hemi_pelvis_pubis_tubercle(gf, side='left'):
    if side == 'left':
        return _make_evaluator(HEMI_PELVIS_LEFT_PT_NODE)
    elif side == 'right':
        return _make_evaluator(HEMI_PELVIS_RIGHT_PT_NODE)


def make_evaluator_hemi_pelvis_ilial_spine(gf, side='left'):
    if side == 'left':
        return _make_evaluator(HEMI_PELVIS_LEFT_IS_NODE)
    elif side == 'right':
        return _make_evaluator(HEMI_PELVIS_RIGHT_IS_NODE)


def make_evaluator_hemi_pelvis_ischial_tuberosity(gf, side='left'):
    if side == 'left':
        return _make_evaluator(HEMI_PELVIS_LEFT_IT_NODE)
    elif side == 'right':
        return _make_evaluator(HEMI_PELVIS_RIGHT_IT_NODE)


def make_evaluator_hemi_pelvis_acetabular_notch(gf, side='left'):
    if side == 'left':
        return _make_evaluator(HEMI_PELVIS_LEFT_AN_NODE)
    elif side == 'right':
        return _make_evaluator(HEMI_PELVIS_RIGHT_AN_NODE)


# ===========================================================================#
# whole pelvis landmarks

PELVIS_LASIS_NODE = 1004
PELVIS_RASIS_NODE = 466
PELVIS_LPSIS_NODE = 924
PELVIS_RPSIS_NODE = 384
PELVIS_LPT_NODE = 642
PELVIS_RPT_NODE = 103
PELVIS_LPS_NODE = 630
PELVIS_RPS_NODE = 92
PELVIS_LIS_NODE = 773
PELVIS_RIS_NODE = 233
PELVIS_LIT_NODE = 566
PELVIS_RIT_NODE = 25
PELVIS_SAC_PLAT_NODE = 1300  # centre-posterior-most point on the vertebral plateau on the sacrum

PELVIS_LHJC_ELEMENTS = [109, 111, 112, 113, 114, 115]   # in flattened mesh
PELVIS_RHJC_ELEMENTS = [36, 38, 39, 40, 41, 42]     # in flattened mesh
PELVIS_RIGHT_HEMI_ELEMENTS = list(range(0, 73))
PELVIS_LEFT_HEMI_ELEMENTS = list(range(73, 146))
PELVIS_SACRUM_ELEMENTS = list(range(146, 260))


def make_evaluator_pelvis_lasis(gf, **kwargs):
    return _make_evaluator(PELVIS_LASIS_NODE)


def make_evaluator_pelvis_rasis(gf, **kwargs):
    return _make_evaluator(PELVIS_RASIS_NODE)


def make_evaluator_pelvis_lpsis(gf, **kwargs):
    return _make_evaluator(PELVIS_LPSIS_NODE)


def make_evaluator_pelvis_rpsis(gf, **kwargs):
    return _make_evaluator(PELVIS_RPSIS_NODE)


def make_evaluator_pelvis_lpt(gf, **kwargs):
    return _make_evaluator(PELVIS_LPT_NODE)


def make_evaluator_pelvis_rpt(gf, **kwargs):
    return _make_evaluator(PELVIS_RPT_NODE)


def make_evaluator_pelvis_lps(gf, **kwargs):
    return _make_evaluator(PELVIS_LPS_NODE)


def make_evaluator_pelvis_rps(gf, **kwargs):
    return _make_evaluator(PELVIS_RPS_NODE)


def make_evaluator_pelvis_lis(gf, **kwargs):
    return _make_evaluator(PELVIS_LIS_NODE)


def make_evaluator_pelvis_ris(gf, **kwargs):
    return _make_evaluator(PELVIS_RIS_NODE)


def make_evaluator_pelvis_lit(gf, **kwargs):
    return _make_evaluator(PELVIS_LIT_NODE)


def make_evaluator_pelvis_rit(gf, **kwargs):
    return _make_evaluator(PELVIS_RIT_NODE)


def make_evaluator_pelvis_sac_plat(gf, **kwargs):
    return _make_evaluator(PELVIS_SAC_PLAT_NODE)


def make_evaluator_pelvis_sacral(gf, **kwargs):
    """
    Mid-point of PSISes
    """

    def eval_pelvis_sacral(mesh_params):
        s = 0.5 * (mesh_params[:, PELVIS_LPSIS_NODE].squeeze() +
                   mesh_params[:, PELVIS_RPSIS_NODE].squeeze()
                   )
        return s

    return eval_pelvis_sacral


def make_evaluator_pelvis_lhjc(gf, disc=5.0, radius=False, side=None):
    # Make evaluator for left acetabulum elements
    acetabElemEval = geometric_field.makeGeometricFieldElementsEvaluatorSparse(
        gf, PELVIS_LHJC_ELEMENTS, disc
    )
    if radius:
        def eval_pelvis_lhjc(mesh_params):
            acetabPoints = acetabElemEval(mesh_params).T
            return geoprimitives.fitSphereAnalytic(acetabPoints)
    else:
        def eval_pelvis_lhjc(mesh_params):
            acetabPoints = acetabElemEval(mesh_params).T
            return geoprimitives.fitSphereAnalytic(acetabPoints)[0]
    return eval_pelvis_lhjc


def make_evaluator_pelvis_rhjc(gf, disc=5.0, radius=False, side=None):
    # Make evaluator for right acetabulum elements
    acetabElemEval = geometric_field.makeGeometricFieldElementsEvaluatorSparse(
        gf, PELVIS_RHJC_ELEMENTS, disc
    )
    if radius:
        def eval_pelvis_rhjc(mesh_params):
            acetabPoints = acetabElemEval(mesh_params).T
            return geoprimitives.fitSphereAnalytic(acetabPoints)
    else:
        def eval_pelvis_rhjc(mesh_params):
            acetabPoints = acetabElemEval(mesh_params).T
            return geoprimitives.fitSphereAnalytic(acetabPoints)[0]
    return eval_pelvis_rhjc


# ===========================================================================#
# tibia fibula combined landmarks

# aligned with opensim tibia model
TIBIA_FIBULA_LEFT_LC_NODE = 257     # 256
TIBIA_FIBULA_RIGHT_LC_NODE = 258
TIBIA_FIBULA_LEFT_MC_NODE = 236     # 235
TIBIA_FIBULA_RIGHT_MC_NODE = 235    # 237
TIBIA_FIBULA_LEFT_LM_NODE = 528
TIBIA_FIBULA_RIGHT_LM_NODE = 527    # 535
TIBIA_FIBULA_LEFT_MM_NODE = 150
TIBIA_FIBULA_RIGHT_MM_NODE = 151
TIBIA_FIBULA_LEFT_TT_NODE = 203
TIBIA_FIBULA_RIGHT_TT_NODE = 203

TIBIA_ELEMENTS = list(range(0, 46))
FIBULA_ELEMENTS = list(range(46, 88))

TIBIA_FIBULA_KNEE_CENTRE_OFFSET = 50.0  # 389.55 from KC to ankle average of 45 subjects from Mousa K.


def make_evaluator_tibia_fibula_lc(gf, side='left'):
    if side == 'left':
        return _make_evaluator(TIBIA_FIBULA_LEFT_LC_NODE)
    elif side == 'right':
        return _make_evaluator(TIBIA_FIBULA_RIGHT_LC_NODE)


def make_evaluator_tibia_fibula_mc(gf, side='left'):
    if side == 'left':
        return _make_evaluator(TIBIA_FIBULA_LEFT_MC_NODE)
    elif side == 'right':
        return _make_evaluator(TIBIA_FIBULA_RIGHT_MC_NODE)


def make_evaluator_tibia_fibula_mm(gf, side='left'):
    if side == 'left':
        return _make_evaluator(TIBIA_FIBULA_LEFT_MM_NODE)
    elif side == 'right':
        return _make_evaluator(TIBIA_FIBULA_RIGHT_MM_NODE)


def make_evaluator_tibia_fibula_lm(gf, side='left'):
    if side == 'left':
        return _make_evaluator(TIBIA_FIBULA_LEFT_LM_NODE)
    elif side == 'right':
        return _make_evaluator(TIBIA_FIBULA_RIGHT_LM_NODE)


def make_evaluator_tibia_fibula_tt(gf, side='left'):
    if side == 'left':
        return _make_evaluator(TIBIA_FIBULA_LEFT_TT_NODE)
    elif side == 'right':
        return _make_evaluator(TIBIA_FIBULA_RIGHT_TT_NODE)


def make_evaluator_tibia_fibula_knee_centre(gf, side='left'):
    evalLC = make_evaluator_tibia_fibula_lc(gf, side=side)
    evalMC = make_evaluator_tibia_fibula_mc(gf, side=side)
    evalLM = make_evaluator_tibia_fibula_lm(gf, side=side)
    evalMM = make_evaluator_tibia_fibula_mm(gf, side=side)

    def eval_tibia_fibula_knee_centre(mesh_params):
        lc = evalLC(mesh_params)
        mc = evalMC(mesh_params)
        lm = evalLM(mesh_params)
        mm = evalMM(mesh_params)

        # calc tibfib ACS
        ic = (mc + lc) / 2.0
        im = (mm + lm) / 2.0  # origin
        # superiorly, IM to IC
        y = geoprimitives.norm(ic - im)
        # anteriorly, normal to plane of IM, LC and MC
        x = geoprimitives.norm(np.cross(lc - im, mc - im))
        # right
        z = geoprimitives.norm(np.cross(x, y))

        # estimate knee centre
        kc = ic + y * TIBIA_FIBULA_KNEE_CENTRE_OFFSET
        return kc

    return eval_tibia_fibula_knee_centre


# ===========================================================================#
# patella landmarks
PATELLA_INF_NODE_LEFT = 29
PATELLA_INF_NODE_RIGHT = 29
PATELLA_SUP_NODE_LEFT = 59
PATELLA_SUP_NODE_RIGHT = 58
PATELLA_LAT_NODE_LEFT = 72
PATELLA_LAT_NODE_RIGHT = 72


def make_evaluator_patella_inf(gf, side='left'):
    if side == 'left':
        return _make_evaluator(PATELLA_INF_NODE_LEFT)
    elif side == 'right':
        return _make_evaluator(PATELLA_INF_NODE_RIGHT)


def make_evaluator_patella_sup(gf, side='left'):
    if side == 'left':
        return _make_evaluator(PATELLA_SUP_NODE_LEFT)
    elif side == 'right':
        return _make_evaluator(PATELLA_SUP_NODE_RIGHT)


def make_evaluator_patella_lat(gf, side='left'):
    if side == 'left':
        return _make_evaluator(PATELLA_LAT_NODE_LEFT)
    elif side == 'right':
        return _make_evaluator(PATELLA_LAT_NODE_RIGHT)


# ===========================================================================#


_landmarkEvaluators = {
    'femur-HC': make_evaluator_femur_head_centre,
    'femur-MEC': make_evaluator_femur_medial_epicondyle,
    'femur-LEC': make_evaluator_femur_lateral_epicondyle,
    'femur-GT': make_evaluator_femur_greater_trochanter,
    'femur-LT': make_evaluator_femur_lesser_trochanter,
    'femur-kneecentre': make_evaluator_femur_knee_centre,
    'hpelvis-ASIS': make_evaluator_hemi_pelvis_asis,
    'hpelvis-PSIS': make_evaluator_hemi_pelvis_psis,
    'hpelvis-PS': make_evaluator_hemi_pelvis_pubis_symphysis,
    'hpelvis-PT': make_evaluator_hemi_pelvis_pubis_tubercle,
    'hpelvis-IS': make_evaluator_hemi_pelvis_ilial_spine,
    'hpelvis-IT': make_evaluator_hemi_pelvis_ischial_tuberosity,
    'hpelvis-AN': make_evaluator_hemi_pelvis_acetabular_notch,
    'pelvis-LASIS': make_evaluator_pelvis_lasis,
    'pelvis-RASIS': make_evaluator_pelvis_rasis,
    'pelvis-LPSIS': make_evaluator_pelvis_lpsis,
    'pelvis-RPSIS': make_evaluator_pelvis_rpsis,
    'pelvis-Sacral': make_evaluator_pelvis_sacral,
    'pelvis-LPS': make_evaluator_pelvis_lps,
    'pelvis-RPS': make_evaluator_pelvis_rps,
    'pelvis-LIS': make_evaluator_pelvis_lis,
    'pelvis-RIS': make_evaluator_pelvis_ris,
    'pelvis-LIT': make_evaluator_pelvis_lit,
    'pelvis-RIT': make_evaluator_pelvis_rit,
    'pelvis-LHJC': make_evaluator_pelvis_lhjc,
    'pelvis-RHJC': make_evaluator_pelvis_rhjc,
    'pelvis-SacPlat': make_evaluator_pelvis_sac_plat,
    'tibiafibula-LC': make_evaluator_tibia_fibula_lc,
    'tibiafibula-MC': make_evaluator_tibia_fibula_mc,
    'tibiafibula-LM': make_evaluator_tibia_fibula_lm,
    'tibiafibula-MM': make_evaluator_tibia_fibula_mm,
    'tibiafibula-TT': make_evaluator_tibia_fibula_tt,
    'tibiafibula-kneecentre': make_evaluator_tibia_fibula_knee_centre,
    'patella-inf': make_evaluator_patella_inf,
    'patella-sup': make_evaluator_patella_sup,
    'patella-lat': make_evaluator_patella_lat,
}

validLandmarks = sorted(_landmarkEvaluators.keys())


def landmark_names():
    """Return a list of implemented landmarks
    """
    return list(_landmarkEvaluators.keys())


def make_landmark_evaluator(name, gf, **args):
    """
    Generate a function to evaluate the named landmark on the given
    geometric_field. Call landmarkNames to get a list of possible landmark
    names.

    inputs
    ------
    name : str of landmark name.
    gf : geometric_field instance on which to evaluate the landmark

    returns
    -------
    func : function that evaluates the named landmark given the field
        parameters of gf. e.g. ldmk_coords = func(gf.field_parameters)
    """
    if args is None:
        args = {}
    try:
        return _landmarkEvaluators[name](gf, **args)
    except KeyError:
        raise ValueError('Unknown landmark name ' + name)


# Define C++ style aliases for functions.
landmarkNames = landmark_names
makeLandmarkEvaluator = make_landmark_evaluator
