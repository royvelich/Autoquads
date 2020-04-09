import nanoid from '../web_modules/nano-id.js';
import * as Enums from './enums.js';

export const CHANGE_ALGORITHM_TYPE = 'CHANGE_ALGORITHM_TYPE';
export const CHANGE_SPLIT_ORIENTATION = 'CHANGE_SPLIT_ORIENTATION';
export const CHANGE_MODEL_VIEWPORT_COLOR = 'CHANGE_MODEL_VIEWPORT_COLOR';
export const CHANGE_SOUP_VIEWPORT_COLOR = 'CHANGE_SOUP_VIEWPORT_COLOR';
export const CHANGE_MODEL_COLOR = 'CHANGE_MODEL_COLOR';
export const CHANGE_SOUP_COLOR = 'CHANGE_SOUP_COLOR';
export const CHANGE_MODEL_WIREFRAME_VISIBILITY = 'CHANGE_MODEL_WIREFRAME_VISIBILITY';
export const CHANGE_SOUP_WIREFRAME_VISIBILITY = 'CHANGE_SOUP_WIREFRAME_VISIBILITY';
export const CHANGE_MODEL_FAT_WIREFRAME_VISIBILITY = 'CHANGE_MODEL_FAT_WIREFRAME_VISIBILITY';
export const CHANGE_SOUP_FAT_WIREFRAME_VISIBILITY = 'CHANGE_SOUP_FAT_WIREFRAME_VISIBILITY';
export const CHANGE_MODEL_VIEW_VISIBILITY = 'CHANGE_MODEL_VIEW_VISIBILITY';
export const CHANGE_SOUP_VIEW_VISIBILITY = 'CHANGE_SOUP_VIEW_VISIBILITY';
export const CHANGE_AUTOCUTS_WEIGHT = 'CHANGE_AUTOCUTS_WEIGHT';
export const CHANGE_DELTA = 'CHANGE_DELTA';
export const CHANGE_LAMBDA = 'CHANGE_LAMBDA';
export const CHANGE_ZETA = 'CHANGE_ZETA';
export const CHANGE_SEAMLESS_WEIGHT = 'CHANGE_SEAMLESS_WEIGHT';
export const CHANGE_SEAMLESS_ANGLE_WEIGHT = 'CHANGE_SEAMLESS_ANGLE_WEIGHT';
export const CHANGE_SEAMLESS_LENGTH_WEIGHT = 'CHANGE_SEAMLESS_LENGTH_WEIGHT';
export const CHANGE_SEAMLESS_TRANSLATION_WEIGHT = 'CHANGE_SEAMLESS_TRANSLATION_WEIGHT';
export const CHANGE_SELECTED_EDGE_SEAMLESS_WEIGHT = 'CHANGE_SELECTED_EDGE_SEAMLESS_WEIGHT';
export const CHANGE_SELECTED_EDGE_SEAMLESS_ANGLE_WEIGHT = 'CHANGE_SELECTED_EDGE_SEAMLESS_ANGLE_WEIGHT';
export const CHANGE_SELECTED_EDGE_SEAMLESS_LENGTH_WEIGHT = 'CHANGE_SELECTED_EDGE_SEAMLESS_LENGTH_WEIGHT';
export const CHANGE_SINGULARITY_WEIGHT = 'CHANGE_SINGULARITY_WEIGHT';
export const CHANGE_SINGULARITY_INTERVAL = 'CHANGE_SINGULARITY_INTERVAL';
export const CHANGE_POSITION_WEIGHT = 'CHANGE_POSITION_WEIGHT';
export const CHANGE_GRID_HORIZONTAL_COLOR = 'CHANGE_GRID_HORIZONTAL_COLOR';
export const CHANGE_GRID_VERTICAL_COLOR = 'CHANGE_GRID_VERTICAL_COLOR';
export const CHANGE_GRID_BACKGROUND_COLOR1 = 'CHANGE_GRID_BACKGROUND_COLOR1';
export const CHANGE_GRID_BACKGROUND_COLOR2 = 'CHANGE_GRID_BACKGROUND_COLOR2';
export const CHANGE_HIGHLIGHTED_FACE_COLOR = 'CHANGE_HIGHLIGHTED_FACE_COLOR';
export const CHANGE_DRAGGED_FACE_COLOR = 'CHANGE_DRAGGED_FACE_COLOR';
export const CHANGE_FIXED_FACE_COLOR = 'CHANGE_FIXED_FACE_COLOR';
export const CHANGE_EDITED_EDGE_COLOR = 'CHANGE_EDITED_EDGE_COLOR';
export const CHANGE_HIGHLIGHTED_EDGE_COLOR = 'CHANGE_HIGHLIGHTED_EDGE_COLOR';
export const CHANGE_GRID_SIZE = 'CHANGE_GRID_SIZE';
export const CHANGE_GRID_TEXTURE_SIZE = 'CHANGE_GRID_TEXTURE_SIZE';
export const CHANGE_GRID_LINE_WIDTH = 'CHANGE_GRID_LINE_WIDTH';
export const CHANGE_UNIT_GRID_VISIBILITY = 'CHANGE_UNIT_GRID_VISIBILITY';
export const CHANGE_SOUP_VIEW_GRID_TEXTURE_VISIBILITY = 'CHANGE_SOUP_VIEW_GRID_TEXTURE_VISIBILITY';
export const CHANGE_OPTIMIZATION_DATA_MONITOR_VISIBILITY = 'CHANGE_OPTIMIZATION_DATA_MONITOR_VISIBILITY';
export const CHANGE_SOLVER_STATE = 'CHANGE_SOLVER_STATE';
export const CHANGE_MODEL_FILENAME = 'CHANGE_MODEL_FILENAME';
export const CHANGE_MODULE_FILENAME = 'CHANGE_MODULE_FILENAME';
export const CHANGE_MODEL_STATE = 'CHANGE_MODEL_STATE';
export const CHANGE_MODULE_STATE = 'CHANGE_MODULE_STATE';
export const CHANGE_OBJECTIVE_FUNCTION_PROPERTY_COLOR = 'CHANGE_OBJECTIVE_FUNCTION_PROPERTY_COLOR';
export const CHANGE_OBJECTIVE_FUNCTION_PROPERTY_VISIBILITY = 'CHANGE_OBJECTIVE_FUNCTION_PROPERTY_VISIBILITY';

export const setAlgorithmType = algorithmType => {
    return {
        type: CHANGE_ALGORITHM_TYPE,
        algorithmType: algorithmType
    };
};

export const setSplitOrientation = orientation => {
    return {
        type: CHANGE_SPLIT_ORIENTATION,
        orientation: orientation
    };
};

export const setModelViewportColor = color => {
    return {
        type: CHANGE_MODEL_VIEWPORT_COLOR,
        color: color
    };
};

export const setSoupViewportColor = color => {
    return {
        type: CHANGE_SOUP_VIEWPORT_COLOR,
        color: color
    };
};

export const setModelColor = color => {
    return {
        type: CHANGE_MODEL_COLOR,
        color: color
    };
};

export const setSoupColor = color => {
    return {
        type: CHANGE_SOUP_COLOR,
        color: color
    };
};

export const showModelWireframe = () => {
    return {
        type: CHANGE_MODEL_WIREFRAME_VISIBILITY,
        visibility: Enums.Visibility.VISIBLE
    };
};

export const hideModelWireframe = () => {
    return {
        type: CHANGE_MODEL_WIREFRAME_VISIBILITY,
        visibility: Enums.Visibility.HIDDEN
    };
};

export const showSoupWireframe = () => {
    return {
        type: CHANGE_SOUP_WIREFRAME_VISIBILITY,
        visibility: Enums.Visibility.VISIBLE
    };
};

export const hideSoupWireframe = () => {
    return {
        type: CHANGE_SOUP_WIREFRAME_VISIBILITY,
        visibility: Enums.Visibility.HIDDEN
    };
};

export const showModelFatWireframe = () => {
    return {
        type: CHANGE_MODEL_FAT_WIREFRAME_VISIBILITY,
        visibility: Enums.Visibility.VISIBLE
    };
};

export const hideModelFatWireframe = () => {
    return {
        type: CHANGE_MODEL_FAT_WIREFRAME_VISIBILITY,
        visibility: Enums.Visibility.HIDDEN
    };
};

export const showSoupFatWireframe = () => {
    return {
        type: CHANGE_SOUP_FAT_WIREFRAME_VISIBILITY,
        visibility: Enums.Visibility.VISIBLE
    };
};

export const hideSoupFatWireframe = () => {
    return {
        type: CHANGE_SOUP_FAT_WIREFRAME_VISIBILITY,
        visibility: Enums.Visibility.HIDDEN
    };
};

export const showModelView = () => {
    return {
        type: CHANGE_MODEL_VIEW_VISIBILITY,
        visibility: Enums.Visibility.VISIBLE
    };
};

export const hideModelView = () => {
    return {
        type: CHANGE_MODEL_VIEW_VISIBILITY,
        visibility: Enums.Visibility.HIDDEN
    };
};

export const showSoupView = () => {
    return {
        type: CHANGE_SOUP_VIEW_VISIBILITY,
        visibility: Enums.Visibility.VISIBLE
    };
};

export const hideSoupView = () => {
    return {
        type: CHANGE_SOUP_VIEW_VISIBILITY,
        visibility: Enums.Visibility.HIDDEN
    };
};

export const setAutocutsWeight = weight => {
    return {
        type: CHANGE_AUTOCUTS_WEIGHT,
        weight: weight
    };
};

export const setDelta = delta => {
    return {
        type: CHANGE_DELTA,
        delta: delta
    };
};

export const setLambda = lambda => {
    return {
        type: CHANGE_LAMBDA,
        lambda: lambda
    };
};

export const setZeta = zeta => {
    return {
        type: CHANGE_ZETA,
        zeta: zeta
    };
};

export const setSeamlessWeight = weight => {
    return {
        type: CHANGE_SEAMLESS_WEIGHT,
        weight: weight
    };
};

export const setSeamlessAngleWeight = weight => {
    return {
        type: CHANGE_SEAMLESS_ANGLE_WEIGHT,
        weight: weight
    };
};

export const setSeamlessLengthWeight = weight => {
    return {
        type: CHANGE_SEAMLESS_LENGTH_WEIGHT,
        weight: weight
    };
};

export const setSeamlessTrnalsationWeight = weight => {
    return {
        type: CHANGE_SEAMLESS_TRANSLATION_WEIGHT,
        weight: weight
    };
};

export const setSelectedEdgeSeamlessAngleWeight = weight => {
    return {
        type: CHANGE_SELECTED_EDGE_SEAMLESS_ANGLE_WEIGHT,
        weight: weight
    };
};

export const setSelectedEdgeSeamlessLengthWeight = weight => {
    return {
        type: CHANGE_SELECTED_EDGE_SEAMLESS_LENGTH_WEIGHT,
        weight: weight
    };
};

export const setSingularityWeight = weight => {
    return {
        type: CHANGE_SINGULARITY_WEIGHT,
        weight: weight
    };
};

export const setSingularityInterval = interval => {
    return {
        type: CHANGE_SINGULARITY_INTERVAL,
        interval: interval
    };
};

export const setPositionWeight = weight => {
    return {
        type: CHANGE_POSITION_WEIGHT,
        weight: weight
    };
};

export const setGridHorizontalColor = color => {
    return {
        type: CHANGE_GRID_HORIZONTAL_COLOR,
        color: color
    };
};

export const setGridVerticalColor = color => {
    return {
        type: CHANGE_GRID_VERTICAL_COLOR,
        color: color
    };
};

export const setGridBackgroundColor1 = color => {
    return {
        type: CHANGE_GRID_BACKGROUND_COLOR1,
        color: color
    };
};

export const setGridBackgroundColor2 = color => {
    return {
        type: CHANGE_GRID_BACKGROUND_COLOR2,
        color: color
    };
};

export const setHighlightedFaceColor = color => {
    return {
        type: CHANGE_HIGHLIGHTED_FACE_COLOR,
        color: color
    };
};

export const setHighlightedEdgeColor = color => {
    return {
        type: CHANGE_HIGHLIGHTED_EDGE_COLOR,
        color: color
    };
};

export const setDraggedFaceColor = color => {
    return {
        type: CHANGE_DRAGGED_FACE_COLOR,
        color: color
    };
};

export const setFixedFaceColor = color => {
    return {
        type: CHANGE_FIXED_FACE_COLOR,
        color: color
    };
};

export const setEditedEdgeColor = color => {
    return {
        type: CHANGE_EDITED_EDGE_COLOR,
        color: color
    };
};

export const setGridSize = size => {
    return {
        type: CHANGE_GRID_SIZE,
        size: size
    };
};

export const setGridTextureSize = size => {
    return {
        type: CHANGE_GRID_TEXTURE_SIZE,
        size: size
    };
};

export const setGridLineWidth = width => {
    return {
        type: CHANGE_GRID_LINE_WIDTH,
        width: width
    };
};

export const showUnitGrid = () => {
    return {
        type: CHANGE_UNIT_GRID_VISIBILITY,
        visibility: Enums.Visibility.VISIBLE
    };
};

export const hideUnitGrid = () => {
    return {
        type: CHANGE_UNIT_GRID_VISIBILITY,
        visibility: Enums.Visibility.HIDDEN
    };
};

export const showSoupViewGridTexture = () => {
    return {
        type: CHANGE_SOUP_VIEW_GRID_TEXTURE_VISIBILITY,
        visibility: Enums.Visibility.VISIBLE
    };
};

export const hideSoupViewGridTexture = () => {
    return {
        type: CHANGE_SOUP_VIEW_GRID_TEXTURE_VISIBILITY,
        visibility: Enums.Visibility.HIDDEN
    };
};

export const showOptimizationDataMonitor = () => {
    return {
        type: CHANGE_OPTIMIZATION_DATA_MONITOR_VISIBILITY,
        visibility: Enums.Visibility.VISIBLE
    };
};

export const hideOptimizationDataMonitor = () => {
    return {
        type: CHANGE_OPTIMIZATION_DATA_MONITOR_VISIBILITY,
        visibility: Enums.Visibility.HIDDEN
    };
};

export const resumeSolver = () => {
    return {
        type: CHANGE_SOLVER_STATE,
        state: Enums.SolverState.ON
    };
};

export const pauseSolver = () => {
    return {
        type: CHANGE_SOLVER_STATE,
        state: Enums.SolverState.OFF
    };
};

export const loadModelFile = filename => {
    return {
        type: CHANGE_MODEL_FILENAME,
        filename: filename
    };
};

export const loadModuleFile = filename => {
    return {
        type: CHANGE_MODULE_FILENAME,
        filename: filename
    };
};

export const setModelState = state => {
    return {
        type: CHANGE_MODEL_STATE,
        state: state
    };
};

export const setModuleState = state => {
    return {
        type: CHANGE_MODULE_STATE,
        state: state
    };
};

export const setObjectiveFunctionPropertyColor = (objectiveFunctionId, propertyId, color) => {
    return {
        type: CHANGE_OBJECTIVE_FUNCTION_PROPERTY_COLOR,
        objectiveFunctionId: objectiveFunctionId,
        propertyId: propertyId,
        color: color
    };
};

export const setObjectiveFunctionPropertyVisibility = (objectiveFunctionId, propertyId, visibility) => {
    return {
        type: CHANGE_OBJECTIVE_FUNCTION_PROPERTY_VISIBILITY,
        objectiveFunctionId: objectiveFunctionId,
        propertyId: propertyId,
        visibility: visibility
    };
};