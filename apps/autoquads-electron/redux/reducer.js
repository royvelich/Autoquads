import {
    CHANGE_SPLIT_ORIENTATION,
    CHANGE_MODEL_VIEWPORT_COLOR,
    CHANGE_SOUP_VIEWPORT_COLOR,
    CHANGE_MODEL_COLOR,
    CHANGE_SOUP_COLOR,
    CHANGE_MODEL_WIREFRAME_VISIBILITY,
    CHANGE_SOUP_WIREFRAME_VISIBILITY, 
    CHANGE_MODEL_VIEW_VISIBILITY,
    CHANGE_SOUP_VIEW_VISIBILITY,
    CHANGE_AUTOCUTS_WEIGHT,    
    CHANGE_DELTA,
    CHANGE_LAMBDA,
    CHANGE_ZETA,
    CHANGE_SEAMLESS_WEIGHT,
    CHANGE_SINGULARITY_WEIGHT,
    CHANGE_SINGULARITY_INTERVAL,
    CHANGE_POSITION_WEIGHT,
    CHANGE_GRID_HORIZONTAL_COLOR,
    CHANGE_GRID_VERTICAL_COLOR,
    CHANGE_GRID_BACKGROUND_COLOR1,
    CHANGE_GRID_BACKGROUND_COLOR2,
    CHANGE_HIGHLIGHTED_FACE_COLOR,
    CHANGE_DRAGGED_FACE_COLOR,
    CHANGE_FIXED_FACE_COLOR,
    CHANGE_GRID_SIZE,
    CHANGE_GRID_TEXTURE_SIZE,
    CHANGE_GRID_LINE_WIDTH,
    CHANGE_UNIT_GRID_VISIBILITY,
    CHANGE_SOUP_VIEW_GRID_TEXTURE_VISIBILITY,
    CHANGE_OPTIMIZATION_DATA_MONITOR_VISIBILITY,
    CHANGE_SOLVER_STATE,
    CHANGE_MODEL_FILENAME,
    CHANGE_MODULE_FILENAME,
    CHANGE_MODEL_STATE,
    CHANGE_MODULE_STATE,
    CHANGE_OBJECTIVE_FUNCTION_PROPERTY_COLOR,
    CHANGE_OBJECTIVE_FUNCTION_PROPERTY_VISIBILITY
} from './actions.js';

import * as EnumsExports from './enums.js';

const INITIAL_STATE = {
    splitOrientation: EnumsExports.SplitOrientation.HORIZONTAL,
    modelViewportColor: 'rgb(200,200,200)',
    soupViewportColor: 'rgb(200,200,200)',
    modelColor: 'rgb(255, 255, 255)',
    soupColor: 'rgb(255, 255, 255)',
    modelWireframeVisibility: EnumsExports.Visibility.HIDDEN,
    soupWireframeVisibility: EnumsExports.Visibility.VISIBLE, 
    modelViewVisibility: EnumsExports.Visibility.VISIBLE,
    soupViewVisibility: EnumsExports.Visibility.VISIBLE,
    autocutsWeight: 1,
    delta: 0.5,
    lambda: 0.2,
    zeta: 0.95,
    seamlessWeight: 0,
    singularityWeight: 0,
    singularityInterval: 1,
    positionWeight: 100,
    gridHorizontalColor: 'rgb(0,0,0)',
    gridVerticalColor: 'rgb(0,0,0)',
    gridBackgroundColor1: 'rgb(170,170,170)',
    gridBackgroundColor2: 'rgb(170,170,170)',
    highlightedFaceColor: 'rgb(0,255,0)',
    draggedFaceColor: 'rgb(0,0,255)',
    fixedFaceColor: 'rgb(255,0,0)',
    gridSize: 7,
    gridTextureSize: 8,
    gridLineWidth: 4,
    unitGridVisibility: EnumsExports.Visibility.HIDDEN,
    soupViewGridTextureVisibility: EnumsExports.Visibility.HIDDEN,
    optimizationDataMonitorVisibility: EnumsExports.Visibility.VISIBLE,
    solverState: EnumsExports.SolverState.OFF,
    modelFilename: '',
    moduleFilename: require('path').join(appRoot, 'node-addon.node'),
    modelState: EnumsExports.LoadState.UNLOADED,
    moduleState: EnumsExports.LoadState.UNLOADED,
    objectiveFunctionsProperties: [
        {
            objectiveFunctionId: 'Separation',
            propertyId: 'value_per_vertex',
            propertyEffectType: 'vertex_color',
            objectiveFunctionName: 'Separation',
            propertyName: 'Value Per Vertex',
            color: 'rgb(255,0,0)',
            visibility: EnumsExports.Visibility.HIDDEN
        },
        {
            objectiveFunctionId: 'Symmetric Dirichlet',
            propertyId: 'value_per_vertex',
            propertyEffectType: 'vertex_color',
            objectiveFunctionName: 'Symmetric Dirichlet',
            propertyName: 'Value Per Vertex',
            color: 'rgb(255,0,0)',
            visibility: EnumsExports.Visibility.HIDDEN            
        },
        {
            objectiveFunctionId: 'Seamless',
            propertyId: 'value_per_vertex',
            propertyEffectType: 'vertex_color',
            objectiveFunctionName: 'Seamless',
            propertyName: 'Value Per Vertex',
            color: 'rgb(255,0,0)',
            visibility: EnumsExports.Visibility.HIDDEN            
        },
        {
            objectiveFunctionId: 'Singular Points',
            propertyId: 'singularity_weight_per_vertex',
            propertyEffectType: 'vertex_color',
            objectiveFunctionName: 'Singular Points',
            propertyName: 'Singularity Weight Per Vertex',
            color: 'rgb(255,0,0)',
            visibility: EnumsExports.Visibility.HIDDEN            
        },
        {
            objectiveFunctionId: 'Singular Points',
            propertyId: 'positive_angular_defect_singularities_indices',
            propertyEffectType: 'scene_objects',
            objectiveFunctionName: 'Singular Points',
            propertyName: 'Positive Angular Defect',
            color: 'rgb(0,0,255)',
            visibility: EnumsExports.Visibility.HIDDEN            
        },
        {
            objectiveFunctionId: 'Singular Points',
            propertyId: 'negative_angular_defect_singularities_indices',
            propertyEffectType: 'scene_objects',
            objectiveFunctionName: 'Singular Points',
            propertyName: 'Negative Angular Defect',
            color: 'rgb(255,0,0)',
            visibility: EnumsExports.Visibility.HIDDEN            
        }             
    ]
};

function matchObjectiveFunctionProperty(objectivePropertyVisualData, action) {
    if(objectivePropertyVisualData.objectiveFunctionId === action.objectiveFunctionId && 
        objectivePropertyVisualData.propertyId === action.propertyId)
    {
        return true;
    }

    return false;
}

export const reducer = (state = INITIAL_STATE, action) => {
    switch (action.type) {
        case CHANGE_SPLIT_ORIENTATION:
            return {
                ...state,
                splitOrientation: action.orientation
            };
        case CHANGE_MODEL_VIEWPORT_COLOR:
            return {
                ...state,
                modelViewportColor: action.color
            };
        case CHANGE_SOUP_VIEWPORT_COLOR:
            return {
                ...state,
                soupViewportColor: action.color
            };
        case CHANGE_MODEL_COLOR:
            return {
                ...state,
                modelColor: action.color
            };
        case CHANGE_SOUP_COLOR:
            return {
                ...state,
                soupColor: action.color
            };
        case CHANGE_MODEL_WIREFRAME_VISIBILITY:
            return {
                ...state,
                modelWireframeVisibility: action.visibility
            };
        case CHANGE_SOUP_WIREFRAME_VISIBILITY:
            return {
                ...state,
                soupWireframeVisibility: action.visibility
            };
        case CHANGE_MODEL_VIEW_VISIBILITY:
            return {
                ...state,
                modelViewVisibility: action.visibility
            };
        case CHANGE_SOUP_VIEW_VISIBILITY:
            return {
                ...state,
                soupViewVisibility: action.visibility
            };
        case CHANGE_AUTOCUTS_WEIGHT:
            return {
                ...state,
                autocutsWeight: action.weight
            };            
        case CHANGE_DELTA:
            return {
                ...state,
                delta: action.delta
            };
        case CHANGE_LAMBDA:
            return {
                ...state,
                lambda: action.lambda
            };
        case CHANGE_ZETA:
            return {
                ...state,
                zeta: action.zeta
            };            
        case CHANGE_SEAMLESS_WEIGHT:
            return {
                ...state,
                seamlessWeight: action.weight
            };
        case CHANGE_SINGULARITY_WEIGHT:
            return {
                ...state,
                singularityWeight: action.weight
            };
        case CHANGE_SINGULARITY_INTERVAL:
            return {
                ...state,
                singularityInterval: action.interval
            };
        case CHANGE_POSITION_WEIGHT:
            return {
                ...state,
                postionWeight: action.weight
            };
        case CHANGE_GRID_HORIZONTAL_COLOR:
            return {
                ...state,
                gridHorizontalColor: action.color
            };
        case CHANGE_GRID_VERTICAL_COLOR:
            return {
                ...state,
                gridVerticalColor: action.color
            };
        case CHANGE_GRID_BACKGROUND_COLOR1:
            return {
                ...state,
                gridBackgroundColor1: action.color
            };
        case CHANGE_GRID_BACKGROUND_COLOR2:
            return {
                ...state,
                gridBackgroundColor2: action.color
            };
        case CHANGE_HIGHLIGHTED_FACE_COLOR:
            return {
                ...state,
                highlightedFaceColor: action.color
            };
        case CHANGE_DRAGGED_FACE_COLOR:
            return {
                ...state,
                draggedFaceColor: action.color
            };
        case CHANGE_FIXED_FACE_COLOR:
            return {
                ...state,
                fixedFaceColor: action.color
            };
        case CHANGE_GRID_SIZE:
            return {
                ...state,
                gridSize: action.size
            };
        case CHANGE_GRID_TEXTURE_SIZE:
            return {
                ...state,
                gridTextureSize: action.size
            };

        case CHANGE_GRID_LINE_WIDTH:
            return {
                ...state,
                gridLineWidth: action.width
            };            
        case CHANGE_UNIT_GRID_VISIBILITY:
            return {
                ...state,
                unitGridVisibility: action.visibility
            };    
        case CHANGE_SOUP_VIEW_GRID_TEXTURE_VISIBILITY:
            return {
                ...state,
                soupViewGridTextureVisibility: action.visibility
            }; 
        case CHANGE_OPTIMIZATION_DATA_MONITOR_VISIBILITY:
            return {
                ...state,
                optimizationDataMonitorVisibility: action.visibility
            }; 
        case CHANGE_SOLVER_STATE:
            return {
                ...state,
                solverState: action.state
            };
        case CHANGE_MODEL_FILENAME:
            return {
                ...state,
                modelFilename: action.filename
            };
        case CHANGE_MODULE_FILENAME:
            return {
                ...state,
                moduleFilename: action.filename
            };
        case CHANGE_MODEL_STATE:
            return {
                ...state,
                modelState: action.state
            };
        case CHANGE_MODULE_STATE:
            return {
                ...state,
                moduleState: action.state
            };
        case CHANGE_OBJECTIVE_FUNCTION_PROPERTY_COLOR:
            var newState = {
                ...state
            }

            newState.objectiveFunctionsProperties = [...newState.objectiveFunctionsProperties];
            for(let objectiveFunctionProperty of newState.objectiveFunctionsProperties)
            {
                if(matchObjectiveFunctionProperty(objectiveFunctionProperty, action)) {
                    objectiveFunctionProperty.color = action.color;
                    break;
                }
            }
            
            return newState;
        case CHANGE_OBJECTIVE_FUNCTION_PROPERTY_VISIBILITY:
            var newState = {
                ...state
            }

            newState.objectiveFunctionsProperties = [...newState.objectiveFunctionsProperties];
            for(let objectiveFunctionProperty of newState.objectiveFunctionsProperties)
            {
                if(matchObjectiveFunctionProperty(objectiveFunctionProperty, action)) {
                    objectiveFunctionProperty.visibility = action.visibility;
                    break;
                }
            }

            return newState;
        default:
            return state;
    }
};