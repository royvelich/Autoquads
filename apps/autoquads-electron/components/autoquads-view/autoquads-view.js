// Web Modules Imports
import { LitElement, html, css } from '../../web_modules/lit-element.js';
import { classMap } from '../../web_modules/lit-html/directives/class-map.js';
import '../../web_modules/@vaadin/vaadin-button.js';
import '../../web_modules/@vaadin/vaadin-split-layout.js';
import '../../web_modules/@vaadin/vaadin-dialog.js';
import { PubSub } from '../../web_modules/pubsub-js.js';
import * as THREE from '../../web_modules/three.js';
import { connect } from '../../web_modules/pwa-helpers.js';

// Components Imports
import '../mesh-view/mesh-view.js';
import '../autoquads-side-bar/autoquads-side-bar.js';
import { MeshProvider } from '../mesh-provider/mesh-provider.js';
import { AutoquadsModelMeshProvider } from '../autoquads-model-mesh-provider/autoquads-model-mesh-provider.js';
import { AutoquadsSoupMeshProvider } from '../autoquads-soup-mesh-provider/autoquads-soup-mesh-provider.js';
import * as ReducerExports from '../../redux/reducer.js';
import * as ActionsExports from '../../redux/actions.js';
import * as EnumsExports from '../../redux/enums.js';
import * as HelpersExports from '../../redux/helpers.js';
import { store } from '../../redux/store.js';

export class AutoquadsView extends connect(store)(LitElement) {
    static get styles() {
        return [css`
            :host {
                display: flex;
                width: 100%;
                height: 100%;
            }

            .outer-split {
                flex: 1;
            }

            .inner-split {
                width: 100%;
                height: 100%;
            }

            /* https://stackoverflow.com/questions/41800823/hide-polymer-element */
            /* [hidden] {
                display: none;
            } */

            .hidden {
                display: none;
            }
        `];
    }

    render() {
        return html`  
            <vaadin-split-layout class="outer-split">
                <autoquads-side-bar></autoquads-side-bar>
                <vaadin-split-layout orientation="horizontal" class="inner-split">
                    <mesh-view 
                        id="model-mesh-view"
                        class="${classMap(this._modelMeshViewClasses)}"
                        use-lights
                        enable-mesh-rotation
                        enable-vertex-selection
                        caption="Mesh View"
                        grid-horizontal-color="${this.gridHorizontalColor}"
                        grid-vertical-color="${this.gridVerticalColor}"
                        grid-background-color1="${this.gridBackgroundColor1}"
                        grid-background-color2="${this.gridBackgroundColor2}"
                        grid-size="${this.gridSize}"
                        grid-texture-size="${this.gridTextureSize}"
                        grid-line-width="${this.gridLineWidth}"
                        ?show-wireframe="${HelpersExports.isVisible(this.modelWireframeVisibility)}"
                        ?show-fat-wireframe="${HelpersExports.isVisible(this.modelFatWireframeVisibility)}"
                        background-color="${this.modelViewportColor}"
                        mesh-color="${this.modelColor}"
                        .meshProvider="${this._modelMeshProvider}"
                        mesh-interaction="${this.meshInteraction}"
                        highlighted-face-color="${this.highlightedFaceColor}"
                        highlighted-edge-color="${this.highlightedEdgeColor}"
                        dragged-face-color="${this.draggedFaceColor}"
                        selected-face-color="${this.fixedFaceColor}"
                        selected-edge-color="${this.editedEdgeColor}"
                        show-grid-texture>
                    </mesh-view>
                    <mesh-view 
                        id="soup-mesh-view"
                        class="${classMap(this._soupMeshViewClasses)}"
                        enable-face-dragging caption="Soup View"
                        ?show-unit-grid="${HelpersExports.isVisible(this.unitGridVisibility)}"
                        grid-horizontal-color="${this.gridHorizontalColor}"
                        grid-vertical-color="${this.gridVerticalColor}"
                        grid-background-color1="${this.gridBackgroundColor1}"
                        grid-background-color2="${this.gridBackgroundColor2}"
                        grid-size="${this.gridSize}"
                        grid-texture-size="${this.gridTextureSize}"
                        grid-line-width="${this.gridLineWidth}"
                        ?show-wireframe="${HelpersExports.isVisible(this.soupWireframeVisibility)}"
                        ?show-fat-wireframe="${HelpersExports.isVisible(this.soupFatWireframeVisibility)}"
                        background-color="${this.soupViewportColor}"
                        mesh-color="${this.soupColor}"
                        .meshProvider="${this._soupMeshProvider}"
                        mesh-interaction="${this.meshInteraction}"
                        highlighted-face-color="${this.highlightedFaceColor}"
                        highlighted-edge-color="${this.highlightedEdgeColor}"
                        dragged-face-color="${this.draggedFaceColor}"
                        selected-face-color="${this.fixedFaceColor}"
                        selected-edge-color="${this.editedEdgeColor}"
                        ?show-debug-data="${HelpersExports.isVisible(this.optimizationDataMonitorVisibility)}"
                        show-grid-texture="${this.showGridTextureInSoupView}">
                    </mesh-view>
                </vaadin-split-layout>
            </vaadin-split-layout>
        `;
    }

    /**
     * Properties static declaration
     */

    static get properties() {
        return {
            splitOrientation: {
                type: String,
                attribute: 'split-orientation'
            },
            modelViewportColor: {
                type: String,
                attribute: 'model-viewport-color'
            },
            soupViewportColor: {
                type: String,
                attribute: 'soup-viewport-color'
            },
            modelColor: {
                type: String,
                attribute: 'model-color'
            },
            soupColor: {
                type: String,
                attribute: 'soup-color'
            },
            modelWireframeVisibility: {
                type: String,
                attribute: 'model-wireframe-visibility'
            },
            soupWireframeVisibility: {
                type: String,
                attribute: 'soup-wireframe-visibility'
            },
            modelFatWireframeVisibility: {
                type: String,
                attribute: 'model-fat-wireframe-visibility'
            },
            soupFatWireframeVisibility: {
                type: String,
                attribute: 'soup-fat-wireframe-visibility'
            },            
            modelViewVisibility: {
                type: String,
                attribute: 'model-view-visibility'
            },
            soupViewVisibility: {
                type: String,
                attribute: 'soup-view-visibility'
            },
            autocutsWeight: {
                type: Number,
                attribute: 'autocuts-weight'
            },            
            delta: {
                type: Number,
                attribute: 'delta'
            },
            lambda: {
                type: Number,
                attribute: 'lambda'
            },
            zeta: {
                type: Number,
                attribute: 'zeta'
            },            
            seamlessWeight: {
                type: Number,
                attribute: 'seamless-weight'
            },
            seamlessAngleWeight: {
                type: Number,
                attribute: 'seamless-angle-weight'
            },
            seamlessLengthWeight: {
                type: Number,
                attribute: 'seamless-length-weight'
            },
            seamlessTranslationWeight: {
                type: Number,
                attribute: 'seamless-translation-weight'
            },
            selectedEdgeSeamlessAngleWeight: {
                type: Number,
                attribute: 'selected-edge-seamless-angle-weight'
            },
            selectedEdgeSeamlessLengthWeight: {
                type: Number,
                attribute: 'selected-edge-seamless-length-weight'
            },
            singularityInterval: {
                type: Number,
                attribute: 'singularity-interval'
            },
            positionWeight: {
                type: Number,
                attribute: 'position-weight'
            },
            gridHorizontalColor: {
                type: String,
                attribute: 'grid-horizontal-color'
            },
            gridVerticalColor: {
                type: String,
                attribute: 'grid-vertical-color'
            },
            gridBackgroundColor1: {
                type: String,
                attribute: 'grid-background-color1'
            },
            gridBackgroundColor2: {
                type: String,
                attribute: 'grid-background-color2'
            },
            highlightedFaceColor: {
                type: String,
                attribute: 'highlighted-face-color'
            },
            highlightedEdgeColor: {
                type: String,
                attribute: 'highlighted-edge-color'
            },
            draggedFaceColor: {
                type: String,
                attribute: 'dragged-face-color'
            },
            fixedFaceColor: {
                type: String,
                attribute: 'fixed-face-color'
            },
            editedEdgeColor: {
                type: String,
                attribute: 'edited-edge-color'
            },
            gridSize: {
                type: Number,
                attribute: 'grid-size'
            },
            gridTextureSize: {
                type: Number,
                attribute: 'grid-texture-size'
            },
            gridLineWidth: {
                type: Number,
                attribute: 'grid-line-width'
            },
            unitGridVisibility: {
                type: String,
                attribute: 'unit-grid-visibility'
            },
            soupViewGridTextureVisibility: {
                type: String,
                attribute: 'soup-view-grid-texture-visibility'
            },
            optimizationDataMonitorVisibility: {
                type: String,
                attribute: 'optimization-data-monitor-visibility'
            },  
            solverState: {
                type: String,
                attribute: 'solver-state'
            },  
            model: {
                type: Object,
                attribute: 'model'
            },  
            module: {
                type: Object,
                attribute: 'module'
            },
            objectiveFunctionsProperties: {
                type: Object,
                attribute: 'objective-functions-properties'
            },
            algorithmType: {
                type: String,
                attribute: 'algorithm-type'
            }
        };
    }

    /**
     * State handling
     */
    stateChanged(state) {
        let localState = {};
        for (let [key, value] of Object.entries(state)) {
            localState[key] = this[key];
        }

        for (let [key, value] of Object.entries(state)) {
            let currentState = store.getState();
            if((state[key] === currentState[key]) && (localState[key] !== state[key])) {
                this[key] = state[key];
            }
        }    
    }

    /**
     * Constructor
     */
    constructor() {
        super();
        this._modelMeshProvider = new MeshProvider();
        this._soupMeshProvider = new MeshProvider();
        this._modelMeshViewClasses = { hidden: false };
        this._soupMeshViewClasses = { hidden: false };
        this._selectedEdge = null;      
    }

    /**
     * Properties
     */
    set splitOrientation(value) {
        const oldValue = this._splitOrientation;
        this._splitOrientation = value;
        this.requestUpdate('splitOrientation', oldValue);
    }

    get splitOrientation() {
        return this._splitOrientation;
    }
    
    set modelViewportColor(value) {
        const oldValue = this._modelViewportColor;
        this._modelViewportColor = value;
        this.requestUpdate('modelViewportColor', oldValue);
    }

    get modelViewportColor() {
        return this._modelViewportColor;
    }

    set soupViewportColor(value) {
        const oldValue = this._soupViewportColor;
        this._soupViewportColor = value;
        this.requestUpdate('soupViewportColor', oldValue);
    }

    get soupViewportColor() {
        return this._soupViewportColor;
    }

    set modelColor(value) {
        const oldValue = this._modelColor;
        this._modelColor = value;
        this._modelMeshProvider.meshColor = this._modelColor;
        this.requestUpdate('modelColor', oldValue);
    }

    get modelColor() {
        return this._modelColor;
    }

    set soupColor(value) {
        const oldValue = this._soupColor;
        this._soupColor = value;
        this._soupMeshProvider.meshColor = this._soupColor;
        this.requestUpdate('soupColor', oldValue);
    }

    get soupColor() {
        return this._soupColor;
    }

    set modelWireframeVisibility(value) {
        const oldValue = this._modelWireframeVisibility;
        this._modelWireframeVisibility = value;
        this.requestUpdate('modelWireframeVisibility', oldValue);
    }

    get modelWireframeVisibility() {
        return this._modelWireframeVisibility;        
    }

    set soupWireframeVisibility(value) {
        const oldValue = this._soupWireframeVisibility;
        this._soupWireframeVisibility = value;
        this.requestUpdate('soupWireframeVisibility', oldValue);
    }

    get soupWireframeVisibility() {
        return this._soupWireframeVisibility;        
    }

    set optimizationDataMonitorVisibility(value) {
        const oldValue = this._optimizationDataMonitorVisibility;
        this._optimizationDataMonitorVisibility = value;
        this.requestUpdate('optimizationDataMonitorVisibility', oldValue);
    }

    get optimizationDataMonitorVisibility() {
        return this._optimizationDataMonitorVisibility;        
    }

    set modelViewVisibility(value) {
        const oldValue = this._modelViewVisibility;
        this._modelViewVisibility = value;
        this._modelMeshViewClasses = {
            hidden: !HelpersExports.isVisible(value)
        }
        this.requestUpdate('modelViewVisibility', oldValue);
    }

    get modelViewVisibility() {
        return this._modelViewVisibility;
    }

    set soupViewVisibility(value) {
        const oldValue = this._soupViewVisibility;
        this._soupViewVisibility = value;
        this._soupMeshViewClasses = {
            hidden: !HelpersExports.isVisible(value)
        }
        this.requestUpdate('soupViewVisibility', oldValue);        
    }

    get soupViewVisibility() {
        return this._soupViewVisibility;
    }

    set autocutsWeight(value) {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            const oldValue = this._autocutsWeight;
            this._autocutsWeight = value;
            this._engine.setObjectiveFunctionProperty('Separation', 'weight', '', this.autocutsWeight * this.lambda);
            this._engine.setObjectiveFunctionProperty('Symmetric Dirichlet', 'weight', '', this.autocutsWeight * (1 - this.lambda));
            this.requestUpdate('autocutsWeight', oldValue);
        }
    }

    get autocutsWeight() {
        return this._autocutsWeight;
    }   

    set delta(value) {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            const oldValue = this._delta;
            this._delta = value;
            this._engine.setObjectiveFunctionProperty('Separation', 'delta', '', value);
            this.requestUpdate('delta', oldValue);
        }
    }

    get delta() {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            return this._engine.getObjectiveFunctionProperty('Separation', 'delta', 'none', '');
        }
    }

    set lambda(value) {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            const oldValue = this._lambda;
            this._lambda = value;
            this._engine.setObjectiveFunctionProperty('Separation', 'weight', '', this.autocutsWeight * value);
            this._engine.setObjectiveFunctionProperty('Symmetric Dirichlet', 'weight', '', this.autocutsWeight * (1 - value));
            this.requestUpdate('lambda', oldValue);
        }
    }

    get lambda() {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            return this._engine.getObjectiveFunctionProperty('Separation', 'weight', 'none', '');
        }
    }

    set zeta(value) {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            const oldValue = this._zeta;
            this._zeta = value;
            this._engine.setObjectiveFunctionProperty('Seamless', 'zeta', '', value);
            this.requestUpdate('zeta', oldValue);
        }
    }

    get zeta() {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            return this._engine.getObjectiveFunctionProperty('Seamless', 'zeta', 'none', '');
        }
    }    

    set seamlessWeight(value) {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            const oldValue = this._seamlessWeight;
            this._seamlessWeight = value;
            this._engine.setObjectiveFunctionProperty('Seamless', 'weight', '', value);
            this.requestUpdate('seamlessWeight', oldValue);
        }
    }

    get seamlessWeight() {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            return this._engine.getObjectiveFunctionProperty('Seamless', 'weight', 'none', '');
        }
    }

    set seamlessAngleWeight(value) {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            const oldValue = this._seamlessAngleWeight;
            this._seamlessAngleWeight = value;
            this._engine.setObjectiveFunctionProperty('Seamless', 'angle_weight', '', value);
            this.requestUpdate('seamlessAngleWeight', oldValue);
        }
    }

    get seamlessAngleWeight() {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            return this._engine.getObjectiveFunctionProperty('Seamless', 'angle_weight', 'none', '');
        }
    }

    set seamlessLengthWeight(value) {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            const oldValue = this._seamlessLengthWeight;
            this._seamlessLengthWeight = value;
            this._engine.setObjectiveFunctionProperty('Seamless', 'length_weight', '', value);
            this.requestUpdate('seamlessLengthWeight', oldValue);
        }
    }

    get seamlessLengthWeight() {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            return this._engine.getObjectiveFunctionProperty('Seamless', 'length_weight', 'none', '');
        }
    }

    set seamlessTranslationWeight(value) {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            const oldValue = this._seamlessTranslationWeight;
            this._seamlessTranslationWeight = value;
            this._engine.setObjectiveFunctionProperty('Seamless', 'translation_weight', '', value);
            this.requestUpdate('seamlessTranslationWeight', oldValue);
        }
    }

    get seamlessTranslationWeight() {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            return this._engine.getObjectiveFunctionProperty('Seamless', 'translation_weight', 'none', '');
        }
    }

    set selectedEdgeSeamlessAngleWeight(value) {
        if(HelpersExports.isModuleLoaded(this.moduleState) && this._selectedEdge !== null) {
            const oldValue = this._selectedEdgeSeamlessAngleWeight;
            this._selectedEdgeSeamlessAngleWeight = value;
            this._engine.setObjectiveFunctionProperty('Seamless', 'edge_angle_weight', this._selectedEdge.id, value);
            this.requestUpdate('selectedEdgeSeamlessAngleWeight', oldValue);
        }
    }

    get selectedEdgeSeamlessAngleWeight() {
        if(HelpersExports.isModuleLoaded(this.moduleState) && this._selectedEdge !== null) {
            return this._engine.getObjectiveFunctionProperty('Seamless', 'edge_angle_weight', 'none', this._selectedEdge.id);
        }
    }

    set selectedEdgeSeamlessLengthWeight(value) {
        if(HelpersExports.isModuleLoaded(this.moduleState) && this._selectedEdge !== null) {
            const oldValue = this._selectedEdgeSeamlessLengthWeight;
            this._selectedEdgeSeamlessLengthWeight = value;
            this._engine.setObjectiveFunctionProperty('Seamless', 'edge_length_weight', this._selectedEdge.id, value);
            this.requestUpdate('selectedEdgeSeamlessLengthWeight', oldValue);
        }
    }

    get selectedEdgeSeamlessLengthWeight() {
        if(HelpersExports.isModuleLoaded(this.moduleState) && this._selectedEdge !== null) {
            return this._engine.getObjectiveFunctionProperty('Seamless', 'edge_length_weight', 'none', this._selectedEdge.id);
        }
    }

    set singularityWeight(value) {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            const oldValue = this._singularityWeight;
            this._singularityWeight = value;
            this._engine.setObjectiveFunctionProperty('Singular Points Position', 'weight', '', value / 2);
            this.requestUpdate('singularityWeight', oldValue);
        }
    }

    get singularityWeight() {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            return this._engine.getObjectiveFunctionProperty('Singular Points Position', 'weight', 'none', '');
        }
    }

    set singularityInterval(value) {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            const oldValue = this._singularityInterval;
            this._singularityInterval = value;
            this._engine.setObjectiveFunctionProperty('Singular Points Position', 'interval', '', value);
            this._engine.setObjectiveFunctionProperty('Seamless', 'translation_interval', '', value);
            this.requestUpdate('singularityInterval', oldValue);
        }
    }

    get singularityInterval() {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            return this._engine.getObjectiveFunctionProperty('Singular Points Position', 'interval', 'none', '');
        }
    }

    set positionWeight(value) {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            const oldValue = this._positionWeight;
            this._positionWeight = value;
            this._engine.positionWeight = value;
            this.requestUpdate('positionWeight', oldValue);
        }
    }

    get positionWeight() {
        return this._positionWeight;
    }

    set gridHorizontalColor(value) {
        const oldValue = this._gridHorizontalColor;
        this._gridHorizontalColor = value;
        this.requestUpdate('gridHorizontalColor', oldValue);
    }

    get gridHorizontalColor() {
        return this._gridHorizontalColor;
    }

    set gridVerticalColor(value) {
        const oldValue = this._gridVerticalColor;
        this._gridVerticalColor = value;
        this.requestUpdate('gridVerticalColor', oldValue);
    }

    get gridVerticalColor() {
        return this._gridVerticalColor;
    } 

    set gridBackgroundColor1(value) {
        const oldValue = this._gridBackgroundColor1;
        this._gridBackgroundColor1 = value;
        this.requestUpdate('gridBackgroundColor1', oldValue);
    }

    get gridBackgroundColor1() {
        return this._gridBackgroundColor1;
    }

    set gridBackgroundColor2(value) {
        const oldValue = this._gridBackgroundColor2;
        this._gridBackgroundColor2 = value;
        this.requestUpdate('gridBackgroundColor2', oldValue);
    }

    get gridBackgroundColor2() {
        return this._gridBackgroundColor2;
    } 

    set highlightedFaceColor(value) {
        const oldValue = this._highlightedFaceColor;
        this._highlightedFaceColor = value;
        this.requestUpdate('highlightedFaceColor', oldValue);
    }

    get highlightedFaceColor() {
        return this._highlightedFaceColor;
    }

    set highlightedEdgeColor(value) {
        const oldValue = this._highlightedEdgeColor;
        this._highlightedEdgeColor = value;
        this.requestUpdate('highlightedEdgeColor', oldValue);
    }

    get highlightedEdgeColor() {
        return this._highlightedEdgeColor;
    } 

    set draggedFaceColor(value) {
        const oldValue = this._draggedFaceColor;
        this._draggedFaceColor = value;
        this.requestUpdate('draggedFaceColor', oldValue);
    }

    get draggedFaceColor() {
        return this._draggedFaceColor;
    }  
    
    set fixedFaceColor(value) {
        const oldValue = this._fixedFaceColor;
        this._fixedFaceColor = value;
        this.requestUpdate('fixedFaceColor', oldValue);
    }

    get fixedFaceColor() {
        return this._fixedFaceColor;
    }

    set editedEdgeColor(value) {
        const oldValue = this._editedEdgeColor;
        this._editedEdgeColor = value;
        this.requestUpdate('editedEdgeColor', oldValue);
    }

    get editedEdgeColor() {
        return this._editedEdgeColor;
    }    
    
    set gridSize(value) {
        const oldValue = this._gridSize;
        this._gridSize = value;
        this.requestUpdate('gridSize', oldValue);
    }

    get gridSize() {
        return this._gridSize;
    }     

    set gridTextureSize(value) {
        const oldValue = this._gridTextureSize;
        this._gridTextureSize = value;
        this.requestUpdate('gridTextureSize', oldValue);
    }

    get gridTextureSize() {
        return this._gridTextureSize;
    }     

    set gridLineWidth(value) {
        const oldValue = this._gridLineWidth;
        this._gridLineWidth = value;
        this.requestUpdate('gridLineWidth', oldValue);
    }

    get gridLineWidth() {
        return this._gridLineWidth;
    } 
   
    set unitGridVisibility(value) {
        const oldValue = this._unitGridVisibility;
        this._unitGridVisibility = value;
        this.requestUpdate('unitGridVisibility', oldValue);
    }

    get unitGridVisibility() {
        return this._unitGridVisibility;
    } 
    
    set soupViewGridTextureVisibility(value) {
        const oldValue = this._soupViewGridTextureVisibility;
        this._soupViewGridTextureVisibility = value;
        this.requestUpdate('soupViewGridTextureVisibility', oldValue);
    }

    get soupViewGridTextureVisibility() {
        return this._soupViewGridTextureVisibility;
    } 

    set optimizationDataMonitorVisibility(value) {
        const oldValue = this._optimizationDataMonitorVisibility;
        this._optimizationDataMonitorVisibility = value;
        this.requestUpdate('optimizationDataMonitorVisibility', oldValue);
    }

    get optimizationDataMonitorVisibility() {
        return this._optimizationDataMonitorVisibility;
    } 

    set solverState(value) {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            const oldValue = this._solverState;
            this._solverState = value;
            switch(this._solverState) {
                case EnumsExports.SolverState.ON:
                    this._engine.resumeSolver();
                    break;
                case EnumsExports.SolverState.OFF:
                    this._engine.pauseSolver();
                    break;                
            }
            this.requestUpdate('solverState', oldValue);
        }
    }

    get solverState() {
        return this._solverState;
    }

    set modelFilename(value) {
        if(value) {
            const oldValue = this._modelFilename;
            this._modelFilename = value;
            this._loadModel(this._modelFilename);
            this.requestUpdate('modelFilename', oldValue);
        }
    }

    get modelFilename() {
        return this._modelFilename;
    }

    set modelState(value) {
        const oldValue = this._modelState;
        this._modelState = value;
        this.requestUpdate('modelState', oldValue);
    }

    get modelState() {
        return this._modelState;
    }

    set moduleFilename(value) {
        if(value) {
            const oldValue = this._moduleFilename; 
            this._moduleFilename = value;
            this._loadModule(this._moduleFilename);
            this.requestUpdate('moduleFilename', oldValue);
        }
    }    

    get moduleFilename() {
        return this._moduleFilename;
    }

    set moduleState(value) {
        const oldValue = this._moduleState;
        this._moduleState = value;
        this.requestUpdate('moduleState', oldValue);
    }

    get moduleState() {
        return this._moduleState;
    }

    set objectiveFunctionsProperties(value) {
        const oldValue = this._objectiveFunctionsProperties;
        this._objectiveFunctionsProperties = value;
        this._modelMeshProvider.objectiveFunctionsProperties = value;
        this._soupMeshProvider.objectiveFunctionsProperties = value;
        this.requestUpdate('objectiveFunctionsProperties', oldValue);
    }

    get objectiveFunctionsProperties() {
        return this._objectiveFunctionsProperties;
    }

    set algorithmType(value) {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            const oldValue = this._algorithmType;
            this._algorithmType = value;
            this._engine.setAlgorithmType(this._algorithmType);
            this.requestUpdate('algorithmType', oldValue);
        }
    }

    get algorithmType() {
        return this._algorithmType;
    }

    /**
     * Element life-cycle callbacks
     */
    connectedCallback() {
        super.connectedCallback();

        /**
         * Face highlighting
         */
        this._meshViewFaceHighlightedSubscriptionToken = PubSub.subscribe('mesh-view-face-highlighted', (name, payload) => {
            PubSub.publish('mesh-view-highlight-face', payload);
        });

        this._meshViewFaceUnhighlightedSubscriptionToken = PubSub.subscribe('mesh-view-face-unhighlighted', (name, payload) => {
            PubSub.publish('mesh-view-unhighlight-face', payload);
        });

        /**
         * Face dragging
         */
        this._meshViewFaceDraggingBeginSubscriptionToken = PubSub.subscribe('mesh-view-face-dragging-begin', (name, payload) => {
            if(!payload.faceSelected) {
                this._engine.constrainFacePosition(payload.face.id);
            }
            PubSub.publish('mesh-view-set-dragged-face', payload);
        });

        this._meshViewFaceDraggingSubscriptionToken = PubSub.subscribe('mesh-view-face-dragging', (name, payload) => {
            this._engine.updateConstrainedFacePosition(payload.face.id, payload.offset.x, payload.offset.y);
        });

        this._meshViewFaceDraggingEndSubscriptionToken = PubSub.subscribe('mesh-view-face-dragging-end', (name, payload) => {
            if(!payload.faceSelected) {
                this._engine.unconstrainFacePosition(payload.face.id);
            }
            else {
                this._engine.reconstrainFacePosition(payload.face.id);
            }
            PubSub.publish('mesh-view-reset-dragged-face', payload);
        });

        /**
         * Face selection
         */
        this._meshViewFaceSelectedSubscriptionToken = PubSub.subscribe('mesh-view-face-selected', (name, payload) => {
            this._engine.constrainFacePosition(payload.face.id);
            PubSub.publish('mesh-view-select-face', payload);
        });

        this._meshViewFaceUnselectedSubscriptionToken = PubSub.subscribe('mesh-view-face-unselected', (name, payload) => {
            this._engine.unconstrainFacePosition(payload.face.id);
            PubSub.publish('mesh-view-unselect-face', payload);
        });
        
        /**
         * Edge selection
         */
        this._meshViewEdgeSelectedSubscriptionToken = PubSub.subscribe('mesh-view-edge-selected', (name, payload) => {
            let previousSelectedEdge = this._selectedEdge;
            if(previousSelectedEdge) {
                payload.edge = previousSelectedEdge;
                PubSub.publish('mesh-view-unselect-edge', payload);
            }

            this._selectedEdge = payload.edge;   

            let edgeAngleWeight = this._engine.getObjectiveFunctionProperty('Seamless', 'edge_angle_weight', 'none', this._selectedEdge.id);
            let edgeLengthWeight = this._engine.getObjectiveFunctionProperty('Seamless', 'edge_length_weight', 'none', this._selectedEdge.id);
            store.dispatch(ActionsExports.setSelectedEdgeSeamlessAngleWeight(edgeAngleWeight));
            store.dispatch(ActionsExports.setSelectedEdgeSeamlessLengthWeight(edgeLengthWeight));
        });

        this._meshViewEdgeUnselectedSubscriptionToken = PubSub.subscribe('mesh-view-edge-unselected', (name, payload) => {
            this._selectedEdge = null;
        });
    }
    
    disconnectedCallback() {
        // TODO: Remove event listeners
        super.disconnectedCallback();
    }

    /**
     * Private Methods
     */

    _loadModel(modelFilename) {
        if(HelpersExports.isModuleLoaded(this.moduleState)) {
            store.dispatch(ActionsExports.setModelState(EnumsExports.LoadState.LOADING));
            try {
                this._engine.loadModel(modelFilename);
                this._modelMeshProvider = new AutoquadsModelMeshProvider(this._engine, this.modelColor, this.objectiveFunctionsProperties);
                this._soupMeshProvider = new AutoquadsSoupMeshProvider(this._engine, this.soupColor, this.objectiveFunctionsProperties);  
                console.log("Model loaded: " + modelFilename);
                store.dispatch(ActionsExports.setModelState(EnumsExports.LoadState.LOADED));
            }
            catch(e) {
                this._modelMeshProvider = new MeshProvider();
                this._soupMeshProvider = new MeshProvider();   
                alert("Couldn't load model file");
                console.error("Failed to load model: " + modelFilename);
                store.dispatch(ActionsExports.setModelState(EnumsExports.LoadState.UNLOADED));
            }
        }
    }

    _loadModule(moduleFilename) {
        store.dispatch(ActionsExports.setModuleState(EnumsExports.LoadState.LOADING));
        try {
            let RDSModule = require(moduleFilename);
            this._engine = new RDSModule.Engine();
            PubSub.publish('module-loaded');
            console.log("Node module loaded: " + moduleFilename);
            store.dispatch(ActionsExports.setModuleState(EnumsExports.LoadState.LOADED));
        }
        catch(e) {
            alert("Couldn't load node module");
            console.error("Failed to node module: " + moduleFilename);
            store.dispatch(ActionsExports.setModuleState(EnumsExports.LoadState.UNLOADED));
        }
        finally {
            this._modelMeshProvider = new MeshProvider();
            this._soupMeshProvider = new MeshProvider();
        }
    }

    _reloadModule() {
        this._loadModule();
        this._loadModel(this._modelFilename);
    }
}

customElements.define('autoquads-view', AutoquadsView);