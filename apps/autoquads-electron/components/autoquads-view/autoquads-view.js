// Web Modules Imports
import { LitElement, html, css } from '../../web_modules/lit-element.js';
import '../../web_modules/@vaadin/vaadin-button.js';
import '../../web_modules/@vaadin/vaadin-split-layout.js';
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
        `];
    }

    render() {
        return html`  
            <vaadin-split-layout class="outer-split">
                <autoquads-side-bar></autoquads-side-bar>
                <vaadin-split-layout orientation="horizontal" class="inner-split">
                    <mesh-view 
                        id="model-mesh-view"
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
                        show-wireframe="${this.showWireframe}"
                        background-color="${this.modelViewportColor}"
                        mesh-color="${this.modelColor}"
                        .meshProvider="${this.modelMeshProvider}"
                        mesh-interaction="${this.meshInteraction}"
                        highlighted-face-color="${this.highlightedFaceColor}"
                        dragged-face-color="${this.draggedFaceColor}"
                        selected-face-color="${this.fixedFaceColor}"
                        show-grid-texture>
                    </mesh-view>
                    <mesh-view 
                        id="soup-mesh-view"
                        enable-face-dragging caption="Soup View"
                        show-grid="${this.showUnitGrid}"
                        grid-horizontal-color="${this.gridHorizontalColor}"
                        grid-vertical-color="${this.gridVerticalColor}"
                        grid-background-color1="${this.gridBackgroundColor1}"
                        grid-background-color2="${this.gridBackgroundColor2}"
                        grid-size="${this.gridSize}"
                        grid-texture-size="${this.gridTextureSize}"
                        grid-line-width="${this.gridLineWidth}"
                        show-wireframe="${this.showWireframe}"
                        background-color="${this.soupViewportColor}"
                        mesh-color="${this.soupColor}"
                        .meshProvider="${this.soupMeshProvider}"
                        mesh-interaction="${this.meshInteraction}"
                        highlighted-face-color="${this.highlightedFaceColor}"
                        dragged-face-color="${this.draggedFaceColor}"
                        selected-face-color="${this.fixedFaceColor}"
                        show-debug-data="${this.showOptimizationDataMonitor}"
                        show-grid-texture="${this.showGridTextureInSoupView}">
                    </mesh-view>
                </vaadin-split-layout>
            </vaadin-split-layout>
        `;
    }

    static get properties() {
        return {
            modelMeshProvider: {
                type: Object,
                attribute: 'model-mesh-provider'
            },
            soupMeshProvider: {
                type: Object,
                attribute: 'soup-mesh-provider'
            },            
            soupViewportColor: {
                type: String,
                attribute: 'soup-viewport-color'
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
            showWireframe: {
                type: Boolean,
                attribute: 'show-wireframe'
            },
            showMeshView: {
                type: Boolean,
                attribute: 'show-mesh-view'
            },
            showSoupView: {
                type: Boolean,
                attribute: 'show-soup-view'
            },
            showSoupGrid: {
                type: Boolean,
                attribute: 'show-soup-grid'
            },
            delta: {
                type: Number,
                attribute: 'delta'
            },
            lambda: {
                type: Number,
                attribute: 'lambda'
            },
            seamlessWeight: {
                type: Number,
                attribute: 'seamless-weight'
            },
            positionWeight: {
                type: Number,
                attribute: 'position-weight'
            },
            splitOrientation: {
                type: String,
                attribute: 'split-orientation'
            },
            solver: {
                type: Boolean,
                attribute: 'solver'
            },
            editingTool: {
                type: String,
                attribute: 'editing-tool'
            },
            meshInteraction: {
                type: String,
                attribute: 'mesh-interaction'
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
            draggedFaceColor: {
                type: String,
                attribute: 'dragged-face-color'
            },
            fixedFaceColor: {
                type: String,
                attribute: 'fixed-face-color'
            },
            vertexEnergyColor: {
                type: Boolean,
                attribute: 'vertex-energy-color'
            },
            vertexEnergyType: {
                type: String,
                attribute: 'vertex-energy-type'
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
            showOptimizationDataMonitor: {
                type: Boolean,
                attribute: 'show-optimization-data-monitor'
            },
            showUnitGrid: {
                type: Boolean,
                attribute: 'show-unit-grid'
            },
            showGridTextureInSoupView: {
                type: Boolean,
                attribute: 'show-grid-texture-in-soup-view'
            },
            modelFilename: {
                type: String,
                attribute: 'model-filename'
            }
        };
    }

    stateChanged(state) {
        this.splitOrientation = state.splitOrientation;
        this.modelViewportColor = state.modelViewportColor;
        this.soupViewportColor = state.soupViewportColor;
        this.modelColor = state.modelColor;
        this.soupColor = state.soupColor;
        this.showWireframe = ReducerExports.isVisible(state.wireframeVisibility);
        this.showMeshView = ReducerExports.isVisible(state.modelViewVisibility);
        this.showSoupView = ReducerExports.isVisible(state.soupViewVisibility);
        this.gridHorizontalColor = state.gridHorizontalColor;
        this.gridVerticalColor = state.gridVerticalColor;
        this.gridBackgroundColor1 = state.gridBackgroundColor1;
        this.gridBackgroundColor2 = state.gridBackgroundColor2;
        this.highlightedFaceColor = state.highlightedFaceColor;
        this.draggedFaceColor = state.draggedFaceColor;
        this.fixedFaceColor = state.fixedFaceColor;
        this.gridSize = state.gridSize;
        this.gridTextureSize = state.gridTextureSize;
        this.gridLineWidth = state.gridLineWidth;
        this.modelFilename = state.modelFilename;
    }    

    constructor() {
        super();
        this.modelMeshProvider = new MeshProvider();
        this.soupMeshProvider = new MeshProvider();           
    }

    firstUpdated() {
        this._loadModule();        
    }

    connectedCallback() {
        super.connectedCallback();
    }
    
    disconnectedCallback() {
        // TODO: Remove event listeners
        super.disconnectedCallback();
    }

    /**
     * Properties
     */

    set vertexEnergyType(value) {
        const oldValue = this._vertexEnergyType;
        this._vertexEnergyType = value;
        
        if(this.modelMeshProvider) {
            this.modelMeshProvider.vertexEnergyType = vertexEnergyType;
        }

        if(this.soupMeshProvider) {
            this.soupMeshProvider.vertexEnergyType = vertexEnergyType;
        }

        this.requestUpdate('vertexEnergyType', oldValue);
    }

    get vertexEnergyType() {
        return this._vertexEnergyType;
    }  
    
    set vertexEnergyColor(value) {
        const oldValue = this._vertexEnergyColor;
        this._vertexEnergyColor = value;
        
        if(this.modelMeshProvider) {
            this.modelMeshProvider.vertexEnergyColor = vertexEnergyColor;
        }

        if(this.soupMeshProvider) {
            this.soupMeshProvider.vertexEnergyColor = vertexEnergyColor;
        }

        this.requestUpdate('vertexEnergyColor', oldValue);
    }

    get vertexEnergyColor() {
        return this._vertexEnergyColor;
    }  
    
    set modelColor(value) {
        const oldValue = this._modelColor;
        this._modelColor = value;
        
        if(this.modelMeshProvider) {
            this.modelMeshProvider.meshColor = this._modelColor;
        }

        this.requestUpdate('modelColor', oldValue);
    }

    get modelColor() {
        return this._modelColor;
    }  
    
    set soupColor(value) {
        const oldValue = this._soupColor;
        this._soupColor = value;
        
        if(this.soupMeshProvider) {
            this.soupMeshProvider.meshColor = this._soupColor;
        }

        this.requestUpdate('soupColor', oldValue);
    }

    get soupColor() {
        return this._soupColor;
    }

    set modelFilename(value) {
        const oldValue = this._modelFilename;
        this._modelFilename = value;
        if(this._modelFilename !== oldValue) {
            this._reloadModel(this._modelFilename);
        }
        this.requestUpdate('modelFilename', oldValue);
    }

    get modelFilename() {
        return this._modelFilename;
    }      
    
    /**
     * Private Methods
     */

    _reloadModel(modelFilename) {
        if(this._engine && this.modelFilename) {
            this._engine.loadModel(modelFilename);   
            this.modelMeshProvider = new AutoquadsModelMeshProvider(this._engine, this.vertexEnergyType, this.vertexEnergyColor, this.modelColor);
            this.soupMeshProvider = new AutoquadsSoupMeshProvider(this._engine, this.vertexEnergyType, this.vertexEnergyColor, this.soupColor);     
        }
    }

    _loadModule() {
        const { join } = require('path');
        let RDSModule = require(join(appRoot, 'node-addon.node'));
        this._engine = new RDSModule.Engine();
        this.modelMeshProvider = new MeshProvider();
        this.soupMeshProvider = new MeshProvider();      
    }

    _reloadModule() {
        this._loadModule();
        this._reloadModel(this._modelFilename);
    }

    _lambdaChanged(lambda) {
        if (!isNaN(lambda)) {
            // this.autoquads.lambda = lambda;
        }
    }

    _deltaChanged(delta) {
        if (!isNaN(delta)) {
            // this.autoquads.delta = delta;
        }
    }

    _integerWeightChanged(integerWeight) {
        if (!isNaN(integerWeight)) {
            // this.autoquads.integerWeight = integerWeight;
        }
    }

    _integerSpacingChanged(integerSpacing) {
        if (!isNaN(integerSpacing)) {
            // this.autoquads.integerSpacing = integerSpacing;
        }
    }

    _seamlessWeightChanged(seamlessWeight) {
        if (!isNaN(seamlessWeight)) {
            // this.autoquads.seamlessWeight = seamlessWeight;
        }
    }

    _positionWeightChanged(positionWeight) {
        if (!isNaN(positionWeight)) {
            // this.autoquads.positionWeight = positionWeight;
        }
    }

    _editingToolChanged(editingTool) {

    }

    _vertexEnergyColorChanged(vertexEnergyColor) {
        this._vertexEnergyColor = new THREE.Color(vertexEnergyColor);
        this._modelMeshProvider.energyColor = this._vertexEnergyColor;
        this._soupMeshProvider.energyColor = this._vertexEnergyColor;
    }

    _vertexEnergyChanged(vertexEnergy) {
        this._modelMeshProvider.vertexEnergy = vertexEnergy;
        this._soupMeshProvider.vertexEnergy = vertexEnergy;
    }

    _modelColorChanged(modelColor) {
        this._modelColor = new THREE.Color(modelColor);
        // this.modelMeshProvider.meshColor = this._modelColor;
    }

    _solverColorChanged(solverColor) {
        this._solverColor = new THREE.Color(solverColor);
        // this.solverMeshProvider.meshColor = this._solverColor;
    }  
}

customElements.define('autoquads-view', AutoquadsView);