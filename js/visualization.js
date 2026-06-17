class MolecularViewer {
    constructor(containerId) {
        this.container = document.getElementById(containerId);
        this.viewer = null;
        this.currentStyle = 'cartoon';
        this.currentColorScheme = 'chain';
        this.glviewer = null;
        this.model = null;
        this.init();
    }

    init() {
        if (typeof $3Dmol !== 'undefined') {
            const config = {
                backgroundColor: 'black',
                antialias: true,
                quality: 'high'
            };
            
            this.glviewer = $3Dmol.createViewer(this.container, config);
            this.viewer = this.glviewer;
            
            this.glviewer.setBackgroundColor(0x0f172a);
            this.glviewer.zoomTo();
            this.glviewer.render();
        }
    }

    loadPDB(pdbText) {
        if (!this.glviewer) return;
        
        this.glviewer.clear();
        this.glviewer.addModel(pdbText, 'pdb');
        this.applyStyle();
        this.glviewer.zoomTo();
        this.glviewer.render();
    }

    loadXYZ(xyzText) {
        if (!this.glviewer) return;
        
        this.glviewer.clear();
        this.glviewer.addModel(xyzText, 'xyz');
        this.applyStyle();
        this.glviewer.zoomTo();
        this.glviewer.render();
    }

    setStyle(style) {
        this.currentStyle = style;
        this.applyStyle();
    }

    setColorScheme(scheme) {
        this.currentColorScheme = scheme;
        this.applyStyle();
    }

    applyStyle() {
        if (!this.glviewer) return;

        const colorConfig = this.getColorConfig();
        
        this.glviewer.setStyle({}, {});

        switch (this.currentStyle) {
            case 'cartoon':
                this.glviewer.setStyle({}, {
                    cartoon: {
                        color: colorConfig,
                        thickness: 0.4,
                        arrows: true,
                        arrowThickness: 0.3
                    }
                });
                break;
                
            case 'stick':
                this.glviewer.setStyle({}, {
                    stick: {
                        colorscheme: colorConfig,
                        radius: 0.15,
                        opacity: 1.0
                    }
                });
                break;
                
            case 'sphere':
                this.glviewer.setStyle({}, {
                    sphere: {
                        colorscheme: colorConfig,
                        radius: 1.2,
                        opacity: 1.0
                    }
                });
                break;
                
            case 'surface':
                this.glviewer.setStyle({}, {
                    cartoon: {
                        color: colorConfig,
                        opacity: 0.8
                    }
                });
                this.glviewer.addSurface($3Dmol.SurfaceType.VDW, {
                    opacity: 0.5,
                    colorscheme: colorConfig
                }, {});
                break;
                
            case 'line':
                this.glviewer.setStyle({}, {
                    line: {
                        colorscheme: colorConfig,
                        linewidth: 1,
                        opacity: 1.0
                    }
                });
                break;
                
            case 'cross':
                this.glviewer.setStyle({}, {
                    cross: {
                        colorscheme: colorConfig,
                        radius: 0.1,
                        scale: 0.3
                    }
                });
                break;
                
            default:
                this.glviewer.setStyle({}, {
                    cartoon: {
                        color: colorConfig,
                        thickness: 0.4
                    }
                });
        }

        this.glviewer.render();
    }

    getColorConfig() {
        switch (this.currentColorScheme) {
            case 'chain':
                return 'chain';
            case 'residue':
                return 'residue';
            case 'atom':
                return 'default';
            case 'secondary':
                return 'ss';
            default:
                return 'chain';
        }
    }

    updateFrame(coordinates) {
        if (!this.glviewer || !coordinates) return;
        
        const model = this.glviewer.getModel();
        if (!model) return;

        const atoms = model.selectedAtoms({});
        for (let i = 0; i < atoms.length && i < coordinates.length; i++) {
            atoms[i].x = coordinates[i][0];
            atoms[i].y = coordinates[i][1];
            atoms[i].z = coordinates[i][2];
        }

        this.glviewer.updateStyle();
        this.glviewer.render();
    }

    zoomTo() {
        if (this.glviewer) {
            this.glviewer.zoomTo();
            this.glviewer.render();
        }
    }

    setBackgroundColor(color) {
        if (this.glviewer) {
            this.glviewer.setBackgroundColor(color);
            this.glviewer.render();
        }
    }

    addSphere(position, radius, color, opacity = 0.8) {
        if (!this.glviewer) return;
        
        const sphere = this.glviewer.addSphere({
            center: { x: position[0], y: position[1], z: position[2] },
            radius: radius,
            color: color,
            alpha: opacity
        });
        
        this.glviewer.render();
        return sphere;
    }

    addArrow(start, end, color, radius = 0.1) {
        if (!this.glviewer) return;
        
        const arrow = this.glviewer.addArrow({
            start: { x: start[0], y: start[1], z: start[2] },
            end: { x: end[0], y: end[1], z: end[2] },
            color: color,
            radius: radius
        });
        
        this.glviewer.render();
        return arrow;
    }

    clearShapes() {
        if (this.glviewer) {
            this.glviewer.clearShapes();
            this.glviewer.render();
        }
    }

    clear() {
        if (this.glviewer) {
            this.glviewer.clear();
            this.glviewer.render();
        }
    }

    resize() {
        if (this.glviewer) {
            this.glviewer.resize();
            this.glviewer.render();
        }
    }

    getViewer() {
        return this.glviewer;
    }

    addHBondVisualization(donorIdx, acceptorIdx, frameCoords) {
        if (!this.glviewer || !frameCoords) return;
        
        const donor = frameCoords[donorIdx];
        const acceptor = frameCoords[acceptorIdx];
        
        if (!donor || !acceptor) return;
        
        const midpoint = [
            (donor[0] + acceptor[0]) / 2,
            (donor[1] + acceptor[1]) / 2,
            (donor[2] + acceptor[2]) / 2
        ];
        
        this.glviewer.addDashedLine({
            start: { x: donor[0], y: donor[1], z: donor[2] },
            end: { x: acceptor[0], y: acceptor[1], z: acceptor[2] },
            color: '#ff6b6b',
            dashed: true,
            radius: 0.1
        });
        
        this.glviewer.render();
    }

    highlightSelection(atomIndices) {
        if (!this.glviewer || !atomIndices || atomIndices.length === 0) return;
        
        this.glviewer.setStyle({}, {
            [this.currentStyle]: {
                opacity: 0.3
            }
        });
        
        atomIndices.forEach(idx => {
            this.glviewer.setStyle({ serial: idx + 1 }, {
                [this.currentStyle]: {
                    opacity: 1.0,
                    color: 'yellow'
                }
            });
        });
        
        this.glviewer.render();
    }

    resetHighlight() {
        this.applyStyle();
    }

    render() {
        if (this.glviewer) {
            this.glviewer.render();
        }
    }

    startPCModeAnimation(pcFrames, atomIndices) {
        if (!this.glviewer || !pcFrames || pcFrames.length === 0) return;

        this.pcModeActive = true;
        this.pcFrames = pcFrames;
        this.pcAtomIndices = atomIndices;
        this.pcFrameIndex = 0;
        this.pcDirection = 1;
        this.pcAnimationId = null;

        this.originalCoords = [];
        const model = this.glviewer.getModel();
        if (model) {
            const atoms = model.selectedAtoms({});
            for (let i = 0; i < atoms.length; i++) {
                this.originalCoords.push([atoms[i].x, atoms[i].y, atoms[i].z]);
            }
        }

        this.animatePCMode();
    }

    animatePCMode() {
        if (!this.pcModeActive) return;

        const frame = this.pcFrames[this.pcFrameIndex];
        if (frame && this.pcAtomIndices) {
            const model = this.glviewer.getModel();
            if (model) {
                const atoms = model.selectedAtoms({});
                
                for (let i = 0; i < this.pcAtomIndices.length; i++) {
                    const atomIdx = this.pcAtomIndices[i];
                    if (atomIdx < atoms.length && i < frame.length) {
                        atoms[atomIdx].x = frame[i][0];
                        atoms[atomIdx].y = frame[i][1];
                        atoms[atomIdx].z = frame[i][2];
                    }
                }
                
                this.glviewer.updateStyle();
                this.glviewer.render();
            }
        }

        this.pcFrameIndex += this.pcDirection;
        if (this.pcFrameIndex >= this.pcFrames.length - 1) {
            this.pcDirection = -1;
        } else if (this.pcFrameIndex <= 0) {
            this.pcDirection = 1;
        }

        this.pcAnimationId = setTimeout(() => this.animatePCMode(), 80);
    }

    stopPCModeAnimation() {
        this.pcModeActive = false;
        
        if (this.pcAnimationId) {
            clearTimeout(this.pcAnimationId);
            this.pcAnimationId = null;
        }

        if (this.originalCoords && this.glviewer) {
            const model = this.glviewer.getModel();
            if (model) {
                const atoms = model.selectedAtoms({});
                for (let i = 0; i < this.originalCoords.length && i < atoms.length; i++) {
                    atoms[i].x = this.originalCoords[i][0];
                    atoms[i].y = this.originalCoords[i][1];
                    atoms[i].z = this.originalCoords[i][2];
                }
                this.glviewer.updateStyle();
                this.glviewer.render();
            }
        }

        this.originalCoords = null;
        this.pcFrames = null;
    }

    addDisplacementArrows(displacements, color = '#ff6b6b') {
        if (!this.glviewer || !displacements) return;

        const step = Math.max(1, Math.floor(displacements.length / 30));
        
        for (let i = 0; i < displacements.length; i += step) {
            const d = displacements[i];
            const dist = this.distance3D(d.start, d.end);
            if (dist > 0.1) {
                this.glviewer.addArrow({
                    start: { x: d.start[0], y: d.start[1], z: d.start[2] },
                    end: { x: d.end[0], y: d.end[1], z: d.end[2] },
                    color: color,
                    radius: 0.1,
                    radiusRatio: 1.5
                });
            }
        }
        
        this.glviewer.render();
    }

    distance3D(a, b) {
        const dx = a[0] - b[0];
        const dy = a[1] - b[1];
        const dz = a[2] - b[2];
        return Math.sqrt(dx * dx + dy * dy + dz * dz);
    }
}

class AnimationController {
    constructor(viewer, trajectory) {
        this.viewer = viewer;
        this.trajectory = trajectory;
        this.isPlaying = false;
        this.currentFrame = 0;
        this.playSpeed = 1;
        this.animationId = null;
        this.onFrameChange = null;
    }

    setTrajectory(trajectory) {
        this.trajectory = trajectory;
        this.currentFrame = 0;
        this.updateView();
    }

    play() {
        if (this.isPlaying || !this.trajectory) return;
        
        this.isPlaying = true;
        this.animate();
    }

    pause() {
        this.isPlaying = false;
        if (this.animationId) {
            cancelAnimationFrame(this.animationId);
            this.animationId = null;
        }
    }

    toggle() {
        if (this.isPlaying) {
            this.pause();
        } else {
            this.play();
        }
    }

    animate() {
        if (!this.isPlaying) return;
        
        this.nextFrame();
        
        const delay = 100 / this.playSpeed;
        this.animationId = setTimeout(() => {
            if (this.isPlaying) {
                this.animate();
            }
        }, delay);
    }

    nextFrame() {
        if (!this.trajectory) return;
        
        if (this.currentFrame < this.trajectory.getNumFrames() - 1) {
            this.currentFrame++;
        } else {
            this.currentFrame = 0;
        }
        this.updateView();
    }

    previousFrame() {
        if (!this.trajectory) return;
        
        if (this.currentFrame > 0) {
            this.currentFrame--;
        } else {
            this.currentFrame = this.trajectory.getNumFrames() - 1;
        }
        this.updateView();
    }

    seekFrame(frameIndex) {
        if (!this.trajectory) return;
        
        frameIndex = Math.max(0, Math.min(this.trajectory.getNumFrames() - 1, parseInt(frameIndex)));
        this.currentFrame = frameIndex;
        this.updateView();
    }

    updateView() {
        if (!this.viewer || !this.trajectory) return;
        
        const coords = this.trajectory.getFrame(this.currentFrame);
        if (coords) {
            this.viewer.updateFrame(coords);
        }
        
        if (this.onFrameChange) {
            this.onFrameChange(this.currentFrame);
        }
    }

    setSpeed(speed) {
        this.playSpeed = parseFloat(speed);
    }

    getCurrentFrame() {
        return this.currentFrame;
    }

    getNumFrames() {
        return this.trajectory ? this.trajectory.getNumFrames() : 0;
    }

    setOnFrameChange(callback) {
        this.onFrameChange = callback;
    }
}
