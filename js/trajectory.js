class Trajectory {
    constructor() {
        this.frames = [];
        this.topology = null;
        this.currentFrame = 0;
        this.isPlaying = false;
        this.playSpeed = 1;
        this.playInterval = null;
        this.timeStep = 0.001;
    }

    setTopology(atoms) {
        this.topology = {
            atoms: atoms,
            numAtoms: atoms.length,
            residues: this.extractResidues(atoms),
            chains: this.extractChains(atoms)
        };
    }

    extractResidues(atoms) {
        const residues = [];
        let currentRes = null;
        atoms.forEach((atom, index) => {
            const resKey = `${atom.resi}_${atom.resn}`;
            if (resKey !== currentRes) {
                currentRes = resKey;
                residues.push({
                    resi: atom.resi,
                    resn: atom.resn,
                    chain: atom.chain,
                    startIndex: index,
                    endIndex: index,
                    atoms: []
                });
            }
            residues[residues.length - 1].endIndex = index;
            residues[residues.length - 1].atoms.push(index);
        });
        return residues;
    }

    extractChains(atoms) {
        const chains = {};
        atoms.forEach((atom, index) => {
            const chainId = atom.chain || 'A';
            if (!chains[chainId]) {
                chains[chainId] = {
                    id: chainId,
                    startIndex: index,
                    endIndex: index,
                    atoms: []
                };
            }
            chains[chainId].endIndex = index;
            chains[chainId].atoms.push(index);
        });
        return Object.values(chains);
    }

    addFrame(coordinates) {
        this.frames.push(coordinates);
    }

    getFrame(index) {
        if (index < 0 || index >= this.frames.length) {
            return null;
        }
        return this.frames[index];
    }

    getCurrentFrame() {
        return this.getFrame(this.currentFrame);
    }

    nextFrame() {
        if (this.currentFrame < this.frames.length - 1) {
            this.currentFrame++;
            return true;
        }
        return false;
    }

    previousFrame() {
        if (this.currentFrame > 0) {
            this.currentFrame--;
            return true;
        }
        return false;
    }

    seekFrame(index) {
        index = Math.max(0, Math.min(this.frames.length - 1, parseInt(index)));
        this.currentFrame = index;
    }

    getNumFrames() {
        return this.frames.length;
    }

    getNumAtoms() {
        return this.topology ? this.topology.numAtoms : 0;
    }

    getNumResidues() {
        return this.topology ? this.topology.residues.length : 0;
    }

    getSimulationTime() {
        return (this.frames.length * this.timeStep).toFixed(3);
    }

    getAtomPositions(frameIndex, atomIndices) {
        const frame = this.getFrame(frameIndex);
        if (!frame) return null;
        
        if (atomIndices) {
            return atomIndices.map(i => frame[i]);
        }
        return frame;
    }

    getResiduePositions(frameIndex, residueIndex) {
        if (!this.topology || residueIndex >= this.topology.residues.length) {
            return null;
        }
        const residue = this.topology.residues[residueIndex];
        return this.getAtomPositions(frameIndex, residue.atoms);
    }

    getAlphaCarbonIndices() {
        if (!this.topology) return [];
        return this.topology.atoms
            .map((atom, index) => atom.atom === 'CA' ? index : -1)
            .filter(i => i >= 0);
    }

    getBackboneIndices() {
        if (!this.topology) return [];
        const backboneAtoms = ['CA', 'C', 'N', 'O'];
        return this.topology.atoms
            .map((atom, index) => backboneAtoms.includes(atom.atom) ? index : -1)
            .filter(i => i >= 0);
    }

    parsePDB(pdbText) {
        const atoms = [];
        const lines = pdbText.split('\n');
        
        lines.forEach(line => {
            if (line.startsWith('ATOM') || line.startsWith('HETATM')) {
                const atom = {
                    record: line.substring(0, 6).trim(),
                    serial: parseInt(line.substring(6, 11)),
                    atom: line.substring(12, 16).trim(),
                    altLoc: line.substring(16, 17).trim(),
                    resn: line.substring(17, 20).trim(),
                    chain: line.substring(21, 22).trim() || 'A',
                    resi: parseInt(line.substring(22, 26)),
                    icode: line.substring(26, 27).trim(),
                    x: parseFloat(line.substring(30, 38)),
                    y: parseFloat(line.substring(38, 46)),
                    z: parseFloat(line.substring(46, 54)),
                    occupancy: parseFloat(line.substring(54, 60)) || 1.0,
                    tempFactor: parseFloat(line.substring(60, 66)) || 0.0,
                    element: line.substring(76, 78).trim() || this.guessElement(line.substring(12, 16).trim()),
                    charge: line.substring(78, 80).trim()
                };
                atoms.push(atom);
            }
        });

        this.setTopology(atoms);
        const coords = atoms.map(a => [a.x, a.y, a.z]);
        this.addFrame(coords);
        
        return atoms.length;
    }

    guessElement(atomName) {
        const firstChar = atomName.charAt(0).toUpperCase();
        const elementMap = {
            'C': 'C', 'N': 'N', 'O': 'O', 'S': 'S', 'P': 'P',
            'H': 'H', 'F': 'F', 'Cl': 'Cl', 'Br': 'Br', 'I': 'I',
            'Fe': 'Fe', 'Zn': 'Zn', 'Mg': 'Mg', 'Ca': 'Ca',
            'Na': 'Na', 'K': 'K', 'Cl': 'Cl'
        };
        return elementMap[firstChar] || 'C';
    }

    parseXYZ(xyzText) {
        const lines = xyzText.trim().split('\n');
        const numAtoms = parseInt(lines[0]);
        const atoms = [];
        const coords = [];

        for (let i = 2; i < 2 + numAtoms && i < lines.length; i++) {
            const parts = lines[i].trim().split(/\s+/);
            if (parts.length >= 4) {
                const element = parts[0];
                const x = parseFloat(parts[1]);
                const y = parseFloat(parts[2]);
                const z = parseFloat(parts[3]);
                atoms.push({
                    serial: i - 1,
                    atom: element,
                    element: element,
                    resn: 'UNK',
                    chain: 'A',
                    resi: 1,
                    x, y, z
                });
                coords.push([x, y, z]);
            }
        }

        this.setTopology(atoms);
        this.addFrame(coords);
        return atoms.length;
    }

    generateSampleTrajectory(numFrames = 50) {
        const numAtoms = 304;
        const atoms = this.generateProteinTopology(numAtoms);
        this.setTopology(atoms);
        
        const baseCoords = atoms.map(a => [a.x, a.y, a.z]);
        this.addFrame([...baseCoords]);

        for (let f = 1; f < numFrames; f++) {
            const t = f / numFrames;
            const frameCoords = baseCoords.map((coord, i) => {
                const atom = atoms[i];
                const fluctuation = Math.sin(t * Math.PI * 2 + i * 0.1) * 0.5;
                const drift = t * 0.3;
                const domainMotion = atom.resi > 10 ? Math.sin(t * Math.PI) * 1.5 : 0;
                
                return [
                    coord[0] + fluctuation + (Math.random() - 0.5) * 0.3,
                    coord[1] + fluctuation * 0.7 + drift + domainMotion + (Math.random() - 0.5) * 0.3,
                    coord[2] + (Math.random() - 0.5) * 0.3
                ];
            });
            this.addFrame(frameCoords);
        }

        return numFrames;
    }

    generateProteinTopology(numAtoms) {
        const atoms = [];
        const residues = [
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
            'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
            'THR', 'TRP', 'TYR', 'VAL'
        ];
        
        let atomIndex = 1;
        let resIndex = 1;
        
        const numResidues = Math.floor(numAtoms / 8);
        
        for (let r = 0; r < numResidues; r++) {
            const resn = residues[r % residues.length];
            const chain = r < numResidues / 2 ? 'A' : 'B';
            
            const phi = (r / numResidues) * Math.PI * 4;
            const radius = 10 + Math.sin(r * 0.3) * 3;
            const x = Math.cos(phi) * radius;
            const y = Math.sin(phi) * radius;
            const z = r * 1.5 - numResidues * 0.75;
            
            atoms.push({
                serial: atomIndex++, atom: 'N', element: 'N',
                resn, chain, resi: resIndex,
                x: x + Math.cos(phi + 0.5) * 0.5,
                y: y + Math.sin(phi + 0.5) * 0.5,
                z: z + 0.5
            });
            
            atoms.push({
                serial: atomIndex++, atom: 'CA', element: 'C',
                resn, chain, resi: resIndex,
                x, y, z
            });
            
            atoms.push({
                serial: atomIndex++, atom: 'C', element: 'C',
                resn, chain, resi: resIndex,
                x: x + Math.cos(phi - 0.5) * 0.8,
                y: y + Math.sin(phi - 0.5) * 0.8,
                z: z - 0.2
            });
            
            atoms.push({
                serial: atomIndex++, atom: 'O', element: 'O',
                resn, chain, resi: resIndex,
                x: x + Math.cos(phi - 1) * 1.2,
                y: y + Math.sin(phi - 1) * 1.2,
                z: z - 0.5
            });
            
            const sideChainAtoms = ['CB', 'CG', 'CD'];
            for (let s = 0; s < Math.min(sideChainAtoms.length, 4); s++) {
                const angle = phi + Math.PI + s * 0.8;
                const dist = 1.5 + s * 1.2;
                atoms.push({
                    serial: atomIndex++, atom: sideChainAtoms[s], element: 'C',
                    resn, chain, resi: resIndex,
                    x: x + Math.cos(angle) * dist,
                    y: y + Math.sin(angle) * dist,
                    z: z + 0.3 + s * 0.2
                });
            }
            
            resIndex++;
        }
        
        return atoms;
    }

    clear() {
        this.frames = [];
        this.topology = null;
        this.currentFrame = 0;
        this.stop();
    }

    play(onFrameChange) {
        if (this.isPlaying) return;
        this.isPlaying = true;
        
        const baseInterval = 100;
        const interval = baseInterval / this.playSpeed;
        
        this.playInterval = setInterval(() => {
            if (!this.nextFrame()) {
                this.currentFrame = 0;
            }
            if (onFrameChange) {
                onFrameChange(this.currentFrame);
            }
        }, interval);
    }

    pause() {
        this.isPlaying = false;
        if (this.playInterval) {
            clearInterval(this.playInterval);
            this.playInterval = null;
        }
    }

    stop() {
        this.pause();
        this.currentFrame = 0;
    }

    setSpeed(speed) {
        this.playSpeed = parseFloat(speed);
        if (this.isPlaying) {
            const wasPlaying = this.isPlaying;
            this.pause();
            if (wasPlaying) {
                this.play(null);
            }
        }
    }

    getSelectionIndices(selection) {
        if (!this.topology) return [];
        
        if (selection === 'all') {
            return this.topology.atoms.map((_, i) => i);
        }
        
        if (selection === 'protein') {
            return this.topology.atoms
                .map((atom, i) => atom.record === 'ATOM' ? i : -1)
                .filter(i => i >= 0);
        }
        
        if (selection === 'backbone') {
            return this.getBackboneIndices();
        }
        
        if (selection === 'calpha') {
            return this.getAlphaCarbonIndices();
        }
        
        if (selection === 'water') {
            return this.topology.atoms
                .map((atom, i) => atom.resn === 'HOH' || atom.resn === 'WAT' ? i : -1)
                .filter(i => i >= 0);
        }
        
        return [];
    }
}

class TrajectoryComparison {
    constructor() {
        this.trajectories = [];
        this.colors = [
            '#6366f1', '#10b981', '#f59e0b', '#ef4444',
            '#8b5cf6', '#06b6d4', '#ec4899', '#84cc16'
        ];
    }

    addTrajectory(trajectory, name) {
        const color = this.colors[this.trajectories.length % this.colors.length];
        this.trajectories.push({
            trajectory,
            name: name || `Trajectory ${this.trajectories.length + 1}`,
            color
        });
        return this.trajectories.length - 1;
    }

    removeTrajectory(index) {
        this.trajectories.splice(index, 1);
    }

    getNumTrajectories() {
        return this.trajectories.length;
    }

    getTrajectory(index) {
        return this.trajectories[index];
    }

    getAllTrajectories() {
        return this.trajectories;
    }

    clear() {
        this.trajectories = [];
    }
}
