class Trajectory {
    constructor() {
        this.frames = [];
        this.topology = null;
        this.currentFrame = 0;
        this.isPlaying = false;
        this.playSpeed = 1;
        this.playInterval = null;
        this.timeStep = 0.001;
        this.keepTopology = false;
        this.summary = null;
        this.frameTimes = [];
    }

    setTopology(atoms) {
        if (this.keepTopology && this.topology) {
            return;
        }
        this.topology = {
            atoms: atoms,
            numAtoms: atoms.length,
            residues: this.extractResidues(atoms),
            chains: this.extractChains(atoms)
        };
    }

    setFramesOnly(frames) {
        if (!this.topology) {
            throw new Error('请先加载拓扑文件');
        }
        const expectedAtoms = this.topology.numAtoms;
        for (let i = 0; i < frames.length; i++) {
            if (frames[i].length !== expectedAtoms) {
                throw new Error(`第${i}帧原子数(${frames[i].length})与拓扑(${expectedAtoms})不匹配`);
            }
        }
        this.frames = frames;
        this.currentFrame = 0;
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

    getFrameTime(index) {
        const idx = index !== undefined ? index : this.currentFrame;
        if (this.frameTimes && this.frameTimes[idx] !== undefined) {
            return this.frameTimes[idx];
        }
        return idx * (this.summary?.timeStep || this.timeStep || 0.001);
    }

    getSummary() {
        if (this.summary) return this.summary;
        if (!this.topology) return null;
        return {
            format: '内部数据',
            numAtoms: this.topology.numAtoms,
            numFrames: this.frames.length,
            timeStep: this.timeStep,
            hasUnitCell: false
        };
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

    parseMultiFramePDB(pdbText) {
        const lines = pdbText.split('\n');
        const frames = [];
        let currentFrameCoords = [];
        let topologyAtoms = [];
        let inModel = false;
        let firstModel = true;
        let hasModels = false;

        for (const line of lines) {
            if (line.startsWith('MODEL')) {
                hasModels = true;
                inModel = true;
                currentFrameCoords = [];
                continue;
            }
            if (line.startsWith('ENDMDL')) {
                inModel = false;
                if (currentFrameCoords.length > 0) {
                    frames.push(currentFrameCoords);
                }
                firstModel = false;
                continue;
            }
            if ((line.startsWith('ATOM') || line.startsWith('HETATM')) && (inModel || !hasModels)) {
                const x = parseFloat(line.substring(30, 38));
                const y = parseFloat(line.substring(38, 46));
                const z = parseFloat(line.substring(46, 54));
                currentFrameCoords.push([x, y, z]);

                if (firstModel || !hasModels) {
                    const atom = {
                        record: line.substring(0, 6).trim(),
                        serial: parseInt(line.substring(6, 11)),
                        atom: line.substring(12, 16).trim(),
                        altLoc: line.substring(16, 17).trim(),
                        resn: line.substring(17, 20).trim(),
                        chain: line.substring(21, 22).trim() || 'A',
                        resi: parseInt(line.substring(22, 26)),
                        icode: line.substring(26, 27).trim(),
                        x, y, z,
                        occupancy: parseFloat(line.substring(54, 60)) || 1.0,
                        tempFactor: parseFloat(line.substring(60, 66)) || 0.0,
                        element: line.substring(76, 78).trim() || this.guessElement(line.substring(12, 16).trim()),
                        charge: line.substring(78, 80).trim()
                    };
                    topologyAtoms.push(atom);
                }
            }
        }

        if (!hasModels && currentFrameCoords.length > 0) {
            frames.push(currentFrameCoords);
        }

        if (topologyAtoms.length === 0) {
            return 0;
        }

        const frameTimes = frames.map((_, i) => i * 0.001);
        const summary = {
            format: 'PDB多帧',
            numAtoms: topologyAtoms.length,
            numFrames: frames.length,
            timeStep: 0.001,
            hasUnitCell: false
        };

        if (this.topology && this.topology.numAtoms > 0) {
            if (frames[0] && frames[0].length !== this.topology.numAtoms) {
                throw new Error(`轨迹原子数(${frames[0].length})与拓扑(${this.topology.numAtoms})不匹配`);
            }
            this.frames = frames;
            this.frameTimes = frameTimes;
            this.summary = summary;
            this.currentFrame = 0;
            return { numAtoms: this.topology.numAtoms, numFrames: frames.length, summary };
        }

        this.setTopology(topologyAtoms);
        this.frames = frames;
        this.frameTimes = frameTimes;
        this.summary = summary;
        this.currentFrame = 0;

        return { numAtoms: topologyAtoms.length, numFrames: frames.length, summary };
    }

    parseMultiFrameXYZ(xyzText) {
        const lines = xyzText.trim().split('\n');
        const frames = [];
        let topologyAtoms = [];
        let firstFrame = true;
        let i = 0;

        while (i < lines.length) {
            const line = lines[i].trim();
            if (!line) { i++; continue; }

            const numAtoms = parseInt(line);
            if (isNaN(numAtoms) || numAtoms <= 0) { i++; continue; }

            i++;
            const commentLine = lines[i] || '';
            i++;

            const frameCoords = [];
            const frameAtoms = [];

            for (let j = 0; j < numAtoms && i < lines.length; j++, i++) {
                const parts = lines[i].trim().split(/\s+/);
                if (parts.length >= 4) {
                    const element = parts[0];
                    const x = parseFloat(parts[1]);
                    const y = parseFloat(parts[2]);
                    const z = parseFloat(parts[3]);

                    frameCoords.push([x, y, z]);

                    if (firstFrame) {
                        frameAtoms.push({
                            serial: j + 1,
                            atom: element,
                            element: element,
                            resn: 'UNK',
                            chain: 'A',
                            resi: 1,
                            x, y, z
                        });
                    }
                }
            }

            if (frameCoords.length > 0) {
                frames.push(frameCoords);
            }

            if (firstFrame && frameAtoms.length > 0) {
                topologyAtoms = frameAtoms;
            }

            firstFrame = false;
        }

        if (topologyAtoms.length === 0) {
            return 0;
        }

        const frameTimes = frames.map((_, i) => i * 0.001);
        const summary = {
            format: 'XYZ多帧',
            numAtoms: topologyAtoms.length,
            numFrames: frames.length,
            timeStep: 0.001,
            hasUnitCell: false
        };

        if (this.topology && this.topology.numAtoms > 0) {
            if (frames[0] && frames[0].length !== this.topology.numAtoms) {
                throw new Error(`轨迹原子数(${frames[0].length})与拓扑(${this.topology.numAtoms})不匹配`);
            }
            this.frames = frames;
            this.frameTimes = frameTimes;
            this.summary = summary;
            this.currentFrame = 0;
            return { numAtoms: this.topology.numAtoms, numFrames: frames.length, summary };
        }

        this.setTopology(topologyAtoms);
        this.frames = frames;
        this.frameTimes = frameTimes;
        this.summary = summary;
        this.currentFrame = 0;

        return { numAtoms: topologyAtoms.length, numFrames: frames.length, summary };
    }

    async parseDCD(arrayBuffer) {
        const view = new DataView(arrayBuffer);
        const byteLength = view.byteLength;

        if (byteLength < 100) {
            throw new Error('DCD文件过小，不是有效的DCD轨迹');
        }

        const beFirst = view.getInt32(0, false);
        const leFirst = view.getInt32(0, true);
        let littleEndian;
        if (beFirst === 84) {
            littleEndian = false;
        } else if (leFirst === 84) {
            littleEndian = true;
        } else {
            throw new Error('DCD文件格式错误：头部块大小不是84，可能不是有效的DCD文件或字节序无法识别');
        }

        const readInt = (off) => view.getInt32(off, littleEndian);
        const readFloat = (off) => view.getFloat32(off, littleEndian);
        const readDouble = (off) => view.getFloat64(off, littleEndian);

        const blockSize1 = readInt(0);
        const magic = String.fromCharCode(
            view.getUint8(4), view.getUint8(5),
            view.getUint8(6), view.getUint8(7)
        );
        if (magic !== 'CORD') {
            throw new Error(`DCD文件magic标识错误：期望"CORD"，实际"${magic}"`);
        }

        const numFramesHeader = readInt(8);
        const charmmVersion = readInt(20);
        const headerEnd = 4 + blockSize1;
        if (readInt(headerEnd) !== blockSize1) {
            throw new Error('DCD文件头部闭合标记与起始不匹配，文件可能损坏');
        }
        let offset = headerEnd + 4;

        const titleBlockSize = readInt(offset);
        offset = offset + 4 + titleBlockSize + 4;

        const natomsBlockSize = readInt(offset);
        const natoms = readInt(offset + 4);
        offset = offset + 4 + natomsBlockSize + 4;

        if (natoms <= 0 || natoms > 10000000) {
            throw new Error(`DCD文件原子数无效: ${natoms}`);
        }

        const expectedCoordBlock = natoms * 4;
        const frames = [];
        const topologyAtoms = [];
        let parseError = null;
        let hasUnitCell = false;
        let totalBoxes = 0;
        let totalNoBoxes = 0;
        const frameTimes = [];
        let frameIdx = 0;

        while (offset + 4 <= byteLength) {
            try {
                const blockStart = offset;
                const blockSize = readInt(blockStart);

                if (blockSize === 48 || blockSize === 56 || blockSize === 72) {
                    hasUnitCell = true;
                    totalBoxes++;
                    const doubles = blockSize / 8;
                    const boxData = [];
                    for (let b = 0; b < doubles && b < 9; b++) {
                        boxData.push(readDouble(blockStart + 4 + b * 8));
                    }
                    if (readInt(blockStart + 4 + blockSize) !== blockSize) {
                        if (frames.length === 0) {
                            parseError = 'DCD晶胞块闭合标记不匹配，文件可能损坏';
                        }
                        break;
                    }
                    offset = blockStart + 4 + blockSize + 4;

                } else if (blockSize === expectedCoordBlock) {
                    if (frames.length === 0) {
                        totalNoBoxes++;
                    }
                } else {
                    if (blockSize > 1 && blockSize % 4 === 0 &&
                        (blockSize <= natoms * 8) &&
                        frames.length === 0) {
                        const probe = readFloat(blockStart + 4);
                        if (isFinite(probe) && Math.abs(probe) < 100000) {
                            totalNoBoxes++;
                        } else {
                            if (frames.length === 0) {
                                parseError = `DCD帧数据格式不匹配（块大小=${blockSize}, 期望坐标块=${expectedCoordBlock}）`;
                            }
                            break;
                        }
                    } else {
                        if (frames.length === 0) {
                            parseError = `DCD帧数据格式不匹配（块大小=${blockSize}, 期望坐标块=${expectedCoordBlock}）`;
                        }
                        break;
                    }
                }

                const xCoords = [], yCoords = [], zCoords = [];
                let dimOk = true;

                for (let dim = 0; dim < 3; dim++) {
                    if (offset + 4 > byteLength) { dimOk = false; break; }
                    const openSize = readInt(offset);
                    if (openSize !== expectedCoordBlock) {
                        dimOk = false;
                        break;
                    }
                    const arr = dim === 0 ? xCoords : dim === 1 ? yCoords : zCoords;
                    if (offset + 4 + expectedCoordBlock + 4 > byteLength) {
                        dimOk = false;
                        break;
                    }
                    for (let i = 0; i < natoms; i++) {
                        arr.push(readFloat(offset + 4 + i * 4));
                    }
                    const closeSize = readInt(offset + 4 + expectedCoordBlock);
                    if (closeSize !== expectedCoordBlock) {
                        dimOk = false;
                        break;
                    }
                    offset = offset + 4 + expectedCoordBlock + 4;
                }

                if (!dimOk) {
                    if (frames.length === 0) {
                        parseError = 'DCD坐标块大小不匹配，文件可能已损坏或字节序错误';
                    }
                    break;
                }

                const frameCoords = [];
                let coordOk = true;
                for (let i = 0; i < natoms; i++) {
                    const x = xCoords[i], y = yCoords[i], z = zCoords[i];
                    if (!isFinite(x) || !isFinite(y) || !isFinite(z) ||
                        Math.abs(x) > 1e5 || Math.abs(y) > 1e5 || Math.abs(z) > 1e5) {
                        if (frames.length === 0) {
                            parseError = 'DCD坐标包含异常值（非有限数或超范围），可能字节序错误';
                        }
                        coordOk = false;
                        break;
                    }
                    frameCoords.push([x, y, z]);
                }
                if (!coordOk) break;

                frames.push(frameCoords);
                frameTimes.push(frameIdx * 0.001);
                frameIdx++;

                if (frames.length === 1) {
                    for (let i = 0; i < natoms; i++) {
                        topologyAtoms.push({
                            serial: i + 1,
                            atom: 'C',
                            element: 'C',
                            resn: 'UNK',
                            chain: 'A',
                            resi: Math.floor(i / 8) + 1,
                            x: xCoords[i],
                            y: yCoords[i],
                            z: zCoords[i]
                        });
                    }
                }

                if (frames.length >= 10000) break;

            } catch (e) {
                if (frames.length === 0) {
                    parseError = e.message;
                }
                break;
            }
        }

        if (frames.length === 0) {
            throw new Error(parseError || '未能从DCD文件中解析出有效帧');
        }

        const summary = {
            format: 'DCD' + (littleEndian ? ' (小端)' : ' (大端)'),
            numAtoms: natoms,
            numFrames: frames.length,
            numFramesHeader: numFramesHeader,
            timeStep: 0.001,
            hasUnitCell: hasUnitCell,
            charmmVersion: charmmVersion,
            byteOrder: littleEndian ? 'little-endian' : 'big-endian'
        };

        if (this.topology && this.topology.numAtoms > 0) {
            if (frames[0].length !== this.topology.numAtoms) {
                throw new Error(`轨迹原子数(${frames[0].length})与拓扑(${this.topology.numAtoms})不匹配`);
            }
            this.frames = frames;
            this.frameTimes = frameTimes;
            this.summary = summary;
            this.currentFrame = 0;
            return { numAtoms: this.topology.numAtoms, numFrames: frames.length, summary };
        }

        this.setTopology(topologyAtoms);
        this.frames = frames;
        this.frameTimes = frameTimes;
        this.summary = summary;
        this.currentFrame = 0;

        return { numAtoms: topologyAtoms.length, numFrames: frames.length, summary };
    }

    async parseXTC(arrayBuffer) {
        const view = new DataView(arrayBuffer);
        const byteLength = view.byteLength;

        if (byteLength < 52) {
            throw new Error('XTC文件过小，不是有效的XTC轨迹');
        }

        const readInt = (off) => view.getInt32(off, false);
        const readFloat = (off) => view.getFloat32(off, false);

        const magic = readInt(0);
        if (magic !== 1995) {
            throw new Error('XTC文件格式错误：文件头magic number不是1995，该文件可能不是有效的XTC轨迹或已损坏');
        }

        const numAtoms = readInt(4);
        if (numAtoms <= 0 || numAtoms > 1000000) {
            throw new Error(`XTC文件原子数无效: ${numAtoms}`);
        }

        const frames = [];
        let topologyAtoms = [];
        let firstFrame = true;
        let offset = 0;

        const headerSize = 4 + 4 + 4 + 4 + 36;

        while (offset + headerSize <= byteLength) {
            const frameMagic = readInt(offset);
            if (frameMagic !== 1995) break;

            const frameNatoms = readInt(offset + 4);
            if (frameNatoms !== numAtoms) break;

            const boxStart = offset + 16;
            let boxValid = true;
            for (let b = 0; b < 9; b++) {
                const bv = readFloat(boxStart + b * 4);
                if (!isFinite(bv)) { boxValid = false; break; }
            }
            if (!boxValid) break;

            const coordStart = offset + headerSize;
            const coordBytes = numAtoms * 3 * 4;

            if (coordStart + coordBytes > byteLength) {
                if (frames.length === 0) {
                    throw new Error('检测到XTC压缩格式。压缩XTC轨迹（GROMACS默认输出）暂不支持在浏览器端解析，请使用未压缩格式，或先用 gmx trjconv 转换为DCD/PDB多帧格式后再导入。');
                }
                break;
            }

            const frameCoords = [];
            let allFinite = true;
            for (let i = 0; i < numAtoms; i++) {
                const base = coordStart + i * 12;
                const x = readFloat(base);
                const y = readFloat(base + 4);
                const z = readFloat(base + 8);
                if (!isFinite(x) || !isFinite(y) || !isFinite(z) ||
                    Math.abs(x) > 1e5 || Math.abs(y) > 1e5 || Math.abs(z) > 1e5) {
                    allFinite = false;
                    break;
                }
                frameCoords.push([x, y, z]);
            }

            if (!allFinite) {
                if (frames.length === 0) {
                    throw new Error('检测到XTC压缩格式。压缩XTC轨迹（GROMACS默认输出）暂不支持在浏览器端解析，请使用未压缩格式，或先用 gmx trjconv 转换为DCD/PDB多帧格式后再导入。');
                }
                break;
            }

            frames.push(frameCoords);

            if (firstFrame) {
                for (let i = 0; i < numAtoms; i++) {
                    topologyAtoms.push({
                        serial: i + 1,
                        atom: 'C',
                        element: 'C',
                        resn: 'UNK',
                        chain: 'A',
                        resi: Math.floor(i / 8) + 1,
                        x: frameCoords[i][0],
                        y: frameCoords[i][1],
                        z: frameCoords[i][2]
                    });
                }
                firstFrame = false;
            }

            offset = coordStart + coordBytes;
            if (frames.length >= 5000) break;
        }

        if (frames.length === 0) {
            throw new Error('未能从XTC文件中解析出有效帧，该文件可能是压缩格式或已损坏');
        }

        const frameTimes = frames.map((_, i) => i * 0.001);
        const summary = {
            format: 'XTC (未压缩)',
            numAtoms: numAtoms,
            numFrames: frames.length,
            timeStep: 0.001,
            hasUnitCell: true
        };

        if (this.topology && this.topology.numAtoms > 0) {
            if (frames[0].length !== this.topology.numAtoms) {
                throw new Error(`轨迹原子数(${frames[0].length})与拓扑(${this.topology.numAtoms})不匹配`);
            }
            this.frames = frames;
            this.frameTimes = frameTimes;
            this.summary = summary;
            this.currentFrame = 0;
            return { numAtoms: this.topology.numAtoms, numFrames: frames.length, summary };
        }

        this.setTopology(topologyAtoms);
        this.frames = frames;
        this.frameTimes = frameTimes;
        this.summary = summary;
        this.currentFrame = 0;

        return { numAtoms: topologyAtoms.length, numFrames: frames.length, summary };
    }

    async parseTRR(arrayBuffer) {
        const view = new DataView(arrayBuffer);
        let offset = 0;

        const readInt = () => {
            const val = view.getInt32(offset, false);
            offset += 4;
            return val;
        };

        const readFloat = () => {
            const val = view.getFloat32(offset, false);
            offset += 4;
            return val;
        };

        const frames = [];
        let topologyAtoms = [];
        let firstFrame = true;
        let frameCount = 0;

        while (offset < view.byteLength - 64) {
            try {
                const startOffset = offset;
                const magic = readInt();
                
                if (magic !== 1995 && magic !== 1993) {
                    offset = startOffset + 1;
                    continue;
                }

                const version = readInt();
                const numAtoms = readInt();
                
                if (numAtoms <= 0 || numAtoms > 1000000) {
                    offset = startOffset + 1;
                    continue;
                }

                const step = readInt();
                const time = readFloat();
                const lambda = readFloat();

                const hasBox = readInt();
                if (hasBox) offset += 9 * 4;

                const hasX = readInt();
                let frameCoords = [];

                if (hasX && numAtoms > 0) {
                    if (offset + numAtoms * 12 > view.byteLength) {
                        throw new Error('TRR文件坐标数据不完整');
                    }
                    for (let i = 0; i < numAtoms; i++) {
                        const x = readFloat() * 10;
                        const y = readFloat() * 10;
                        const z = readFloat() * 10;
                        frameCoords.push([x, y, z]);
                    }
                    frames.push(frameCoords);
                    frameCount++;

                    if (firstFrame) {
                        for (let i = 0; i < numAtoms; i++) {
                            topologyAtoms.push({
                                serial: i + 1,
                                atom: 'C',
                                element: 'C',
                                resn: 'UNK',
                                chain: 'A',
                                resi: Math.floor(i / 8) + 1,
                                x: frameCoords[i][0],
                                y: frameCoords[i][1],
                                z: frameCoords[i][2]
                            });
                        }
                        firstFrame = false;
                    }
                } else {
                    break;
                }

                const hasV = readInt();
                if (hasV) {
                    if (offset + numAtoms * 12 > view.byteLength) break;
                    offset += numAtoms * 12;
                }

                const hasF = readInt();
                if (hasF) {
                    if (offset + numAtoms * 12 > view.byteLength) break;
                    offset += numAtoms * 12;
                }

                if (frameCount >= 5000) break;

            } catch (e) {
                break;
            }
        }

        if (frames.length === 0) {
            throw new Error('未能从TRR文件中解析出有效帧');
        }

        const frameTimes = frames.map((_, i) => i * 0.001);
        const summary = {
            format: 'TRR',
            numAtoms: frames[0].length,
            numFrames: frames.length,
            timeStep: 0.001,
            hasUnitCell: true
        };

        if (this.topology && this.topology.numAtoms > 0) {
            if (frames[0] && frames[0].length !== this.topology.numAtoms) {
                throw new Error(`轨迹原子数(${frames[0].length})与拓扑(${this.topology.numAtoms})不匹配`);
            }
            this.frames = frames;
            this.frameTimes = frameTimes;
            this.summary = summary;
            this.currentFrame = 0;
            return { numAtoms: this.topology.numAtoms, numFrames: frames.length, summary };
        }

        this.setTopology(topologyAtoms);
        this.frames = frames;
        this.frameTimes = frameTimes;
        this.summary = summary;
        this.currentFrame = 0;

        return { numAtoms: topologyAtoms.length, numFrames: frames.length, summary };
    }

    generatePDBString(frameIndex = 0) {
        if (!this.topology) return '';
        
        const atoms = this.topology.atoms;
        const coords = this.getFrame(frameIndex) || [];
        let pdb = '';
        
        atoms.forEach((atom, i) => {
            const coord = coords[i] || [atom.x || 0, atom.y || 0, atom.z || 0];
            const record = (atom.record || 'ATOM').padEnd(6);
            const serial = ((atom.serial || i + 1) + '').padStart(5);
            const atomName = (atom.atom || 'C').padEnd(4);
            const resn = (atom.resn || 'UNK').padEnd(3);
            const chain = (atom.chain || 'A').padStart(1);
            const resi = ((atom.resi || 1) + '').padStart(4);
            const x = coord[0].toFixed(3).padStart(8);
            const y = coord[1].toFixed(3).padStart(8);
            const z = coord[2].toFixed(3).padStart(8);
            const occ = ((atom.occupancy != null) ? atom.occupancy : 1.0).toFixed(2).padStart(6);
            const temp = ((atom.tempFactor != null) ? atom.tempFactor : 0.0).toFixed(2).padStart(6);
            const element = (atom.element || 'C').padStart(2);
            
            pdb += `${record}${serial} ${atomName}${resn} ${chain}${resi}    ${x}${y}${z}${occ}${temp}          ${element}\n`;
        });
        
        pdb += 'END\n';
        return pdb;
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
