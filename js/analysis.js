class MDAnalysis {
    constructor(trajectory) {
        this.trajectory = trajectory;
        this.results = {};
    }

    setTrajectory(trajectory) {
        this.trajectory = trajectory;
        this.results = {};
    }

    calculateRMSD(selection = 'backbone', referenceFrame = 0) {
        if (!this.trajectory || this.trajectory.getNumFrames() === 0) {
            return [];
        }

        const atomIndices = this.trajectory.getSelectionIndices(selection);
        if (atomIndices.length === 0) {
            return [];
        }

        const refCoords = this.trajectory.getAtomPositions(referenceFrame, atomIndices);
        const refCentroid = this.calculateCentroid(refCoords);
        const refAligned = this.translateCoords(refCoords, this.negateVector(refCentroid));

        const rmsdValues = [];
        const numFrames = this.trajectory.getNumFrames();

        for (let f = 0; f < numFrames; f++) {
            const frameCoords = this.trajectory.getAtomPositions(f, atomIndices);
            const centroid = this.calculateCentroid(frameCoords);
            const translated = this.translateCoords(frameCoords, this.negateVector(centroid));
            
            const aligned = this.kabschAlign(translated, refAligned);
            const rmsd = this.calculateRMSDValue(aligned, refAligned);
            rmsdValues.push(rmsd);
        }

        this.results.rmsd = {
            values: rmsdValues,
            selection,
            referenceFrame,
            atomCount: atomIndices.length
        };

        return rmsdValues;
    }

    calculateRMSF(selection = 'backbone', byResidue = true) {
        if (!this.trajectory || this.trajectory.getNumFrames() === 0) {
            return [];
        }

        const atomIndices = this.trajectory.getSelectionIndices(selection);
        if (atomIndices.length === 0) {
            return [];
        }

        const numFrames = this.trajectory.getNumFrames();
        const numAtoms = atomIndices.length;

        const avgPositions = new Array(numAtoms).fill(null).map(() => [0, 0, 0]);

        for (let f = 0; f < numFrames; f++) {
            const coords = this.trajectory.getAtomPositions(f, atomIndices);
            for (let i = 0; i < numAtoms; i++) {
                avgPositions[i][0] += coords[i][0];
                avgPositions[i][1] += coords[i][1];
                avgPositions[i][2] += coords[i][2];
            }
        }

        for (let i = 0; i < numAtoms; i++) {
            avgPositions[i][0] /= numFrames;
            avgPositions[i][1] /= numFrames;
            avgPositions[i][2] /= numFrames;
        }

        const rmsfValues = new Array(numAtoms).fill(0);
        for (let f = 0; f < numFrames; f++) {
            const coords = this.trajectory.getAtomPositions(f, atomIndices);
            for (let i = 0; i < numAtoms; i++) {
                const dx = coords[i][0] - avgPositions[i][0];
                const dy = coords[i][1] - avgPositions[i][1];
                const dz = coords[i][2] - avgPositions[i][2];
                rmsfValues[i] += dx * dx + dy * dy + dz * dz;
            }
        }

        for (let i = 0; i < numAtoms; i++) {
            rmsfValues[i] = Math.sqrt(rmsfValues[i] / numFrames);
        }

        if (byResidue && this.trajectory.topology) {
            const residues = this.trajectory.topology.residues;
            const resRMSF = [];
            
            for (let r = 0; r < residues.length; r++) {
                const res = residues[r];
                let sum = 0;
                let count = 0;
                
                for (let j = 0; j < atomIndices.length; j++) {
                    if (atomIndices[j] >= res.startIndex && atomIndices[j] <= res.endIndex) {
                        sum += rmsfValues[j];
                        count++;
                    }
                }
                
                if (count > 0) {
                    resRMSF.push({
                        resi: res.resi,
                        resn: res.resn,
                        value: sum / count
                    });
                }
            }

            this.results.rmsf = {
                values: resRMSF,
                byResidue: true,
                selection,
                residueCount: resRMSF.length
            };

            return resRMSF;
        }

        this.results.rmsf = {
            values: rmsfValues,
            byResidue: false,
            selection,
            atomCount: numAtoms
        };

        return rmsfValues;
    }

    calculateRDF(referenceSelection = 'all', targetSelection = 'all', maxDistance = 10, binWidth = 0.1) {
        if (!this.trajectory || this.trajectory.getNumFrames() === 0) {
            return { distances: [], rdf: [] };
        }

        const refIndices = this.trajectory.getSelectionIndices(referenceSelection);
        const targetIndices = this.trajectory.getSelectionIndices(targetSelection);

        if (refIndices.length === 0 || targetIndices.length === 0) {
            return { distances: [], rdf: [] };
        }

        const numBins = Math.ceil(maxDistance / binWidth);
        const histogram = new Array(numBins).fill(0);
        const distances = [];
        
        for (let i = 0; i < numBins; i++) {
            distances.push((i + 0.5) * binWidth);
        }

        const numFrames = Math.min(this.trajectory.getNumFrames(), 20);

        for (let f = 0; f < numFrames; f++) {
            const refCoords = this.trajectory.getAtomPositions(f, refIndices);
            const targetCoords = this.trajectory.getAtomPositions(f, targetIndices);

            for (let i = 0; i < refCoords.length; i++) {
                for (let j = 0; j < targetCoords.length; j++) {
                    if (refIndices[i] === targetIndices[j]) continue;
                    
                    const dist = this.distance3D(refCoords[i], targetCoords[j]);
                    const binIndex = Math.floor(dist / binWidth);
                    
                    if (binIndex >= 0 && binIndex < numBins) {
                        histogram[binIndex]++;
                    }
                }
            }
        }

        const rdf = histogram.map((count, i) => {
            const r1 = i * binWidth;
            const r2 = (i + 1) * binWidth;
            const volume = (4/3) * Math.PI * (r2 * r2 * r2 - r1 * r1 * r1);
            const density = targetIndices.length / ((4/3) * Math.PI * Math.pow(maxDistance, 3));
            const expected = density * volume * numFrames * refIndices.length;
            
            return expected > 0 ? count / expected : 0;
        });

        this.results.rdf = {
            distances,
            rdf,
            maxDistance,
            binWidth,
            referenceSelection,
            targetSelection
        };

        return { distances, rdf };
    }

    analyzeHydrogenBonds(donorResidues = null, acceptorResidues = null, maxDistance = 3.5, minAngle = 120) {
        if (!this.trajectory || !this.trajectory.topology) {
            return { pairs: [], perFrameCounts: [] };
        }

        const atoms = this.trajectory.topology.atoms;
        const numFrames = this.trajectory.getNumFrames();

        const donorAtomNames = ['N', 'O', 'S', 'ND1', 'ND2', 'NE1', 'NE2', 'NZ', 'OG', 'OG1', 'SG', 'OH'];
        const acceptorAtomNames = ['O', 'N', 'S', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'SG', 'OH'];

        const donorIndices = [];
        const acceptorIndices = [];

        atoms.forEach((atom, index) => {
            const atomName = (atom.atom || '').toUpperCase();
            const element = (atom.element || atomName.charAt(0)).toUpperCase();
            const resn = (atom.resn || '').toUpperCase();

            const isDonorElement = ['N', 'O', 'S'].includes(element.charAt(0));
            const isAcceptorElement = ['O', 'N', 'S'].includes(element.charAt(0));

            const passesDonorFilter = !donorResidues || donorResidues.length === 0 || donorResidues.includes(resn);
            const passesAcceptorFilter = !acceptorResidues || acceptorResidues.length === 0 || acceptorResidues.includes(resn);

            if (isDonorElement && passesDonorFilter) {
                donorIndices.push(index);
            }

            if (isAcceptorElement && passesAcceptorFilter) {
                acceptorIndices.push(index);
            }
        });

        const perFrameCounts = new Array(numFrames).fill(0);
        const hbondPairs = {};

        for (let f = 0; f < numFrames; f++) {
            const coords = this.trajectory.getFrame(f);
            const frameHBonds = new Set();

            for (let di = 0; di < donorIndices.length; di++) {
                const donorIdx = donorIndices[di];
                const donorCoord = coords[donorIdx];

                for (let ai = 0; ai < acceptorIndices.length; ai++) {
                    const acceptorIdx = acceptorIndices[ai];
                    if (donorIdx === acceptorIdx) continue;

                    const donorAtom = atoms[donorIdx];
                    const acceptorAtom = atoms[acceptorIdx];

                    if (donorAtom.resi === acceptorAtom.resi && donorAtom.chain === acceptorAtom.chain) {
                        continue;
                    }

                    const acceptorCoord = coords[acceptorIdx];
                    const dist = this.distance3D(donorCoord, acceptorCoord);

                    if (dist <= maxDistance) {
                        const pairKey = `${Math.min(donorIdx, acceptorIdx)}-${Math.max(donorIdx, acceptorIdx)}`;

                        if (!frameHBonds.has(pairKey)) {
                            frameHBonds.add(pairKey);
                            perFrameCounts[f]++;

                            if (!hbondPairs[pairKey]) {
                                hbondPairs[pairKey] = {
                                    donor: donorAtom,
                                    acceptor: acceptorAtom,
                                    donorIndex: donorIdx,
                                    acceptorIndex: acceptorIdx,
                                    count: 0
                                };
                            }
                            hbondPairs[pairKey].count++;
                        }
                    }
                }
            }
        }

        const results = Object.values(hbondPairs)
            .map(pair => ({
                ...pair,
                frequency: ((pair.count / numFrames) * 100).toFixed(1),
                pairLabel: `${pair.donor.resn}${pair.donor.resi}:${pair.donor.atom} - ${pair.acceptor.resn}${pair.acceptor.resi}:${pair.acceptor.atom}`
            }))
            .sort((a, b) => parseFloat(b.frequency) - parseFloat(a.frequency))
            .slice(0, 20);

        this.results.hbond = {
            pairs: results,
            perFrameCounts,
            maxDistance,
            minAngle,
            numFrames,
            donorResidues,
            acceptorResidues
        };

        return { pairs: results, perFrameCounts };
    }

    performPCA(selection = 'backbone') {
        if (!this.trajectory || this.trajectory.getNumFrames() < 3) {
            return { eigenvalues: [], eigenvectors: [], projections: [], varianceExplained: [] };
        }

        const atomIndices = this.trajectory.getSelectionIndices(selection);
        if (atomIndices.length < 3) {
            return { eigenvalues: [], eigenvectors: [], projections: [], varianceExplained: [] };
        }

        const numFrames = this.trajectory.getNumFrames();
        const numAtoms = atomIndices.length;
        const numCoords = numAtoms * 3;

        const trajectories = [];
        for (let f = 0; f < numFrames; f++) {
            const coords = this.trajectory.getAtomPositions(f, atomIndices);
            const flat = [];
            for (let i = 0; i < numAtoms; i++) {
                flat.push(coords[i][0], coords[i][1], coords[i][2]);
            }
            trajectories.push(flat);
        }

        const meanCoords = new Array(numCoords).fill(0);
        for (let f = 0; f < numFrames; f++) {
            for (let i = 0; i < numCoords; i++) {
                meanCoords[i] += trajectories[f][i];
            }
        }
        for (let i = 0; i < numCoords; i++) {
            meanCoords[i] /= numFrames;
        }

        const centered = trajectories.map(frame => 
            frame.map((val, i) => val - meanCoords[i])
        );

        let eigenvalues, eigenvectors;
        const maxFullPCA = 300;

        if (numCoords <= maxFullPCA) {
            const covarianceMatrix = [];
            for (let i = 0; i < numCoords; i++) {
                covarianceMatrix[i] = [];
                for (let j = 0; j < numCoords; j++) {
                    let sum = 0;
                    for (let f = 0; f < numFrames; f++) {
                        sum += centered[f][i] * centered[f][j];
                    }
                    covarianceMatrix[i][j] = sum / (numFrames - 1);
                }
            }
            const pcaResult = this.powerIterationPCA(covarianceMatrix, Math.min(10, numCoords));
            eigenvalues = pcaResult.eigenvalues;
            eigenvectors = pcaResult.eigenvectors;
        } else {
            const innerProduct = [];
            for (let i = 0; i < numFrames; i++) {
                innerProduct[i] = [];
                for (let j = 0; j < numFrames; j++) {
                    let sum = 0;
                    for (let k = 0; k < numCoords; k++) {
                        sum += centered[i][k] * centered[j][k];
                    }
                    innerProduct[i][j] = sum / (numFrames - 1);
                }
            }
            const innerPCA = this.powerIterationPCA(innerProduct, Math.min(10, numFrames));
            
            eigenvalues = innerPCA.eigenvalues.map(v => v);
            eigenvectors = [];
            
            const numPCs = innerPCA.eigenvectors.length;
            for (let p = 0; p < numPCs; p++) {
                const ev = innerPCA.eigenvectors[p];
                const eigVal = Math.max(eigenvalues[p], 1e-10);
                const factor = 1.0 / Math.sqrt((numFrames - 1) * eigVal);
                const fullEV = new Array(numCoords).fill(0);
                for (let f = 0; f < numFrames; f++) {
                    const coef = ev[f] * factor;
                    for (let k = 0; k < numCoords; k++) {
                        fullEV[k] += coef * centered[f][k];
                    }
                }
                let norm = 0;
                for (let k = 0; k < numCoords; k++) norm += fullEV[k] * fullEV[k];
                norm = Math.sqrt(norm);
                if (norm > 0) {
                    for (let k = 0; k < numCoords; k++) fullEV[k] /= norm;
                }
                eigenvectors.push(fullEV);
            }
        }

        const projections = [];
        for (let f = 0; f < numFrames; f++) {
            const proj = [];
            for (let p = 0; p < eigenvectors.length; p++) {
                let sum = 0;
                const ev = eigenvectors[p];
                const minLen = Math.min(ev.length, centered[f].length);
                for (let i = 0; i < minLen; i++) {
                    sum += centered[f][i] * ev[i];
                }
                proj.push(sum);
            }
            projections.push(proj);
        }

        const totalVariance = eigenvalues.reduce((a, b) => a + b, 0);
        const varianceExplained = eigenvalues.map(ev => totalVariance > 1e-12 ? (ev / totalVariance * 100).toFixed(1) : '0.0');

        let coordScale = 0;
        for (let i = 0; i < meanCoords.length; i++) {
            coordScale += meanCoords[i] * meanCoords[i];
        }
        coordScale = Math.sqrt(coordScale / Math.max(meanCoords.length, 1));
        const lowVariance = totalVariance < 1e-6 || (coordScale > 0 && Math.sqrt(totalVariance) < coordScale * 1e-4);

        this.results.pca = {
            eigenvalues,
            eigenvectors,
            projections,
            varianceExplained,
            meanCoords,
            selection,
            atomIndices,
            numAtoms,
            totalVariance,
            lowVariance
        };

        return { eigenvalues, eigenvectors, projections, varianceExplained };
    }

    generatePCModeFrames(pcIndex, numFrames = 20, amplitude = 2.0) {
        if (!this.results.pca) {
            this.performPCA('backbone');
        }

        const pca = this.results.pca;
        if (!pca || pcIndex >= pca.eigenvectors.length) {
            return null;
        }

        const { meanCoords, eigenvectors, atomIndices, numAtoms, lowVariance } = pca;
        const eigenvalue = pca.eigenvalues[pcIndex];
        const eigenvector = eigenvectors[pcIndex];

        if (!meanCoords || !eigenvector || meanCoords.length < numAtoms * 3 || eigenvector.length < numAtoms * 3) {
            return null;
        }

        let maxEVComp = 0;
        for (let k = 0; k < numAtoms * 3; k++) {
            const v = Math.abs(eigenvector[k]);
            if (isFinite(v) && v > maxEVComp) maxEVComp = v;
        }

        if (maxEVComp < 1e-10 || !isFinite(maxEVComp)) {
            return {
                frames: this.generateStaticFrames(meanCoords, numAtoms, numFrames),
                atomIndices,
                numAtoms,
                eigenvalue: 0,
                varianceExplained: '0.0',
                lowVariance: true,
                staticMode: true
            };
        }

        const effectiveAmplitude = lowVariance ? Math.min(amplitude, 0.4) : amplitude;
        const scale = effectiveAmplitude / maxEVComp;
        const frames = [];

        for (let f = 0; f < numFrames; f++) {
            const t = Math.sin((f / (numFrames - 1)) * Math.PI * 2) * 0.5;
            const factor = t * scale;

            const frameCoords = [];
            for (let i = 0; i < numAtoms; i++) {
                const baseIdx = i * 3;
                let x = meanCoords[baseIdx] + factor * eigenvector[baseIdx];
                let y = meanCoords[baseIdx + 1] + factor * eigenvector[baseIdx + 1];
                let z = meanCoords[baseIdx + 2] + factor * eigenvector[baseIdx + 2];

                if (!isFinite(x)) x = meanCoords[baseIdx];
                if (!isFinite(y)) y = meanCoords[baseIdx + 1];
                if (!isFinite(z)) z = meanCoords[baseIdx + 2];

                frameCoords.push([x, y, z]);
            }
            frames.push(frameCoords);
        }

        return {
            frames,
            atomIndices,
            numAtoms,
            eigenvalue,
            varianceExplained: pca.varianceExplained[pcIndex],
            lowVariance: !!lowVariance
        };
    }

    generateStaticFrames(meanCoords, numAtoms, numFrames) {
        const frames = [];
        for (let f = 0; f < numFrames; f++) {
            const frameCoords = [];
            for (let i = 0; i < numAtoms; i++) {
                const baseIdx = i * 3;
                frameCoords.push([
                    isFinite(meanCoords[baseIdx]) ? meanCoords[baseIdx] : 0,
                    isFinite(meanCoords[baseIdx + 1]) ? meanCoords[baseIdx + 1] : 0,
                    isFinite(meanCoords[baseIdx + 2]) ? meanCoords[baseIdx + 2] : 0
                ]);
            }
            frames.push(frameCoords);
        }
        return frames;
    }

    getPCDisplacementVectors(pcIndex, amplitude = 2.0) {
        if (!this.results.pca) {
            this.performPCA('backbone');
        }

        const pca = this.results.pca;
        if (!pca || pcIndex >= pca.eigenvectors.length) {
            return null;
        }

        const { meanCoords, eigenvectors, atomIndices, numAtoms } = pca;
        const eigenvalue = pca.eigenvalues[pcIndex];
        const eigenvector = eigenvectors[pcIndex];

        const scale = amplitude * Math.sqrt(Math.max(eigenvalue, 0.001));

        const displacements = [];
        for (let i = 0; i < numAtoms; i++) {
            const baseIdx = i * 3;
            const start = [
                meanCoords[baseIdx],
                meanCoords[baseIdx + 1],
                meanCoords[baseIdx + 2]
            ];
            const end = [
                meanCoords[baseIdx] + scale * eigenvector[baseIdx],
                meanCoords[baseIdx + 1] + scale * eigenvector[baseIdx + 1],
                meanCoords[baseIdx + 2] + scale * eigenvector[baseIdx + 2]
            ];
            displacements.push({ start, end, atomIndex: atomIndices[i] });
        }

        return {
            displacements,
            atomIndices,
            eigenvalue,
            varianceExplained: pca.varianceExplained[pcIndex]
        };
    }

    powerIterationPCA(matrix, numComponents) {
        const n = matrix.length;
        const eigenvalues = [];
        const eigenvectors = [];
        const residual = matrix.map(row => [...row]);

        for (let k = 0; k < numComponents; k++) {
            let vector = new Array(n).fill(0).map(() => Math.random() - 0.5);
            let norm = Math.sqrt(vector.reduce((s, v) => s + v * v, 0));
            if (norm < 1e-12) {
                vector = new Array(n).fill(0);
                vector[0] = 1;
                norm = 1;
            }
            vector = vector.map(v => v / norm);

            let eigenvalue = 0;

            for (let iter = 0; iter < 100; iter++) {
                const newVector = [];
                for (let i = 0; i < n; i++) {
                    let sum = 0;
                    for (let j = 0; j < n; j++) {
                        sum += residual[i][j] * vector[j];
                    }
                    newVector.push(sum);
                }

                const newNorm = Math.sqrt(newVector.reduce((s, v) => s + v * v, 0));

                if (newNorm < 1e-12) {
                    eigenvalue = 0;
                    vector = new Array(n).fill(0);
                    vector[k % n] = 1;
                    break;
                }

                const newEigenvalue = newNorm;

                if (Math.abs(newEigenvalue - eigenvalue) < 1e-10) {
                    eigenvalue = newEigenvalue;
                    vector = newVector.map(v => v / newNorm);
                    break;
                }

                eigenvalue = newEigenvalue;
                vector = newVector.map(v => v / newNorm);
            }

            eigenvalues.push(eigenvalue);
            eigenvectors.push(vector);

            for (let i = 0; i < n; i++) {
                for (let j = 0; j < n; j++) {
                    residual[i][j] -= eigenvalue * vector[i] * vector[j];
                }
            }
        }

        return { eigenvalues, eigenvectors };
    }

    calculateRadiusOfGyration(selection = 'backbone') {
        if (!this.trajectory || this.trajectory.getNumFrames() === 0) {
            return [];
        }

        const atomIndices = this.trajectory.getSelectionIndices(selection);
        if (atomIndices.length < 2) {
            return [];
        }

        const numFrames = this.trajectory.getNumFrames();
        const rgValues = [];

        for (let f = 0; f < numFrames; f++) {
            const coords = this.trajectory.getAtomPositions(f, atomIndices);
            
            const centroid = this.calculateCentroid(coords);
            
            let sumDistSq = 0;
            for (let i = 0; i < coords.length; i++) {
                const dx = coords[i][0] - centroid[0];
                const dy = coords[i][1] - centroid[1];
                const dz = coords[i][2] - centroid[2];
                sumDistSq += dx * dx + dy * dy + dz * dz;
            }
            
            const rg = Math.sqrt(sumDistSq / coords.length);
            rgValues.push(rg);
        }

        return rgValues;
    }

    getReactionCoordinate(type) {
        if (!this.trajectory || this.trajectory.getNumFrames() === 0) {
            return [];
        }

        switch (type) {
            case 'pc1':
            case 'pc2':
            case 'pc3':
                if (!this.results.pca) {
                    this.performPCA('backbone');
                }
                if (this.results.pca && this.results.pca.projections) {
                    const pcIdx = parseInt(type.charAt(2)) - 1;
                    return this.results.pca.projections.map(p => p[pcIdx] || 0);
                }
                return [];
                
            case 'rmsd':
                if (!this.results.rmsd) {
                    this.calculateRMSD('backbone', 0);
                }
                return this.results.rmsd?.values || [];
                
            case 'rg':
                return this.calculateRadiusOfGyration('backbone');
                
            case 'rmsf':
                if (!this.results.rmsf) {
                    this.calculateRMSF('backbone', false);
                }
                const rmsf = this.results.rmsf?.values;
                if (rmsf && Array.isArray(rmsf)) {
                    const result = [];
                    for (let i = 0; i < this.trajectory.getNumFrames(); i++) {
                        result.push(rmsf[i % rmsf.length]);
                    }
                    return result;
                }
                return [];
                
            default:
                return [];
        }
    }

    calculateFreeEnergySurface(xData, yData, temperature = 300) {
        if (!xData || !yData || xData.length === 0) {
            return { x: [], y: [], fes: [], minFes: 0 };
        }

        const R = 0.0019872;

        const xMin = Math.min(...xData);
        const xMax = Math.max(...xData);
        const yMin = Math.min(...yData);
        const yMax = Math.max(...yData);

        const numBins = 50;
        const xBinWidth = (xMax - xMin) / numBins;
        const yBinWidth = (yMax - yMin) / numBins;

        const histogram = [];
        const xCenters = [];
        const yCenters = [];

        for (let i = 0; i < numBins; i++) {
            histogram[i] = new Array(numBins).fill(0);
            xCenters.push(xMin + (i + 0.5) * xBinWidth);
            yCenters.push(yMin + (i + 0.5) * yBinWidth);
        }

        let maxCount = 0;
        for (let k = 0; k < xData.length; k++) {
            const xBin = Math.floor((xData[k] - xMin) / xBinWidth);
            const yBin = Math.floor((yData[k] - yMin) / yBinWidth);
            
            if (xBin >= 0 && xBin < numBins && yBin >= 0 && yBin < numBins) {
                histogram[xBin][yBin]++;
                if (histogram[xBin][yBin] > maxCount) {
                    maxCount = histogram[xBin][yBin];
                }
            }
        }

        const fes = [];
        let minFes = Infinity;

        for (let i = 0; i < numBins; i++) {
            fes[i] = [];
            for (let j = 0; j < numBins; j++) {
                if (histogram[i][j] > 0) {
                    const energy = -R * temperature * Math.log(histogram[i][j] / maxCount);
                    fes[i][j] = energy;
                    if (energy < minFes) {
                        minFes = energy;
                    }
                } else {
                    fes[i][j] = null;
                }
            }
        }

        this.results.fes = {
            x: xCenters,
            y: yCenters,
            fes,
            minFes,
            temperature,
            numBins
        };

        return { x: xCenters, y: yCenters, fes, minFes };
    }

    calculateRMSDValue(coords1, coords2) {
        if (coords1.length !== coords2.length) return 0;
        
        let sum = 0;
        for (let i = 0; i < coords1.length; i++) {
            const dx = coords1[i][0] - coords2[i][0];
            const dy = coords1[i][1] - coords2[i][1];
            const dz = coords1[i][2] - coords2[i][2];
            sum += dx * dx + dy * dy + dz * dz;
        }
        
        return Math.sqrt(sum / coords1.length);
    }

    calculateCentroid(coords) {
        const centroid = [0, 0, 0];
        for (const coord of coords) {
            centroid[0] += coord[0];
            centroid[1] += coord[1];
            centroid[2] += coord[2];
        }
        centroid[0] /= coords.length;
        centroid[1] /= coords.length;
        centroid[2] /= coords.length;
        return centroid;
    }

    translateCoords(coords, translation) {
        return coords.map(coord => [
            coord[0] + translation[0],
            coord[1] + translation[1],
            coord[2] + translation[2]
        ]);
    }

    negateVector(v) {
        return [-v[0], -v[1], -v[2]];
    }

    distance3D(a, b) {
        const dx = a[0] - b[0];
        const dy = a[1] - b[1];
        const dz = a[2] - b[2];
        return Math.sqrt(dx * dx + dy * dy + dz * dz);
    }

    kabschAlign(P, Q) {
        const { rotationMatrix } = this.kabschRotation(P, Q);
        return this.rotateCoords(P, rotationMatrix);
    }

    kabschRotation(P, Q) {
        const n = P.length;
        if (n < 3 || Q.length !== n) {
            return { rotationMatrix: this.identityMatrix(), rmsd: 0 };
        }

        const H = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < 3; j++) {
                for (let k = 0; k < 3; k++) {
                    H[j][k] += P[i][j] * Q[i][k];
                }
            }
        }

        const { U, S, V } = this.svd3x3(H);

        const d = this.det3x3(U) * this.det3x3(V);
        const sign = d >= 0 ? 1 : -1;

        const I = this.identityMatrix();
        I[2][2] = sign;

        const R = this.matMul(V, this.matMul(I, this.transpose3x3(U)));

        let rmsd = 0;
        for (let i = 0; i < n; i++) {
            const rotated = this.rotateVector(P[i], R);
            const dx = rotated[0] - Q[i][0];
            const dy = rotated[1] - Q[i][1];
            const dz = rotated[2] - Q[i][2];
            rmsd += dx * dx + dy * dy + dz * dz;
        }
        rmsd = Math.sqrt(rmsd / n);

        return { rotationMatrix: R, rmsd };
    }

    svd3x3(A) {
        const AtA = this.matMul(this.transpose3x3(A), A);
        const { eigenvalues, eigenvectors } = this.symmetricEigen(AtA);

        const S = eigenvalues.map(v => Math.sqrt(Math.max(v, 0)));

        const V = eigenvectors;

        const AV = this.matMul(A, V);
        const U = [];
        for (let j = 0; j < 3; j++) {
            U.push([0, 0, 0]);
        }
        for (let j = 0; j < 3; j++) {
            if (S[j] > 1e-10) {
                for (let i = 0; i < 3; i++) {
                    U[i][j] = AV[i][j] / S[j];
                }
            }
        }

        if (S[2] > 1e-10) {
            const u0 = [U[0][0], U[1][0], U[2][0]];
            const u1 = [U[0][1], U[1][1], U[2][1]];
            const u2 = this.cross(u0, u1);
            for (let i = 0; i < 3; i++) {
                U[i][2] = u2[i];
            }
        }

        return { U, S, V };
    }

    symmetricEigen(matrix) {
        const A = matrix.map(row => [...row]);
        const n = 3;
        const V = this.identityMatrix();

        for (let iter = 0; iter < 100; iter++) {
            let max = 0;
            let p = 0, q = 0;

            for (let i = 0; i < n; i++) {
                for (let j = i + 1; j < n; j++) {
                    if (Math.abs(A[i][j]) > max) {
                        max = Math.abs(A[i][j]);
                        p = i;
                        q = j;
                    }
                }
            }

            if (max < 1e-10) break;

            let theta = 0.5 * Math.atan2(2 * A[p][q], A[q][q] - A[p][p]);
            const c = Math.cos(theta);
            const s = Math.sin(theta);

            const app = A[p][p];
            const aqq = A[q][q];
            const apq = A[p][q];

            A[p][p] = c * c * app - 2 * s * c * apq + s * s * aqq;
            A[q][q] = s * s * app + 2 * s * c * apq + c * c * aqq;
            A[p][q] = 0;
            A[q][p] = 0;

            for (let i = 0; i < n; i++) {
                if (i !== p && i !== q) {
                    const aip = A[i][p];
                    const aiq = A[i][q];
                    A[i][p] = c * aip - s * aiq;
                    A[p][i] = A[i][p];
                    A[i][q] = s * aip + c * aiq;
                    A[q][i] = A[i][q];
                }
            }

            for (let i = 0; i < n; i++) {
                const vip = V[i][p];
                const viq = V[i][q];
                V[i][p] = c * vip - s * viq;
                V[i][q] = s * vip + c * viq;
            }
        }

        const eigenvalues = [A[0][0], A[1][1], A[2][2]];
        const eigenvectors = [];

        const indices = [0, 1, 2];
        indices.sort((a, b) => eigenvalues[b] - eigenvalues[a]);

        const sortedEigenvalues = indices.map(i => eigenvalues[i]);
        const sortedEigenvectors = [];
        for (const idx of indices) {
            sortedEigenvectors.push([V[0][idx], V[1][idx], V[2][idx]]);
        }

        const VMatrix = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
        for (let j = 0; j < 3; j++) {
            for (let i = 0; i < 3; i++) {
                VMatrix[i][j] = sortedEigenvectors[j][i];
            }
        }
        
        return { eigenvalues: sortedEigenvalues, eigenvectors: VMatrix };
    }

    identityMatrix() {
        return [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    }

    transpose3x3(M) {
        return [
            [M[0][0], M[1][0], M[2][0]],
            [M[0][1], M[1][1], M[2][1]],
            [M[0][2], M[1][2], M[2][2]]
        ];
    }

    matMul(A, B) {
        const result = [];
        for (let i = 0; i < 3; i++) {
            result.push([0, 0, 0]);
            for (let j = 0; j < 3; j++) {
                for (let k = 0; k < 3; k++) {
                    result[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return result;
    }

    det3x3(M) {
        return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
             - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
             + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
    }

    cross(a, b) {
        return [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]
        ];
    }

    rotateCoords(coords, rotationMatrix) {
        return coords.map(coord => this.rotateVector(coord, rotationMatrix));
    }

    rotateVector(v, R) {
        return [
            R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2],
            R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2],
            R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2]
        ];
    }

    getStatistics(values) {
        if (!values || values.length === 0) {
            return { avg: 0, max: 0, min: 0, std: 0, final: 0 };
        }

        const sum = values.reduce((a, b) => a + b, 0);
        const avg = sum / values.length;
        const max = Math.max(...values);
        const min = Math.min(...values);
        const final = values[values.length - 1];

        const variance = values.reduce((s, v) => s + (v - avg) * (v - avg), 0) / values.length;
        const std = Math.sqrt(variance);

        return { avg, max, min, std, final };
    }
}
