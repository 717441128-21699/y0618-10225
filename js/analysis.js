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
            return { pairs: [], frequencies: [] };
        }

        const atoms = this.trajectory.topology.atoms;
        const numFrames = this.trajectory.getNumFrames();

        const donorIndices = [];
        const acceptorIndices = [];

        const hbondDonorAtoms = ['N', 'O', 'S'];
        const hbondAcceptorAtoms = ['O', 'N', 'S'];

        atoms.forEach((atom, index) => {
            const element = (atom.element || atom.atom).toUpperCase();
            
            if (hbondDonorAtoms.includes(element.charAt(0))) {
                if (!donorResidues || donorResidues.includes(atom.resn.toUpperCase())) {
                    donorIndices.push(index);
                }
            }
            
            if (hbondAcceptorAtoms.includes(element.charAt(0))) {
                if (!acceptorResidues || acceptorResidues.includes(atom.resn.toUpperCase())) {
                    acceptorIndices.push(index);
                }
            }
        });

        const hbondPairs = {};

        const sampleFrames = Math.min(numFrames, 50);
        const frameStep = Math.max(1, Math.floor(numFrames / sampleFrames));

        for (let f = 0; f < numFrames; f += frameStep) {
            const coords = this.trajectory.getFrame(f);
            
            for (const donorIdx of donorIndices) {
                for (const acceptorIdx of acceptorIndices) {
                    if (donorIdx === acceptorIdx) continue;
                    
                    const dist = this.distance3D(coords[donorIdx], coords[acceptorIdx]);
                    
                    if (dist <= maxDistance) {
                        const pairKey = `${Math.min(donorIdx, acceptorIdx)}-${Math.max(donorIdx, acceptorIdx)}`;
                        
                        if (!hbondPairs[pairKey]) {
                            hbondPairs[pairKey] = {
                                donor: atoms[donorIdx],
                                acceptor: atoms[acceptorIdx],
                                donorIndex: donorIdx,
                                acceptorIndex: acceptorIdx,
                                count: 0,
                                totalFrames: 0
                            };
                        }
                        hbondPairs[pairKey].count++;
                        hbondPairs[pairKey].totalFrames++;
                    }
                }
            }
        }

        const actualFrames = Math.ceil(numFrames / frameStep);
        const results = Object.values(hbondPairs)
            .map(pair => ({
                ...pair,
                frequency: (pair.count / actualFrames * 100).toFixed(1),
                pairLabel: `${pair.donor.resn}${pair.donor.resi}:${pair.donor.atom} - ${pair.acceptor.resn}${pair.acceptor.resi}:${pair.acceptor.atom}`
            }))
            .sort((a, b) => b.frequency - a.frequency)
            .slice(0, 20);

        this.results.hbond = {
            pairs: results,
            maxDistance,
            minAngle,
            numFrames: actualFrames
        };

        return results;
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

        const covarianceMatrix = [];
        const sampleSize = Math.min(numCoords, 100);
        
        for (let i = 0; i < sampleSize; i++) {
            covarianceMatrix[i] = [];
            for (let j = 0; j < sampleSize; j++) {
                let sum = 0;
                for (let f = 0; f < numFrames; f++) {
                    sum += centered[f][i] * centered[f][j];
                }
                covarianceMatrix[i][j] = sum / (numFrames - 1);
            }
        }

        const { eigenvalues, eigenvectors } = this.powerIterationPCA(covarianceMatrix, Math.min(10, sampleSize));

        const projections = [];
        for (let f = 0; f < numFrames; f++) {
            const proj = [];
            for (let p = 0; p < eigenvalues.length; p++) {
                let sum = 0;
                for (let i = 0; i < sampleSize; i++) {
                    sum += centered[f][i] * eigenvectors[p][i];
                }
                proj.push(sum);
            }
            projections.push(proj);
        }

        const totalVariance = eigenvalues.reduce((a, b) => a + b, 0);
        const varianceExplained = eigenvalues.map(ev => (ev / totalVariance * 100).toFixed(1));

        this.results.pca = {
            eigenvalues,
            eigenvectors,
            projections,
            varianceExplained,
            meanCoords,
            selection,
            atomIndices,
            numAtoms
        };

        return { eigenvalues, eigenvectors, projections, varianceExplained };
    }

    powerIterationPCA(matrix, numComponents) {
        const n = matrix.length;
        const eigenvalues = [];
        const eigenvectors = [];
        const residual = matrix.map(row => [...row]);

        for (let k = 0; k < numComponents; k++) {
            let vector = new Array(n).fill(0).map(() => Math.random() - 0.5);
            let norm = Math.sqrt(vector.reduce((s, v) => s + v * v, 0));
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
        return P;
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
