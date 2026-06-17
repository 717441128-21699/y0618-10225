let viewer = null;
let trajectory = null;
let analysis = null;
let chartManager = null;
let animationController = null;
let trajectoryComparison = null;

document.addEventListener('DOMContentLoaded', function() {
    initApp();
});

function initApp() {
    trajectory = new Trajectory();
    analysis = new MDAnalysis(trajectory);
    chartManager = new ChartManager();
    trajectoryComparison = new TrajectoryComparison();
    
    viewer = new MolecularViewer('viewer');
    animationController = new AnimationController(viewer, trajectory);
    
    animationController.setOnFrameChange(onFrameChange);
    
    chartManager.createRMSDChart('rmsdChart');
    chartManager.createRMSFChart('rmsfChart');
    chartManager.createRDFChart('rdfChart');
    chartManager.createHBondChart('hbondChart');
    chartManager.createPCAChart('pcaChart');
    chartManager.createFESChart('fesChart');
    
    setupFileUpload();
    setupDragAndDrop();
    
    setTimeout(() => {
        if (viewer) {
            viewer.resize();
        }
    }, 100);
}

function setupFileUpload() {
    const topologyInput = document.getElementById('topologyFile');
    const trajectoryInput = document.getElementById('trajectoryFile');
    
    topologyInput.addEventListener('change', (e) => {
        const file = e.target.files[0];
        if (file) {
            loadTopologyFile(file);
        }
    });
    
    trajectoryInput.addEventListener('change', (e) => {
        const file = e.target.files[0];
        if (file) {
            loadTrajectoryFile(file);
        }
    });
}

function setupDragAndDrop() {
    const topologyArea = document.getElementById('topologyUploadArea');
    const trajectoryArea = document.getElementById('trajectoryUploadArea');
    const allAreas = [topologyArea, trajectoryArea];
    
    allAreas.forEach(area => {
        if (!area) return;
        
        ['dragenter', 'dragover', 'dragleave', 'drop'].forEach(eventName => {
            area.addEventListener(eventName, preventDefaults, false);
        });
        
        ['dragenter', 'dragover'].forEach(eventName => {
            area.addEventListener(eventName, () => area.classList.add('drag-over'), false);
        });
        
        ['dragleave', 'drop'].forEach(eventName => {
            area.addEventListener(eventName, () => area.classList.remove('drag-over'), false);
        });
    });
    
    function preventDefaults(e) {
        e.preventDefault();
        e.stopPropagation();
    }
    
    topologyArea.addEventListener('drop', (e) => {
        const file = e.dataTransfer.files[0];
        if (file) {
            const ext = file.name.split('.').pop().toLowerCase();
            if (['pdb', 'gro', 'xyz'].includes(ext)) {
                loadTopologyFile(file);
            } else {
                showToast('请将拓扑文件拖到拓扑区域，轨迹文件拖到轨迹区域', 'info');
            }
        }
    }, false);
    
    trajectoryArea.addEventListener('drop', (e) => {
        const file = e.dataTransfer.files[0];
        if (file) {
            loadTrajectoryFile(file);
        }
    }, false);
}

function loadTopologyFile(file) {
    const reader = new FileReader();
    const ext = file.name.split('.').pop().toLowerCase();
    
    showLoading('加载拓扑文件...');
    
    reader.onload = function(e) {
        const content = e.target.result;
        
        let result = null;
        let pdbContent = content;
        
        try {
            if (ext === 'pdb') {
                result = trajectory.parseMultiFramePDB(content);
                pdbContent = trajectory.generatePDBString(0);
            } else if (ext === 'xyz') {
                result = trajectory.parseMultiFrameXYZ(content);
                pdbContent = trajectory.generatePDBString(0);
            } else if (ext === 'gro') {
                result = parseGRO(content);
                pdbContent = trajectory.generatePDBString(0);
            }
            
            if (result) {
                viewer.loadPDB(pdbContent);
                analysis.setTrajectory(trajectory);
                animationController.setTrajectory(trajectory);
                
                updateTopologyUI(file.name);
                updateInfoBadges();
                
                if (result.numFrames > 1) {
                    updateTrajectoryUI(file.name, result.numFrames);
                    showToast(`${file.name} 加载成功 (${result.numAtoms}个原子, ${result.numFrames}帧)`, 'success');
                } else {
                    showToast(`${file.name} 加载成功 (${result.numAtoms}个原子)`, 'success');
                }
            }
        } catch (err) {
            showToast(`加载失败: ${err.message}`, 'error');
        }
        
        hideLoading();
    };
    
    reader.onerror = function() {
        hideLoading();
        showToast('文件读取失败', 'error');
    };
    
    reader.readAsText(file);
}

function loadTrajectoryFile(file) {
    const ext = file.name.split('.').pop().toLowerCase();
    const isBinary = ['xtc', 'dcd', 'trr'].includes(ext);
    
    if (!trajectory.topology) {
        showToast('请先加载拓扑文件', 'warning');
        return;
    }
    
    showLoading('加载轨迹文件...');
    
    if (isBinary) {
        const reader = new FileReader();
        
        reader.onload = async function(e) {
            try {
                const buffer = e.target.result;
                let result;
                
                if (ext === 'dcd') {
                    result = await trajectory.parseDCD(buffer);
                } else if (ext === 'xtc') {
                    result = await trajectory.parseXTC(buffer);
                } else if (ext === 'trr') {
                    result = await trajectory.parseTRR(buffer);
                }
                
                if (result && result.numFrames > 0) {
                    const pdbContent = trajectory.generatePDBString(0);
                    viewer.loadPDB(pdbContent);
                    analysis.setTrajectory(trajectory);
                    animationController.setTrajectory(trajectory);
                    
                    updateTrajectoryUI(file.name, result.numFrames);
                    updateInfoBadges();
                    
                    showToast(`轨迹加载成功 (${result.numFrames}帧)`, 'success');
                }
            } catch (err) {
                showToast(`轨迹加载失败: ${err.message}`, 'error');
            }
            
            hideLoading();
        };
        
        reader.onerror = function() {
            hideLoading();
            showToast('文件读取失败', 'error');
        };
        
        reader.readAsArrayBuffer(file);
    } else {
        const reader = new FileReader();
        
        reader.onload = function(e) {
            try {
                const content = e.target.result;
                let result;
                
                if (ext === 'pdb') {
                    result = trajectory.parseMultiFramePDB(content);
                } else if (ext === 'xyz') {
                    result = trajectory.parseMultiFrameXYZ(content);
                }
                
                if (result && result.numFrames > 0) {
                    const pdbContent = trajectory.generatePDBString(0);
                    viewer.loadPDB(pdbContent);
                    analysis.setTrajectory(trajectory);
                    animationController.setTrajectory(trajectory);
                    
                    updateTrajectoryUI(file.name, result.numFrames);
                    updateInfoBadges();
                    
                    showToast(`轨迹加载成功 (${result.numFrames}帧)`, 'success');
                }
            } catch (err) {
                showToast(`轨迹加载失败: ${err.message}`, 'error');
            }
            
            hideLoading();
        };
        
        reader.onerror = function() {
            hideLoading();
            showToast('文件读取失败', 'error');
        };
        
        reader.readAsText(file);
    }
}

function parseGRO(groText) {
    const lines = groText.trim().split('\n');
    if (lines.length < 3) return null;
    
    const numAtoms = parseInt(lines[1].trim());
    const atoms = [];
    const coords = [];
    
    for (let i = 2; i < 2 + numAtoms && i < lines.length; i++) {
        const line = lines[i];
        const resi = parseInt(line.substring(0, 5));
        const resn = line.substring(5, 10).trim();
        const atom = line.substring(10, 15).trim();
        const serial = parseInt(line.substring(15, 20));
        const x = parseFloat(line.substring(20, 28)) * 10;
        const y = parseFloat(line.substring(28, 36)) * 10;
        const z = parseFloat(line.substring(36, 44)) * 10;
        
        atoms.push({
            serial, atom, resn, resi,
            chain: 'A',
            element: atom.charAt(0).toUpperCase(),
            x, y, z
        });
        coords.push([x, y, z]);
    }
    
    trajectory.setTopology(atoms);
    trajectory.addFrame(coords);
    
    return { numAtoms: atoms.length, numFrames: 1 };
}

function updateTopologyUI(fileName) {
    const placeholder = document.querySelector('#topologyUploadArea .upload-placeholder');
    const selected = document.getElementById('topologySelected');
    const fileNameEl = document.getElementById('topologyFileName');
    
    if (placeholder) placeholder.style.display = 'none';
    if (selected) selected.style.display = 'flex';
    if (fileNameEl) fileNameEl.textContent = fileName;
}

function updateTrajectoryUI(fileName, numFrames) {
    const placeholder = document.querySelector('#trajectoryUploadArea .upload-placeholder');
    const selected = document.getElementById('trajectorySelected');
    const fileNameEl = document.getElementById('trajectoryFileName');
    const frameCountEl = document.getElementById('trajectoryFrameCount');
    
    if (placeholder) placeholder.style.display = 'none';
    if (selected) selected.style.display = 'flex';
    if (fileNameEl) fileNameEl.textContent = fileName;
    if (frameCountEl) frameCountEl.textContent = `${numFrames} 帧`;
}

function loadSampleData() {
    trajectory.clear();
    const numFrames = trajectory.generateSampleTrajectory(50);
    
    const pdbContent = generatePDBFromTrajectory();
    viewer.loadPDB(pdbContent);
    
    analysis.setTrajectory(trajectory);
    animationController.setTrajectory(trajectory);
    
    updateTopologyUI('sample_protein.pdb');
    updateTrajectoryUI('sample_trajectory.xtc', numFrames);
    updateInfoBadges();
    
    showToast(`示例数据加载成功 (${numFrames} 帧)`, 'success');
}

function generatePDBFromTrajectory() {
    if (!trajectory.topology) return '';
    
    const atoms = trajectory.topology.atoms;
    const coords = trajectory.getFrame(0);
    let pdb = '';
    
    atoms.forEach((atom, i) => {
        const coord = coords[i] || [0, 0, 0];
        const record = atom.record || 'ATOM';
        const serial = (atom.serial || i + 1).toString().padStart(5);
        const atomName = (atom.atom || 'C').padEnd(4);
        const resn = (atom.resn || 'UNK').padEnd(3);
        const chain = (atom.chain || 'A').padStart(1);
        const resi = (atom.resi || 1).toString().padStart(4);
        const x = coord[0].toFixed(3).padStart(8);
        const y = coord[1].toFixed(3).padStart(8);
        const z = coord[2].toFixed(3).padStart(8);
        const occ = (atom.occupancy || 1.0).toFixed(2).padStart(6);
        const temp = (atom.tempFactor || 0.0).toFixed(2).padStart(6);
        const element = (atom.element || 'C').padStart(2);
        
        pdb += `${record.padEnd(6)}${serial} ${atomName}${resn} ${chain}${resi}    ${x}${y}${z}${occ}${temp}          ${element}\n`;
    });
    
    pdb += 'END\n';
    return pdb;
}

function updateFileList(name, type) {
    const fileList = document.getElementById('fileList');
    const emptyItem = fileList.querySelector('.empty');
    if (emptyItem) {
        emptyItem.remove();
    }
    
    const fileItem = document.createElement('div');
    fileItem.className = 'file-item';
    fileItem.innerHTML = `
        <span class="file-name">${name}</span>
        <span class="file-type">${type === 'topology' ? '拓扑' : '轨迹'}</span>
        <span class="file-remove" onclick="removeFile(this)">×</span>
    `;
    fileList.appendChild(fileItem);
}

function removeFile(element) {
    element.parentElement.remove();
    const fileList = document.getElementById('fileList');
    if (fileList.children.length === 0) {
        fileList.innerHTML = '<div class="file-item empty"><span>暂无文件</span></div>';
    }
}

function updateInfoBadges() {
    document.getElementById('atomCount').textContent = `原子: ${trajectory.getNumAtoms()}`;
    document.getElementById('residueCount').textContent = `残基: ${trajectory.getNumResidues()}`;
    document.getElementById('frameCount').textContent = `帧数: ${trajectory.getNumFrames()}`;
    document.getElementById('simulationTime').textContent = `模拟时间: ${trajectory.getSimulationTime()} ns`;
    
    const frameSlider = document.getElementById('frameSlider');
    frameSlider.max = Math.max(0, trajectory.getNumFrames() - 1);
    frameSlider.value = 0;
    
    updateFrameInfo();
}

function onFrameChange(frameIndex) {
    document.getElementById('frameSlider').value = frameIndex;
    updateFrameInfo();
}

function updateFrameInfo() {
    const current = trajectory ? trajectory.currentFrame + 1 : 0;
    const total = trajectory ? trajectory.getNumFrames() : 0;
    document.getElementById('frameInfo').textContent = `帧: ${current} / ${total}`;
}

function togglePlay() {
    if (!trajectory || trajectory.getNumFrames() < 2) {
        showToast('请先加载轨迹数据', 'info');
        return;
    }
    
    animationController.toggle();
    
    const playIcon = document.getElementById('playIcon');
    const pauseIcon = document.getElementById('pauseIcon');
    
    if (animationController.isPlaying) {
        playIcon.style.display = 'none';
        pauseIcon.style.display = 'block';
    } else {
        playIcon.style.display = 'block';
        pauseIcon.style.display = 'none';
    }
}

function nextFrame() {
    if (animationController) {
        animationController.nextFrame();
    }
}

function previousFrame() {
    if (animationController) {
        animationController.previousFrame();
    }
}

function seekFrame(value) {
    if (animationController) {
        animationController.seekFrame(value);
    }
}

function setPlaySpeed(value) {
    if (animationController) {
        animationController.setSpeed(value);
    }
}

function setDisplayMode(mode) {
    if (viewer) {
        viewer.setStyle(mode);
    }
    
    document.querySelectorAll('.mode-btn').forEach(btn => {
        btn.classList.remove('active');
    });
    document.querySelector(`.mode-btn[data-mode="${mode}"]`)?.classList.add('active');
}

function setColorScheme(scheme) {
    if (viewer) {
        viewer.setColorScheme(scheme);
    }
    
    document.querySelectorAll('.color-btn').forEach(btn => {
        btn.classList.remove('active');
    });
    document.querySelector(`.color-btn[data-scheme="${scheme}"]`)?.classList.add('active');
}

function runAnalysis(type) {
    if (!trajectory || trajectory.getNumFrames() < 2) {
        showToast('请先加载轨迹数据（至少2帧）', 'info');
        return;
    }
    
    document.getElementById('analysisPanel').classList.remove('collapsed');
    switchTab(type);
    
    switch (type) {
        case 'rmsd':
            runRMSDAnalysis();
            break;
        case 'rmsf':
            runRMSFAnalysis();
            break;
        case 'rdf':
            runRDFAnalysis();
            break;
        case 'hbond':
            runHBondAnalysis();
            break;
        case 'pca':
            runPCAnalysis();
            break;
        case 'fes':
            runFESAnalysis();
            break;
    }
}

function runRMSDAnalysis() {
    showLoading('正在计算 RMSD...');
    
    setTimeout(() => {
        const rmsdValues = analysis.calculateRMSD('backbone', 0);
        chartManager.updateRMSDChart('rmsdChart', rmsdValues);
        
        const stats = analysis.getStatistics(rmsdValues);
        document.getElementById('rmsdAvg').textContent = stats.avg.toFixed(3) + ' Å';
        document.getElementById('rmsdMax').textContent = stats.max.toFixed(3) + ' Å';
        document.getElementById('rmsdFinal').textContent = stats.final.toFixed(3) + ' Å';
        document.getElementById('rmsdFluct').textContent = stats.std.toFixed(3) + ' Å';
        
        hideLoading();
        showToast('RMSD 分析完成', 'success');
    }, 100);
}

function runRMSFAnalysis() {
    showLoading('正在计算 RMSF...');
    
    setTimeout(() => {
        const byResidue = document.getElementById('rmsfByResidue').checked;
        const rmsfValues = analysis.calculateRMSF('backbone', byResidue);
        chartManager.updateRMSFChart('rmsfChart', rmsfValues);
        
        hideLoading();
        showToast('RMSF 分析完成', 'success');
    }, 100);
}

function updateRMSFChart() {
    if (!trajectory || trajectory.getNumFrames() < 2) return;
    
    const byResidue = document.getElementById('rmsfByResidue').checked;
    const rmsfValues = analysis.calculateRMSF('backbone', byResidue);
    chartManager.updateRMSFChart('rmsfChart', rmsfValues);
}

function runRDFAnalysis() {
    showLoading('正在计算 RDF...');
    
    setTimeout(() => {
        const maxDist = parseFloat(document.getElementById('rdfMaxDist').value) || 10;
        const reference = document.getElementById('rdfReference').value;
        
        const { distances, rdf } = analysis.calculateRDF(reference, reference, maxDist, 0.1);
        chartManager.updateRDFChart('rdfChart', distances, rdf);
        
        hideLoading();
        showToast('RDF 分析完成', 'success');
    }, 100);
}

function recalculateRDF() {
    runRDFAnalysis();
}

function runHBondAnalysis() {
    showLoading('正在分析氢键...');
    
    setTimeout(() => {
        const donorRes = document.getElementById('hbondDonors').value;
        const acceptorRes = document.getElementById('hbondAcceptors').value;
        
        const donorList = donorRes ? donorRes.split(',').map(r => r.trim().toUpperCase()) : null;
        const acceptorList = acceptorRes ? acceptorRes.split(',').map(r => r.trim().toUpperCase()) : null;
        
        const { pairs, perFrameCounts } = analysis.analyzeHydrogenBonds(donorList, acceptorList, 3.5, 120);
        
        chartManager.updateHBondChart('hbondChart', perFrameCounts);
        
        const hbondList = document.getElementById('hbondList');
        hbondList.innerHTML = '';
        
        if (pairs.length > 0) {
            pairs.slice(0, 10).forEach(pair => {
                const item = document.createElement('div');
                item.className = 'hbond-item';
                item.innerHTML = `
                    <span class="hbond-pair">${pair.pairLabel}</span>
                    <span class="hbond-frequency">${pair.frequency}%</span>
                `;
                hbondList.appendChild(item);
            });
        } else {
            hbondList.innerHTML = '<div class="hbond-item"><span style="color: var(--text-muted)">未找到氢键</span></div>';
        }
        
        hideLoading();
        showToast(`氢键分析完成，发现 ${pairs.length} 对氢键`, 'success');
    }, 100);
}

function analyzeHBonds() {
    runHBondAnalysis();
}

function runPCAnalysis() {
    showLoading('正在执行 PCA 分析...');
    
    setTimeout(() => {
        const { eigenvalues, eigenvectors, projections, varianceExplained } = analysis.performPCA('backbone');
        
        if (projections.length > 0) {
            chartManager.updatePCAChart('pcaChart', projections);
            
            if (varianceExplained.length > 0) {
                document.getElementById('pc1Value').textContent = varianceExplained[0] + '%';
                document.getElementById('pc1Progress').style.width = varianceExplained[0] + '%';
            }
            if (varianceExplained.length > 1) {
                document.getElementById('pc2Value').textContent = varianceExplained[1] + '%';
                document.getElementById('pc2Progress').style.width = varianceExplained[1] + '%';
            }
            if (varianceExplained.length > 2) {
                document.getElementById('pc3Value').textContent = varianceExplained[2] + '%';
                document.getElementById('pc3Progress').style.width = varianceExplained[2] + '%';
            }
        }
        
        hideLoading();
        showToast('PCA 分析完成', 'success');
    }, 100);
}

let currentPCMode = null;

function visualizePCMode(pcIndex) {
    if (!trajectory || trajectory.getNumFrames() < 3) {
        showToast('请先加载轨迹数据（至少3帧）', 'info');
        return;
    }

    if (currentPCMode === pcIndex) {
        stopPCModeAnimation();
        return;
    }

    if (currentPCMode !== null) {
        stopPCModeAnimation();
    }

    showLoading(`正在生成 PC${pcIndex} 运动模式...`);

    setTimeout(() => {
        if (!analysis.results.pca) {
            analysis.performPCA('backbone');
        }

        const pcModeData = analysis.generatePCModeFrames(pcIndex - 1, 30, 2.5);
        
        if (pcModeData) {
            viewer.startPCModeAnimation(pcModeData.frames, pcModeData.atomIndices);
            currentPCMode = pcIndex;
            
            document.querySelectorAll('.pca-controls .btn').forEach(btn => {
                btn.classList.remove('active');
            });
            event?.target?.classList.add('active');
            
            showToast(`PC${pcIndex} 运动模式已启动，点击再次按钮停止`, 'success');
        } else {
            showToast('PC模式生成失败', 'error');
        }
        
        hideLoading();
    }, 100);
}

function stopPCModeAnimation() {
    if (viewer) {
        viewer.stopPCModeAnimation();
    }
    currentPCMode = null;
    
    document.querySelectorAll('.pca-controls .btn').forEach(btn => {
        btn.classList.remove('active');
    });
    
    showToast('PC模式动画已停止', 'info');
}

let customFESData = null;

function runFESAnalysis() {
    showLoading('正在计算自由能面...');
    
    setTimeout(() => {
        const xAxis = document.getElementById('fesXAxis').value;
        const yAxis = document.getElementById('fesYAxis').value;
        const temperature = parseFloat(document.getElementById('fesTemperature').value) || 300;
        
        let xData, yData;
        
        if (xAxis === 'custom' || yAxis === 'custom') {
            if (!customFESData) {
                showToast('请先导入自定义数据文件', 'warning');
                hideLoading();
                return;
            }
            xData = customFESData.x;
            yData = customFESData.y;
        } else {
            xData = analysis.getReactionCoordinate(xAxis);
            yData = analysis.getReactionCoordinate(yAxis);
        }
        
        if (xData.length === 0 || yData.length === 0) {
            showToast('无法获取反应坐标数据', 'error');
            hideLoading();
            return;
        }
        
        const minLen = Math.min(xData.length, yData.length);
        xData = xData.slice(0, minLen);
        yData = yData.slice(0, minLen);
        
        const { x, y, fes, minFes } = analysis.calculateFreeEnergySurface(xData, yData, temperature);
        chartManager.updateFESChart('fesChart', x, y, fes);
        
        const fesChart = chartManager.getChart('fesChart');
        if (fesChart) {
            fesChart.options.scales.x.title.text = getAxisLabel(xAxis);
            fesChart.options.scales.y.title.text = getAxisLabel(yAxis);
            fesChart.update('none');
        }
        
        hideLoading();
        showToast('自由能面计算完成', 'success');
    }, 100);
}

function getAxisLabel(axis) {
    const labels = {
        'pc1': 'PC1',
        'pc2': 'PC2',
        'pc3': 'PC3',
        'rmsd': 'RMSD (Å)',
        'rg': '回转半径 (Å)',
        'rmsf': 'RMSF (Å)',
        'custom': '自定义'
    };
    return labels[axis] || axis;
}

function updateFESAxisLabel() {
}

function loadFESData(event) {
    const file = event.target.files[0];
    if (!file) return;
    
    const reader = new FileReader();
    
    reader.onload = function(e) {
        const content = e.target.result;
        const lines = content.trim().split('\n');
        const xData = [];
        const yData = [];
        
        for (const line of lines) {
            const trimmed = line.trim();
            if (trimmed === '' || trimmed.startsWith('#') || trimmed.startsWith('@')) {
                continue;
            }
            
            const parts = trimmed.split(/[\s,]+/);
            if (parts.length >= 2) {
                const x = parseFloat(parts[0]);
                const y = parseFloat(parts[1]);
                if (!isNaN(x) && !isNaN(y)) {
                    xData.push(x);
                    yData.push(y);
                }
            }
        }
        
        if (xData.length > 0) {
            customFESData = { x: xData, y: yData, fileName: file.name };
            document.getElementById('fesDataStatus').textContent = `${file.name} (${xData.length} 点)`;
            document.getElementById('fesDataStatus').style.color = 'var(--success-color)';
            
            document.getElementById('fesXAxis').value = 'custom';
            document.getElementById('fesYAxis').value = 'custom';
            
            showToast(`成功导入 ${xData.length} 个数据点`, 'success');
        } else {
            showToast('未找到有效数据', 'error');
        }
    };
    
    reader.readAsText(file);
}

function calculateFES() {
    runFESAnalysis();
}

function switchTab(tabName) {
    document.querySelectorAll('.tab-btn').forEach(btn => {
        btn.classList.remove('active');
    });
    document.querySelector(`.tab-btn[data-tab="${tabName}"]`)?.classList.add('active');
    
    document.querySelectorAll('.tab-content').forEach(content => {
        content.classList.remove('active');
    });
    document.getElementById(`tab-${tabName}`)?.classList.add('active');
}

function closeAnalysisPanel() {
    document.getElementById('analysisPanel').classList.add('collapsed');
}

function addTrajectoryForComparison() {
    if (!trajectory || trajectory.getNumFrames() === 0) {
        showToast('请先加载主轨迹', 'info');
        return;
    }
    
    const index = trajectoryComparison.getNumTrajectories();
    const newTraj = new Trajectory();
    newTraj.generateSampleTrajectory(50);
    
    trajectoryComparison.addTrajectory(newTraj, `Trajectory ${index + 1}`);
    updateCompareList();
    
    showToast('对比轨迹已添加', 'success');
    
    if (index >= 1) {
        compareTrajectories();
    }
}

function updateCompareList() {
    const compareList = document.getElementById('compareList');
    const trajs = trajectoryComparison.getAllTrajectories();
    
    if (trajs.length === 0) {
        compareList.innerHTML = '<div class="compare-item empty"><span>添加多条轨迹进行对比</span></div>';
        return;
    }
    
    compareList.innerHTML = '';
    trajs.forEach((traj, index) => {
        const item = document.createElement('div');
        item.className = 'compare-item';
        item.innerHTML = `
            <span class="color-indicator" style="background: ${traj.color}"></span>
            <span class="compare-name">${traj.name}</span>
        `;
        compareList.appendChild(item);
    });
}

function compareTrajectories() {
    const trajs = trajectoryComparison.getAllTrajectories();
    
    const datasets = [];
    
    const mainRMSD = analysis.calculateRMSD('backbone', 0);
    datasets.push({ name: '主轨迹', data: mainRMSD });
    
    trajs.forEach(traj => {
        const tempAnalysis = new MDAnalysis(traj.trajectory);
        const rmsd = tempAnalysis.calculateRMSD('backbone', 0);
        datasets.push({ name: traj.name, data: rmsd });
    });
    
    switchTab('rmsd');
    chartManager.updateRMSDComparison('rmsdChart', datasets);
    
    showToast('轨迹对比完成', 'success');
}

function showHelp() {
    document.getElementById('helpModal').style.display = 'flex';
}

function closeHelp() {
    document.getElementById('helpModal').style.display = 'none';
}

function showToast(message, type = 'info') {
    const toast = document.createElement('div');
    toast.className = `toast ${type}`;
    toast.textContent = message;
    document.body.appendChild(toast);
    
    setTimeout(() => {
        toast.style.animation = 'slideIn 0.3s ease reverse';
        setTimeout(() => toast.remove(), 300);
    }, 3000);
}

function showLoading(text = '加载中...') {
    let overlay = document.querySelector('.loading-overlay');
    if (!overlay) {
        overlay = document.createElement('div');
        overlay.className = 'loading-overlay';
        overlay.innerHTML = `
            <div style="text-align: center; color: var(--text-secondary);">
                <div class="loading-spinner" style="margin: 0 auto 10px;"></div>
                <span class="loading-text">${text}</span>
            </div>
        `;
        document.querySelector('.viewer-container')?.appendChild(overlay);
    } else {
        overlay.querySelector('.loading-text').textContent = text;
    }
}

function hideLoading() {
    const overlay = document.querySelector('.loading-overlay');
    if (overlay) {
        overlay.remove();
    }
}

window.addEventListener('resize', function() {
    if (viewer) {
        viewer.resize();
    }
});
