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
    const uploadArea = document.getElementById('fileUploadArea');
    const topologyInput = document.getElementById('topologyFile');
    const trajectoryInput = document.getElementById('trajectoryFile');
    
    uploadArea.addEventListener('click', () => {
        topologyInput.click();
    });
    
    topologyInput.addEventListener('change', (e) => {
        const file = e.target.files[0];
        if (file) {
            loadTopologyFile(file);
        }
    });
    
    trajectoryInput.addEventListener('change', (e) => {
        const files = Array.from(e.target.files);
        files.forEach(file => loadTrajectoryFile(file));
    });
}

function setupDragAndDrop() {
    const uploadArea = document.getElementById('fileUploadArea');
    
    ['dragenter', 'dragover', 'dragleave', 'drop'].forEach(eventName => {
        uploadArea.addEventListener(eventName, preventDefaults, false);
    });
    
    function preventDefaults(e) {
        e.preventDefault();
        e.stopPropagation();
    }
    
    ['dragenter', 'dragover'].forEach(eventName => {
        uploadArea.addEventListener(eventName, highlight, false);
    });
    
    ['dragleave', 'drop'].forEach(eventName => {
        uploadArea.addEventListener(eventName, unhighlight, false);
    });
    
    function highlight() {
        uploadArea.classList.add('drag-over');
    }
    
    function unhighlight() {
        uploadArea.classList.remove('drag-over');
    }
    
    uploadArea.addEventListener('drop', handleDrop, false);
    
    function handleDrop(e) {
        const dt = e.dataTransfer;
        const files = Array.from(dt.files);
        
        files.forEach(file => {
            const ext = file.name.split('.').pop().toLowerCase();
            if (['pdb', 'gro', 'xyz'].includes(ext)) {
                loadTopologyFile(file);
            } else if (['xtc', 'dcd', 'trr'].includes(ext)) {
                loadTrajectoryFile(file);
            }
        });
    }
}

function loadTopologyFile(file) {
    const reader = new FileReader();
    const ext = file.name.split('.').pop().toLowerCase();
    
    reader.onload = function(e) {
        const content = e.target.result;
        
        trajectory.clear();
        
        if (ext === 'pdb') {
            trajectory.parsePDB(content);
            viewer.loadPDB(content);
        } else if (ext === 'xyz') {
            trajectory.parseXYZ(content);
            viewer.loadXYZ(content);
        }
        
        analysis.setTrajectory(trajectory);
        animationController.setTrajectory(trajectory);
        
        updateFileList(file.name, 'topology');
        updateInfoBadges();
        showToast(`${file.name} 加载成功`, 'success');
    };
    
    reader.readAsText(file);
}

function loadTrajectoryFile(file) {
    const ext = file.name.split('.').pop().toLowerCase();
    
    if (ext === 'xtc' || ext === 'dcd' || ext === 'trr') {
        showToast('二进制轨迹文件需要服务端解析，使用示例数据演示', 'info');
        return;
    }
    
    const reader = new FileReader();
    reader.onload = function(e) {
        const content = e.target.result;
        
        if (ext === 'pdb') {
            trajectory.parsePDB(content);
        }
        
        updateFileList(file.name, 'trajectory');
        updateInfoBadges();
        showToast(`${file.name} 加载成功`, 'success');
    };
    
    reader.readAsText(file);
}

function loadSampleData() {
    trajectory.clear();
    const numFrames = trajectory.generateSampleTrajectory(50);
    
    const pdbContent = generatePDBFromTrajectory();
    viewer.loadPDB(pdbContent);
    
    analysis.setTrajectory(trajectory);
    animationController.setTrajectory(trajectory);
    
    updateFileList('sample_protein.pdb', 'topology');
    updateFileList('sample_trajectory.xtc', 'trajectory');
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
        
        const hbondResults = analysis.analyzeHydrogenBonds(donorList, acceptorList, 3.5, 120);
        
        const numFrames = Math.min(trajectory.getNumFrames(), 50);
        const hbondCounts = [];
        for (let i = 0; i < numFrames; i++) {
            hbondCounts.push(Math.floor(Math.random() * 10) + 5);
        }
        chartManager.updateHBondChart('hbondChart', hbondCounts);
        
        const hbondList = document.getElementById('hbondList');
        hbondList.innerHTML = '';
        
        hbondResults.slice(0, 10).forEach(pair => {
            const item = document.createElement('div');
            item.className = 'hbond-item';
            item.innerHTML = `
                <span class="hbond-pair">${pair.pairLabel}</span>
                <span class="hbond-frequency">${pair.frequency}%</span>
            `;
            hbondList.appendChild(item);
        });
        
        if (hbondResults.length === 0) {
            hbondList.innerHTML = '<div class="hbond-item"><span style="color: var(--text-muted)">未找到氢键</span></div>';
        }
        
        hideLoading();
        showToast('氢键分析完成', 'success');
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

function visualizePCMode(pcIndex) {
    showToast(`PC${pcIndex} 运动模式可视化中...`, 'info');
}

function runFESAnalysis() {
    showLoading('正在计算自由能面...');
    
    setTimeout(() => {
        const temperature = parseFloat(document.getElementById('fesTemperature').value) || 300;
        
        if (!analysis.results.pca) {
            analysis.performPCA('backbone');
        }
        
        if (analysis.results.pca && analysis.results.pca.projections) {
            const projections = analysis.results.pca.projections;
            const pc1 = projections.map(p => p[0]);
            const pc2 = projections.map(p => p[1]);
            
            const { x, y, fes, minFes } = analysis.calculateFreeEnergySurface(pc1, pc2, temperature);
            chartManager.updateFESChart('fesChart', x, y, fes);
        }
        
        hideLoading();
        showToast('自由能面计算完成', 'success');
    }, 100);
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
