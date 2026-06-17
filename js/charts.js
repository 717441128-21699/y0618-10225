class ChartManager {
    constructor() {
        this.charts = {};
        this.colors = [
            '#6366f1', '#10b981', '#f59e0b', '#ef4444',
            '#8b5cf6', '#06b6d4', '#ec4899', '#84cc16'
        ];
    }

    createRMSDChart(canvasId) {
        const ctx = document.getElementById(canvasId);
        if (!ctx) return null;

        const chart = new Chart(ctx, {
            type: 'line',
            data: {
                labels: [],
                datasets: [{
                    label: 'RMSD (Å)',
                    data: [],
                    borderColor: '#6366f1',
                    backgroundColor: 'rgba(99, 102, 241, 0.1)',
                    fill: true,
                    tension: 0.3,
                    pointRadius: 0,
                    pointHoverRadius: 4,
                    borderWidth: 2
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                interaction: {
                    intersect: false,
                    mode: 'index'
                },
                onClick: (event, elements) => {
                    if (elements && elements.length > 0 && typeof window.onChartFrameClick === 'function') {
                        window.onChartFrameClick(canvasId, elements[0].index);
                    }
                },
                plugins: {
                    legend: {
                        display: false
                    },
                    tooltip: {
                        backgroundColor: 'rgba(30, 41, 59, 0.9)',
                        titleColor: '#f1f5f9',
                        bodyColor: '#94a3b8',
                        borderColor: '#334155',
                        borderWidth: 1,
                        padding: 10,
                        cornerRadius: 6,
                        callbacks: {
                            title: (items) => `帧: ${items[0].label}`,
                            label: (item) => `RMSD: ${item.raw.toFixed(3)} Å`
                        }
                    }
                },
                scales: {
                    x: {
                        title: {
                            display: true,
                            text: '帧 / 模拟时间 (ns)',
                            color: '#94a3b8',
                            font: { size: 11 }
                        },
                        grid: {
                            color: 'rgba(51, 65, 85, 0.5)',
                            drawBorder: false
                        },
                        ticks: {
                            color: '#64748b',
                            font: { size: 10 }
                        }
                    },
                    y: {
                        title: {
                            display: true,
                            text: 'RMSD (Å)',
                            color: '#94a3b8',
                            font: { size: 11 }
                        },
                        grid: {
                            color: 'rgba(51, 65, 85, 0.5)',
                            drawBorder: false
                        },
                        ticks: {
                            color: '#64748b',
                            font: { size: 10 }
                        },
                        beginAtZero: true
                    }
                }
            }
        });

        this.charts[canvasId] = chart;
        return chart;
    }

    createRMSFChart(canvasId) {
        const ctx = document.getElementById(canvasId);
        if (!ctx) return null;

        const chart = new Chart(ctx, {
            type: 'bar',
            data: {
                labels: [],
                datasets: [{
                    label: 'RMSF (Å)',
                    data: [],
                    backgroundColor: function(context) {
                        const value = context.raw;
                        if (value === undefined) return '#6366f1';
                        if (value > 3) return '#ef4444';
                        if (value > 2) return '#f59e0b';
                        return '#10b981';
                    },
                    borderRadius: 2,
                    borderSkipped: false
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                plugins: {
                    legend: {
                        display: false
                    },
                    tooltip: {
                        backgroundColor: 'rgba(30, 41, 59, 0.9)',
                        titleColor: '#f1f5f9',
                        bodyColor: '#94a3b8',
                        borderColor: '#334155',
                        borderWidth: 1,
                        padding: 10,
                        cornerRadius: 6
                    }
                },
                scales: {
                    x: {
                        title: {
                            display: true,
                            text: '残基',
                            color: '#94a3b8',
                            font: { size: 11 }
                        },
                        grid: {
                            display: false
                        },
                        ticks: {
                            color: '#64748b',
                            font: { size: 10 },
                            maxTicksLimit: 20
                        }
                    },
                    y: {
                        title: {
                            display: true,
                            text: 'RMSF (Å)',
                            color: '#94a3b8',
                            font: { size: 11 }
                        },
                        grid: {
                            color: 'rgba(51, 65, 85, 0.5)',
                            drawBorder: false
                        },
                        ticks: {
                            color: '#64748b',
                            font: { size: 10 }
                        },
                        beginAtZero: true
                    }
                }
            }
        });

        this.charts[canvasId] = chart;
        return chart;
    }

    createRDFChart(canvasId) {
        const ctx = document.getElementById(canvasId);
        if (!ctx) return null;

        const chart = new Chart(ctx, {
            type: 'line',
            data: {
                labels: [],
                datasets: [{
                    label: 'g(r)',
                    data: [],
                    borderColor: '#10b981',
                    backgroundColor: 'rgba(16, 185, 129, 0.1)',
                    fill: true,
                    tension: 0.1,
                    pointRadius: 0,
                    borderWidth: 2
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                interaction: {
                    intersect: false,
                    mode: 'index'
                },
                plugins: {
                    legend: {
                        display: false
                    },
                    tooltip: {
                        backgroundColor: 'rgba(30, 41, 59, 0.9)',
                        titleColor: '#f1f5f9',
                        bodyColor: '#94a3b8',
                        borderColor: '#334155',
                        borderWidth: 1,
                        padding: 10,
                        cornerRadius: 6,
                        callbacks: {
                            title: (items) => `距离: ${items[0].label} Å`,
                            label: (item) => `g(r): ${item.raw.toFixed(3)}`
                        }
                    }
                },
                scales: {
                    x: {
                        title: {
                            display: true,
                            text: '距离 r (Å)',
                            color: '#94a3b8',
                            font: { size: 11 }
                        },
                        grid: {
                            color: 'rgba(51, 65, 85, 0.5)',
                            drawBorder: false
                        },
                        ticks: {
                            color: '#64748b',
                            font: { size: 10 }
                        }
                    },
                    y: {
                        title: {
                            display: true,
                            text: 'g(r)',
                            color: '#94a3b8',
                            font: { size: 11 }
                        },
                        grid: {
                            color: 'rgba(51, 65, 85, 0.5)',
                            drawBorder: false
                        },
                        ticks: {
                            color: '#64748b',
                            font: { size: 10 }
                        },
                        beginAtZero: true
                    }
                }
            }
        });

        this.charts[canvasId] = chart;
        return chart;
    }

    createHBondChart(canvasId) {
        const ctx = document.getElementById(canvasId);
        if (!ctx) return null;

        const chart = new Chart(ctx, {
            type: 'line',
            data: {
                labels: [],
                datasets: [{
                    label: '氢键数量',
                    data: [],
                    borderColor: '#f59e0b',
                    backgroundColor: 'rgba(245, 158, 11, 0.1)',
                    fill: true,
                    tension: 0.3,
                    pointRadius: 0,
                    borderWidth: 2
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                interaction: {
                    intersect: false,
                    mode: 'index'
                },
                plugins: {
                    legend: {
                        display: false
                    },
                    tooltip: {
                        backgroundColor: 'rgba(30, 41, 59, 0.9)',
                        titleColor: '#f1f5f9',
                        bodyColor: '#94a3b8',
                        borderColor: '#334155',
                        borderWidth: 1,
                        padding: 10,
                        cornerRadius: 6
                    }
                },
                scales: {
                    x: {
                        title: {
                            display: true,
                            text: '帧',
                            color: '#94a3b8',
                            font: { size: 11 }
                        },
                        grid: {
                            color: 'rgba(51, 65, 85, 0.5)',
                            drawBorder: false
                        },
                        ticks: {
                            color: '#64748b',
                            font: { size: 10 }
                        }
                    },
                    y: {
                        title: {
                            display: true,
                            text: '氢键数量',
                            color: '#94a3b8',
                            font: { size: 11 }
                        },
                        grid: {
                            color: 'rgba(51, 65, 85, 0.5)',
                            drawBorder: false
                        },
                        ticks: {
                            color: '#64748b',
                            font: { size: 10 }
                        },
                        beginAtZero: true
                    }
                }
            }
        });

        this.charts[canvasId] = chart;
        return chart;
    }

    createPCAChart(canvasId) {
        const ctx = document.getElementById(canvasId);
        if (!ctx) return null;

        const chart = new Chart(ctx, {
            type: 'scatter',
            data: {
                datasets: [{
                    label: '构象',
                    data: [],
                    backgroundColor: function(context) {
                        const index = context.dataIndex;
                        const total = context.dataset.data.length;
                        const hue = (index / total) * 360;
                        return `hsla(${hue}, 70%, 60%, 0.8)`;
                    },
                    pointRadius: 4,
                    pointHoverRadius: 6
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                onClick: (event, elements) => {
                    if (elements && elements.length > 0 && typeof window.onChartFrameClick === 'function') {
                        window.onChartFrameClick(canvasId, elements[0].index);
                    }
                },
                plugins: {
                    legend: {
                        display: false
                    },
                    tooltip: {
                        backgroundColor: 'rgba(30, 41, 59, 0.9)',
                        titleColor: '#f1f5f9',
                        bodyColor: '#94a3b8',
                        borderColor: '#334155',
                        borderWidth: 1,
                        padding: 10,
                        cornerRadius: 6,
                        callbacks: {
                            title: (items) => `帧: ${items[0].dataIndex + 1}`,
                            label: (item) => [
                                `PC1: ${item.raw.x.toFixed(2)}`,
                                `PC2: ${item.raw.y.toFixed(2)}`
                            ]
                        }
                    }
                },
                scales: {
                    x: {
                        title: {
                            display: true,
                            text: 'PC1',
                            color: '#94a3b8',
                            font: { size: 11 }
                        },
                        grid: {
                            color: 'rgba(51, 65, 85, 0.5)',
                            drawBorder: false
                        },
                        ticks: {
                            color: '#64748b',
                            font: { size: 10 }
                        }
                    },
                    y: {
                        title: {
                            display: true,
                            text: 'PC2',
                            color: '#94a3b8',
                            font: { size: 11 }
                        },
                        grid: {
                            color: 'rgba(51, 65, 85, 0.5)',
                            drawBorder: false
                        },
                        ticks: {
                            color: '#64748b',
                            font: { size: 10 }
                        }
                    }
                }
            }
        });

        this.charts[canvasId] = chart;
        return chart;
    }

    createFESChart(canvasId) {
        const ctx = document.getElementById(canvasId);
        if (!ctx) return null;

        const chart = new Chart(ctx, {
            type: 'scatter',
            data: {
                datasets: [{
                    label: '自由能 (kcal/mol)',
                    data: [],
                    backgroundColor: function(context) {
                        const value = context.raw?.v;
                        if (value === undefined) return '#6366f1';
                        const maxV = 5;
                        const ratio = Math.min(Math.max(value / maxV, 0), 1);
                        const r = Math.round(16 + ratio * 239);
                        const g = Math.round(185 * (1 - ratio) + 68 * ratio);
                        const b = Math.round(253 * (1 - ratio));
                        return `rgba(${r}, ${g}, ${b}, 0.8)`;
                    },
                    pointRadius: 6,
                    pointHoverRadius: 8,
                    pointStyle: 'rectRot'
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                plugins: {
                    legend: {
                        display: false
                    },
                    tooltip: {
                        backgroundColor: 'rgba(30, 41, 59, 0.9)',
                        titleColor: '#f1f5f9',
                        bodyColor: '#94a3b8',
                        borderColor: '#334155',
                        borderWidth: 1,
                        padding: 10,
                        cornerRadius: 6,
                        callbacks: {
                            label: function(context) {
                                const raw = context.raw;
                                if (raw && raw.v !== undefined) {
                                    return `自由能: ${raw.v.toFixed(2)} kcal/mol`;
                                }
                                return '';
                            }
                        }
                    }
                },
                scales: {
                    x: {
                        title: {
                            display: true,
                            text: 'PC1',
                            color: '#94a3b8',
                            font: { size: 11 }
                        },
                        grid: {
                            color: 'rgba(51, 65, 85, 0.5)',
                            drawBorder: false
                        },
                        ticks: {
                            color: '#64748b',
                            font: { size: 10 }
                        }
                    },
                    y: {
                        title: {
                            display: true,
                            text: 'PC2',
                            color: '#94a3b8',
                            font: { size: 11 }
                        },
                        grid: {
                            color: 'rgba(51, 65, 85, 0.5)',
                            drawBorder: false
                        },
                        ticks: {
                            color: '#64748b',
                            font: { size: 10 }
                        }
                    }
                }
            }
        });

        this.charts[canvasId] = chart;
        return chart;
    }

    updateRMSDChart(canvasId, values, timeStep = 0.001) {
        const chart = this.charts[canvasId];
        if (!chart) return;

        const labels = values.map((_, i) => i);
        
        chart.data.labels = labels;
        chart.data.datasets[0].data = values;
        chart.update('none');
    }

    updateRMSDComparison(canvasId, datasets) {
        const chart = this.charts[canvasId];
        if (!chart) return;

        const maxLen = Math.max(...datasets.map(d => d.data.length));
        const labels = Array.from({ length: maxLen }, (_, i) => i);

        chart.data.labels = labels;
        chart.data.datasets = datasets.map((ds, i) => ({
            label: ds.name,
            data: ds.data,
            borderColor: this.colors[i % this.colors.length],
            backgroundColor: this.colors[i % this.colors.length] + '20',
            fill: false,
            tension: 0.3,
            pointRadius: 0,
            borderWidth: 2
        }));
        chart.update('none');
    }

    updateRMSFChart(canvasId, data) {
        const chart = this.charts[canvasId];
        if (!chart) return;

        if (Array.isArray(data) && data.length > 0 && typeof data[0] === 'object') {
            chart.data.labels = data.map(d => d.resi);
            chart.data.datasets[0].data = data.map(d => d.value);
        } else {
            chart.data.labels = data.map((_, i) => i + 1);
            chart.data.datasets[0].data = data;
        }
        chart.update('none');
    }

    updateRDFChart(canvasId, distances, rdf) {
        const chart = this.charts[canvasId];
        if (!chart) return;

        chart.data.labels = distances.map(d => d.toFixed(1));
        chart.data.datasets[0].data = rdf;
        chart.update('none');
    }

    updateHBondChart(canvasId, hbondCount) {
        const chart = this.charts[canvasId];
        if (!chart) return;

        const labels = hbondCount.map((_, i) => i);
        chart.data.labels = labels;
        chart.data.datasets[0].data = hbondCount;
        chart.update('none');
    }

    updatePCAChart(canvasId, projections) {
        const chart = this.charts[canvasId];
        if (!chart || projections.length === 0) return;

        const data = projections.map(p => ({ x: p[0], y: p[1] }));
        chart.data.datasets[0].data = data;
        chart.update('none');
    }

    updateFESChart(canvasId, xData, yData, fesData) {
        const chart = this.charts[canvasId];
        if (!chart) return;

        const scatterData = [];
        for (let i = 0; i < fesData.length; i++) {
            for (let j = 0; j < fesData[i].length; j++) {
                if (fesData[i][j] !== null && fesData[i][j] !== undefined) {
                    scatterData.push({
                        x: xData[i],
                        y: yData[j],
                        v: fesData[i][j]
                    });
                }
            }
        }

        chart.data.labels = xData;
        chart.data.datasets = [{
            label: '自由能 (kcal/mol)',
            data: scatterData,
            backgroundColor: function(context) {
                const value = context.raw?.v;
                if (value === undefined) return '#6366f1';
                const maxV = 5;
                const ratio = Math.min(value / maxV, 1);
                const r = Math.round(239 + ratio * 16);
                const g = Math.round(68 * (1 - ratio) + 153 * ratio);
                const b = Math.round(68 * (1 - ratio));
                return `rgba(${r}, ${g}, ${b}, 0.7)`;
            },
            pointRadius: 5,
            pointHoverRadius: 7
        }];
        chart.update('none');
    }

    getChart(canvasId) {
        return this.charts[canvasId];
    }

    destroyChart(canvasId) {
        if (this.charts[canvasId]) {
            this.charts[canvasId].destroy();
            delete this.charts[canvasId];
        }
    }

    destroyAll() {
        Object.keys(this.charts).forEach(id => {
            this.charts[id].destroy();
        });
        this.charts = {};
    }
}
