#!/usr/bin/env python3
"""
MutationScan Web GUI - Flask Backend
Provides a local web interface for running the MutationScan pipeline.
"""

import os
import json
import time
import threading
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional
from queue import Queue

from flask import Flask, render_template, request, jsonify, Response, send_file, session
from werkzeug.utils import secure_filename
import secrets

# Initialize Flask app
app = Flask(__name__)
app.secret_key = secrets.token_hex(16)
app.config['UPLOAD_FOLDER'] = Path('uploads')
app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024  # 50MB max file size

# Ensure upload directory exists
app.config['UPLOAD_FOLDER'].mkdir(exist_ok=True)

# Global state for pipeline execution
pipeline_state = {
    'status': 'idle',  # idle, running, completed, error
    'current_stage': None,
    'current_genome': None,
    'progress': 0,
    'total_genomes': 0,
    'stages': [
        {'id': 1, 'name': 'Genome Harvester', 'status': 'pending', 'progress': 0},
        {'id': 2, 'name': 'Gene Annotator', 'status': 'pending', 'progress': 0},
        {'id': 3, 'name': 'Sequence Extractor', 'status': 'pending', 'progress': 0},
        {'id': 4, 'name': 'Wild-type Aligner', 'status': 'pending', 'progress': 0},
        {'id': 5, 'name': 'Mutation Analyzer', 'status': 'pending', 'progress': 0},
        {'id': 6, 'name': 'Co-occurrence Analyzer', 'status': 'pending', 'progress': 0},
        {'id': 7, 'name': 'Report Generator', 'status': 'pending', 'progress': 0},
    ],
    'start_time': None,
    'end_time': None,
    'error_message': None,
    'log_messages': [],
    'output_dir': None,
    'report_path': None
}

# Queue for server-sent events
event_queue = Queue()


# ============================================================================
# ROUTES - Main Pages
# ============================================================================

@app.route('/')
def index():
    """Redirect to setup page"""
    return render_template('setup.html')


@app.route('/setup')
def setup():
    """Setup screen for configuring pipeline parameters"""
    return render_template('setup.html')


@app.route('/pipeline')
def pipeline():
    """Pipeline execution screen with live progress tracking"""
    return render_template('pipeline.html')


@app.route('/results')
def results():
    """Results screen for viewing and downloading reports"""
    return render_template('results.html')


# ============================================================================
# API ENDPOINTS - Configuration & File Upload
# ============================================================================

@app.route('/api/health', methods=['GET'])
def health_check():
    """Check system health and dependencies"""
    health_status = {
        'python': True,
        'emboss': check_command_exists('needle'),
        'abricate': check_command_exists('abricate'),
        'disk_space_gb': get_disk_space(),
        'internet': check_internet_connection(),
        'timestamp': datetime.now().isoformat()
    }
    
    return jsonify(health_status)


@app.route('/api/upload/accessions', methods=['POST'])
def upload_accessions():
    """Upload accession list file"""
    if 'file' not in request.files:
        return jsonify({'error': 'No file provided'}), 400
    
    file = request.files['file']
    if file.filename == '':
        return jsonify({'error': 'No file selected'}), 400
    
    if file and allowed_file(file.filename, ['txt', 'csv']):
        filename = secure_filename(file.filename)
        filepath = app.config['UPLOAD_FOLDER'] / filename
        file.save(filepath)
        
        # Read and validate accessions
        try:
            with open(filepath, 'r') as f:
                accessions = [line.strip() for line in f if line.strip()]
            
            return jsonify({
                'success': True,
                'filepath': str(filepath),
                'count': len(accessions),
                'preview': accessions[:5]
            })
        except Exception as e:
            return jsonify({'error': f'Failed to read file: {str(e)}'}), 400
    
    return jsonify({'error': 'Invalid file type'}), 400


@app.route('/api/upload/genes', methods=['POST'])
def upload_genes():
    """Upload gene list file"""
    if 'file' not in request.files:
        return jsonify({'error': 'No file provided'}), 400
    
    file = request.files['file']
    if file.filename == '':
        return jsonify({'error': 'No file selected'}), 400
    
    if file and allowed_file(file.filename, ['txt', 'csv']):
        filename = secure_filename(file.filename)
        filepath = app.config['UPLOAD_FOLDER'] / filename
        file.save(filepath)
        
        # Read and validate genes
        try:
            with open(filepath, 'r') as f:
                genes = [line.strip() for line in f if line.strip()]
            
            return jsonify({
                'success': True,
                'filepath': str(filepath),
                'count': len(genes),
                'preview': genes[:10]
            })
        except Exception as e:
            return jsonify({'error': f'Failed to read file: {str(e)}'}), 400
    
    return jsonify({'error': 'Invalid file type'}), 400


@app.route('/api/config/templates', methods=['GET'])
def get_templates():
    """Get available pipeline templates"""
    templates = [
        {
            'id': 'quick-start',
            'name': 'Quick Start',
            'description': 'Basic AMR analysis with default settings',
            'genes': ['gyrA', 'parC', 'blaKPC'],
            'database': 'ncbi'
        },
        {
            'id': 'fluoroquinolone',
            'name': 'Fluoroquinolone Resistance',
            'description': 'Pre-configured for gyrA, parC, gyrB analysis',
            'genes': ['gyrA', 'parC', 'gyrB', 'parE'],
            'database': 'ncbi'
        },
        {
            'id': 'beta-lactamase',
            'name': 'Beta-Lactamase Analysis',
            'description': 'Focused on bla genes',
            'genes': ['blaKPC', 'blaNDM', 'blaOXA', 'blaCTX'],
            'database': 'ncbi'
        }
    ]
    return jsonify(templates)


# ============================================================================
# API ENDPOINTS - Pipeline Execution
# ============================================================================

@app.route('/api/pipeline/start', methods=['POST'])
def start_pipeline():
    """Start pipeline execution with provided configuration"""
    global pipeline_state
    
    # Check if pipeline is already running
    if pipeline_state['status'] == 'running':
        return jsonify({'error': 'Pipeline is already running'}), 400
    
    # Get configuration from request
    config = request.json
    
    # Validate configuration
    if not config.get('accessions_file'):
        return jsonify({'error': 'Accessions file is required'}), 400
    if not config.get('genes_file') and not config.get('genes_list'):
        return jsonify({'error': 'Genes must be provided'}), 400
    
    # Reset pipeline state
    reset_pipeline_state()
    pipeline_state['status'] = 'running'
    pipeline_state['start_time'] = datetime.now().isoformat()
    
    # Create output directory with timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir = Path(config.get('output_dir', 'results')) / f'run_{timestamp}'
    output_dir.mkdir(parents=True, exist_ok=True)
    pipeline_state['output_dir'] = str(output_dir)
    
    # Start pipeline in background thread
    thread = threading.Thread(
        target=run_pipeline_worker,
        args=(config, output_dir),
        daemon=True
    )
    thread.start()
    
    add_log_message('info', 'Pipeline started successfully')
    
    return jsonify({
        'success': True,
        'output_dir': str(output_dir),
        'message': 'Pipeline started successfully'
    })


@app.route('/api/pipeline/status', methods=['GET'])
def get_pipeline_status():
    """Get current pipeline status"""
    return jsonify(pipeline_state)


@app.route('/api/pipeline/stop', methods=['POST'])
def stop_pipeline():
    """Stop pipeline execution"""
    # TODO: Implement graceful pipeline termination
    global pipeline_state
    
    if pipeline_state['status'] == 'running':
        pipeline_state['status'] = 'stopped'
        pipeline_state['end_time'] = datetime.now().isoformat()
        add_log_message('warning', 'Pipeline stopped by user')
        return jsonify({'success': True, 'message': 'Pipeline stopped'})
    
    return jsonify({'error': 'No pipeline is running'}), 400


# ============================================================================
# API ENDPOINTS - Results & Downloads
# ============================================================================

@app.route('/api/results/list', methods=['GET'])
def list_results():
    """List all available result directories"""
    results_dir = Path('results')
    if not results_dir.exists():
        return jsonify([])
    
    results = []
    for run_dir in sorted(results_dir.iterdir(), reverse=True):
        if run_dir.is_dir() and run_dir.name.startswith('run_'):
            report_path = run_dir / 'final_report.html'
            results.append({
                'directory': run_dir.name,
                'timestamp': run_dir.stat().st_mtime,
                'has_report': report_path.exists(),
                'report_path': str(report_path) if report_path.exists() else None
            })
    
    return jsonify(results)


@app.route('/api/results/report/<path:report_path>')
def view_report(report_path):
    """Serve HTML report"""
    full_path = Path(report_path)
    if full_path.exists() and full_path.suffix == '.html':
        return send_file(full_path)
    return jsonify({'error': 'Report not found'}), 404


@app.route('/api/results/download/<path:file_path>')
def download_file(file_path):
    """Download result file"""
    full_path = Path(file_path)
    if full_path.exists():
        return send_file(full_path, as_attachment=True)
    return jsonify({'error': 'File not found'}), 404


# ============================================================================
# SERVER-SENT EVENTS (SSE) - Live Progress Updates
# ============================================================================

@app.route('/api/stream/progress')
def stream_progress():
    """Server-sent events endpoint for real-time progress updates"""
    def generate():
        while True:
            # Send current pipeline state every second
            data = json.dumps(pipeline_state)
            yield f"data: {data}\n\n"
            time.sleep(1)
            
            # Stop streaming if pipeline is completed or in error state
            if pipeline_state['status'] in ['completed', 'error', 'stopped']:
                break
    
    return Response(generate(), mimetype='text/event-stream')


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def allowed_file(filename: str, extensions: List[str]) -> bool:
    """Check if file extension is allowed"""
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in extensions


def check_command_exists(command: str) -> bool:
    """Check if a command-line tool is available"""
    try:
        subprocess.run(
            [command, '--version'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=5
        )
        return True
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.CalledProcessError):
        return False


def get_disk_space() -> float:
    """Get available disk space in GB"""
    import shutil
    total, used, free = shutil.disk_usage('.')
    return round(free / (1024**3), 2)


def check_internet_connection() -> bool:
    """Check if internet connection is available"""
    try:
        import urllib.request
        urllib.request.urlopen('https://www.ncbi.nlm.nih.gov', timeout=3)
        return True
    except:
        return False


def reset_pipeline_state():
    """Reset pipeline state to initial values"""
    global pipeline_state
    pipeline_state['status'] = 'idle'
    pipeline_state['current_stage'] = None
    pipeline_state['current_genome'] = None
    pipeline_state['progress'] = 0
    pipeline_state['total_genomes'] = 0
    pipeline_state['start_time'] = None
    pipeline_state['end_time'] = None
    pipeline_state['error_message'] = None
    pipeline_state['log_messages'] = []
    pipeline_state['output_dir'] = None
    pipeline_state['report_path'] = None
    
    for stage in pipeline_state['stages']:
        stage['status'] = 'pending'
        stage['progress'] = 0


def add_log_message(level: str, message: str):
    """Add a log message to the pipeline state"""
    pipeline_state['log_messages'].append({
        'timestamp': datetime.now().isoformat(),
        'level': level,
        'message': message
    })


def update_stage_progress(stage_id: int, progress: int, status: str = 'running'):
    """Update progress for a specific pipeline stage"""
    for stage in pipeline_state['stages']:
        if stage['id'] == stage_id:
            stage['progress'] = progress
            stage['status'] = status
            pipeline_state['current_stage'] = stage['name']
            break


def run_pipeline_worker(config: Dict, output_dir: Path):
    """
    Background worker that executes the pipeline
    This is a placeholder - will be implemented in TODO 4
    """
    global pipeline_state
    
    try:
        add_log_message('info', f'Output directory: {output_dir}')
        
        # Simulate pipeline execution (will be replaced with actual pipeline calls)
        stages = [
            ('Genome Harvester', 'run_harvester.py'),
            ('Gene Annotator', 'run_annotator.py'),
            ('Sequence Extractor', 'run_extractor.py'),
            ('Wild-type Aligner', 'run_aligner.py'),
            ('Mutation Analyzer', 'run_analyzer.py'),
            ('Co-occurrence Analyzer', 'run_cooccurrence_analyzer.py'),
            ('Report Generator', 'run_reporter.py'),
        ]
        
        for idx, (stage_name, script) in enumerate(stages, 1):
            add_log_message('info', f'Starting {stage_name}...')
            update_stage_progress(idx, 0, 'running')
            
            # Simulate work (will be replaced with actual subprocess calls in TODO 4)
            for progress in range(0, 101, 20):
                time.sleep(0.5)  # Simulate processing time
                update_stage_progress(idx, progress, 'running')
            
            update_stage_progress(idx, 100, 'completed')
            add_log_message('success', f'{stage_name} completed')
        
        # Mark pipeline as completed
        pipeline_state['status'] = 'completed'
        pipeline_state['end_time'] = datetime.now().isoformat()
        pipeline_state['report_path'] = str(output_dir / 'final_report.html')
        add_log_message('success', 'Pipeline completed successfully!')
        
    except Exception as e:
        pipeline_state['status'] = 'error'
        pipeline_state['end_time'] = datetime.now().isoformat()
        pipeline_state['error_message'] = str(e)
        add_log_message('error', f'Pipeline failed: {str(e)}')


# ============================================================================
# MAIN
# ============================================================================

if __name__ == '__main__':
    app.run(debug=True, host='127.0.0.1', port=5000, threaded=True)
