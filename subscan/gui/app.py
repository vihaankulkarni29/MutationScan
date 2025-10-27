#!/usr/bin/env python3
"""
MutationScan Web GUI - Flask Backend
Provides a local web interface for running the MutationScan pipeline.
"""

import os
import sys
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
    """Check system health and dependencies using health_check.py tool"""
    try:
        # Call health_check.py with --json flag
        health_check_path = Path(__file__).parent.parent / 'tools' / 'health_check.py'
        
        result = subprocess.run(
            [sys.executable, str(health_check_path), '--json'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=15,
            text=True
        )
        
        if result.returncode in [0, 1]:  # 0 = good/warning, 1 = error but still valid JSON
            health_status = json.loads(result.stdout)
            health_status['timestamp'] = datetime.now().isoformat()
            return jsonify(health_status)
        else:
            raise Exception(f"Health check failed: {result.stderr}")
            
    except Exception as e:
        # Fallback to basic checks if health_check.py fails
        return jsonify({
            'python': {'status': 'good', 'installed': True, 'version': f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"},
            'emboss': {'status': 'good' if check_command_exists('needle') else 'error', 'installed': check_command_exists('needle')},
            'abricate': {'status': 'good' if check_command_exists('abricate') else 'error', 'installed': check_command_exists('abricate')},
            'disk': {'free_gb': get_disk_space(), 'status': 'good' if get_disk_space() >= 10 else 'warning'},
            'internet': {'connected': check_internet_connection(), 'status': 'good' if check_internet_connection() else 'warning'},
            'overall_status': 'warning',
            'error': str(e),
            'timestamp': datetime.now().isoformat()
        })


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


def run_subprocess_with_progress(cmd: List[str], stage_id: int, stage_name: str) -> tuple:
    """
    Run a subprocess command and capture output
    Returns: (success: bool, stdout: str, stderr: str)
    """
    try:
        # Run the command and capture output
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            universal_newlines=True
        )
        
        stdout_lines = []
        stderr_lines = []
        
        # Track progress (incremental updates)
        progress = 10
        progress_increment = 15  # Increment by 15% for each output line (max 90%)
        
        # Read stdout in real-time
        for line in process.stdout:
            line = line.strip()
            if line:
                stdout_lines.append(line)
                add_log_message('info', f'[{stage_name}] {line}')
                
                # Update progress incrementally
                progress = min(progress + progress_increment, 90)
                update_stage_progress(stage_id, progress, 'running')
        
        # Wait for completion
        process.wait()
        
        # Read any remaining stderr
        stderr_output = process.stderr.read()
        if stderr_output:
            stderr_lines = stderr_output.strip().split('\n')
            for line in stderr_lines:
                if line.strip():
                    add_log_message('warning', f'[{stage_name}] {line.strip()}')
        
        stdout = '\n'.join(stdout_lines)
        stderr = '\n'.join(stderr_lines)
        
        success = process.returncode == 0
        
        if not success:
            add_log_message('error', f'{stage_name} failed with exit code {process.returncode}')
        
        return success, stdout, stderr
        
    except Exception as e:
        add_log_message('error', f'Failed to run {stage_name}: {str(e)}')
        return False, '', str(e)


def run_pipeline_worker(config: Dict, output_dir: Path):
    """
    Background worker that executes the pipeline
    Calls actual domino tools as subprocesses
    """
    global pipeline_state
    
    try:
        add_log_message('info', f'Output directory: {output_dir}')
        add_log_message('info', 'Preparing pipeline execution...')
        
        # Get tools directory
        tools_dir = Path(__file__).parent.parent / 'tools'
        
        # Prepare gene list (from file, template, or manual entry)
        gene_list_path = None
        if config.get('genes_file'):
            gene_list_path = config['genes_file']
        elif config.get('genes_list'):
            # Create temporary gene list file
            gene_list_path = output_dir / 'genes.txt'
            with open(gene_list_path, 'w') as f:
                for gene in config['genes_list']:
                    f.write(f"{gene.strip()}\n")
            add_log_message('info', f'Created gene list with {len(config["genes_list"])} genes')
        
        # Stage 1: Genome Harvester
        stage_1_dir = output_dir / '01_harvester_results'
        stage_1_dir.mkdir(exist_ok=True)
        
        add_log_message('info', 'Stage 1/7: Downloading genomes...')
        update_stage_progress(1, 10, 'running')
        
        harvester_cmd = [
            sys.executable,
            str(tools_dir / 'run_harvester.py'),
            '--accessions', config['accessions_file'],
            '--email', config.get('email', 'user@example.com'),
            '--output-dir', str(stage_1_dir),
            '--database', config.get('database', 'ncbi')
        ]
        
        success, stdout, stderr = run_subprocess_with_progress(harvester_cmd, 1, 'Genome Harvester')
        if not success:
            raise Exception(f"Harvester failed: {stderr}")
        
        update_stage_progress(1, 100, 'completed')
        add_log_message('success', 'Genome download completed')
        
        genome_manifest = stage_1_dir / 'genome_manifest.json'
        if not genome_manifest.exists():
            raise Exception('Genome manifest not created')
        
        # Stage 2: Gene Annotator
        stage_2_dir = output_dir / '02_annotator_results'
        stage_2_dir.mkdir(exist_ok=True)
        
        add_log_message('info', 'Stage 2/7: Identifying AMR genes...')
        update_stage_progress(2, 10, 'running')
        
        annotator_cmd = [
            sys.executable,
            str(tools_dir / 'run_annotator.py'),
            '--manifest', str(genome_manifest),
            '--output-dir', str(stage_2_dir),
            '--threads', str(config.get('threads', 'auto'))
        ]
        
        success, stdout, stderr = run_subprocess_with_progress(annotator_cmd, 2, 'Gene Annotator')
        if not success:
            raise Exception(f"Annotator failed: {stderr}")
        
        update_stage_progress(2, 100, 'completed')
        add_log_message('success', 'AMR gene annotation completed')
        
        annotation_manifest = stage_2_dir / 'annotation_manifest.json'
        if not annotation_manifest.exists():
            raise Exception('Annotation manifest not created')
        
        # Stage 3: Sequence Extractor
        stage_3_dir = output_dir / '03_extractor_results'
        stage_3_dir.mkdir(exist_ok=True)
        
        add_log_message('info', 'Stage 3/7: Extracting protein sequences...')
        update_stage_progress(3, 10, 'running')
        
        extractor_cmd = [
            sys.executable,
            str(tools_dir / 'run_extractor.py'),
            '--manifest', str(annotation_manifest),
            '--gene-list', str(gene_list_path),
            '--output-dir', str(stage_3_dir),
            '--threads', str(config.get('threads', 'auto'))
        ]
        
        success, stdout, stderr = run_subprocess_with_progress(extractor_cmd, 3, 'Sequence Extractor')
        if not success:
            raise Exception(f"Extractor failed: {stderr}")
        
        update_stage_progress(3, 100, 'completed')
        add_log_message('success', 'Sequence extraction completed')
        
        protein_manifest = stage_3_dir / 'protein_manifest.json'
        if not protein_manifest.exists():
            raise Exception('Protein manifest not created')
        
        # Stage 4: Wild-type Aligner
        stage_4_dir = output_dir / '04_aligner_results'
        stage_4_dir.mkdir(exist_ok=True)
        
        add_log_message('info', 'Stage 4/7: Aligning to reference sequences...')
        update_stage_progress(4, 10, 'running')
        
        aligner_cmd = [
            sys.executable,
            str(tools_dir / 'run_aligner.py'),
            '--manifest', str(protein_manifest),
            '--output-dir', str(stage_4_dir),
            '--threads', str(config.get('threads', 'auto')),
            '--sepi-species', config.get('species', 'Escherichia coli')
        ]
        
        success, stdout, stderr = run_subprocess_with_progress(aligner_cmd, 4, 'Wild-type Aligner')
        if not success:
            raise Exception(f"Aligner failed: {stderr}")
        
        update_stage_progress(4, 100, 'completed')
        add_log_message('success', 'Alignment completed')
        
        alignment_manifest = stage_4_dir / 'alignment_manifest.json'
        if not alignment_manifest.exists():
            raise Exception('Alignment manifest not created')
        
        # Stage 5: Mutation Analyzer
        stage_5_dir = output_dir / '05_analyzer_results'
        stage_5_dir.mkdir(exist_ok=True)
        
        add_log_message('info', 'Stage 5/7: Analyzing mutations...')
        update_stage_progress(5, 10, 'running')
        
        analyzer_cmd = [
            sys.executable,
            str(tools_dir / 'run_analyzer.py'),
            '--manifest', str(alignment_manifest),
            '--output-dir', str(stage_5_dir),
            '--threads', str(config.get('threads', 'auto'))
        ]
        
        success, stdout, stderr = run_subprocess_with_progress(analyzer_cmd, 5, 'Mutation Analyzer')
        if not success:
            raise Exception(f"Analyzer failed: {stderr}")
        
        update_stage_progress(5, 100, 'completed')
        add_log_message('success', 'Mutation analysis completed')
        
        analysis_manifest = stage_5_dir / 'analysis_manifest.json'
        if not analysis_manifest.exists():
            raise Exception('Analysis manifest not created')
        
        # Stage 6: Co-occurrence Analyzer
        stage_6_dir = output_dir / '06_cooccurrence_results'
        stage_6_dir.mkdir(exist_ok=True)
        
        add_log_message('info', 'Stage 6/7: Finding mutation patterns...')
        update_stage_progress(6, 10, 'running')
        
        cooccurrence_cmd = [
            sys.executable,
            str(tools_dir / 'run_cooccurrence_analyzer.py'),
            '--manifest', str(analysis_manifest),
            '--output-dir', str(stage_6_dir)
        ]
        
        success, stdout, stderr = run_subprocess_with_progress(cooccurrence_cmd, 6, 'Co-occurrence Analyzer')
        if not success:
            raise Exception(f"Co-occurrence analyzer failed: {stderr}")
        
        update_stage_progress(6, 100, 'completed')
        add_log_message('success', 'Co-occurrence analysis completed')
        
        cooccurrence_manifest = stage_6_dir / 'cooccurrence_manifest.json'
        if not cooccurrence_manifest.exists():
            raise Exception('Co-occurrence manifest not created')
        
        # Stage 7: Report Generator
        stage_7_dir = output_dir / '07_reporter_results'
        stage_7_dir.mkdir(exist_ok=True)
        
        add_log_message('info', 'Stage 7/7: Generating final report...')
        update_stage_progress(7, 10, 'running')
        
        reporter_cmd = [
            sys.executable,
            str(tools_dir / 'run_reporter.py'),
            '--manifest', str(cooccurrence_manifest),
            '--output-dir', str(stage_7_dir)
        ]
        
        success, stdout, stderr = run_subprocess_with_progress(reporter_cmd, 7, 'Report Generator')
        if not success:
            raise Exception(f"Reporter failed: {stderr}")
        
        update_stage_progress(7, 100, 'completed')
        add_log_message('success', 'Report generation completed')
        
        # Check for final report
        final_report = stage_7_dir / 'final_report.html'
        if not final_report.exists():
            # Try alternative name
            final_report = stage_7_dir / 'mutation_analysis_report.html'
        
        # Mark pipeline as completed
        pipeline_state['status'] = 'completed'
        pipeline_state['end_time'] = datetime.now().isoformat()
        pipeline_state['report_path'] = str(final_report) if final_report.exists() else None
        add_log_message('success', 'Pipeline completed successfully!')
        add_log_message('success', f'Report available at: {final_report}')
        
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
