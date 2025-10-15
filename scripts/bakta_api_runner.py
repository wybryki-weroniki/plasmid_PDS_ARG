#!/usr/bin/env python3
"""
Master script to run Bakta API on multiple FASTA files with extracted parameters.

This script:
1. Loads extracted parameters from extract_bakta_parameters.py
2. Submits jobs to Bakta API with optimized batching
3. Monitors job progress with error handling
4. Downloads completed results (JSON.gz format)
5. Provides logging and progress tracking

Usage:
    python bakta_api_runner.py [--dry-run] [--max-concurrent 5] [--batch-size 10]
"""

import os
import json
import time
import requests
import logging
import argparse
from pathlib import Path
from typing import Dict, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime


# Set up logging
def setup_logging(log_dir: Path):
    """Setup comprehensive logging."""
    log_dir.mkdir(exist_ok=True)

    # Create logger
    logger = logging.getLogger('bakta_api_runner')
    logger.setLevel(logging.INFO)

    # Create formatters
    detailed_formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(message)s'
    )
    simple_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # File handler for detailed logs
    file_handler = logging.FileHandler(log_dir / 'bakta_api_runner.log')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(detailed_formatter)

    # Console handler for user feedback
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(simple_formatter)

    # Add handlers
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger


@dataclass
class BaktaJob:
    """Class to track Bakta API job status."""
    file_id: str
    fasta_path: str
    job_id: Optional[str] = None
    secret: Optional[str] = None
    status: str = "pending"
    submit_time: Optional[datetime] = None
    complete_time: Optional[datetime] = None
    error_message: Optional[str] = None
    result_path: Optional[str] = None
    retry_count: int = 0


class BaktaAPIClient:
    """Client for interacting with Bakta API."""

    def __init__(self, base_url: str = "https://api.bakta.computational.bio", timeout: int = 30):
        self.base_url = base_url.rstrip('/')
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'BaktaAPIRunner/1.0'
        })
        self.logger = logging.getLogger('bakta_api_runner.client')

    def init_job(self, name: str, replicon_table_type: str = "CSV") -> Tuple[Optional[str], Optional[str], Optional[Dict]]:
        """
        Initialize a new Bakta job.

        Returns:
            Tuple of (job_id, secret, upload_links) or (None, None, None) on failure
        """
        try:
            response = self.session.post(
                f"{self.base_url}/api/v1/job/init",
                json={
                    "name": name,
                    "repliconTableType": replicon_table_type
                },
                timeout=self.timeout
            )

            if response.status_code == 200:
                data = response.json()
                # Extract job info from nested 'job' object
                job_info = data.get('job', {})
                job_id = job_info.get('jobID')  # Note: API uses 'jobID' not 'jobId'
                secret = job_info.get('secret')
                upload_links = {
                    'fasta': data.get('uploadLinkFasta'),
                    'prodigal': data.get('uploadLinkProdigal'),
                    'replicons': data.get('uploadLinkReplicons')
                }
                self.logger.info(f"Job initialized successfully: {job_id}")
                return job_id, secret, upload_links
            else:
                self.logger.error(f"Job initialization failed: {response.status_code} - {response.text}")
                return None, None, None

        except Exception as e:
            self.logger.error(f"Error initializing job: {str(e)}")
            return None, None, None

    def upload_file(self, upload_url: str, file_path: str) -> bool:
        """Upload file to pre-signed S3 URL."""
        try:
            # Check if file exists
            if not os.path.exists(file_path):
                self.logger.error(f"File does not exist: {file_path}")
                return False

            # Get file size for logging
            file_size = os.path.getsize(file_path)
            self.logger.debug(f"Uploading file: {file_path} ({file_size} bytes)")

            with open(file_path, 'rb') as f:
                response = self.session.put(upload_url, data=f, timeout=300)  # 5 min timeout for uploads

            self.logger.debug(f"Upload response status: {response.status_code}")
            if response.status_code == 200:
                self.logger.info(f"File uploaded successfully: {file_path} ({file_size} bytes)")
                return True
            else:
                self.logger.error(f"File upload failed: {response.status_code} - {file_path}")
                self.logger.error(f"Upload response: {response.text}")
                return False

        except Exception as e:
            self.logger.error(f"Error uploading file {file_path}: {str(e)}")
            return False

    def start_job(self, job_id: str, secret: str, genus: str = "", species: str = "", strain: str = "", complete_genome: bool = False, locus_tag: str = "", locus: str = "") -> bool:
        """Start the Bakta annotation job."""
        try:
            # Generate default locus_tag and locus if not provided
            if not locus_tag and genus and species:
                # Create locus tag from genus/species (e.g., "ECOLI" for Escherichia coli)
                locus_tag = f"{genus[:2].upper()}{species[:3].upper()}" if len(genus) >= 2 and len(species) >= 3 else "BAKTA"
            elif not locus_tag:
                locus_tag = "BAKTA"

            if not locus:
                locus = f"{locus_tag}_{job_id[:8]}"  # Use first 8 chars of job ID

            config = {
                "translationTable": 11,
                "completeGenome": complete_genome,
                "keepContigHeaders": False,
                "minContigLength": 1,  # Minimum contig length (1 bp = include all contigs)
                "compliant": True,  
                "genus": genus,
                "species": species,
                "strain": strain if strain else "",  # Ensure strain is never None
                "locusTag": locus_tag,         # Required: locus tag prefix
                "locus": locus,                # Required: locus identifier
                "hasReplicons": False,         # Required: replicon flag
                "dermType": None,              # Dermatophyte type (optional)
                "plasmid": None,               # Plasmid name (optional)
                "prodigalTrainingFile": None   # Training file (optional)
            }

            response = self.session.post(
                f"{self.base_url}/api/v1/job/start",
                json={
                    "job": {"jobID": job_id, "secret": secret},
                    "config": config
                },
                timeout=self.timeout
            )

            if response.status_code == 200:
                self.logger.info(f"Job started successfully: {job_id}")
                return True
            else:
                self.logger.error(f"Job start failed: {response.status_code} - {response.text}")
                return False

        except Exception as e:
            self.logger.error(f"Error starting job {job_id}: {str(e)}")
            return False

    def get_job_status(self, job_id: str, secret: str) -> Optional[Dict]:
        """Get job status and results."""
        try:
            # Bakta API uses /job/list endpoint with POST request
            response = self.session.post(
                f"{self.base_url}/api/v1/job/list",
                json={"jobs": [{"jobID": job_id, "secret": secret}]},
                timeout=self.timeout
            )

            if response.status_code == 200:
                data = response.json()
                # Response contains jobs array, get the first one
                jobs = data.get('jobs', [])
                if jobs and len(jobs) > 0:
                    job_status = jobs[0]
                    return {
                        'status': job_status.get('jobStatus', 'UNKNOWN'),
                        'result': job_status.get('result', {}),
                        'error': job_status.get('error', None)
                    }
                else:
                    self.logger.error(f"No job data returned for {job_id}")
                    return None
            else:
                self.logger.error(f"Status check failed: {response.status_code} - {response.text}")
                return None

        except Exception as e:
            self.logger.error(f"Error checking status for job {job_id}: {str(e)}")
            return None

    def get_job_logs(self, job_id: str, secret: str) -> Optional[str]:
        """Get job logs for debugging."""
        try:
            response = self.session.get(
                f"{self.base_url}/api/v1/job/{job_id}/logs",
                params={"secret": secret},
                timeout=self.timeout
            )

            if response.status_code == 200:
                return response.text
            else:
                self.logger.debug(f"Logs not available: {response.status_code}")
                return None

        except Exception as e:
            self.logger.debug(f"Error getting logs for job {job_id}: {str(e)}")
            return None

    def get_job_results(self, job_id: str, secret: str) -> Optional[Dict]:
        """Get result files URLs for a completed job using the correct endpoint."""
        try:
            # Use the result endpoint: POST /api/v1/job/result
            response = self.session.post(
                f"{self.base_url}/api/v1/job/result",
                json={"jobID": job_id, "secret": secret},
                timeout=self.timeout
            )

            if response.status_code == 200:
                result_data = response.json()
                self.logger.info(f"Retrieved result URLs for job {job_id}")
                return result_data
            else:
                self.logger.error(f"Failed to get results for job {job_id}: {response.status_code} - {response.text}")
                return None

        except Exception as e:
            self.logger.error(f"Error getting results for job {job_id}: {str(e)}")
            return None

    def download_result(self, download_url: str, output_path: str) -> bool:
        """Download result file."""
        try:
            response = self.session.get(download_url, timeout=600)  # 10 min timeout for downloads

            if response.status_code == 200:
                with open(output_path, 'wb') as f:
                    f.write(response.content)
                self.logger.info(f"Result downloaded: {output_path}")
                return True
            else:
                self.logger.error(f"Download failed: {response.status_code}")
                return False

        except Exception as e:
            self.logger.error(f"Error downloading result: {str(e)}")
            return False


class BaktaAPIRunner:
    """Main class for running Bakta API jobs."""

    def __init__(self, max_concurrent: int = 5, max_retries: int = 3):
        self.client = BaktaAPIClient()
        self.max_concurrent = max_concurrent
        self.max_retries = max_retries
        self.logger = logging.getLogger('bakta_api_runner.runner')
        self.jobs: Dict[str, BaktaJob] = {}
        self.results_dir = Path("bakta_results")
        self.results_dir.mkdir(exist_ok=True)

    def _generate_locus_tag(self, genus: str, species: str, file_id: str) -> str:
        """Generate a locus tag from genus, species, and file ID."""
        if genus and species:
            # Create locus tag from genus/species (e.g., "ECOLI" for Escherichia coli)
            if len(genus) >= 2 and len(species) >= 3:
                return f"{genus[:2].upper()}{species[:3].upper()}"
            else:
                return f"{genus.upper()}{species.upper()}"[:5]
        else:
            # Fall back to using file ID
            clean_file_id = ''.join(c.upper() for c in file_id if c.isalnum())[:5]
            return clean_file_id if clean_file_id else "BAKTA"

    def load_parameters(self, params_file: str) -> Dict:
        """Load parameters from JSON file created by extract_bakta_parameters.py"""
        try:
            with open(params_file, 'r') as f:
                params = json.load(f)
            self.logger.info(f"Loaded parameters for {len(params)} files")
            return params
        except Exception as e:
            self.logger.error(f"Error loading parameters: {str(e)}")
            return {}

    def create_replicon_table_for_file(self, file_id: str, params: Dict, replicon_df) -> Optional[str]:
        """Create a replicon table file for a specific FASTA file."""
        try:
            # Filter replicon table for this specific file
            sequence_id = params.get('sequence_id', file_id)
            file_replicons = replicon_df[replicon_df['locus'] == sequence_id]

            if file_replicons.empty:
                self.logger.warning(f"No replicon data found for {file_id}")
                return None

            # Create temporary replicon file
            replicon_file = self.results_dir / f"{file_id}_replicons.csv"
            file_replicons.to_csv(replicon_file, index=False)
            return str(replicon_file)

        except Exception as e:
            self.logger.error(f"Error creating replicon table for {file_id}: {str(e)}")
            return None

    def submit_job(self, file_id: str, params: Dict, replicon_df=None) -> BaktaJob:
        """Submit a single job to Bakta API."""
        job = BaktaJob(file_id=file_id, fasta_path=params['fasta_file'])

        try:
            # Initialize job
            job_id, secret, upload_links = self.client.init_job(f"bakta_job_{file_id}")

            if not job_id:
                job.status = "failed"
                job.error_message = "Job initialization failed"
                return job

            job.job_id = job_id
            job.secret = secret
            job.submit_time = datetime.now()

            # Upload FASTA file
            self.logger.info(f"Uploading FASTA file: {job.fasta_path}")
            self.logger.debug(f"Upload URL: {upload_links['fasta']}")

            if not self.client.upload_file(upload_links['fasta'], job.fasta_path):
                job.status = "failed"
                job.error_message = "FASTA upload failed"
                return job
            else:
                self.logger.info(f"FASTA file uploaded successfully: {job.fasta_path}")

            # Upload replicon table if available
            if replicon_df is not None:
                replicon_file = self.create_replicon_table_for_file(file_id, params, replicon_df)
                if replicon_file:
                    if not self.client.upload_file(upload_links['replicons'], replicon_file):
                        self.logger.warning(f"Replicon table upload failed for {file_id}, continuing without it")

            # Start job with taxonomic info
            genus = params.get('genus', '')
            species = params.get('species', '')
            strain = params.get('strain', '')
            complete_genome = params.get('circular', False)  # Use circular info as complete genome indicator

            # Generate locus tag and locus identifiers
            locus_tag = self._generate_locus_tag(genus, species, file_id)
            locus = f"{locus_tag}_{job_id[:8]}"

            if not self.client.start_job(job_id, secret, genus, species, strain, complete_genome, locus_tag, locus):
                job.status = "failed"
                job.error_message = "Job start failed"
                return job

            job.status = "running"
            self.logger.info(f"Job submitted successfully: {file_id} -> {job_id}")

        except Exception as e:
            job.status = "failed"
            job.error_message = str(e)
            self.logger.error(f"Error submitting job for {file_id}: {str(e)}")

        return job

    def monitor_job(self, job: BaktaJob) -> bool:
        """Monitor a single job until completion."""
        if not job.job_id or not job.secret:
            return False

        max_wait_time = 3600  # 1 hour timeout
        check_interval = 30   # Check every 30 seconds
        elapsed_time = 0

        while elapsed_time < max_wait_time:
            status_data = self.client.get_job_status(job.job_id, job.secret)

            if not status_data:
                self.logger.error(f"Could not get status for job {job.job_id}")
                return False

            current_status = status_data.get('status', 'unknown')
            self.logger.debug(f"Job {job.file_id} status: {current_status}")

            if current_status == 'SUCCESSFUL':
                job.status = "completed"
                job.complete_time = datetime.now()

                # DOWNLOAD RESULTS using the API endpoint
                self.logger.info(f"Job {job.file_id} completed successfully, retrieving results...")

                # Use the documented result endpoint: POST /api/v1/job/result
                result_data = self.client.get_job_results(job.job_id, job.secret)

                if result_data and 'ResultFiles' in result_data:
                    result_files = result_data['ResultFiles']
                    self.logger.info(f"Retrieved {len(result_files)} result files for job {job.file_id}")

                    # Priority order: JSON (main results), then JSON.gz if available
                    download_candidates = [
                        ('JSON', '.json'),
                        ('JSONGZ', '.json.gz'),  # In case they add compressed version later
                    ]

                    downloaded = False
                    for file_type, extension in download_candidates:
                        if file_type in result_files:
                            download_url = result_files[file_type]
                            output_file = self.results_dir / f"{job.file_id}_bakta_results{extension}"

                            self.logger.info(f"Downloading {file_type} results for {job.file_id}...")
                            if self.client.download_result(download_url, str(output_file)):
                                job.result_path = str(output_file)
                                job.status = "completed"
                                self.logger.info(f"Successfully downloaded {file_type} results: {job.file_id}")
                                downloaded = True
                                break
                            else:
                                self.logger.warning(f"Failed to download {file_type} results for {job.file_id}")

                    if downloaded:
                        return True
                    else:
                        job.status = "download_failed"
                        self.logger.error(f"Failed to download any result files for job {job.file_id}")
                        return False

                else:
                    job.status = "no_results"
                    self.logger.error(f"Could not retrieve result URLs for job {job.file_id}")
                    if result_data:
                        self.logger.debug(f"Result data received: {result_data}")
                    return False

            elif current_status == 'ERROR':
                job.status = "failed"
                job.error_message = status_data.get('error', 'Unknown error')
                self.logger.error(f"Job failed: {job.file_id} - {job.error_message}")
                return False

            elif current_status in ['PENDING', 'SUBMITTED', 'RUNNING', 'INIT', 'UPLOADING']:
                time.sleep(check_interval)
                elapsed_time += check_interval
            else:
                self.logger.warning(f"Unknown job status: {current_status}")
                time.sleep(check_interval)
                elapsed_time += check_interval

        # Timeout
        job.status = "timeout"
        self.logger.error(f"Job timed out: {job.file_id}")
        return False

    def run_batch(self, parameters: Dict, replicon_df=None, batch_size: int = 10, dry_run: bool = False):
        """Run Bakta API for all files in batches."""
        file_ids = list(parameters.keys())
        total_files = len(file_ids)

        self.logger.info(f"Starting Bakta API run for {total_files} files")
        self.logger.info(f"Batch size: {batch_size}, Max concurrent: {self.max_concurrent}")

        if dry_run:
            self.logger.info("DRY RUN MODE - No actual API calls will be made")
            for file_id in file_ids:
                params = parameters[file_id]
                self.logger.info(f"Would process: {file_id} -> {params['fasta_file']}")
            return

        # Process in batches
        for batch_start in range(0, total_files, batch_size):
            batch_end = min(batch_start + batch_size, total_files)
            batch_files = file_ids[batch_start:batch_end]

            self.logger.info(f"Processing batch {batch_start//batch_size + 1}: files {batch_start+1}-{batch_end}")

            # Submit jobs in parallel
            with ThreadPoolExecutor(max_workers=self.max_concurrent) as executor:
                # Submit jobs
                future_to_file = {}
                for file_id in batch_files:
                    params = parameters[file_id]
                    future = executor.submit(self.submit_job, file_id, params, replicon_df)
                    future_to_file[future] = file_id

                # Collect submitted jobs
                submitted_jobs = []
                for future in as_completed(future_to_file):
                    file_id = future_to_file[future]
                    try:
                        job = future.result()
                        self.jobs[file_id] = job
                        if job.status == "running":
                            submitted_jobs.append(job)
                        else:
                            self.logger.error(f"Job submission failed for {file_id}: {job.error_message}")
                    except Exception as e:
                        self.logger.error(f"Exception during job submission for {file_id}: {str(e)}")

                self.logger.info(f"Submitted {len(submitted_jobs)} jobs in batch")

                # Monitor jobs
                if submitted_jobs:
                    monitor_futures = {}
                    for job in submitted_jobs:
                        future = executor.submit(self.monitor_job, job)
                        monitor_futures[future] = job.file_id

                    # Wait for all jobs to complete
                    completed_count = 0
                    for future in as_completed(monitor_futures):
                        file_id = monitor_futures[future]
                        try:
                            success = future.result()
                            completed_count += 1
                            if success:
                                self.logger.info(f"Batch job completed successfully: {file_id} ({completed_count}/{len(submitted_jobs)})")
                            else:
                                self.logger.error(f"Batch job failed: {file_id}")
                        except Exception as e:
                            self.logger.error(f"Exception during job monitoring for {file_id}: {str(e)}")

            # Brief pause between batches
            if batch_end < total_files:
                self.logger.info("Pausing 10 seconds between batches...")
                time.sleep(10)

        self.generate_final_report()

    def generate_final_report(self):
        """Generate final summary report."""
        report_path = self.results_dir / "bakta_run_report.txt"

        successful = [j for j in self.jobs.values() if j.status == "completed"]
        failed = [j for j in self.jobs.values() if j.status in ["failed", "timeout", "download_failed"]]

        with open(report_path, 'w') as f:
            f.write("BAKTA API RUN SUMMARY REPORT\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Total files processed: {len(self.jobs)}\n")
            f.write(f"Successful: {len(successful)}\n")
            f.write(f"Failed: {len(failed)}\n\n")

            if successful:
                f.write("SUCCESSFUL JOBS:\n")
                f.write("-" * 20 + "\n")
                for job in successful:
                    duration = (job.complete_time - job.submit_time).total_seconds() / 60
                    f.write(f"{job.file_id}: {job.result_path} ({duration:.1f} minutes)\n")
                f.write("\n")

            if failed:
                f.write("FAILED JOBS:\n")
                f.write("-" * 20 + "\n")
                for job in failed:
                    f.write(f"{job.file_id}: {job.status} - {job.error_message}\n")

        self.logger.info(f"Final report generated: {report_path}")
        self.logger.info(f"SUMMARY: {len(successful)} successful, {len(failed)} failed out of {len(self.jobs)} total")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Run Bakta API on multiple FASTA files")
    parser.add_argument("--dry-run", action="store_true", help="Perform dry run without actual API calls")
    parser.add_argument("--max-concurrent", type=int, default=5, help="Maximum concurrent jobs (default: 5)")
    parser.add_argument("--batch-size", type=int, default=10, help="Batch size for processing (default: 10)")
    parser.add_argument("--params-file", default="bakta_params/bakta_parameters.json", help="Parameters JSON file")
    parser.add_argument("--replicon-file", default="bakta_params/replicons.csv", help="Replicon table CSV file")

    args = parser.parse_args()

    # Setup logging
    log_dir = Path("logs")
    logger = setup_logging(log_dir)

    # Check if parameter files exist
    script_dir = Path(__file__).parent
    params_file = script_dir / args.params_file
    replicon_file = script_dir / args.replicon_file

    if not params_file.exists():
        logger.error(f"Parameters file not found: {params_file}")
        logger.error("Please run extract_bakta_parameters.py first!")
        return

    # Initialize runner
    runner = BaktaAPIRunner(max_concurrent=args.max_concurrent)

    # Load parameters
    parameters = runner.load_parameters(str(params_file))
    if not parameters:
        logger.error("No parameters loaded, exiting")
        return

    # Load replicon table if available
    replicon_df = None
    if replicon_file.exists():
        try:
            import pandas as pd
            replicon_df = pd.read_csv(replicon_file)
            logger.info(f"Loaded replicon table with {len(replicon_df)} entries")
        except Exception as e:
            logger.warning(f"Could not load replicon table: {str(e)}")

    # Run the batch
    runner.run_batch(
        parameters=parameters,
        replicon_df=replicon_df,
        batch_size=args.batch_size,
        dry_run=args.dry_run
    )


if __name__ == "__main__":
    main()