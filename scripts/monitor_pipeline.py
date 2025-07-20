#!/usr/bin/env python3
"""
Pipeline Performance Monitor
============================

This script monitors Snakemake pipeline execution and collects performance metrics
including CPU usage, memory consumption, and rule execution times.
"""

import psutil
import time
import json
import subprocess
import threading
import logging
from pathlib import Path
from datetime import datetime
import sys

class PipelineMonitor:
    def __init__(self, output_file="pipeline_performance.json"):
        self.output_file = output_file
        self.start_time = None
        self.end_time = None
        self.monitoring = False
        self.metrics = {
            "start_time": None,
            "end_time": None,
            "total_duration_seconds": None,
            "cpu_usage": [],
            "memory_usage": [],
            "disk_usage": [],
            "snakemake_process": None,
            "peak_memory_mb": 0,
            "average_cpu_percent": 0,
            "rule_times": {},
            "system_info": {}
        }
        
        # Set up logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler('pipeline_monitor.log'),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)

    def collect_system_info(self):
        """Collect basic system information"""
        self.metrics["system_info"] = {
            "cpu_count": psutil.cpu_count(),
            "cpu_count_logical": psutil.cpu_count(logical=True),
            "total_memory_gb": round(psutil.virtual_memory().total / (1024**3), 2),
            "available_memory_gb": round(psutil.virtual_memory().available / (1024**3), 2),
            "disk_total_gb": round(psutil.disk_usage('.').total / (1024**3), 2),
            "disk_free_gb": round(psutil.disk_usage('.').free / (1024**3), 2),
            "python_version": sys.version,
            "platform": psutil.WINDOWS if psutil.WINDOWS else "Linux/Unix"
        }

    def monitor_resources(self):
        """Monitor system resources while pipeline runs"""
        self.logger.info("Starting resource monitoring...")
        
        while self.monitoring:
            try:
                # CPU usage
                cpu_percent = psutil.cpu_percent(interval=1)
                
                # Memory usage
                memory = psutil.virtual_memory()
                memory_mb = memory.used / (1024**2)
                
                # Disk usage
                disk = psutil.disk_usage('.')
                disk_usage_percent = (disk.used / disk.total) * 100
                
                # Update peak memory
                if memory_mb > self.metrics["peak_memory_mb"]:
                    self.metrics["peak_memory_mb"] = memory_mb
                
                # Store metrics
                timestamp = time.time()
                self.metrics["cpu_usage"].append({
                    "timestamp": timestamp,
                    "cpu_percent": cpu_percent
                })
                
                self.metrics["memory_usage"].append({
                    "timestamp": timestamp,
                    "memory_mb": memory_mb,
                    "memory_percent": memory.percent
                })
                
                self.metrics["disk_usage"].append({
                    "timestamp": timestamp,
                    "disk_percent": disk_usage_percent,
                    "disk_free_gb": disk.free / (1024**3)
                })
                
                # Log every 30 seconds
                if len(self.metrics["cpu_usage"]) % 30 == 0:
                    self.logger.info(f"CPU: {cpu_percent:.1f}%, Memory: {memory_mb:.1f}MB ({memory.percent:.1f}%)")
                
            except Exception as e:
                self.logger.warning(f"Error collecting metrics: {e}")
            
            time.sleep(1)

    def parse_snakemake_log(self, log_file="snakemake.log"):
        """Parse Snakemake log to extract rule execution times"""
        if not Path(log_file).exists():
            self.logger.warning(f"Snakemake log file not found: {log_file}")
            return
        
        try:
            with open(log_file, 'r') as f:
                lines = f.readlines()
            
            current_rule = None
            rule_start_time = None
            
            for line in lines:
                # Look for rule start patterns
                if "rule " in line and ":" in line:
                    rule_match = line.split("rule ")[1].split(":")[0].strip()
                    if rule_match:
                        current_rule = rule_match
                        # Try to extract timestamp (this is approximate)
                        rule_start_time = time.time()
                
                # Look for rule completion patterns
                if current_rule and ("Finished job" in line or "completed" in line):
                    if rule_start_time:
                        duration = time.time() - rule_start_time
                        if current_rule not in self.metrics["rule_times"]:
                            self.metrics["rule_times"][current_rule] = []
                        self.metrics["rule_times"][current_rule].append(duration)
                    current_rule = None
                    rule_start_time = None
                    
        except Exception as e:
            self.logger.warning(f"Error parsing Snakemake log: {e}")

    def start_monitoring(self):
        """Start monitoring the pipeline"""
        self.start_time = time.time()
        self.metrics["start_time"] = datetime.now().isoformat()
        self.monitoring = True
        
        self.logger.info("Pipeline monitoring started")
        self.collect_system_info()
        
        # Start resource monitoring in a separate thread
        self.monitor_thread = threading.Thread(target=self.monitor_resources)
        self.monitor_thread.daemon = True
        self.monitor_thread.start()

    def stop_monitoring(self):
        """Stop monitoring and generate report"""
        self.monitoring = False
        self.end_time = time.time()
        self.metrics["end_time"] = datetime.now().isoformat()
        self.metrics["total_duration_seconds"] = self.end_time - self.start_time
        
        # Calculate averages
        if self.metrics["cpu_usage"]:
            self.metrics["average_cpu_percent"] = sum(
                m["cpu_percent"] for m in self.metrics["cpu_usage"]
            ) / len(self.metrics["cpu_usage"])
        
        # Parse Snakemake log if available
        self.parse_snakemake_log()
        
        # Save metrics to file
        with open(self.output_file, 'w') as f:
            json.dump(self.metrics, f, indent=2)
        
        self.logger.info(f"Monitoring stopped. Report saved to {self.output_file}")
        self.print_summary()

    def print_summary(self):
        """Print a summary of the performance metrics"""
        duration_minutes = self.metrics["total_duration_seconds"] / 60
        
        print("\n" + "="*50)
        print("PIPELINE PERFORMANCE SUMMARY")
        print("="*50)
        print(f"Total Duration: {duration_minutes:.1f} minutes ({self.metrics['total_duration_seconds']:.1f} seconds)")
        print(f"Average CPU Usage: {self.metrics['average_cpu_percent']:.1f}%")
        print(f"Peak Memory Usage: {self.metrics['peak_memory_mb']:.1f} MB")
        print(f"System CPU Cores: {self.metrics['system_info']['cpu_count']} physical, {self.metrics['system_info']['cpu_count_logical']} logical")
        print(f"System Memory: {self.metrics['system_info']['total_memory_gb']:.1f} GB total")
        
        if self.metrics["rule_times"]:
            print(f"\nRule Execution Times:")
            for rule, times in self.metrics["rule_times"].items():
                avg_time = sum(times) / len(times)
                print(f"  {rule}: {avg_time:.1f}s average ({len(times)} executions)")
        
        print(f"\nDetailed metrics saved to: {self.output_file}")
        print("="*50)


def run_snakemake_with_monitoring(snakemake_args=None, output_file="pipeline_performance.json"):
    """Run Snakemake with performance monitoring"""
    
    if snakemake_args is None:
        snakemake_args = ["snakemake", "--cores", "8", "all_final_reports"]
    
    monitor = PipelineMonitor(output_file)
    
    try:
        # Start monitoring
        monitor.start_monitoring()
        
        # Run Snakemake
        monitor.logger.info(f"Starting Snakemake with args: {' '.join(snakemake_args)}")
        
        # Capture Snakemake output to log file
        with open("snakemake.log", "w") as log_file:
            process = subprocess.Popen(
                snakemake_args,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                text=True
            )
            
            # Wait for completion
            return_code = process.wait()
            
            if return_code == 0:
                monitor.logger.info("Snakemake completed successfully")
            else:
                monitor.logger.error(f"Snakemake failed with return code {return_code}")
    
    except KeyboardInterrupt:
        monitor.logger.info("Pipeline interrupted by user")
        if 'process' in locals():
            process.terminate()
    
    except Exception as e:
        monitor.logger.error(f"Error running pipeline: {e}")
    
    finally:
        # Stop monitoring
        monitor.stop_monitoring()
    
    return monitor


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Monitor Snakemake pipeline performance")
    parser.add_argument("--cores", type=int, default=8, help="Number of cores for Snakemake")
    parser.add_argument("--target", default="all_final_reports", help="Snakemake target")
    parser.add_argument("--output", default="pipeline_performance.json", help="Output file for metrics")
    parser.add_argument("--dry-run", action="store_true", help="Perform dry run")
    
    args = parser.parse_args()
    
    snakemake_cmd = ["snakemake", "--cores", str(args.cores)]
    
    if args.dry_run:
        snakemake_cmd.append("--dry-run")
    
    snakemake_cmd.append(args.target)
    
    print(f"Running: {' '.join(snakemake_cmd)}")
    
    monitor = run_snakemake_with_monitoring(snakemake_cmd, args.output)