"""
3-1_Alignment.py: STAR two-pass alignment를 자동화하는 스크립트.

Usage:
  3-1_Alignment.py [--work_dir=<dir>] [--in_dir=<dir>] [--out_dir=<dir>] [--genome_dir=<dir>] [--threads=<n>] --picard_path=<dir> [--sort_out_dir=<dir>]

Options:
  --work_dir=<dir>       작업 디렉토리.
  --in_dir=<dir>         입력 FASTQ 파일들이 위치한 디렉토리. (trimming 결과나 raw fastq)
  --out_dir=<dir>        결과가 저장될 출력 디렉토리.
  --genome_dir=<dir>     STAR genome index가 위치한 디렉토리.
  --threads=<n>          사용할 스레드 수. [default: 20]
  --picard_path=<dir>    picard 저장된 pathway
  --sort_out_dir=<>      sort out한 결과가 저장될 출력디렉토리 aligned_queryname_sorted.bam


  각 샘플의 파일명은 <sample>_1.fastq.gz와 <sample>_2.fastq.gz 형태여야 하며,
  trimming 결과물이 존재하면 <sample>_1_paired.fastq.gz, <sample>_2_paired.fastq.gz 파일을 사용합니다.
"""


import os
import sys
import glob
import shutil
import logging
import subprocess
from pytz import timezone
from docopt import docopt
from datetime import datetime

from MakeGenomeref import MakeGenomeref as mg
from TrimmingFastq import TrimmingFastq as tf



def make_logger(log_dir, name="app", consoleset=True, streamset=True) -> logging.Logger:
    """
    Logger 생성 함수.
    
    - 콘솔 핸들러와 파일 핸들러를 추가하여, 로그 메시지를 콘솔 및 파일에 기록합니다.
    - 로그 시간은 Asia/Seoul 시간대를 사용합니다.
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    loggerformat = "%(asctime)s - %(name)s - %(module)s - %(levelname)s - %(message)s"
    formatter = logging.Formatter(loggerformat)
        
    #  change time zone 시간대 설정
    tz = timezone("Asia/Seoul")
    def timetz(*args):
        return datetime.now(tz).timetuple()
    logging.Formatter.converter = timetz

    if consoleset:
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        console.setFormatter(formatter)
        logger.addHandler(console)
    if streamset:
        loggerfile = os.path.join(log_dir, name)
        file_handler = logging.FileHandler(filename=f"{loggerfile}.log")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    return logger




class AlignFastq:
    @staticmethod
    def check_genome_dir(genome_dir: str, work_dir: str, logger) -> str:
        """
        STAR genome index 디렉토리가 존재하는지 확인합니다.
        genome_dir이 지정되지 않으면, work_dir/genome_ref를 기본값으로 사용합니다.
        
        Returns:
            실제 사용된 genome_dir 경로
        """
        if not genome_dir:
            genome_dir = os.path.join(work_dir, "genome_ref")
            logger.info(f"genome_dir이 지정되지 않아 기본값으로 {genome_dir}를 사용합니다.")
            
        if not os.path.isdir(genome_dir):
            logger.error(f"STAR genome index 디렉토리가 존재하지 않습니다: {genome_dir}")
            sys.exit(1)
            
        return genome_dir

    @staticmethod
    def make_dir(dir_path: str) -> None:
        os.makedirs(dir_path, exist_ok=True)



    @staticmethod
    def run_alignment(work_dir: str, in_dir: str, out_dir: str, genome_dir: str, threads: int,logger) -> None:
        """
        STAR Alignment 실행.
        
        *주의*
         - trimming을 거쳤다면, in_dir은 trimming 결과물이 저장된 디렉토리여야 합니다.
           (예: work_dir/trimmed, 여기에는 <sample>_1_paired.fastq.gz 파일이 있음)
         - trimming을 거치지 않았다면, in_dir은 원본 FASTQ 파일이 있는 디렉토리여야 합니다.
           (예: work_dir/fastq, 여기에는 <sample>_1.fastq.gz 파일이 있음)
         
        이 함수는 in_dir 내에서 _1_paired.fastq.gz 파일이 있으면 trimming 결과물을 사용하고,
        없으면 원본 FASTQ 파일을 사용합니다.
        """
        # STAR 설치 확인 (MakeGenomeref의 함수를 사용)
        mg.check_star_installed(logger)

        # genome_dir 존재 여부 확인
        genome_dir = AlignFastq.check_genome_dir(genome_dir, work_dir, logger)

        # 작업 디렉토리 변경
        try:
            os.chdir(work_dir)
        except Exception as exc:
            sys.stderr.write(f"작업 디렉토리 변경 에러: {exc}\n")
            sys.exit(1)


        # out_dir이 지정되어 있으면 사용, 아니면 기본값 지정
        if out_dir:
            final_out_dir = out_dir
            logger.info(f"사용자 지정 out_dir: {final_out_dir}")
        else:
            final_out_dir = os.path.join(work_dir, "star_alignment")
            logger.info(f"사용자가 out_dir을 지정하지 않았으므로 기본값을 사용합니다: {final_out_dir}")

        # final_out_dir이 없으면 생성
        if not os.path.isdir(final_out_dir):
            os.makedirs(final_out_dir, exist_ok=True)
            logger.info(f"{final_out_dir} 이(가) 존재하지 않아 새로 생성되었습니다.")
        else:
            logger.info(f"{final_out_dir} 이(가) 이미 존재합니다.")
        
        # 작업 디렉토리를 final_out_dir로 변경
        os.chdir(final_out_dir)


        
        sample_files = [f for f in os.listdir(in_dir) if f.endswith("_1.fastq.gz") or f.endswith("_1_paired.fastq.gz")]
        if not sample_files:
            logger.error("FASTQ 파일이 존재하지 않습니다.")
            sys.exit(1)
        
        samples = []
        for f in sample_files:
            if f.endswith("_1_paired.fastq.gz"):
                sample = f.replace("_1_paired.fastq.gz", "")
            else:
                sample = f.replace("_1.fastq.gz", "")
            if sample not in samples:
                samples.append(sample)
        
        for sample in samples:
            paired_f1 = os.path.join(in_dir, f"{sample}_1_paired.fastq.gz")
            paired_f2 = os.path.join(in_dir, f"{sample}_2_paired.fastq.gz")
            if os.path.isfile(paired_f1) and os.path.isfile(paired_f2):
                f1 = paired_f1
                f2 = paired_f2
            else:
                f1 = os.path.join(in_dir, f"{sample}_1.fastq.gz")
                f2 = os.path.join(in_dir, f"{sample}_2.fastq.gz")
            
            sample_out_dir = os.path.join(final_out_dir, sample)
            AlignFastq.make_dir(sample_out_dir)
            command = (
                f"STAR --runThreadN {threads} "
                f"--runMode alignReads "
                f"--genomeDir {genome_dir} "
                f"--readFilesIn {f1} {f2} "
                f"--readFilesCommand zcat "
                f"--twopassMode Basic "
                f"--outFileNamePrefix {sample_out_dir}/ "
                f"--sjdbOverhang 100 "
                f"--outSAMtype BAM SortedByCoordinate"
            )
            logger.info(f"STAR alignment for sample {sample}...")
            logger.info(command)
            try:
                subprocess.run(command, shell=True, check=True)
                logger.info(f"{sample} aligned.")
            except subprocess.CalledProcessError as error:
                logger.error(f"STAR 에러 ({sample}): {error}")
                continue
            
            bam_file = os.path.join(sample_out_dir, "Aligned.sortedByCoord.out.bam")
            if os.path.isfile(bam_file):
                index_command = f"samtools index {bam_file}"
                try:
                    subprocess.run(index_command, shell=True, check=True)
                    logger.info(f"{sample} BAM index created.")
                except subprocess.CalledProcessError as error:
                    logger.error(f"samtools index 에러 ({sample}): {error}")
            else:
                logger.error(f"{sample} BAM 파일 생성 실패: {bam_file}")
        logger.info("STAR alignment step completed.")
        return samples


class SortBam:
    @staticmethod
    def check_picard(picard_path: str, logger: logging.Logger) -> None:
        if not os.path.exists(picard_path):
            logger.error(f"Picard 경로가 존재하지 않습니다: {picard_path}")
            sys.exit(1)

    @staticmethod
    def run_sort(alignment_dir: str, sort_out_dir: str, picard_path: str, samples: list, logger: logging.Logger) -> None:
        """
        Picard SortSam을 사용하여 alignment 결과 BAM 파일을 query name 순으로 정렬하여 sort_out_dir에 저장합니다.
        
        alignment_dir: STAR alignment 결과가 저장된 상위 디렉토리 (각 샘플 폴더 내에 Aligned.sortedByCoord.out.bam 존재)
        sort_out_dir: 정렬된 BAM 파일이 저장될 디렉토리.
        """
        AlignFastq.make_dir(sort_out_dir)
        for sample in samples:
            input_bam = os.path.join(alignment_dir, sample, "Aligned.sortedByCoord.out.bam")
            output_bam = os.path.join(sort_out_dir, f"{sample}_aligned_queryname_sorted.bam")
            if not os.path.isfile(input_bam):
                logger.error(f"{sample}의 입력 BAM 파일이 존재하지 않습니다: {input_bam}")
                continue
            command = (
                f"java -jar {picard_path} SortSam "
                f"I={input_bam} "
                f"O={output_bam} "
                f"SORT_ORDER=queryname"
            )
            logger.info(f"{sample}의 BAM query name 정렬 시작...")
            try:
                subprocess.run(command, shell=True, check=True)
                logger.info(f"{sample} 정렬 완료: {output_bam}")
            except subprocess.CalledProcessError as error:
                logger.error(f"{sample} 정렬 에러: {error}")
                continue
        logger.info("BAM query name 정렬 단계 완료.")

def main() -> None:
    arguments = docopt(__doc__)

    # work_dir가 제공되지 않으면 현재 디렉토리를 기본값으로 사용합니다.
    work_dir = arguments.get("--work_dir")
    if not work_dir:
        work_dir = os.getcwd()

    # --in_dir 옵션이 지정되지 않으면 원본 FASTQ 디렉토리 (예: work_dir/trimmed) 사용
    in_dir = arguments.get("--in_dir")
    if not in_dir:
        trimmed_dir = os.path.join(work_dir, "trimmed")
        if os.path.isdir(trimmed_dir):
            in_dir = trimmed_dir
        else:
            in_dir = os.path.join(work_dir, "fastq")

    out_dir = arguments["--out_dir"]
    if not out_dir:
        out_dir = os.path.join(work_dir, "star_alignment")
    genome_dir = arguments["--genome_dir"]
    threads = int(arguments["--threads"])
    picard_path = arguments["--picard_path"]

    sort_out_dir = arguments.get("--sort_out_dir")
    if not sort_out_dir:
        sort_out_dir = os.path.join(work_dir, "star_alignment_sorted_bam")
    
    # 로거 생성 (work_dir에 로그 파일 생성)
    logger = make_logger(work_dir, name="pipeline_logger")

    # TrimmingFastq 모듈의 디렉토리 체크 함수에 인자 전달
    tf.check_in_dir(in_dir)
    tf.check_workdir(work_dir)

    # STAR Alignment 실행 후, 처리된 샘플 리스트 반환
    samples = AlignFastq.run_alignment(work_dir, in_dir, out_dir, genome_dir, threads, logger)
    
    # Picard SortSam을 이용한 BAM 정렬 단계 실행
    SortBam.check_picard(picard_path, logger)
    SortBam.run_sort(out_dir, sort_out_dir, picard_path, samples, logger)
    

if __name__ == "__main__":
    main()
