"""
MergeBamUnbam.py: Star로 alignment한 bam file 과 unaligned 된 bam을 merge하는 스크립트.

Usage:
   MergeBamUnbam.py [--work_dir=<dir>] [--fastq_dir=<dir>] [--aligned_dir=<dir>] [--out_dir=<dir>] [--ref_genome=<dir>] [--picard_path=<dir>]

Options:
  --work_dir=<dir>       작업 디렉토리. (지정하지 않으면 현재 작업 디렉토리)
  --fastq_dir=<dir>      FASTQ 파일들이 위치한 디렉토리. (지정하지 않으면 work_dir/trimmed 있으면 trimmed, 없으면 work_dir/fastq)
  --aligned_dir=<dir>    정렬된 BAM 파일들이 위치한 디렉토리. (지정하지 않으면 ./star_alignment)
  --out_dir=<dir>        결과 파일들이 저장될 출력 디렉토리. (지정하지 않으면 work_dir/merged_bam)
  --ref_genome=<dir>     fasta 레퍼런스 파일이 있는 디렉토리. (반드시 ref_genome/genome.fa 파일이 존재해야 함)
  --picard_path=<dir>    Picard 실행 파일 (jar)의 경로.
"""

import os
import sys
import logging
import subprocess
import TrimmingFastq as tf
from docopt import docopt
from pytz import timezone
from datetime import datetime




def make_logger(log_dir, name="app", consoleset=True, streamset=True) -> logging.Logger:
    """
    Logger 생성 함수.
      - 콘솔과 파일에 로그를 남깁니다.
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


class MergeBamUnbam:
    logger = None

    @staticmethod
    def check_dir(dir_path: str, desc: str) -> None:
        """디렉토리 존재 여부 확인"""
        if not os.path.isdir(dir_path):
            MergeBamUnbam.logger.error(f"{desc} 디렉토리가 존재하지 않습니다: {dir_path}")
            sys.exit(1)
        MergeBamUnbam.logger.info(f"{desc} 디렉토리 확인: {dir_path}")

    @staticmethod
    def check_file(file_path: str, desc: str) -> None:
        """파일 존재 여부 확인"""
        if not os.path.isfile(file_path):
            MergeBamUnbam.logger.error(f"{desc} 파일이 존재하지 않습니다: {file_path}")
            sys.exit(1)
        MergeBamUnbam.logger.info(f"{desc} 파일 확인: {file_path}")

    @staticmethod
    def get_samples(aligned_dir: str) -> list:
        """
        aligned_dir 내의 파일 이름이 {SAMPLE}_aligned_queryname_sorted.bam 인 파일들에서
        샘플 이름({SAMPLE})을 추출합니다.
        """
        bam_files = [f for f in os.listdir(aligned_dir) if f.endswith("_aligned_queryname_sorted.bam")]
        if not bam_files:
            MergeBamUnbam.logger.error("aligned_dir에서 정렬된 BAM 파일을 찾을 수 없습니다.")
            sys.exit(1)
        samples = [f.replace("_aligned_queryname_sorted.bam", "") for f in bam_files]
        MergeBamUnbam.logger.info(f"추출된 샘플 리스트: {samples}")
        return samples

    @staticmethod
    def generate_unaligned_bam(sample: str, fastq_dir: str, unaligned_dir: str, picard_path: str) -> str:
        """
        Picard FastqToSam 명령어를 이용해 FASTQ 파일로부터 unaligned BAM 파일을 생성합니다.
        trimmed 결과가 있다면 _1_paired.fastq.gz, _2_paired.fastq.gz 파일을 사용하고,
        없으면 원시 FASTQ 파일을 사용합니다.
        출력 파일: {unaligned_dir}/{sample}_unaligned.bam
        """
        # 먼저 trimmed FASTQ 파일이 존재하는지 확인
        paired_f1 = os.path.join(fastq_dir, f"{sample}_1_paired.fastq.gz")
        paired_f2 = os.path.join(fastq_dir, f"{sample}_2_paired.fastq.gz")
        if os.path.isfile(paired_f1) and os.path.isfile(paired_f2):
            f1 = paired_f1
            f2 = paired_f2
            MergeBamUnbam.logger.info(f"{sample}: trimmed FASTQ 파일 사용.")
        else:
            f1 = os.path.join(fastq_dir, f"{sample}_1.fastq.gz")
            f2 = os.path.join(fastq_dir, f"{sample}_2.fastq.gz")
            MergeBamUnbam.logger.info(f"{sample}: 원시 FASTQ 파일 사용.")
        
        MergeBamUnbam.check_file(f1, f"FASTQ (forward) for {sample}")
        MergeBamUnbam.check_file(f2, f"FASTQ (reverse) for {sample}")
        
        os.makedirs(unaligned_dir, exist_ok=True)
        out_bam = os.path.join(unaligned_dir, f"{sample}_unaligned.bam")
        cmd = (
            f"java -jar {picard_path} FastqToSam "
            f"F1={f1} F2={f2} O={out_bam} SM={sample}"
        )
        MergeBamUnbam.logger.info(f"{sample}: FastqToSam 실행 중...")
        try:
            subprocess.run(cmd, shell=True, check=True)
            MergeBamUnbam.logger.info(f"{sample}: unaligned BAM 생성 완료: {out_bam}")
        except subprocess.CalledProcessError as e:
            MergeBamUnbam.logger.error(f"{sample}: FastqToSam 실행 실패: {e}")
            sys.exit(1)
        return out_bam

    @staticmethod
    def merge_bams(sample: str, aligned_dir: str, unaligned_bam: str, out_dir: str, ref_fa: str, picard_path: str) -> None:
        """
        Picard MergeBamAlignment 명령어를 이용해 aligned BAM과 unaligned BAM을 병합합니다.
        출력 파일: {out_dir}/{sample}_merged.bam
        """
        input_aligned = os.path.join(aligned_dir, f"{sample}_aligned_queryname_sorted.bam")
        MergeBamUnbam.check_file(input_aligned, f"Aligned BAM for {sample}")
        MergeBamUnbam.check_file(unaligned_bam, f"Unaligned BAM for {sample}")
        
        os.makedirs(out_dir, exist_ok=True)
        output_merged = os.path.join(out_dir, f"{sample}_merged.bam")
        
        cmd = (
            f"java -jar {picard_path} MergeBamAlignment "
            f"ALIGNED={input_aligned} "
            f"UNMAPPED={unaligned_bam} "
            f"O={output_merged} "
            f"R={ref_fa} "
            "CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT"
        )
        MergeBamUnbam.logger.info(f"{sample}: MergeBamAlignment 실행 중...")
        try:
            subprocess.run(cmd, shell=True, check=True)
            MergeBamUnbam.logger.info(f"{sample}: MergeBamAlignment 완료: {output_merged}")
        except subprocess.CalledProcessError as e:
            MergeBamUnbam.logger.error(f"{sample}: MergeBamAlignment 실행 실패: {e}")

    @staticmethod
    def main() -> None:
        arguments = docopt(__doc__)
        # 작업 디렉토리 (없으면 현재 작업 디렉토리)
        work_dir = arguments.get("--work_dir") or os.getcwd()
        # FASTQ 파일들이 위치한 디렉토리; 사용자가 지정하지 않으면,
        # 먼저 work_dir/trimmed 폴더가 있으면 그곳을, 없으면 work_dir/fastq 폴더를 사용합니다.
        fastq_dir = arguments.get("--fastq_dir")
        if not fastq_dir:
            trimmed_dir = os.path.join(work_dir, "trimmed")
            if os.path.isdir(trimmed_dir):
                fastq_dir = trimmed_dir
            else:
                fastq_dir = os.path.join(work_dir, "fastq")

        ref_genome_dir = arguments.get("--ref_genome")
        if not ref_genome_dir:
            ref_genome_dir = "/data/reference/human/fasta"
            # 또는 필수 옵션으로 처리 가능
            # print("Error: --ref_genome 옵션이 필요합니다.")
            # sys.exit(1)
        ref_fa = os.path.join(ref_genome_dir, "genome.fa")
        
        # 정렬된 BAM 파일들이 위치한 디렉토리 (없으면 work_dir/star_alignment_sorted_bam)
        aligned_dir = arguments.get("--aligned_dir") or os.path.join(work_dir, "star_alignment_sorted_bam")
                
        # 출력 디렉토리 (없으면 work_dir/merged_bam)
        out_dir = arguments.get("--out_dir") or os.path.join(work_dir, "merged_bam")
        # unaligned BAM을 저장할 임시 디렉토리 (out_dir 내부에 생성)
        unaligned_dir = os.path.join(out_dir, "unaligned_bam")   
        os.makedirs(unaligned_dir, exist_ok=True)    
        
        picard_path = arguments.get("--picard_path")
        if not picard_path:
            print("Error: --picard_path 옵션은 필수입니다.")
            sys.exit(1)
        
        # 각 디렉토리 존재 여부 확인
        MergeBamUnbam.check_dir(fastq_dir, "FASTQ")
        MergeBamUnbam.check_dir(aligned_dir, "Aligned BAM")
        MergeBamUnbam.check_dir(work_dir, "작업")
        MergeBamUnbam.check_dir(ref_genome_dir, "Reference genome")
        MergeBamUnbam.check_file(ref_fa, "Reference FASTA")
        MergeBamUnbam.check_file(picard_path, "Picard")
        tf.TrimmingFastq.check_java(MergeBamUnbam.logger)
    

        # logger 초기화 후 메시지 출력
        MergeBamUnbam.logger.info("MergeBamUnbam 프로세스 시작.")
        if os.path.isdir(os.path.join(work_dir, "trimmed")):
            MergeBamUnbam.logger.info("trimmed 폴더가 존재하여 trimmed FASTQ 파일을 사용합니다.")
        else:
            MergeBamUnbam.logger.info("trimmed 폴더가 없어 원시 FASTQ 파일을 사용합니다.")
        
        # aligned_dir에서 샘플 리스트 자동 추출
        samples = MergeBamUnbam.get_samples(aligned_dir)
        for sample in samples:
            # FASTQ 파일로부터 unaligned BAM 생성 (덮어쓰기)
            unaligned_bam = MergeBamUnbam.generate_unaligned_bam(sample, fastq_dir, unaligned_dir, picard_path)
            # aligned BAM과 unaligned BAM을 merge하여 최종 merged BAM 생성
            MergeBamUnbam.merge_bams(sample, aligned_dir, unaligned_bam, out_dir, ref_fa, picard_path)
        
        MergeBamUnbam.logger.info("모든 샘플 병합 완료.")

if __name__ == "__main__":
    MergeBamUnbam.logger = make_logger(os.getcwd(), name="MergeBamUnbam")
    MergeBamUnbam.main()


