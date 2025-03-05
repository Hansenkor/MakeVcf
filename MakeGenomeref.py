#gtf file was downloaded from GENCODE; same release to 10x providing gtf, but that of 10x is modifed
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/file_name from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/
#primary assembly annotation should be used, as GATK4 does not support patched annotation: https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19

#fasta file from Broad v0 google cloud: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false

#GRCh37.p13 and p14 difference; release year; p13 used as reference: https://www.ncbi.nlm.nih.gov/grc/human
#genome reference including decoy should be used, but hg19 is too old and hg38 is seem to contain decoy


#conda install bioconda::star
"""
    date: 2025.02.28
    object: 3_make_genomeref.py (STAR genome index 생성 스크립트)
    author: hansen

This script requires that Java is installed and that Trimmomatic is accessible via the provided jar file.
Usage:
  3_MakeGenomeref.py [--ref_dir=<dir>] [--work_dir=<dir>] [--out_dir=<dir>] [--threads=<n>] [--sjdb_overhang=<n>]

Options:
  --ref_dir=<dir>       작업 디렉토리. [default: /data/reference/human]
  --work_dir=<dir>      작업 디렉토리. 지정되지 않으면 현재 디렉토리가 사용됩니다.
  --out_dir=<dir>       결과가 저장될 출력 디렉토리. 지정되지 않으면 ./genome_ref  생성.
  --threads=<n>         사용할 스레드 수. [default: 20]
  --sjdb_overhang=<n>   STAR의 sjdbOverhang 옵션. [default: 100]
"""


import os
import subprocess
import sys
import shutil
import logging
from docopt import docopt


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

class MakeGenomeref:
    @staticmethod
    def check_star_installed(logger) -> None:
        """
        STAR가 설치되어 있는지 확인합니다.
        STAR가 PATH에 없거나 버전 확인에 실패하면 에러 메시지를 출력하고 종료합니다.
        """
        if shutil.which("STAR") is None:
            if logger:
                logger.error("STAR가 설치되어 있지 않거나 PATH에 없습니다. STAR를 설치해주세요.")
            else:
                sys.stderr.write("STAR가 설치되어 있지 않거나 PATH에 없습니다. STAR를 설치해주세요.\n")
            sys.exit(1)
        else:
            try:
                subprocess.run(
                    ["STAR", "--version"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    check=True
                )
            except Exception:
                if logger:
                    logger.error("STAR 설치 확인 중 에러가 발생했습니다.")
                else:
                    sys.stderr.write("STAR 설치 확인 중 에러가 발생했습니다.\n")
                sys.exit(1)
                



    @staticmethod
    def check_ref_dir(ref_dir: str) -> None:
        """REF 디렉토리가 존재하는지 확인합니다."""
        if not os.path.isdir(ref_dir):
            sys.stderr.write(f"REF 디렉토리가 존재하지 않습니다: {ref_dir}\n"
                            "ref_dir 경로를 확인해 주세요.\n")
            sys.exit(1)

    @staticmethod
    def make_genome_index(ref_dir: str, work_dir: str, out_dir:str, threads: int, sjdb_overhang: int) -> None:
        """
        STAR의 genomeGenerate 모드를 실행하여 genome index를 생성합니다.
        
        Parameters:
          ref_dir (str): 참조 유전체 파일들이 위치한 디렉토리 (디폴트 "/data/reference/human")
          out_dir (str): STAR genome index가 저장될 디렉토리. (디폴트 "./genome_reference")
          threads (int): STAR 실행 시 사용할 스레드 수. (디폴트 20)
          sjdb_overhang (int): STAR의 sjdbOverhang 옵션 값. (디폴트 100)
        """

        # 만약 out_dir이 지정되어 있다면 그 디렉토리를 사용하고,
        # 지정되어 있지 않다면 work_dir/genome_ref를 기본값으로 사용합니다.
        if out_dir:
            final_out_dir = out_dir
            print(f"사용자 지정 out_dir: {final_out_dir}")
        else:
            final_out_dir = os.path.join(work_dir, "genome_ref")
            print("사용자가 out_dir을 지정하지 않았으므로 기본값을 사용합니다:", final_out_dir)
        
        # final_out_dir이 존재하지 않으면 생성
        if not os.path.isdir(final_out_dir):
            os.makedirs(final_out_dir, exist_ok=True)
            print(f"{final_out_dir} 이(가) 존재하지 않아 새로 생성되었습니다.")
        else:
            print(f"{final_out_dir} 이(가) 이미 존재합니다.")
        
        # 작업 디렉토리를 final_out_dir로 변경
        os.chdir(final_out_dir)
        
        
        # STAR genomeGenerate 명령어 구성
        
        command = (
            f"STAR --runMode genomeGenerate --genomeDir {final_out_dir} --runThreadN {threads} "
            f"--genomeFastaFiles {ref_dir}/fasta/genome.fa "
            f"--sjdbGTFfile {ref_dir}/genes/genes.gtf "
            f"--sjdbOverhang {sjdb_overhang}"
        )
        print("Running STAR genomeGenerate with the following command:")

        print(command)
        
        try:
            subprocess.run(command, shell=True, check=True)
            print("STAR genome index generation completed successfully.")
        except subprocess.CalledProcessError:
            sys.stderr.write("Error occurred while running STAR genomeGenerate.\n")
            sys.exit(1)


    @staticmethod
    def main() -> None:
        arguments = docopt(__doc__)
        ref_dir = arguments["--ref_dir"]
        work_dir = arguments["--work_dir"]
        out_dir = arguments["--out_dir"]
        threads = int(arguments["--threads"])
        sjdb_overhang = int(arguments["--sjdb_overhang"])

        # work_dir 기본값 처리: 지정되지 않으면 현재 작업 디렉토리 사용
        if not work_dir:
            work_dir = os.getcwd()
        
        # 필수 도구 및 파일 존재 여부 확인
        MakeGenomeref.check_star_installed()
        MakeGenomeref.check_ref_dir(ref_dir)
        
        MakeGenomeref.make_genome_index(ref_dir, work_dir, out_dir, threads, sjdb_overhang)


if __name__ == "__main__":
    MakeGenomeref.main()
