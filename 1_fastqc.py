# 필요한 패키지 설치 (conda 환경에서 미리 설치되지 않은 경우 자동 설치 시도)
# conda install -c bioconda gatk4 
# conda install -c bioconda fastqc 
# conda install -c bioconda samtools bcftools
#
"""
1_fastqc.py: FASTQ 파일들에 대해 fastqc 분석을 자동화하는 스크립트.

This script requires that fastqc is installed via conda:
    conda install -c bioconda fastqc

Usage:
  fastqc_run.py --workdir=<dir> [--threads=<n>] [--outdir=<dir>]

Options:
  --workdir=<dir>       FASTQ 파일들이 위치한 작업 디렉토리.
  --threads=<n>         사용할 스레드 수. [default: 16]
  --outdir=<dir>        결과가 저장될 출력 디렉토리. 지정하지 않으면 workdir/fastqc가 사용됨.
"""


import glob
import os
import shutil
import subprocess
import sys
from docopt import docopt



def check_and_install_fastqc() -> None:
    """fastqc 명령어가 없으면 conda를 통해 설치를 시도합니다."""
    if shutil.which("fastqc") is None:
        sys.stderr.write("fastqc 명령어를 찾을 수 없습니다. 설치를 시도합니다...\n")
        try:
            subprocess.run("conda install -c bioconda fastqc -y", shell=True, check=True)
        except subprocess.CalledProcessError as error:
            sys.stderr.write(f"fastqc 설치 실패: {error}\n")
            sys.exit(1)
        if shutil.which("fastqc") is None:
            sys.stderr.write("fastqc 설치 후에도 fastqc 명령어를 찾을 수 없습니다.\n")
            sys.exit(1)
        sys.stderr.write("fastqc 설치 완료.\n")


def run_fastqc(workdir: str, threads: int, outdir: str) -> None:
    """하위 폴더를 포함한 FASTQ 파일들에 대해 fastqc를 실행하는 함수.

    Args:
        workdir: FASTQ 파일들이 있는 최상위 작업 디렉토리.
        threads: fastqc 실행 시 사용할 스레드 수.
        outdir: fastqc 결과가 저장될 디렉토리.
    """
    # fastqc가 설치되어 있는지 확인하고, 없으면 설치 시도
    check_and_install_fastqc()

    # outdir가 지정되지 않은 경우 기본값 설정
    if outdir is None:
        outdir = os.path.join(workdir, "fastqc")

    # 작업 디렉토리로 이동
    try:
        os.chdir(workdir)
    except Exception as exc:
        sys.stderr.write(f"작업 디렉토리 변경 중 에러 발생: {exc}\n")
        sys.exit(1)

    # 출력 디렉토리 생성 (존재하지 않으면)
    os.makedirs(outdir, exist_ok=True)

    # workdir 내 하위 폴더들에서 모든 *.fastq.gz 파일들을 재귀적으로 검색
    fastq_files = glob.glob(os.path.join(workdir, "**", "*.fastq.gz"), recursive=True)
    if not fastq_files:
        sys.stderr.write("FASTQ 파일을 찾을 수 없습니다.\n")
        sys.exit(1)

    # fastqc 명령어 구성: fastqc -o {outdir} -t {threads} <파일리스트>
    command = f"fastqc -o {outdir} -t {threads} " + " ".join(fastq_files)

    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as error:
        sys.stderr.write(f"fastqc 실행 중 에러 발생: {error}\n")
        sys.exit(1)

    print("All done!")


def main() -> None:
    """메인 함수: 인자 파싱 및 fastqc 실행."""
    arguments = docopt(__doc__)

    workdir = arguments["--workdir"]
    outdir = arguments["--outdir"]

    try:
        threads = int(arguments["--threads"])
    except ValueError:
        sys.stderr.write("Error: --threads 값은 정수여야 합니다.\n")
        sys.exit(1)

    run_fastqc(workdir, threads, outdir)


if __name__ == "__main__":
    """이 조건문은 현재 모듈이 직접 실행되고 있는지 확인합니다."""
    main()