"""
    date: 2025.02.28
    object: 2_trimming.py (FASTQ 파일들에 대해 Trimmomatic trimming 분석을 자동화하는 스크립트)
    author: hansen

This script requires that Java is installed and that Trimmomatic is accessible via the provided jar file.
Usage:
  2_trimming.py --work_dir=<dir> --in_dir=<dir> [--out_dir=<dir>] [--jar_path=<jar>] [--threads=<n>] [--adapter_fa=<path>] [--smm=<val>] [--pal=<val>] [--thr=<val>] --samples=<list> [--lead=<val>] [--trail=<val>] [--slidewin=<val>] [--minlen=<val>]

Options:
  --work_dir=<dir>      작업 디렉토리.
  --in_dir=<dir>        입력 FASTQ 파일들이 위치한 디렉토리.
  --out_dir=<dir>       결과가 저장될 출력 디렉토리. 지정되지 않으면 workdif/trimmed 생성됨.
  --jar_path=<jar>      Trimmomatic jar 파일의 경로. [default: /data/share/Trimmomatic-0.39/trimmomatic-0.39.jar]
  --threads=<n>         사용할 스레드 수. [default: 30]
  --adapter_fa=<path>   ILLUMINACLIP에 사용할 adapter 파일에 얼마나 서열을  [default: /data/share/Trimmomatic-0.39/adapters/TruSeq3-PE.fa]
  --smm=<val>           seed mismatches로 허용가능한 seed 서열의 mismatch 개수. [default: 2]
  --pal=<val>           palindrome clip threshold로 paired end에 해당하는 경우 target read 앞뒤에 존재하는 adapter 서열을 역상보적으로 삭제 [default: 30] 
  --thr=<val>           simple clip threshold로 설정기준에 부합하는 충분히 정확한 match가 확인되면 적절히 clipping함. [default: 10] 
  --samples=<list>      샘플 리스트 (따옴표(" ") 내에 공백으로 구분된 샘플 이름들). ex: "A-1 A-2 B-1 B-2 B-3 B-4 B-5 C-1 C-2"
  --lead=<val>          LEADING 옵션 값. [default: 3]
  --trail=<val>         TRAILING 옵션 값. [default: 3]
  --slidewin=<val>      SLIDINGWINDOW 옵션 값 (window:quality). [default: 4:15]
  --minlen=<val>        MINLEN 옵션 값. [default: 36]

  <sample>_1.fastq.gz <sample>_2.fastq.gz [ex : A-1_1.fastq.gz , A-1_2.fastq.gz ] 이런식으로 파일명이 되어있어야 --samples의 샘플리스트를 읽을 수 있습니다.
"""


import os
import sys
import subprocess
import shutil
from docopt import docopt# conda install -c conda-forge docopt 해야 

class TrimmingFastq:
    @staticmethod
    def check_java():
        """Java가 설치되어 있는지 확인합니다."""
        try:
            subprocess.run(["java", "-version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        except subprocess.CalledProcessError:
            sys.stderr.write("Java가 설치되어 있지 않거나 실행할 수 없습니다.\n")
            sys.exit(1)

    @staticmethod
    def check_jar(jar_path: str) -> None:
        """Trimmomatic jar 파일이 존재하는지 확인합니다."""
        if not os.path.isfile(jar_path):
            sys.stderr.write(f"Jar 파일이 존재하지 않습니다: {jar_path}\n"
                            "Jar 파일의 경로를 확인하거나, 올바른 경로를 지정해 주세요.\n")
            sys.exit(1)

    
    @staticmethod
    def check_in_dir(in_dir: str) -> None:
        """입력 디렉토리가 존재하는지 확인합니다."""
        if not os.path.isdir(in_dir):
            sys.stderr.write(f"입력 디렉토리가 존재하지 않습니다: {in_dir}\n"
                            "in_dir 경로를 확인해 주세요.\n")
            sys.exit(1)

    @staticmethod
    def check_workdir(work_dir: str) -> None:
        """작업 디렉토리가 존재하는지 확인합니다."""
        if not os.path.isdir(work_dir):
            sys.stderr.write(f"작업 디렉토리가 존재하지 않습니다: {work_dir}\n"
                            "work_dir 경로를 확인해 주세요.\n")
            sys.exit(1)

    @staticmethod
    def check_sample_files(in_dir: str, samples: str) -> None:
        """샘플 목록에 따른 FASTQ 파일들이 존재하는지 미리 확인합니다."""
        missing_samples = []
        for sample in samples.split():
            in_file1 = os.path.join(in_dir, f"{sample}_1.fastq.gz")
            # 단일 파일 입력 여부도 확인할 수 있지만, 여기서는 기본적으로 _1 파일의 존재 여부만 확인
            if not os.path.isfile(in_file1):
                missing_samples.append(sample)
        if missing_samples:
            sys.stderr.write(f"다음 샘플의 FASTQ 파일이 존재하지 않습니다: {', '.join(missing_samples)}\n")
            sys.exit(1)

    @staticmethod
    def check_adapter(adapter_fa: str) -> None:
        """Adapter 파일이 존재하는지 확인합니다."""
        if not os.path.isfile(adapter_fa):
            sys.stderr.write(f"Adapter 파일이 존재하지 않습니다: {adapter_fa}\n"
                            "Adapter 파일 경로를 확인하거나, 올바른 경로를 지정해 주세요.\n")
            sys.exit(1)
            
    @staticmethod
    def run_trimmomatic(work_dir: str, in_dir: str, out_dir: str, jar_path: str, threads: int, adapter_fa: str, smm: int, pal: int, thr: int, samples: str, lead: int, trail: int, slidewin: int, minlen: int) -> None:
        """
        각 샘플에 대해 Trimmomatic을 실행합니다.
        
        각 샘플의 입력 파일은 다음과 같이 사용합니다:
        - Paired End (PE): {in_dir}/{SAMPLE}_1.fastq.gz 와 {in_dir}/{SAMPLE}_2.fastq.gz가 모두 존재할 경우,
            출력 파일은 {out_dir}/{SAMPLE}_1_paired.fastq.gz, {out_dir}/{SAMPLE}_1_unpaired.fastq.gz,
                    {out_dir}/{SAMPLE}_2_paired.fastq.gz, {out_dir}/{SAMPLE}_2_unpaired.fastq.gz 로 생성됩니다.
        - Single End (SE): {in_dir}/{SAMPLE}_1.fastq.gz 만 존재할 경우,
            출력 파일은 {out_dir}/{SAMPLE}_trimmed.fastq.gz 로 생성됩니다.
        """
        # out_dir이 지정되지 않았거나 빈 문자열이면 workdir/trimmed를 사용합니다.
        if not out_dir:
            out_dir = os.path.join(work_dir, "trimmed")
        os.makedirs(out_dir, exist_ok=True)

        # 작업 디렉토리로 이동
        try:
            os.chdir(work_dir)
        except Exception as exc:
            sys.stderr.write(f"작업 디렉토리 변경 중 에러 발생: {exc}\n")
            sys.exit(1)

        sample_list = samples.split()
        for sample in sample_list:
            in_file1 = os.path.join(in_dir, f"{sample}_1.fastq.gz")
            in_file2 = os.path.join(in_dir, f"{sample}_2.fastq.gz")
            
            # SE trimming: _1.fastq.gz만 있을 경우
            if os.path.isfile(in_file1) and not os.path.isfile(in_file2):
                out_file = os.path.join(out_dir, f"{sample}_trimmed.fastq.gz")
                command = (
                    f"java -jar {jar_path} SE -threads {threads} -phred33 "
                    f"{in_file1} {out_file} "
                    f"ILLUMINACLIP:{adapter_fa}:{smm}:0:{thr} LEADING:{lead} TRAILING:{trail} SLIDINGWINDOW:{slidewin} MINLEN:{minlen}"
                )
                sys.stdout.write(f"Running Trimmomatic (SE) for sample {sample}...\n")
                try:
                    subprocess.run(command, shell=True, check=True)
                    sys.stdout.write(f"{sample} (SE) has processed!\n")
                except subprocess.CalledProcessError as error:
                    sys.stderr.write(f"Trimmomatic 실행 중 에러 발생 ({sample}, SE): {error}\n")
                    continue

            # PE trimming: _1.fastq.gz와 _2.fastq.gz 모두 존재할 경우
            elif os.path.isfile(in_file1) and os.path.isfile(in_file2):
                out1_paired = os.path.join(out_dir, f"{sample}_1_paired.fastq.gz")
                out1_unpaired = os.path.join(out_dir, f"{sample}_1_unpaired.fastq.gz")
                out2_paired = os.path.join(out_dir, f"{sample}_2_paired.fastq.gz")
                out2_unpaired = os.path.join(out_dir, f"{sample}_2_unpaired.fastq.gz")
                
                command = (
                    f"java -jar {jar_path} PE -threads {threads} -phred33 "
                    f"{in_file1} {in_file2} "
                    f"{out1_paired} {out1_unpaired} {out2_paired} {out2_unpaired} "
                    f"ILLUMINACLIP:{adapter_fa}:{smm}:{pal}:{thr} LEADING:{lead} TRAILING:{trail} SLIDINGWINDOW:{slidewin} MINLEN:{minlen}"
                )
                sys.stdout.write(f"Running Trimmomatic (PE) for sample {sample}...\n")
                try:
                    subprocess.run(command, shell=True, check=True)
                    sys.stdout.write(f"{sample} (PE) has processed!\n")
                except subprocess.CalledProcessError as error:
                    sys.stderr.write(f"Trimmomatic 실행 중 에러 발생 ({sample}, PE): {error}\n")
                    continue

            else:
                sys.stderr.write(f"{sample}의 입력 파일이 올바르지 않습니다: {in_file1} 혹은 {in_file2}\n")
                continue

        sys.stdout.write("All done!\n")


    """
    --SE                  Single-end read input. Default input choice is single-
                        end if nothing is specified
    --PE                  Paired-end read input. Must have the exact same file
                        name and end with _F for the forward read and _R for
                        the reverse read


    Trimmomatic options:
    The order the options will be run are: ILLUMINACLIP, LEADING,
    TRAILING, SLIDINGWINDOW, MINLEN

    --clip=ILLUMINACLIP
                        ILLUMINACLIP options. MiSeq & HiSeq usually
                        TruSeq3.fa; GAII usually TruSeq2.fa. Default is
                        ILLUMINACLIP:TruSeq3-SE.fa:2:30:10. Usage:
                        --clip=<adapterseqs>:<seed mismatches>:<palindrome
                        clip threshold>:<simple clip threshold>
    --lead=LEADING      Set the minimun quality required to keep a base.
                        Default is LEADING=3. Usage: --lead=<quality>
    --trail=TRAILING    Set the minimum quality required to keep a base.
                        Default is TRAILING=3. Usage: --trail=<quality>
    --slidewin=SLIDINGWINDOW
                        SLIDINGWINDOW options. Default is SLIDINGWINDOW:4:15.
                        Usage: --slidewin=<window_size>:<required_quality>

    HTSeq options:
    --stranded=STRANDED
                        Stranded options: yes, no, reverse. Default is
                        --stranded=reverse. Usage: --stranded=yes/no/reverse
    --order=ORDER       Order options: name, pos. Usage: --order=name/pos.
    --minqual=MINQUAL   Skip all reads with quality lower than the given
                        value. Default is --minqual=10. Usage:
                        --minqual=<value>
    --idattr=IDATTR     Feature ID from the GTF file to identify counts in the
                        output table Default is --idattr=gene_id. Usage:
                        --idattr=<id attribute>
    --mode=MODE         Mode to handle reads overlapping more than one
                        feature. Default is --mode=union. Usage: --mode=union
                        /intersection-strict/intersection-nonempty

    """
    @staticmethod
    def main() -> None:
        arguments = docopt(__doc__)
        
        work_dir = arguments["--work_dir"]
        in_dir = arguments["--in_dir"]
        out_dir = arguments["--out_dir"]
        jar_path = arguments["--jar_path"]
        threads = int(arguments["--threads"])
        adapter_fa = arguments["--adapter_fa"]
        smm = int(arguments["--smm"])
        pal = int(arguments["--pal"])
        thr = int(arguments["--thr"])
        samples = arguments["--samples"]
        lead = int(arguments["--lead"])
        trail = int(arguments["--trail"])
        slidewin = arguments["--slidewin"]
        minlen = int(arguments["--minlen"])
        
        # 필수 도구 및 파일 존재 여부 확인
        TrimmingFastq.check_java()
        TrimmingFastq.check_adapter(adapter_fa)
        TrimmingFastq.check_jar(jar_path)
        TrimmingFastq.check_workdir(work_dir)
        TrimmingFastq.check_in_dir(in_dir)
        TrimmingFastq.check_sample_files(in_dir, samples)
        
        TrimmingFastq.run_trimmomatic(work_dir, in_dir, out_dir, jar_path, threads, adapter_fa, smm, pal, thr, samples, lead, trail, slidewin, minlen)


if __name__ == "__main__":
    TrimmingFastq.main()